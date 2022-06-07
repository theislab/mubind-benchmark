import copy

import numpy as np
import pandas as pd
import torch
import torch.optim as topti
import torch.utils.data as tdata
import time
import multibind as mb


def calculate_enrichment(data, approx=True, cols=[0, 1]):
    data["p0"] = data[cols[0]] / np.sum(data[cols[0]])
    data["p1"] = data[cols[1]] / np.sum(data[cols[1]])
    data["enr"] = data["p1"] / data["p0"]
    if approx:
        data["enr_approx"] = np.where(data["p0"] == 0, data["p1"] / (data["p0"] + 1e-06), data["enr"])
    return data


def create_datasets(data_file):
    # read data and calculate additional columns
    data = pd.read_csv(data_file, sep="\t", header=None)
    data.columns = ["seq", 0, 1]
    data = calculate_enrichment(data)
    # divide in train and test data
    test_dataframe = data.sample(frac=0.001)
    train_dataframe = data.drop(test_dataframe.index)
    test_dataframe.index = range(len(test_dataframe))
    train_dataframe.index = range(len(train_dataframe))
    # create datasets and dataloaders
    train_data = mb.datasets.SelexDataset(data_frame=train_dataframe)
    train_loader = tdata.DataLoader(dataset=train_data, batch_size=256, shuffle=True)
    test_data = mb.datasets.SelexDataset(data_frame=test_dataframe)
    test_loader = tdata.DataLoader(dataset=test_data, batch_size=1, shuffle=False)
    return train_loader, test_loader


def test_network(model, dataloader, device):
    all_seqs = []
    all_targets = np.zeros((len(dataloader.dataset), dataloader.dataset.n_rounds+1), dtype=np.float32)
    all_preds = np.zeros((len(dataloader.dataset), dataloader.dataset.n_rounds+1), dtype=np.float32)
    position = 0
    with torch.no_grad():  # we don't need gradients in the testing phase
        for i, batch in enumerate(dataloader):
            # Get a batch and potentially send it to GPU memory.
            mononuc = batch["mononuc"].to(device)
            b = batch["batch"].to(device) if "batch" in batch else None
            rounds = batch["rounds"].to(device) if "rounds" in batch else None
            seqlen = batch["seqlen"] if "seqlen" in batch else None
            countsum = batch["countsum"].to(device) if "countsum" in batch else None
            seq = batch["seq"] if "seq" in batch else None

            inputs = (mononuc, b, seqlen, countsum)
            output = model(inputs)

            # print(output.squeeze().cpu().detach().numpy().shape, 'vs', len(seq))
            all_preds[position:(position + len(seq)), :] = output.squeeze().cpu().detach().numpy()
            all_targets[position:(position + len(seq)), :] = rounds.cpu().detach().numpy()
            all_seqs.extend(seq)
            position += len(seq)
    return all_seqs, all_targets, all_preds


# loss_history = []


# if early_stopping is positive, training is stopped if over the length of early_stopping no improvement happened or
# num_epochs is reached.
def train_network(
    model,
    train_dataloader,
    device,
    optimiser,
    criterion,
    num_epochs=15,
    early_stopping=-1,
    dirichlet_regularization=0,
    exp_max=40,  # if this value is negative, the exponential barrier will not be used.
    log_each=None,
    verbose=0,
):
    # global loss_history
    loss_history = []
    best_loss = None
    best_epoch = -1
    if verbose != 0:
        print('optimizing using', str(type(optimiser)), 'and', str(type(criterion)))

    for f in ['lr', 'weight_decay']:
        if f in optimiser.param_groups[0]:
            if verbose != 0:
                print('%s=' % f, optimiser.param_groups[0][f], end=', ')

    if verbose != 0:
        print('dir weight=', dirichlet_regularization)

    is_LBFGS = 'LBFGS' in str(optimiser)

    t0 = time.time()
    for epoch in range(num_epochs):
        running_loss = 0
        for i, batch in enumerate(train_dataloader):
            # Get a batch and potentially send it to GPU memory.
            mononuc = batch["mononuc"].to(device)
            b = batch["batch"].to(device) if "batch" in batch else None
            batch["target"].to(device) if "target" in batch else None
            rounds = batch["rounds"].to(device) if "rounds" in batch else None
            batch["is_count_data"] if "is_count_data" in batch else None
            seqlen = batch["seqlen"] if "seqlen" in batch else None
            countsum = batch["countsum"].to(device) if "countsum" in batch else None

            inputs = (mononuc, b, seqlen, countsum)
            # PyTorch calculates gradients by accumulating contributions to them (useful for
            # RNNs).  Hence we must manully set them to zero before calculating them.
            # print('outputs', rounds)
            # print('rounds', rounds)
            loss = None
            if not is_LBFGS:
                optimiser.zero_grad()
                outputs = model(inputs)  # Forward pass through the network.
                loss = criterion(outputs, rounds) + dirichlet_regularization*model.dirichlet_regularization()
                if exp_max >= 0:
                    loss += model.exp_barrier(exp_max)
                loss.backward()  # Calculate gradients.
                optimiser.step()

            else:
                def closure():
                    optimiser.zero_grad()
                    # this statement here is mandatory to
                    outputs = model(inputs)
                    loss = criterion(outputs, rounds) + dirichlet_regularization*model.dirichlet_regularization()
                    if exp_max >= 0:
                        loss += model.exp_barrier(exp_max)
                    loss.backward() # retain_graph=True)
                    return loss
                loss = optimiser.step(closure)  # Step to minimise the loss according to the gradient.
            running_loss += loss.item()


        loss_final = running_loss / len(train_dataloader)
        if log_each != -1 and (epoch % log_each == 0):
            if verbose != 0:
                print("Epoch: %2d, Loss: %.6f" % (epoch + 1, loss_final),
                      ', best epoch: %i' % best_epoch, 'secs per epoch: %.3f s' % ((time.time() - t0) / max(epoch, 1)))

        if best_loss is None or loss_final < best_loss:
            best_loss = loss_final
            best_epoch = epoch
            model.best_model_state = copy.deepcopy(model.state_dict())
            model.best_loss = best_loss

        # print("Epoch: %2d, Loss: %.3f" % (epoch + 1, running_loss / len(train_dataloader)))
        loss_history.append(loss_final)

        if early_stopping > 0 and epoch >= best_epoch + early_stopping:
            if verbose != 0:
                print("Epoch: %2d, Loss: %.4f" % (epoch + 1, loss_final),
                  ', best epoch: %i' % best_epoch, 'secs per epoch: %.3f s' % ((time.time() - t0) / max(epoch, 1)))
            if verbose != 0:
                print('early stop!')
            break

    model.loss_history += loss_history


def train_iterative(
    train,
    device,
    n_kernels=4,
    w=15,
    # min_w=10,
    # max_w=20,
    num_epochs=100,
    early_stopping=15,
    log_each=10,
    optimize_motif_shift=True,
    show_logo=True,
    optimiser=None,
    criterion=None,
    seed=None,
    lr=0.01,
    ignore_kernel=False,
    weight_decay=0.001,
    stop_at_kernel=None,
    dirichlet_regularization=0,
    verbose=2,
    exp_max=40,
    **kwargs
):

    # model_by_k = {}
    res = []
    loss_history = []
    n_rounds = train.dataset.n_rounds
    n_batches = train.dataset.n_batches
    enr_series = train.dataset.enr_series

    if verbose != 0:
        print('# rounds', n_rounds)
        print('# batches', n_batches)
        print('# enr_series', enr_series)

    if criterion is None:
        criterion = mb.tl.PoissonLoss()

    # color for visualization of history
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628"]

    # for w in range(min_w, max_w, 2):
    # step 1) freeze everything before the current binding mode
    if verbose != 0:
        print("next w", w, type(w))
    # assert False
    model = mb.models.DinucSelex(kernels=[0] + [w] * (n_kernels - 1), n_rounds=n_rounds,
                                 n_batches=n_batches, enr_series=enr_series, **kwargs).to(device)

    # this sets up the seed at the first positoin
    if seed is not None:
        model.set_seed(seed, 1)
        model = model.to(device)


    for i in range(0, n_kernels):
        if verbose != 0:
            print("\nKernel to optimize %i" % i)
            print('\nFreezing kernels')
        for ki in range(n_kernels):
            if verbose != 0:
                print("setting grad status of kernel at %i to %i" % (ki, ki == i))
            mb.tl.update_grad(model, ki, ki == i)
        print('\n')

        if show_logo:
            if verbose != 0:
                print("before kernel optimization.")
            mb.pl.plot_activities(model, train)
            mb.pl.conv_mono(model)
            mb.pl.conv_mono(model, flip=True, log=False)


        next_lr = lr if not isinstance(lr, list) else lr[i]
        next_weight_decay = weight_decay if not isinstance(weight_decay, list) else weight_decay[i]

        next_optimiser = topti.Adam(model.parameters(),
                                    lr=next_lr, weight_decay=next_weight_decay) if optimiser is None else optimiser(model.parameters(), lr=next_lr)

        # mask kernels to avoid using weights from further steps into early ones.
        if ignore_kernel:
            model.ignore_kernel = np.array([0 for i in range(i + 1)] + [1 for i in range(i + 1, n_kernels)])

        if verbose != 0:
            print('kernels mask', model.ignore_kernel)

        # assert False
        mb.tl.train_network(
            model,
            train,
            device,
            next_optimiser,
            criterion,
            num_epochs=num_epochs,
            early_stopping=early_stopping,
            log_each=log_each,
            dirichlet_regularization=dirichlet_regularization,
            exp_max=exp_max,
        )
        # print('next color', colors[i])
        model.loss_color += list(np.repeat(colors[i], len(model.loss_history)))
        # probably here load the state of the best epoch and save
        model.load_state_dict(model.best_model_state)
        k_parms = "%i" % w
        # store model parameters and fit for later visualization
        model_by_k[k_parms] = copy.deepcopy(model)
        # optimizer for left / right flanks
        best_loss = model_by_k[k_parms].best_loss

        if show_logo:
            print("\n##After kernel opt / before shift optim.")
            mb.pl.plot_activities(model, train)
            mb.pl.conv_mono(model)
            mb.pl.conv_mono(model, flip=True, log=False)
            mb.pl.plot_loss(model)

        # print(model_by_k[k_parms].loss_color)
        #######
        # optimize the flanks through +1/-1 shifts
        #######
        n_attempts = 0

        if optimize_motif_shift and i != 0:
            next_loss = None
            while next_loss is None or next_loss < best_loss:
                n_attempts += 1

                if verbose != 0:
                    print(
                        "\noptimize_motif_shift (%s)..." % ("first" if next_loss is None else "again"),
                        end="",
                    )
                model = copy.deepcopy(model)
                best_loss = model.best_loss
                next_color = colors[-(1 if n_attempts % 2 == 0 else -2)]


                all_shifts = []
                for shift in [-3, -2, -1, 1, 2, 3]:
                    if verbose != 0:
                        print('next shift:', shift)
                    model_shift = copy.deepcopy(model)
                    model_shift.loss_history = []
                    model_shift.loss_color = []

                    next_optimiser = topti.Adam(model.parameters(),
                                                lr=next_lr, weight_decay=next_weight_decay) if optimiser is None else optimiser(model.parameters(), lr=next_lr)
                    model_left = mb.tl.train_shift(
                        model_shift,
                        train,
                        kernel_i=i,
                        shift=shift,
                        device=device,
                        num_epochs=num_epochs,
                        early_stopping=early_stopping,
                        log_each=log_each,
                        update_grad_i=i,
                        lr=lr, weight_decay=weight_decay,
                        optimiser=next_optimiser,
                        dirichlet_regularization=dirichlet_regularization,
                        exp_max=exp_max,
                        **kwargs,
                    )
                    model_shift.loss_color += list(np.repeat(next_color, len(model_shift.loss_history)))
                    # print('history left', len(model_left.loss_history))
                    all_shifts.append([shift, model_shift, model_shift.best_loss])
                    # print('\n')

                # for shift, model_shift, loss in all_shifts:
                #     print('shift=%i' % shift, 'loss=%.4f' % loss)
                best = sorted(all_shifts + [[0, model, best_loss]],
                    key=lambda x: x[-1],
                )
                if verbose != 0:
                    print('sorted')
                best_df = pd.DataFrame([[shift, loss] for shift, model_shift, loss in best],
                                       columns=['shift', 'loss'])
                if verbose != 0:
                    print(best_df.sort_values('loss'))
                # for shift, model_shift, loss in best:
                #     print('shift=%i' % shift, 'loss=%.4f' % loss)

                # print('\n history len')
                next_position, next_model, next_loss = best[0]
                if verbose != 0:
                    print('action: %s\n' % next_position)

                if next_position != 0:
                    next_model.loss_history = model.loss_history + next_model.loss_history
                    next_model.loss_color = model.loss_color + next_model.loss_color

                model = copy.deepcopy(next_model)

        # model = model_by_k[k_parms]

        n_feat = sum(
            np.prod(layer.kernel_size)
            for conv in [model.conv_mono, model.conv_di]
            for layer in conv
            if layer is not None
        )
        l_best = model.best_loss

        if show_logo:
            if verbose != 0:
                print("after shift optimz model")
            mb.pl.plot_activities(model, train)
            mb.pl.conv_mono(model)
            mb.pl.conv_mono(model, flip=True, log=False)
            mb.pl.plot_loss(model)
            print("")

        if verbose != 0:
            print('\n\nfinal refinement step (after shift)...')
            print('\nunfreezing all layers for final refinement')
        for ki in range(n_kernels):
            if verbose != 0:
                print("kernel grad (%i) = %i" % (ki, True), sep=', ', end='')
            mb.tl.update_grad(model, ki, ki == i)
        if verbose != 0:
            print('')

        next_optimiser = topti.Adam(model.parameters(),
                                    lr=next_lr, weight_decay=next_weight_decay) if optimiser is None else optimiser(model.parameters(), lr=next_lr)

        # mask kernels to avoid using weights from further steps into early ones.
        if ignore_kernel:
            model.ignore_kernel = np.array([0 for i in range(i + 1)] + [1 for i in range(i + 1, n_kernels)])
        if verbose != 0:
            print('kernels mask', model.ignore_kernel)
        # assert False
        mb.tl.train_network(
            model,
            train,
            device,
            next_optimiser,
            criterion,
            num_epochs=num_epochs,
            early_stopping=early_stopping,
            log_each=log_each,
            dirichlet_regularization=dirichlet_regularization,
        )

        if stop_at_kernel is not None and stop_at_kernel == i:
            break

    # r = [k_parms, w, n_feat, l_best]
    # # print(r)
    # res.append(r)

    return model, model.best_loss


def update_grad(model, position, value):
    if model.conv_mono[position] is not None:
        model.conv_mono[position].weight.requires_grad = value
        # print("mono grad", position, model.conv_mono[position].weight.grad)

    if model.conv_di[position] is not None:
        model.conv_di[position].weight.requires_grad = value
        # print("di grad", position, model.conv_mono[position].weight.grad)

    model.log_activities[position].requires_grad = value
    if not value:
        model.log_activities[position].grad = None
        if model.kernels[position] != 0:
            model.conv_mono[position].weight.grad = None
            model.conv_di[position].weight.grad = None

    # padding required?
    # model.padding[position].weight.requires_grad = valueassert False


def train_shift(
    model,
    train,
    shift=0,
    device=None,
    num_epochs=500,
    early_stopping=15,
    log_each=-1,
    update_grad_i=None,
    kernel_i=None,
    lr=0.01,
    weight_decay=0.001,
    optimiser=None,
    criterion=None,
    dirichlet_regularization=0,
    exp_max=40,
    **kwargs,
):

    # shift mono
    for i, m in enumerate(model.conv_mono):
        if kernel_i is not None and kernel_i != i:
            continue
        if m is None:
            continue
        # update the weight
        if shift >= 1:
            m.weight = torch.nn.Parameter(torch.cat([m.weight[:, :, :, shift:], torch.zeros(1, 1, 4, shift).to(device)], dim=3))
        elif shift <= -1:
            m.weight = torch.nn.Parameter(
                torch.cat(
                    [
                        torch.zeros(1, 1, 4, -shift).to(device),
                        m.weight[:, :, :, :-shift],
                    ],
                    dim=3,
                )
            )
    # shift di
    for i, m in enumerate(model.conv_di):
        if kernel_i is not None and kernel_i != i:
            continue
        if m is None:
            continue
        # update the weight
        if shift >= 1:
            m.weight = torch.nn.Parameter(
                torch.cat([m.weight[:, :, :, shift:], torch.zeros(1, 1, 16, shift).to(device)], dim=3)
            )
        elif shift <= -1:
            m.weight = torch.nn.Parameter(
                torch.cat([torch.zeros(1, 1, 16, -shift).to(device), m.weight[:, :, :, :-shift]], dim=3)
            )

    # requires grad update
    n_kernels = len(model.conv_mono)
    for ki in range(n_kernels):
        update_grad(model, ki, ki == update_grad_i)

    if criterion is None:
        criterion = mb.tl.PoissonLoss()

    mb.tl.train_network(
        model,
        train,
        device,
        optimiser,
        criterion,
        num_epochs=num_epochs,
        early_stopping=early_stopping,
        log_each=log_each,
        dirichlet_regularization=dirichlet_regularization,
        exp_max=exp_max,
    )

    return model


# returns the last loss value if it exists
def get_last_loss_value():
    global loss_history
    if len(loss_history) > 0:
        return loss_history[-1]
    return None


def create_simulated_data(motif="GATA", n_batch=None, n_trials=20000, seqlen=100, batch_sizes=10):
    x2, y2 = mb.datasets.simulate_xy(motif, n_trials=n_trials, seqlen=seqlen, max_mismatches=-1, batch=1)

    # print('skip normalizing...')
    # y2 = ((y2 - y2.min()) / (np.max(y2) - np.min(y2))).astype(np.float32)

    # data = pd.DataFrame({'seq': x1, 'enr_approx': y1})
    data = pd.DataFrame({"seq": x2, "target": y2})
    batch = mb.tl.onehot_covar(np.random.random_integers(0, n_batch - 1, len(y2)))

    assert n_batch == len(batch_sizes)
    batch_mult = np.array(batch_sizes)[(np.argmax(batch, axis=1))]
    # print(batch_mult)
    # df = pd.DataFrame()
    # df['y.before'] = np.array(data['target'])
    # print(np.argmax(batch, axis=1))
    # df['batch'] = np.argmax(batch, axis=1)
    # df['mult'] = batch_mult
    # df['y.after'] = np.array(data['target']) * np.array(batch_mult)
    # print('head...')
    # print(df.head(25))
    # print('tail')
    # print(df.tail(25))

    data["target"] = np.array(data["target"]) * np.array(batch_mult)
    data["is_count_data"] = np.repeat(1, data.shape[0])
    # assert False

    # divide in train and test data -- copied from above, organize differently!
    test_dataframe = data.sample(frac=0.01)
    train_dataframe = data.drop(test_dataframe.index)
    batch_test = np.argmax(batch[test_dataframe.index], axis=1)
    batch_train = np.argmax(batch[train_dataframe.index], axis=1)

    test_dataframe.index = range(len(test_dataframe))
    train_dataframe.index = range(len(train_dataframe))

    n_train = len(train_dataframe.index)
    n_test = len(test_dataframe.index)
    n_train + n_test

    # create datasets and dataloaders
    # print(train_dataframe)
    # print(batch, n_train)
    # print(batch[:n_train])
    train_data = mb.datasets.ChipSeqDataset(data_frame=train_dataframe)
    train_data.batch = batch_train
    train_data.seq = train_dataframe.seq

    train_loader = tdata.DataLoader(dataset=train_data, batch_size=256, shuffle=True)
    test_data = mb.datasets.ChipSeqDataset(data_frame=test_dataframe)
    test_data.batch = batch_test
    test_data.seq = test_dataframe.seq

    test_loader = tdata.DataLoader(dataset=test_data, batch_size=1, shuffle=False)
    return train_loader, test_loader


def create_multi_data(n_chip=100, n_selex=100, n_batch_selex=3):
    # chip seq
    x_chip, y_chip = mb.datasets.gata_remap(n_sample=n_chip) if n_chip > 0 else [[], []]
    df_chip = pd.DataFrame({"seq": x_chip, "target": y_chip})
    df_chip["is_count_data"] = np.repeat(0, df_chip.shape[0])
    df_chip["batch"] = "0_chip"
    # selex
    # seqlen = 10
    # x_selex, y_selex = mb.datasets.simulate_xy('GATA', n_trials=n_selex, seqlen=seqlen, max_mismatches=-1, batch=50)
    # df_selex = pd.DataFrame({'seq': x_selex, 'target': y_selex})
    # df_selex['is_count_data'] = np.repeat(1, df_selex.shape[0])
    # n_batch = n_batch_selex
    # batch sizes. Dataset 2 has many times more reads than Dataset 1
    batch_sizes = [int(5 * 10**i) for i in range(n_batch_selex)]
    train1, test1 = mb.tl.create_simulated_data(
        motif="GATA",
        n_batch=n_batch_selex,
        n_trials=n_selex,
        seqlen=10,
        batch_sizes=batch_sizes,
    )  # multiplier=100)
    # print(train1.dataset.seq, train1.shape)
    # assert False
    df_selex = pd.DataFrame({"seq": train1.dataset.seq, "target": train1.dataset.target})
    df_selex["batch"] = pd.Series(train1.dataset.batch).astype(str).values + "_selex"
    df_selex["is_count_data"] = np.repeat(1, df_selex.shape[0])
    datasets = [df_chip, df_selex]
    data = pd.concat(datasets)

    # pad nucleotides for leaking regions

    data["seqlen"] = data["seq"].str.len()
    # print(data.head())
    # print(data.tail())
    # assert False

    max_seqlen = data["seq"].str.len().max()
    print("max seqlen", max_seqlen)
    max_seqlen = int(max_seqlen)
    data["k"] = "_"
    padded_seq = data["seq"] + data["k"].str.repeat(max_seqlen - data["seqlen"].astype(int))
    data["seq"] = np.where(data["seqlen"] < max_seqlen, padded_seq, data["seq"])
    data["seqlen"] = data["seq"].str.len()
    assert len(set(data["seqlen"])) == 1

    # batches
    # n_batch = len(datasets)
    # batch = []
    # for i, df in enumerate(datasets):
    #     batch.append(np.repeat(i, df.shape[0]))
    # batch = np.array(np.concatenate(batch))
    # data['batch'] = batch
    # data['batch'] -= np.min(data['batch'])

    # print(data.batch.value_counts())
    # assert False

    #     print(data['batch'].value_counts())
    #     print(data.head())

    # print(data.columns)

    # divide in train and test data -- copied from above, organize differently!
    test_dataframe = data.sample(frac=0.01)
    train_dataframe = data.drop(test_dataframe.index)
    batch_test = test_dataframe["batch"]
    batch_train = train_dataframe["batch"]

    #     print('train batch labels')
    #     print(batch_train)
    #     print(batch_train.value_counts())

    test_dataframe.index = range(len(test_dataframe))
    train_dataframe.index = range(len(train_dataframe))
    n_train = len(train_dataframe.index)
    n_test = len(test_dataframe.index)
    n_train + n_test

    train_data = mb.datasets.MultiDataset(data_frame=train_dataframe)
    train_data.batch = np.array(batch_train)
    train_data.batch_one_hot = mb.tl.onehot_covar(train_data.batch)
    train_loader = tdata.DataLoader(dataset=train_data, batch_size=256, shuffle=True)
    # print(np.bincount(batch))

    test_data = mb.datasets.MultiDataset(data_frame=test_dataframe)
    test_data.batch = np.array(batch_test)
    test_data.batch_one_hot = mb.tl.onehot_covar(test_data.batch)
    test_loader = tdata.DataLoader(dataset=test_data, batch_size=1, shuffle=False)

    return train_loader, test_loader
