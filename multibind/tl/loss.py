import torch
import torch.nn as tnn



class PoissonLoss(tnn.Module):
    def __init__(self):
        super(PoissonLoss, self).__init__()

    def forward(self, y_pred, y_true, is_count_data):
        return torch.mean(y_pred - y_true*torch.log(y_pred))
    
# (negative) Log-likelihood of the Poisson distribution
class MultiDatasetLoss(tnn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, y_pred, y_true, is_count_data):
        
        a, b = y_pred[is_count_data == 1], y_pred[is_count_data == 0]        
        # print(a.shape, b.shape)
                
        # add poisson loss
        poisson_loss = None
        if a.shape[0] > 0:
            loss = a - y_true[is_count_data == 1]*torch.log(a)
            loss = torch.abs(loss)
            poisson_loss = torch.mean(loss)

        m = torch.nn.Sigmoid()
        bce = torch.nn.BCELoss()
        bce_loss = None
        if b.shape[0] > 0:
            bce_loss = bce(m(b), y_true[is_count_data == 0])
        
        # print(poisson_loss, bce_loss)
        if a.shape[0] != 0 and b.shape[0] != 0:
            return poisson_loss + bce_loss
        elif a.shape[0] != 0:
            return poisson_loss
        elif b.shape[0] != 0:
            return bce_loss
        assert False # problem with the input data


# Custom loss function
class CustomLoss(tnn.Module):
    def __init__(self, weight=None, size_average=True):
        super().__init__()

    def forward(self, inputs, rounds, avoid_zero=True):
        if avoid_zero:
            rounds = rounds + 0.1
        f = inputs*rounds[:, 0]
        return -torch.sum(rounds[:, 1]*torch.log(f+0.0001) - f)