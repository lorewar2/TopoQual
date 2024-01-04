import torch
import util

class quality_model_1_layer(torch.nn.Module):
    # Object Constructor
    def __init__(self, base_context_count, extra_count):
        super().__init__()
        self.tensor_length = pow(5, base_context_count) + extra_count
        self.linear = torch.nn.Linear(self.tensor_length, int(self.tensor_length / 2), bias = True)
        self.linear2 = torch.nn.Linear(int(self.tensor_length / 2), int(self.tensor_length / 4), bias = True)
        self.linear3 = torch.nn.Linear(int(self.tensor_length / 4), 1, bias = True)
        self.sig = torch.nn.Sigmoid()
    # define the forward function for prediction
    def forward(self, x):
        out = self.linear(x)
        out = self.linear2(out)
        out = self.linear3(out)
        out = self.sig(out)
        return out