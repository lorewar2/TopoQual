import torch
import numpy as np
import random
import model
from dataset import QualityDataset
from torch.utils.data import DataLoader
import random
import math

DATA_PATH = "/data1/hifi_consensus/processed_data/chr2_ip_pw_filtered.txt"
RAW_PATH = "/data1/hifi_consensus/processed_data/chr1_ip_pw.txt"
MODEL_PATH = "./result/model/multi_layered_model.pt"
CONTEXT_COUNT = 3
EXTRA_COUNT = 20

def main():
    # set the seed
    torch.manual_seed(1)
    random.seed(3)
    np.random.seed(2)
    #evaluate_model()
    return

# this function will evalute the model and aggregate the results (output of the model for wrong and right)
def evaluate_model():
    tensor_length = pow(5, CONTEXT_COUNT) + EXTRA_COUNT
    # arrays to save the result
    error_counts = [0] * 94
    all_counts = [0] * 94
    batch_size = 1024
    # get the data to test
    eval_dataset = QualityDataset (DATA_PATH, False, CONTEXT_COUNT)
    eval_loader = DataLoader (
        dataset = eval_dataset,
        batch_size = batch_size,
        num_workers = 64,
        shuffle = False,
        drop_last = True
    )
    eval_len = len(eval_loader)
    # load the model
    lr_model = model.quality_model_1_layer(CONTEXT_COUNT, EXTRA_COUNT)
    checkpoint = torch.load(MODEL_PATH)
    lr_model.load_state_dict(checkpoint['model_state_dict'])
    # run the data
    with torch.no_grad():
        lr_model.eval()
        for batch_idx, (batch_inputs, batch_labels) in enumerate(eval_loader):
            pred = lr_model(batch_inputs)
            for i in range(len(batch_inputs)):
                pacbio_qual = batch_inputs[i][0][tensor_length - 5].item()
                position = int(-10 * math.log((1.0 - pred[i].item()) + 0.000000000001, 10))
                if position > 92:
                    position = 93
                all_counts[position] += 1
                if batch_labels[i].item() < 0.5 and pacbio_qual > 0.001:
                    error_counts[position] += 1
                # save the output and stuff
            print("Evaluating {}/{}".format(batch_idx, eval_len))
    print(all_counts)
    print(error_counts)


if __name__ == "__main__":
    main()