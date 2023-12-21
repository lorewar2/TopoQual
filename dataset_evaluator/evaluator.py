import torch
import numpy as np
import random
import model
from dataset import QualityDataset
from torch.utils.data import DataLoader
import random
import math
import os
import pysam

DATA_PATH = "./intermediate/"
MODEL_PATH = "./dataset_evaluator/model/multi_layered_model_new.pt"
CONTEXT_COUNT = 3
EXTRA_COUNT = 20

def main():
    # set the seed
    torch.manual_seed(1)
    random.seed(3)
    np.random.seed(2)
    # get the arguments read bam and output bam paths
    create_modified_bam()
    return

def create_modified_bam():
    infile = pysam.AlignmentFile("./sample_files/test.ccs.bam", "rb", check_sq = False)
    outfile = pysam.AlignmentFile("read_modified.bam", "wb", template = infile)
    for read in infile:
        required_query_name = read.query_name
        # find the query in the intermediate folder
        required_path = "{}{}".format(DATA_PATH, required_query_name.replace("/", "."))
        if os.path.isfile(required_path):
            (polished_seq, polished_qual) = evaluate(required_path)
        else:
            continue
        print(polished_seq)
        print(polished_qual)
        # save the sequence and qualities to the new bam
        read.query_sequence = polished_seq
        read.query_qualities = pysam.qualitystring_to_array(polished_qual)
        outfile.write(read)
    return

# this function will evalute the model and return the sequence and the quality scores
def evaluate(file_path):
    tensor_length = pow(5, CONTEXT_COUNT) + EXTRA_COUNT
    # arrays to save the result
    polished_sequence_arr = []
    polished_quality_arr = []
    # get the data from the file
    eval_dataset = QualityDataset (file_path, False, CONTEXT_COUNT)
    eval_loader = DataLoader (
        dataset = eval_dataset,
        batch_size = 1,
        num_workers = 64,
        shuffle = False,
        drop_last = False
    )
    eval_len = len(eval_loader)
    # load the model
    lr_model = model.quality_model_1_layer(CONTEXT_COUNT, EXTRA_COUNT)
    checkpoint = torch.load(MODEL_PATH)
    lr_model.load_state_dict(checkpoint['model_state_dict'])
    # run the data
    with torch.no_grad():
        lr_model.eval()
        for batch_idx, (batch_inputs, calling_base) in enumerate(eval_loader):
            if calling_base == "X":
                continue
            # get the quality prediction
            pred = lr_model(batch_inputs)
            quality = int(-10 * math.log((1.0 - pred.item()) + 0.000000000001, 10))
            # cut off
            if quality > 52:
                quality = 52
            if quality < 20:
                quality = 20
            # add base and qual to array
            polished_sequence_arr.append(chr(calling_base))
            polished_quality_arr.append(chr(quality))
            print("Evaluating {}/{}".format(batch_idx, eval_len))
    # convert both to string
    polished_sequence_string = "".join([str(i) for i in polished_sequence_arr])
    polisehed_quality_string = "".join([str(i) for i in polished_quality_arr])
    return (polished_sequence_string, polisehed_quality_string)

if __name__ == "__main__":  
    main()