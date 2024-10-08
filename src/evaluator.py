import torch
import sys, getopt
import numpy as np
import random
import model
from dataset import QualityDataset
from torch.utils.data import DataLoader
import random
import math
import pysam
import os, os.path

DATA_PATH = "./intermediate/"
MODEL_PATH = "./model/multi_layered_model_new.pt"
CONTEXT_COUNT = 3
EXTRA_COUNT = 20
INPUT_FILE = "./sample_files/test.ccs.bam"
OUTPUT_FILE = "read_modified.bam"

# def main(argv):
#     # get arguments input file and output file
#     inputfile = ''
#     outputfile = ''
#     opts, args = getopt.getopt(argv, "hi:o:",["ifile=","ofile="])
#     for opt, arg in opts:
#         if opt == '-h':
#             print ('evaluate.py -i <inputbam> -o <outputbam>')
#             sys.exit()
#         elif opt in ("-i", "--ifile"):
#             inputfile = arg
#         elif opt in ("-o", "--ofile"):
#             outputfile = arg
#     print ('Input bam is ', inputfile)
#     print ('Output bam is ', outputfile)
#     # set the seed
#     torch.manual_seed(1)
#     random.seed(3)
#     np.random.seed(2)
#     # get the arguments read bam and output bam paths
#     create_modified_bam(inputfile, outputfile)
#     return
def main():
    evaluate("test.txt")
    return

def create_modified_bam(inputfile, outputfile):
    inbam = pysam.AlignmentFile(inputfile, "rb", check_sq = False)
    outbam = pysam.AlignmentFile(outputfile, "wb", template = inbam)
    total_len = len([name for name in os.listdir(DATA_PATH) if os.path.isfile(os.path.join(DATA_PATH, name))])
    index = 0
    for read in inbam:
        required_query_name = read.query_name
        # find the query in the intermediate folder
        required_path = "{}{}".format(DATA_PATH, required_query_name.replace("/", "."))
        if os.path.isfile(required_path):
            (polished_seq, polished_qual) = evaluate(required_path)
        else:
            continue
        #print(polished_seq)
        #print(polished_qual)
        # save the sequence and qualities to the new bam
        read.query_sequence = polished_seq
        read.query_qualities = pysam.qualitystring_to_array(polished_qual)
        outbam.write(read)
        index += 1
        print("Progress {} / {}".format(index, total_len))
    return

# this function will evalute the model and return the sequence and the quality scores
def evaluate(file_path):
    batch_size = 1024
    # get the data from the file
    eval_dataset = QualityDataset (file_path, False, CONTEXT_COUNT)
    eval_loader = DataLoader (
        dataset = eval_dataset,
        batch_size = batch_size,
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
            # get the quality prediction
            pred = lr_model(batch_inputs)
            for idx in range(len(calling_base)):
                if calling_base[idx] == "X":
                    continue
                # get the quality prediction
                quality = int(-10 * math.log((1.0 - pred[idx].item()) + 0.000000000001, 10))
                # cut off
                print("before ", quality)
                quality = quality + 33
                if quality > 52:
                    quality = 52
                if quality < 20:
                    quality = 20
                print("after ", quality)
                # add base and qual to array
            print("Evaluating {}/{}".format(batch_idx, eval_len))
    return

if __name__ == "__main__":  
    main(sys.argv[1:])