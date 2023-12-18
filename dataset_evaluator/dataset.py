import torch
import os
import numpy as np

class QualityDataset(torch.utils.data.Dataset):
    def __init__(self, file_loc, shuffle_all, base_context_count):
        self.file_loc = file_loc
        self.shuffle_all = shuffle_all
        self.base_context_count = base_context_count
        self.tensor_length = pow(5, base_context_count) + 20
        # get len and save it
        with open(file_loc) as f:
            f.seek(0, 2)
            offset = f.tell()
            self.len = int((offset - 90) / 90) - 1
        # load all data
        if self.shuffle_all:
            # make a shuffled index
            self.index_array = np.arange(0, self.len)
            np.random.shuffle(self.index_array)
            print("initializing complete")

    def __len__(self):
        return self.len

    def reshuffle(self):
        if self.shuffle_all:
            np.random.shuffle(self.index_array)
            print("reshuffling complete")
        return

    def __getitem__(self, index):
        if self.shuffle_all == True:
            (input_tensor, calling_base) = self.retrieve_item_from_disk(self.index_array[index])
        else:
            (input_tensor, calling_base) = self.retrieve_item_from_disk(index)
        return (input_tensor, calling_base)

    def retrieve_item_from_disk(self, index):
        # search the index file to file the location # index offset is 90
        retrieved_line = ""
        with open(self.file_loc) as f1:
            f1.seek(index * 90)
            retrieved_line = f1.readline()
        split_txt = retrieved_line.split(" ")
        # case of corrupted data $dont use this$
        if len(split_txt) != 16:
            print("====================ERROR=========================")
            return (torch.zeros(1, self.tensor_length), "X")
        # get the calling base and the state if opao
        calling_base, poa_state = self.get_state_info(split_txt[5])
        if calling_base == "X":
            calling_base = split_txt[4][3]
        # get the required base context
        if self.base_context_count == 3:
            base_context = [split_txt[4][2], calling_base, split_txt[4][4]]
        elif self.base_context_count == 5:
            base_context = [split_txt[4][1], split_txt[4][2], calling_base, split_txt[4][4], split_txt[4][5]]
        else:
            base_context = [split_txt[4][0], split_txt[4][1], split_txt[4][2], calling_base, split_txt[4][4], split_txt[4][5], split_txt[4][6]]
        # get the number from the base context
        converted_number = self.convert_bases_to_bits(base_context, self.base_context_count)
        hot_encoded = [0.0] * pow(5, self.base_context_count)
        hot_encoded[converted_number] = 1.0
        # get the read length and position information
        base_pos_info = self.read_len_info(int(split_txt[2]), int(split_txt[3]))
        # get quality in float
        quality = float(split_txt[0]) / 100
        # get the num of parallel bases in float
        parallel_vec_f = self.clean_string_get_array([split_txt[12], split_txt[13], split_txt[14], split_txt[15]])
        # retrieve sn
        sn_vec_f = self.clean_string_get_array([split_txt[6], split_txt[7], split_txt[8], split_txt[9]])
        # get the required sn details
        sn_info = self.process_sn_info([split_txt[4][2], split_txt[4][3], split_txt[4][4]], calling_base, sn_vec_f)
        # get pw and ip
        ip = float(split_txt[10])
        pw = float(split_txt[11])
        # rearrange so that the calling base num first and rest in decending order
        sorted_vec = self.rearrange_sort_parallel_bases(parallel_vec_f, calling_base)
        # make and append to the input tensor,
        input_tensor = torch.tensor([hot_encoded + poa_state + base_pos_info + sn_info + [ip, pw, quality] + sorted_vec])
        return (input_tensor, calling_base)

    def process_sn_info(self, three_base_context, calling_base, sn_vec):
        if calling_base == "X":
            calling_base = three_base_context[1]
        if three_base_context[0] == "X":
            three_base_context = [calling_base, calling_base, calling_base]
        if three_base_context[2] == "X":
            three_base_context = [calling_base, calling_base, calling_base]
        # base sn
        base_sn = sn_vec[self.get_base_to_int(calling_base)]
        # left base sn diff
        left_diff = abs(base_sn - sn_vec[self.get_base_to_int(three_base_context[0])])
        # right base sn diff
        right_diff = abs(base_sn - sn_vec[self.get_base_to_int(three_base_context[2])])
        # sn vec with least diff to highest
        sn_vec_diff = [abs(base_sn - sn_vec[0]), abs(base_sn - sn_vec[1]), abs(base_sn - sn_vec[2]), abs(base_sn - sn_vec[3])]
        sn_vec_diff.sort()
        return [base_sn, left_diff, right_diff] + sn_vec_diff

    def convert_bases_to_bits(self, base_array, count):
        converted_number = 0
        for index in range(0, count):
            base_number = self.get_base_to_int(base_array[count - index - 1])
            converted_number += pow(5, index) * base_number
        return converted_number

    def get_base_to_int(self, base):
        result = 4
        if base == "A":
            result = 0
        elif base == "C":
            result = 1
        elif base == "G":
            result = 2
        elif base == "T":
            result = 3
        elif base == "X":
            result = 4
        return result

    def get_state_info(self, status):
        calling_base = "X"
        poa_state = [0.0] * 3
        if status[0] == "O":
            calling_base = status[3]
            poa_state[0] = 1.0
        elif status[0] == "S":
            calling_base = status[3]
            poa_state[1] = 1.0
        else:
            poa_state[2] = 1.0
        return calling_base, poa_state

    def read_len_info (self, read_position, read_len):
        base_pos_info = [0.0] * 3
        if (read_position >= read_len - 3) or (read_position <= 3):
            base_pos_info[0] = 1.0
        if (read_position >= read_len - 10) or (read_position <= 10):
            base_pos_info[1] = 1.0
        if (read_position >= read_len - 100) or (read_position <= 100):
            base_pos_info[2] = 1.0
        return base_pos_info
        
    def clean_string_get_array(self, string_array):
        char_remov = ["]", "[", ",", "\n"]
        for char in char_remov:
            for index_s in range(len(string_array)):
                temp = string_array[index_s].replace(char, "")
                string_array[index_s] = temp
        vec_f = []
        for parallel in string_array:
            vec_f.append(float(parallel))
        return vec_f
    
    def rearrange_sort_parallel_bases(self, parallel_vec, base):
        if base == "A":
            selected_base = parallel_vec[0]
            del parallel_vec[0]
        elif base == "C":
            selected_base = parallel_vec[1]
            del parallel_vec[1]
        elif base == "G":
            selected_base = parallel_vec[2]
            del parallel_vec[2]
        elif base == "T":
            selected_base = parallel_vec[3]
            del parallel_vec[3]
        elif base == "X":
            return [0.0, 0.0, 0.0, 0.0]
        else:
            selected_base = parallel_vec[0]
            del parallel_vec[0]
        parallel_vec.sort(reverse = True)
        parallel_vec = [selected_base] + parallel_vec
        return parallel_vec