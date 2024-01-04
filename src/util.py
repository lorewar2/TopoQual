
import torch
import random
import numpy as np
import scipy.special

CORRECT_PATH = "/data1/hifi_consensus/processed_data/mutation_data/chr2_mutation_5base_all.txt"
ERROR_PATH = "/data1/hifi_consensus/processed_data/mutation_data/chr2_mutation_5base_err.txt"
WRITE_PATH = "/data1/hifi_consensus/processed_data/mutation_data/chr2_mutation_5base_final.txt"
READ_MUTATION_PATH = "/data1/hifi_consensus/processed_data/mutation_data/chr2_mutation_3base_final.txt"

def output_the_base_corrections (data_path, write_folder):
    ttoa = [[0 for i in range(16)] for j in range(5)] # total error poa del parallel
    ttoc = [[0 for i in range(16)] for j in range(5)] # total error poa del parallel
    ttog = [[0 for i in range(16)] for j in range(5)] # total error poa del parallel

    index = 0
    with open(data_path) as f1:
        for line in f1:
            index += 1
            split_txt = line.split(" ")
            if len(split_txt) != 18:
                continue
            else:
                # get the converted number
                base_correction_context = [split_txt[1][2], split_txt[1][4]]
                converted_number = convert_bases_to_2bit(base_correction_context)
                calling_base, poa_state = get_state_info(split_txt[7])
                if split_txt[1][3] == "T":
                    ttoa[0][converted_number] += 1
                    if split_txt[6][3] == "A":
                        ttoa[1][converted_number] += 1
                        if calling_base == split_txt[1][3]:
                            ttoa[2][converted_number] += 1
                        elif poa_state[2] > 0.5:
                            ttoa[3][converted_number] += 1
                        else:
                            parallel_vec_f = clean_string_get_array([split_txt[14], split_txt[15], split_txt[16], split_txt[17]])
                            if (parallel_vec_f [3] >= parallel_vec_f[0]):
                                ttoa[4][converted_number] += 1
                    elif split_txt[6][3] == "C":
                        ttoc[1][converted_number] += 1
                        if calling_base == split_txt[1][3]:
                            ttoc[2][converted_number] += 1
                        elif poa_state[2] > 0.5:
                            ttoc[3][converted_number] += 1
                        else:
                            parallel_vec_f = clean_string_get_array([split_txt[14], split_txt[15], split_txt[16], split_txt[17]])
                            if (parallel_vec_f [3] >= parallel_vec_f[1]):
                                ttoc[4][converted_number] += 1
                    elif split_txt[6][3] == "G":
                        ttog[1][converted_number] += 1
                        if calling_base == split_txt[1][3]:
                            ttog[2][converted_number] += 1
                        elif poa_state[2] > 0.5:
                            ttog[3][converted_number] += 1
                        else:
                            parallel_vec_f = clean_string_get_array([split_txt[14], split_txt[15], split_txt[16], split_txt[17]])
                            if (parallel_vec_f [3] >= parallel_vec_f[2]):
                                ttog[4][converted_number] += 1
            if index % 10000 == 0:
                if index > 0:
                    break
                print("lines {}".format(index))
    with open("{}/ttoa.txt".format(write_folder), 'a+') as fw:
        for i in range(16):
           fw.write("{},{},{},{},{}\n".format(ttoa[0][i], ttoa[1][i], ttoa[2][i], ttoa[3][i], ttoa[4][i]))
    with open("{}/ttoc.txt".format(write_folder), 'a+') as fw:
        for i in range(16):
           fw.write("{},{},{},{},{}\n".format(ttoc[0][i], ttoc[1][i], ttoc[2][i], ttoc[3][i], ttoc[4][i]))
    with open("{}/ttog.txt".format(write_folder), 'a+') as fw:
        for i in range(16):
           fw.write("{},{},{},{},{}\n".format(ttog[0][i], ttog[1][i], ttog[2][i], ttog[3][i], ttog[4][i]))
    return

def convert_bases_to_2bit(base_array):
    converted_number = 0
    count = 2
    for index in range(0, count):
        base_number = get_base_to_int2(base_array[count - index - 1])
        converted_number += pow(4, index) * base_number
    return converted_number

def convert_2bit_to_bases(converted_number):
    base_array = []
    count = 2
    for _ in range(0, count):
        base_array.append(get_int_to_base2(converted_number % 4))
        converted_number = int(converted_number / 4)
    base_array.reverse()
    return base_array

def get_base_to_int2(base):
    result = 0
    if base == "A":
        result = 0
    elif base == "C":
        result = 1
    elif base == "G":
        result = 2
    elif base == "T":
        result = 3
    elif base == "X":
        result = 0
    return result

def get_int_to_base2(number):
    base = 'A'
    if number == 0:
        base = 'A'
    elif number == 1:
        base = 'C'
    elif number == 2:
        base = 'G'
    elif number == 3:
        base = 'T'
    return base

def clean_string_get_array(string_array):
    char_remov = ["]", "[", ",", "\n"]
    for char in char_remov:
        for index_s in range(len(string_array)):
            temp = string_array[index_s].replace(char, "")
            string_array[index_s] = temp
    vec_f = []
    for parallel in string_array:
        vec_f.append(float(parallel))
    return vec_f

def get_state_info(status):
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

def check_line_sizes_in_file(file_loc):
    with open(file_loc) as f:
        for index in range(0, 1000000):
            f.seek(index * 108)
            line = f.readline()
            print(line)
            if (len(line) > 10):
                break
    return

def filter_data_using_confident_germline_indel_depth(chromosone, data_path, filter_path, write_path):
    # ALL DATA IN ORDER
    # read the confident file put relavant chromosone data in array
    confident_regions = []
    path = "{}/confident_data.txt".format(filter_path)
    with open(path, 'r') as cr:
        for index, line in enumerate(cr):
            split_txt = line.split(" ")
            if chromosone == split_txt[0]:
                start = int(split_txt[1])
                end = int(split_txt[3])
                confident_regions.append((start, end))
    # read germline file put relavant chromosone data in array
    germline_locations = []
    path = "{}/germline_data.txt".format(filter_path)
    with open(path, 'r') as gr:
        for index, line in enumerate(gr):
            split_txt = line.split(" ")
            if chromosone == split_txt[0]:
                location = int(split_txt[1])
                count = len(split_txt[2])
                germline_locations.append((location, count))
    # read the data file, go line by line
    modified_lines = []
    read_file = open(data_path, 'r')
    confident_index = 0
    germline_index = 0
    print("confident regions : {} germline locations : {}".format(len(confident_regions), len(germline_locations)))
    with open(write_path, 'a') as fw:
        for index, line in enumerate(read_file):
            split_txt = line.split(" ")
            if len(split_txt) != 18:
                continue
            current_location = int(split_txt[0])
            # iterate to correct area of confident region
            while current_location > confident_regions[confident_index][1]:
                if confident_index + 1 >= len(confident_regions):
                    break
                confident_index += 1
            # iterate to correct area of germline region
            while current_location > germline_locations[germline_index][0]:
                if germline_index + 1 >= len(germline_locations):
                    break
                germline_index += 1
            # check if in confident region if not continue
            if (current_location < confident_regions[confident_index][0]) or (current_location > confident_regions[confident_index][1]):
                #print("Not confident region {} start: {} end: {} ".format(current_location, confident_regions[confident_index][0], confident_regions[confident_index][1]))
                continue
            # check if germline variant
            if (current_location >= germline_locations[germline_index][0]) and (current_location <= (germline_locations[germline_index][0] + germline_locations[germline_index][1])):
                #print("Germline variant location {} == {} +- {}".format(current_location, germline_locations[germline_index][0], germline_locations[germline_index][1]))
                continue
            # this is run if not filtered
            # filter the line
            char_remov = ["]", "[", ",", "\n"]
            stripped_line = line
            for char in char_remov:
                stripped_line = stripped_line.replace(char, "")
            split_txt = stripped_line.split(" ")
            location = split_txt[0].zfill(9)
            seven_base_ref = split_txt[1]
            quality = split_txt[2].zfill(2)
            read_position = split_txt[4].zfill(5)
            read_length = split_txt[5].zfill(5)
            seven_base_cont = split_txt[6]
            base = split_txt[7]
            sn1 = split_txt[8].zfill(7)
            sn2 = split_txt[9].zfill(8)
            sn3 = split_txt[10].zfill(8)
            sn4 = split_txt[11].zfill(8)
            # check the
            pw = split_txt[12].zfill(3)
            ip = split_txt[13].zfill(3)
            parallel1 = split_txt[14].zfill(2)
            parallel2 = split_txt[15].zfill(2)
            parallel3 = split_txt[16].zfill(2)
            parallel4 = split_txt[17].zfill(2)
            modified_line = "{} {} {} : {} {} {} {} [{} {} {} {}] {} {} [{} {} {} {}]\n".format(location, seven_base_ref, quality, read_position, read_length, seven_base_cont, base, sn1, sn2, sn3, sn4, pw, ip, parallel1, parallel2, parallel3, parallel4)
            modified_lines.append(modified_line)
            if index % 1_000_000 == 0:
                for write_line in modified_lines:
                    fw.write(write_line)
                modified_lines.clear()
                print("processed {} records, {}/{}".format(index, germline_index, len(germline_locations)))
    return

def get_summary_ip_pw(read_path):
    error_ip_array = []
    error_ip_sum = 0
    error_pw_array = []
    error_pw_sum = 0
    correct_ip_array = []
    correct_ip_sum = 0
    correct_pw_array = []
    correct_pw_sum = 0
    # arrays to save the result
    file = open(read_path, "r")
    for index, line in enumerate(file):
        split_txt = line.split(" ")
        if len(split_txt) != 18:
            continue
        ip = int(split_txt[12])
        pw = int(split_txt[13])
        ref_base = split_txt[1][3]
        call_base = split_txt[6][3]
        if ref_base != call_base:
            error_ip_array.append(ip)
            error_ip_sum += ip
            error_pw_array.append(pw)
            error_pw_sum += pw
        else:
            correct_ip_array.append(ip)
            correct_ip_sum += ip
            correct_pw_array.append(pw)
            correct_pw_sum += pw
        if index % 100_000 == 0:
            print("current line {}".format(index))
    correct_pw_avg = correct_pw_sum / len(correct_pw_array)
    error_pw_avg = error_pw_sum / len(error_pw_array)
    correct_ip_avg = correct_ip_sum / len(correct_ip_array)
    error_ip_avg = error_ip_sum / len(error_ip_array)
    print("correct: pw avg: {} ip avg: {}".format(correct_pw_avg, correct_ip_avg))
    print("error: pw avg: {} ip avg: {}".format(error_pw_avg, error_ip_avg))
    print("sorting arrays...")
    correct_ip_array.sort()
    error_ip_array.sort()
    correct_pw_array.sort()
    error_pw_array.sort()
    print("Correct summary:")
    print("pw avg:{} max:{} min:{} median:{} 1stQ:{} 3rdQ:{} iqr:{}".format(correct_pw_avg, correct_pw_array[len(correct_pw_array) - 1], correct_pw_array[0], correct_pw_array[int(len(correct_pw_array)/2)], correct_pw_array[int(len(correct_pw_array)/4)], correct_pw_array[int(3 * len(correct_pw_array)/4)], -correct_pw_array[int(len(correct_pw_array)/4)] + correct_pw_array[int(3 * len(correct_pw_array)/4)]))
    print("ip avg:{} max:{} min:{} median:{} 1stQ:{} 3rdQ:{} iqr:{}".format(correct_ip_avg, correct_ip_array[len(correct_ip_array) - 1], correct_ip_array[0], correct_ip_array[int(len(correct_ip_array)/2)], correct_ip_array[int(len(correct_ip_array)/4)], correct_ip_array[int(3 * len(correct_ip_array)/4)], -correct_ip_array[int(len(correct_ip_array)/4)] + correct_ip_array[int(3 * len(correct_ip_array)/4)]))
    print("Error summary:")
    print("pw avg:{} max:{} min:{} median:{} 1stQ:{} 3rdQ:{} iqr:{}".format(error_pw_avg, error_pw_array[len(error_pw_array) - 1], error_pw_array[0], error_pw_array[int(len(error_pw_array)/2)], error_pw_array[int(len(error_pw_array)/4)], error_pw_array[int(3 * len(error_pw_array)/4)], -error_pw_array[int(len(error_pw_array)/4)] + error_pw_array[int(3 * len(error_pw_array)/4)]))
    print("ip avg:{} max:{} min:{} median:{} 1stQ:{} 3rdQ:{} iqr:{}".format(error_ip_avg, error_ip_array[len(error_ip_array) - 1], error_ip_array[0], error_ip_array[int(len(error_ip_array)/2)], error_ip_array[int(len(error_ip_array)/4)], error_ip_array[int(3 * len(error_ip_array)/4)], -error_ip_array[int(len(error_ip_array)/4)] + error_ip_array[int(3 * len(error_ip_array)/4)]))
    print(read_path)
    return

def print_deep_scores(read_path):
    # arrays to save the result
    error_counts = [0] * 94
    all_counts = [0] * 94
    file = open(read_path, "r")
    for index, line in enumerate(file):
        split_txt = line.split(" ")
        if len(split_txt) != 5:
            continue
        base_quality = int(split_txt[2])
        ref_base = split_txt[1]
        call_base = split_txt[4].strip()
        all_counts[base_quality] += 1
        if ref_base != call_base:
            error_counts[base_quality] += 1
        if index % 100000 == 0:
            print("Running line {}".format(index))
            break
    print(error_counts)
    print(all_counts)
    print(read_path)
    return

def filter_data_deep_consensus(chromosone, data_path, filter_path, write_path):
    # ALL DATA IN ORDER
    # read the confident file put relavant chromosone data in array
    confident_regions = []
    path = "{}/confident_data.txt".format(filter_path)
    with open(path, 'r') as cr:
        for index, line in enumerate(cr):
            split_txt = line.split(" ")
            if chromosone == split_txt[0]:
                start = int(split_txt[1])
                end = int(split_txt[3])
                confident_regions.append((start, end))
    # read germline file put relavant chromosone data in array
    germline_locations = []
    path = "{}/germline_data.txt".format(filter_path)
    with open(path, 'r') as gr:
        for index, line in enumerate(gr):
            split_txt = line.split(" ")
            if chromosone == split_txt[0]:
                location = int(split_txt[1])
                count = len(split_txt[2])
                germline_locations.append((location, count))
    # read the data file, go line by line
    modified_lines = []
    read_file = open(data_path, 'r')
    confident_index = 0
    germline_index = 0
    print("confident regions : {} germline locations : {}".format(len(confident_regions), len(germline_locations)))
    with open(write_path, 'a') as fw:
        for index, line in enumerate(read_file):
            split_txt = line.split(" ")
            if len(split_txt) != 18:
                continue
            current_location = int(split_txt[0])
            # iterate to correct area of confident region
            while current_location > confident_regions[confident_index][1]:
                if confident_index + 1 >= len(confident_regions):
                    break
                confident_index += 1
            # iterate to correct area of germline region
            while current_location > germline_locations[germline_index][0]:
                if germline_index + 1 >= len(germline_locations):
                    break
                germline_index += 1
            # check if in confident region if not continue
            if (current_location < confident_regions[confident_index][0]) or (current_location > confident_regions[confident_index][1]):
                #print("Not confident region {} start: {} end: {} ".format(current_location, confident_regions[confident_index][0], confident_regions[confident_index][1]))
                continue
            # check if germline variant
            if (current_location >= germline_locations[germline_index][0]) and (current_location <= (germline_locations[germline_index][0] + germline_locations[germline_index][1])):
                #print("Germline variant location {} == {} +- {}".format(current_location, germline_locations[germline_index][0], germline_locations[germline_index][1]))
                continue
            # this is run if not filtered
            # filter the line
            modified_lines.append(line)
            if index % 1_000_000 == 0:
                for write_line in modified_lines:
                    fw.write(write_line)
                modified_lines.clear()
                print("processed {} records, {}/{}".format(index, germline_index, len(germline_locations)))
    return

def concancate_quality_scores_from_files ():
    correct_count = [0] * 94
    error_count = [0] * 94
    for index in range(0, 64):
        current_path = "/data1/hifi_consensus/deepresult/{}_data.txt".format(index)
        with open(current_path, 'r') as hr:
            for line_index, line in enumerate(hr):
                cleaned_line = line
                char_remov = ["]", "[", ",", "\n"]
                for char in char_remov:
                    cleaned_line = cleaned_line.replace(char, "")
                #print(cleaned_line)
                split_line = cleaned_line.split(" ")
                if line_index == 0:
                    for array_index in range(0, 94):
                        correct_count[array_index] += int(split_line[array_index])
                elif line_index == 1:
                    for array_index in range(0, 94):
                        error_count[array_index] += int(split_line[array_index])
    print(correct_count)
    print(error_count)
    return

def get_mutation_probablility_array_prof (context_count):
    normal_correct_rate = 0.85
    normal_error_rate = 1 - normal_correct_rate
    array_length = pow(5, context_count)
    prob_array = [0.90] * array_length
    error_rates_array = [0.0] * array_length
    # get the average error rate and all error rates for each context
    total_error_count = 0
    total_all_count = 0
    with open(READ_MUTATION_PATH, 'r') as hr:
        for _, line in enumerate(hr):
            split_line = line.split(" ")
            index = int(split_line[0])
            all_count = int(split_line[1])
            error_count = int(split_line[2].strip())
            error_rates_array[index] = error_count / (all_count + 0.000000001)
            total_error_count += error_count
            total_all_count += all_count
    average_error_rate = total_error_count / (total_all_count + 0.000000001)
    # calculate the probablity and save it in array
    for index in range(0, array_length):
        multiplier = error_rates_array[index] / (average_error_rate + 0.000000001)
        context_error_rate = normal_error_rate * multiplier
        context_correct_rate = 1 - context_error_rate
        prob_array[index] = context_correct_rate
    print(prob_array)
    return prob_array
    
def fix_the_mutation_file(context_count):
    array_length = pow(5, context_count)
    error_array = [0] * array_length
    correct_array = [0] * array_length
    with open(CORRECT_PATH, 'r') as hr:
        for _, line in enumerate(hr):
            split_line = line.split(" ")
            index = int(split_line[0])
            correct = int(split_line[2])
            correct_array[index] = correct
    with open(ERROR_PATH, 'r') as hr:
        for _, line in enumerate(hr):
            split_line = line.split(" ")
            index = int(split_line[0])
            error = int(split_line[2])
            error_array[index] = error
    print(error_array)
    print(correct_array)
    with open(WRITE_PATH, 'a') as fw:
        for index in range(0, array_length):
            fw.write("{} {} {}\n".format(index, correct_array[index], error_array[index]))
    return

def get_mutation_probablility_array (context_count):
    array_length = pow(5, context_count)
    prob_array = [0.90] * array_length
    with open(CORRECT_PATH, 'r') as hr:
        for _, line in enumerate(hr):
            split_line = line.split(" ")
            index = int(split_line[0])
            quality_array = [int(split_line[8].strip()), int(split_line[7]), int(split_line[6]), int(split_line[5]), int(split_line[4])]
            temp_prob = 0.70
            for quality in quality_array:
                if quality > 0:
                    prob_array[index] = temp_prob
                    break
                else:
                    temp_prob += 0.05
            print(temp_prob)
    print(prob_array)
    return prob_array

def get_base_context_from_file(data_path, write_path1, write_path2, write_path3, prob):
    # read the correct file and take the indices in
    index_list_3 = []
    index_list_5 = []
    index_list_7 = []
    with open("chr2_mutation_context_correct.txt", 'r') as rf:
        for line in rf:
            split_txt = line.split(" ")
            if len(split_txt[1]) == 3:
                index_list_3.append(int(split_txt[0]))
            if len(split_txt[1]) == 5:
                index_list_5.append(int(split_txt[0]))
            if len(split_txt[1]) == 7:
                index_list_7.append(int(split_txt[0]))
    # initialize the arrays
    three_base_context_info = []
    for _ in range(0, pow(5, 3)):
        three_base_context_info.append([0, 0, 0, 0, 0, 0, 0])
    five_base_context_info = []
    for _ in range(0, pow(5, 5)):
        five_base_context_info.append([0, 0, 0, 0, 0, 0, 0])
    seven_base_context_info = []
    for _ in range(0, pow(5, 7)):
        seven_base_context_info.append([0, 0, 0, 0, 0, 0, 0])
    read_file = open(data_path, 'r')
    for index, line in enumerate(read_file):
        if index % 1_000_000 == 0:
            print("Running line {}".format(index))
        split_txt = line.split(" ")
        if len(split_txt) != 11:
            continue
        # get the quality
        calling_base = split_txt[5]
        ref_base = split_txt[1][1]
        parallel_vec_s = [split_txt[7], split_txt[8], split_txt[9], split_txt[10]]
        char_remov = ["]", "[", ",", "\n"]
        for char in char_remov:
            for index_s in range(len(parallel_vec_s)):
                temp = parallel_vec_s[index_s].replace(char, "")
                parallel_vec_s[index_s] = temp
        parallel_vec_f = []
        for parallel in parallel_vec_s:
            parallel_vec_f.append(float(parallel))
        # three base context
        three_base_num = convert_bases_to_bits([split_txt[3][2], split_txt[3][3], split_txt[3][4]], 3)
        # five base context
        five_base_num = convert_bases_to_bits([split_txt[3][1], split_txt[3][2], split_txt[3][3], split_txt[3][4], split_txt[3][5]], 5)
        # seven base context
        seven_base_num = convert_bases_to_bits([split_txt[3][0], split_txt[3][1], split_txt[3][2], split_txt[3][3], split_txt[3][4], split_txt[3][5], split_txt[3][6]], 7)
        
        three_base_context_info[three_base_num][0] += 1
        five_base_context_info[five_base_num][0] += 1
        seven_base_context_info[seven_base_num][0] += 1
        # if error
        if ref_base != calling_base:
            recalculated_score = int(calculate_topology_score(calling_base, parallel_vec_f[0], parallel_vec_f[1], parallel_vec_f[2], parallel_vec_f[3], (parallel_vec_f[0] + parallel_vec_f[1] + parallel_vec_f[2] + parallel_vec_f[3]), prob))
            three_base_context_info[three_base_num][1] += 1
            five_base_context_info[five_base_num][1] += 1
            seven_base_context_info[seven_base_num][1] += 1
            if recalculated_score > 50:
                three_base_context_info[three_base_num][2] += 1
                five_base_context_info[five_base_num][2] += 1
                seven_base_context_info[seven_base_num][2] += 1
            if recalculated_score > 60:
                three_base_context_info[three_base_num][3] += 1
                five_base_context_info[five_base_num][3] += 1
                seven_base_context_info[seven_base_num][3] += 1
            if recalculated_score > 70:
                three_base_context_info[three_base_num][4] += 1
                five_base_context_info[five_base_num][4] += 1
                seven_base_context_info[seven_base_num][4] += 1
            if recalculated_score > 80:
                three_base_context_info[three_base_num][5] += 1
                five_base_context_info[five_base_num][5] += 1
                seven_base_context_info[seven_base_num][5] += 1
            if recalculated_score > 90:
                three_base_context_info[three_base_num][6] += 1
                five_base_context_info[five_base_num][6] += 1
                seven_base_context_info[seven_base_num][6] += 1
    # write the data to file
    with open(write_path1, 'a') as fw:
        for index, info in enumerate(three_base_context_info):
            bases = convert_bits_to_bases(index, 3)
            bases_str = "{}{}{}".format(bases[0], bases[1], bases[2])
            info_str = "{} {} {} {} {} {} {}".format(info[0], info[1], info[2], info[3], info[4], info[5], info[6])
            #if info[0] == 0 and info [1] == 0:
                #continue
            if index in index_list_3:
                fw.write("{} {} {}\n".format(index, bases_str, info_str))
    with open(write_path2, 'a') as fw:
        for index, info in enumerate(five_base_context_info):
            bases = convert_bits_to_bases(index, 5)
            bases_str = "{}{}{}{}{}".format(bases[0], bases[1], bases[2], bases[3], bases[4])
            info_str = "{} {} {} {} {} {} {}".format(info[0], info[1], info[2], info[3], info[4], info[5], info[6])
            #if info[0] == 0 and info [1] == 0:
                #continue
            if index in index_list_5:
                fw.write("{} {} {}\n".format(index, bases_str, info_str))
    with open(write_path3, 'a') as fw:
        for index, info in enumerate(seven_base_context_info):
            bases = convert_bits_to_bases(index, 7)
            bases_str = "{}{}{}{}{}{}{}".format(bases[0], bases[1], bases[2], bases[3], bases[4], bases[5], bases[6])
            info_str = "{} {} {} {} {} {} {}".format(info[0], info[1], info[2], info[3], info[4], info[5], info[6])
            #if info[0] == 0 and info [1] == 0:
                #continue
            if index in index_list_7:
                fw.write("{} {} {}\n".format(index, bases_str, info_str))
    return

def convert_bases_to_bits(base_array, count):
    converted_number = 0
    for index in range(0, count):
        base_number = get_base_to_int(base_array[count - index - 1])
        converted_number += pow(5, index) * base_number
    return converted_number

def convert_bits_to_bases(converted_number, count):
    base_array = []
    for _ in range(0, count):
        base_array.append(get_int_to_base(converted_number % 5))
        converted_number = int(converted_number / 5)
    base_array.reverse()
    return base_array

def get_base_to_int(base):
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

def get_int_to_base(number):
    base = 'X'
    if number == 0:
        base = 'A'
    elif number == 1:
        base = 'C'
    elif number == 2:
        base = 'G'
    elif number == 3:
        base = 'T'
    elif number == 4:
        base = 'X'
    return base

def write_errors_to_file(read_path, write_path):
    # arrays to save the result
    total_error_count = 0
    file = open(read_path, "r")
    modified_lines = []
    with open(write_path, 'a') as fw:
        for index, line in enumerate(file):
            split_txt = line.split(" ")
            if len(split_txt) != 11:
                continue
            ref_base = split_txt[1][1]
            call_base = split_txt[5]
            if ref_base != call_base:
                total_error_count += 1
                modified_lines.append(line)
            if index % 100000 == 0:
                for write_line in modified_lines:
                    fw.write(write_line)
                modified_lines.clear()
                print("processed {} records, {}".format(index, total_error_count))
        for write_line in modified_lines:
            fw.write(write_line)
    print(read_path)
    return

def calculate_topology_score(calling_base, base_A_count, base_C_count, base_G_count, base_T_count, num_of_reads, prob):
    ln_prob_base_A = np.log(0.25)
    ln_prob_base_C = np.log(0.25)
    ln_prob_base_G = np.log(0.25)
    ln_prob_base_T = np.log(0.25)
    
    ln_prob_data_given_A = np.log(calculate_binomial(num_of_reads, base_A_count, prob))
    ln_prob_data_given_C = np.log(calculate_binomial(num_of_reads, base_C_count, prob))
    ln_prob_data_given_G = np.log(calculate_binomial(num_of_reads, base_G_count, prob))
    ln_prob_data_given_T = np.log(calculate_binomial(num_of_reads, base_T_count, prob))

    ln_sum_of_probabilities = ln_prob_data_given_A + ln_prob_base_A
    ln_sum_of_probabilities = np.logaddexp(ln_sum_of_probabilities, ln_prob_data_given_C + ln_prob_base_C)
    ln_sum_of_probabilities = np.logaddexp(ln_sum_of_probabilities, ln_prob_data_given_G + ln_prob_base_G)
    ln_sum_of_probabilities = np.logaddexp(ln_sum_of_probabilities, ln_prob_data_given_T + ln_prob_base_T)

    if calling_base == "A":
        correct_rate = np.exp(ln_prob_data_given_A + ln_prob_base_A - ln_sum_of_probabilities)
    elif calling_base == "C":
        correct_rate = np.exp(ln_prob_data_given_C + ln_prob_base_C - ln_sum_of_probabilities)
    elif calling_base == "G":
        correct_rate = np.exp(ln_prob_data_given_G + ln_prob_base_G - ln_sum_of_probabilities)
    elif calling_base == "T":
        correct_rate = np.exp(ln_prob_data_given_T + ln_prob_base_T - ln_sum_of_probabilities)
    else:
        correct_rate = np.exp(ln_prob_data_given_A + ln_prob_base_A - ln_sum_of_probabilities)

    error_rate = 1.0 - correct_rate
    quality_score = (-10.00) * np.log10(error_rate + 0.000000000000000000001)
    #print(quality_score)
    return quality_score

def calculate_topology_score_variable_prob (mutation_list, base_context, calling_base, base_A_count, base_C_count, base_G_count, base_T_count, num_of_reads):
    converted_number = convert_bases_to_bits(base_context, 3)
    prob = mutation_list[converted_number]
    # read the file
    ln_prob_base_A = np.log(0.25)
    ln_prob_base_C = np.log(0.25)
    ln_prob_base_G = np.log(0.25)
    ln_prob_base_T = np.log(0.25)
    
    ln_prob_data_given_A = np.log(calculate_binomial(num_of_reads, base_A_count, prob) + 0.000000000000000000001)
    ln_prob_data_given_C = np.log(calculate_binomial(num_of_reads, base_C_count, prob) + 0.000000000000000000001)
    ln_prob_data_given_G = np.log(calculate_binomial(num_of_reads, base_G_count, prob) + 0.000000000000000000001)
    ln_prob_data_given_T = np.log(calculate_binomial(num_of_reads, base_T_count, prob) + 0.000000000000000000001)

    ln_sum_of_probabilities = ln_prob_data_given_A + ln_prob_base_A
    ln_sum_of_probabilities = np.logaddexp(ln_sum_of_probabilities, ln_prob_data_given_C + ln_prob_base_C)
    ln_sum_of_probabilities = np.logaddexp(ln_sum_of_probabilities, ln_prob_data_given_G + ln_prob_base_G)
    ln_sum_of_probabilities = np.logaddexp(ln_sum_of_probabilities, ln_prob_data_given_T + ln_prob_base_T)

    if calling_base == "A":
        correct_rate = np.exp(ln_prob_data_given_A + ln_prob_base_A - ln_sum_of_probabilities)
    elif calling_base == "C":
        correct_rate = np.exp(ln_prob_data_given_C + ln_prob_base_C - ln_sum_of_probabilities)
    elif calling_base == "G":
        correct_rate = np.exp(ln_prob_data_given_G + ln_prob_base_G - ln_sum_of_probabilities)
    elif calling_base == "T":
        correct_rate = np.exp(ln_prob_data_given_T + ln_prob_base_T - ln_sum_of_probabilities)
    else:
        correct_rate = np.exp(ln_prob_data_given_A + ln_prob_base_A - ln_sum_of_probabilities)

    error_rate = 1.0 - correct_rate
    quality_score = (-10.00) * np.log10(error_rate + 0.000000000000000000001)
    if quality_score > 190:
        quality_score = 190
    #print(quality_score)
    return quality_score

def pipeline_calculate_topology_score_with_probability(read_path):
    # get the prob list
    mutation_list = get_mutation_probablility_array (7)
    # arrays to save the result
    error_counts = [0] * 300
    #all_counts = [0] * 300
    file = open(read_path, "r")
    for index, line in enumerate(file):
        split_txt = line.split(" ")
        if len(split_txt) != 11:
            continue
        calling_base = split_txt[5]
        ref_base = split_txt[1][1]
        base_context = [split_txt[3][0], split_txt[3][1], split_txt[3][2], split_txt[3][3], split_txt[3][4], split_txt[3][5], split_txt[3][6]]
        parallel_vec_s = [split_txt[7], split_txt[8], split_txt[9], split_txt[10]]
        char_remov = ["]", "[", ",", "\n"]
        for char in char_remov:
            for index_s in range(len(parallel_vec_s)):
                temp = parallel_vec_s[index_s].replace(char, "")
                parallel_vec_s[index_s] = temp
        parallel_vec_f = []
        for parallel in parallel_vec_s:
            parallel_vec_f.append(float(parallel))
        #all_counts[recalculated_score] += 1
        if ref_base != calling_base:
            recalculated_score = int(calculate_topology_score_variable_prob(mutation_list, base_context, calling_base, parallel_vec_f[0], parallel_vec_f[1], parallel_vec_f[2], parallel_vec_f[3], (parallel_vec_f[0] + parallel_vec_f[1] + parallel_vec_f[2] + parallel_vec_f[3])))
            error_counts[recalculated_score] += 1
        if index % 100 == 0:
            print("Running line {}".format(index))
    print(error_counts)
    #print(all_counts)
    return

def add_corrected_errors_to_file(read_path, write_path):
    file = open(read_path, "r")
    with open(write_path, 'a') as fw:
        for index, line in enumerate(file):
            split_txt = line.split(" ")
            if len(split_txt) != 9:
                continue
            location = split_txt[0].zfill(9)
            three_base = split_txt[1]
            quality = split_txt[2].zfill(2)
            base = split_txt[3]
            count = split_txt[4].zfill(2)
            # get the parallel bases
            parallel_vec_s = [split_txt[5], split_txt[6], split_txt[7], split_txt[8]]
            char_remov = ["]", "[", ",", "\n"]
            for char in char_remov:
                for index_s in range(len(parallel_vec_s)):
                    temp = parallel_vec_s[index_s].replace(char, "")
                    parallel_vec_s[index_s] = temp
            parallel_vec_mod = []
            for parallel in parallel_vec_s:
                parallel_vec_mod.append(parallel.zfill(2))
            modified_line = "{} {} {} {} {} [{} {} {} {}]\n".format(location, three_base, quality, base, count, parallel_vec_mod[0], parallel_vec_mod[1], parallel_vec_mod[2], parallel_vec_mod[3])
            fw.write(modified_line)

def remove_errors_from_file(read_path, write_path):
    # arrays to save the result
    total_error_count = 0
    file = open(read_path, "r")
    modified_lines = []
    with open(write_path, 'a') as fw:
        for index, line in enumerate(file):
            split_txt = line.split(" ")
            if len(split_txt) != 9:
                continue
            ref_base = split_txt[1][1]
            call_base = split_txt[3]
            if ref_base == call_base:
                total_error_count += 1
                modified_lines.append(line)
            if index % 100000 == 0:
                for write_line in modified_lines:
                    fw.write(write_line)
                modified_lines.clear()
                print("processed {} records, {}".format(index, total_error_count))
        for write_line in modified_lines:
            fw.write(write_line)
    print(read_path)
    return

def use_himut_file_to_identify_errors(chromosone, data_path, filter_path, write_path):
    # ALL DATA IN ORDER
    # read the himut file put relavant chromosone data in array
    error_locations = []
    path = "{}/himut_data.txt".format(filter_path)
    with open(path, 'r') as hr:
        for index, line in enumerate(hr):
            split_txt = line.split(" ")
            if chromosone == split_txt[0]:
                location = int(split_txt[1])
                ref = split_txt[2]
                alt = split_txt[4].strip()
                error_locations.append((location, ref, alt))
    # read the data file, go line by line
    modified_lines = []
    read_file = open(data_path, 'r')
    himut_index = 0
    with open(write_path, 'a') as fw:
        for index, line in enumerate(read_file):
            split_txt = line.split(" ")
            if len(split_txt) != 9:
                continue
            current_location = int(split_txt[0])
            # iterate to correct area of confident region
            while current_location > error_locations[himut_index][0]:
                if himut_index + 1 >= len(error_locations):
                    break
                himut_index += 1
            if current_location != error_locations[himut_index][0]:
                # check if error, if error do not append
                ref_base = split_txt[1][1]
                calling_base = split_txt[3]
                if ref_base == calling_base:
                    modified_lines.append(line)
            else:
                # check if correct error, if not ignore
                ref_base = split_txt[1][1]
                calling_base = split_txt[3]
                if ref_base == error_locations[himut_index][1] and calling_base == error_locations[himut_index][2]:
                    modified_lines.append(line)
                elif ref_base == error_locations[himut_index][1] and calling_base == error_locations[himut_index][1]:
                    modified_lines.append(line)
            if index % 1_000_000 == 0:
                for write_line in modified_lines:
                    fw.write(write_line)
                modified_lines.clear()
                print("processed {} records, {}/{}".format(index, himut_index, len(error_locations)))
    return

def go_through_and_get_high_qual_errors(read_path):
    file = open(read_path, "r")
    for index, line in enumerate(file):
        split_txt = line.split(" ")
        if len(split_txt) != 9:
            continue
        base_quality = int(split_txt[2])
        ref_base = split_txt[1][1]
        call_base = split_txt[3]
        if ref_base != call_base and base_quality >= 93:
            if random.choice(False, False, False, False, True):
                print(line)
        if index % 100000 == 0:
            print("Running line {}".format(index))
            break
    return

def print_pacbio_scores(read_path):
    # arrays to save the result
    error_counts = [0] * 94
    all_counts = [0] * 94
    file = open(read_path, "r")
    for index, line in enumerate(file):
        split_txt = line.split(" ")
        if len(split_txt) != 9:
            continue
        base_quality = int(split_txt[2])
        ref_base = split_txt[1][1]
        call_base = split_txt[3]
        all_counts[base_quality] += 1
        if ref_base != call_base:
            error_counts[base_quality] += 1
        if index % 100000 == 0:
            print("Running line {}".format(index))
            break
    print(error_counts)
    print(all_counts)
    print(read_path)
    return

def make_sub_array(error_lines, location):
    range = 100
    sub_error_array = []
    closest_error_value = 1000000
    closest_error_index = 0
    for index, error_line in enumerate(error_lines):
        if closest_error_value > abs(error_line[0] - location):
            closest_error_value = abs(error_line[0] - location)
            closest_error_index = index
    if len(error_lines) < closest_error_index + range:
        sub_error_array = error_lines[closest_error_index - range: len(error_lines)]
    elif closest_error_index < range:
        sub_error_array = error_lines[0: closest_error_index + range]
    else:
        sub_error_array = error_lines[closest_error_index - range: closest_error_index + range]
    sub_array_low = sub_error_array[0][0]
    sub_array_high = sub_error_array[len(sub_error_array) - 1][0]
    return sub_error_array, sub_array_low, sub_array_high

def make_unfiltered(read_path, error_path, write_path):
    # error list save
    error_lines = []
    error_count = 0
    last_error_location = 0
    modified_lines = []
    sub_error_array = []
    sub_array_low = -1
    sub_array_high = -1
    # open error file
    error_file = open(error_path, "r")
    for _, line in enumerate(error_file):
        split_txt = line.split(" ")
        location = int(split_txt[0])
        base = split_txt[1][0]
        error_lines.append((location, base))
    # open filtered data file
    read_file = open(read_path, 'r')
    with open(write_path, 'a') as fw:
        for index, line in enumerate(read_file):
            split_txt = line.split(" ")
            if len(split_txt) != 11:
                continue
            location = int(split_txt[0])
            base_context = split_txt[2]
            base_1 = split_txt[3]
            pac_qual = split_txt[4]
            base_2 = split_txt[5]
            total_count = split_txt[6]
            parallel1 = split_txt[7]
            parallel2 = split_txt[8]
            parallel3 = split_txt[9]
            parallel4 = split_txt[10]
            if not (location >= sub_array_low and location <= sub_array_high):
                sub_array, sub_array_low, sub_array_high = make_sub_array(error_lines, location)
            try:
                required_index = [y[0] for y in sub_array].index(location)
                if base_1 == sub_array[required_index][1]:
                    #if location == 37666995:
                    #    print("{} {}".format(base_1, sub_array[required_index]))
                    result = "true"
                    if last_error_location != location:
                        last_error_location = location
                        error_count += 1
                else:
                    result = "false"
            except ValueError:
                result = "false"
            modified_lines.append("{} {} {} {} {} {} {} {} {} {} {}".format(location, result, base_context, base_1, pac_qual, base_2, total_count, parallel1, parallel2, parallel3, parallel4))
            if index % 1_000_000 == 0:
                for write_line in modified_lines:
                    fw.write(write_line)
                modified_lines.clear()
                print("processed {} records, errors {}/{}".format(index, error_count, len(error_lines)))
    return

def old_format_to_new_format_converter(read_path, write_path):
    # array to save offsets
    modified_lines = []
    # indices to output
    write_index = 1
    read_index = 1
    # open files read and write files
    file = open(read_path, "r")
    with open(write_path, 'a') as fw:
        for index, line in enumerate(file):
            split_txt = line.split(" ")
            if len(split_txt) != 10:
                continue
            location = split_txt[0]
            base_context = split_txt[1]
            pac_qual = split_txt[2]
            base_1 = split_txt[3]
            total_count = split_txt[5]
            parallel1 = split_txt[6]
            parallel2 = split_txt[7]
            parallel3 = split_txt[8]
            parallel4 = split_txt[9]
            modified_lines.append("{} {} {} {} {} {} {} {} {}".format(location, base_context, pac_qual, base_1, total_count, parallel1, parallel2, parallel3, parallel4))
            if index % 1_000_000 == 0:
                for write_line in modified_lines:
                    fw.write(write_line)
                modified_lines.clear()
                print("indexed {} records".format(index))
    return

def check_and_clean_data (path):
    other_bigger_count = 0
    # open the file with ml data
    file = open(path, "r")
    # go line by line
    for index, line in enumerate(file):
        if index % 100000 == 0:
            print("current line {} other bigger count {}".format(index, other_bigger_count))
        split_txt = line.split(" ")
        if len(split_txt) != 11:
            print("line number {} is invalid, line: {}".format(index, line))
            continue
        #encoded_bases = self.one_hot_encoding_bases(split_txt[1][0]) + self.one_hot_encoding_bases(split_txt[1][1]) + self.one_hot_encoding_bases(split_txt[1][2])
        three_context_bases = [split_txt[2][0], split_txt[2][1], split_txt[2][2]]
        # get quality in float
        quality = float(split_txt[4]) / 100
        # get the num of parallel bases in float
        parallel_vec_s = [split_txt[7], split_txt[8], split_txt[9], split_txt[10]]
        char_remov = ["]", "[", ",", "\n"]
        for char in char_remov:
            for index_s in range(len(parallel_vec_s)):
                temp = parallel_vec_s[index_s].replace(char, "")
                parallel_vec_s[index_s] = temp
        parallel_vec_f = []
        for parallel in parallel_vec_s:
            parallel_vec_f.append(float(parallel))
        # rearrange so that the calling base num first and rest in decending order
        sorted_vec = rearrange_sort_parallel_bases(parallel_vec_f, split_txt[2][1])
        if sorted_vec[1] > sorted_vec[0] and split_txt[1] == "false":
            #print("line number {} is invalid, true with parallel higher line: {}".format(index, line))
            other_bigger_count += 1
    return

def rearrange_sort_parallel_bases(parallel_vec, base):
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
    else:
        selected_base = parallel_vec[0]
        del parallel_vec[0]
    parallel_vec.sort(reverse = True)
    parallel_vec = [selected_base] + parallel_vec
    return parallel_vec

def one_hot_encoding_bases(base):
    one_hot_base = [0.0] * 4
    if base == "A":
        one_hot_base[0] = 1.0
    elif base == "C":
        one_hot_base[1] = 1.0
    elif base == "G":
        one_hot_base[2] = 1.0
    elif base == "T":
        one_hot_base[3] = 1.0
    return one_hot_base

def list_corrected_errors_rust_input(read_path):
    result_array = np.zeros((4, 4, 16))
    # open the file with ml data
    file = open(read_path, "r")
    # go line by line
    for index, line in enumerate(file):
        split_txt = line.split(" ")
        ref_base = get_base_to_int(split_txt[0][1])
        alt_base = get_base_to_int(split_txt[2][0])
        first_in_three_base = get_base_to_int(split_txt[0][0])
        third_in_three_base = get_base_to_int(split_txt[0][2])
        one_three_value = first_in_three_base * 4 + third_in_three_base
        result_array[ref_base][alt_base][one_three_value] += 1
    print(result_array)
    return

def list_corrected_errors(read_path, write_path):
    # open the file with ml data
    file = open(read_path, "r")
    # go line by line
    for index, line in enumerate(file):
        split_txt = line.split(" ")
        # check if line is valid
        if len(split_txt) != 12:
            continue
        result = split_txt[1]
        # only check the false ones
        if result == 'false':
            parallel_vec_s = [split_txt[8], split_txt[9], split_txt[10], split_txt[11]]
            char_remov = ["]", "[", ",", "\n"]
            for char in char_remov:
                for index_s in range(len(parallel_vec_s)):
                    temp = parallel_vec_s[index_s].replace(char, "")
                    parallel_vec_s[index_s] = temp
            parallel_vec_f = []
            for parallel in parallel_vec_s:
                parallel_vec_f.append(float(parallel))
            # save required data
            calling_base = split_txt[2][1]
            calling_base_seq_num = 0
            top_base = "A"
            top_base_seq_num = max(parallel_vec_f)
            max_index = parallel_vec_f.index(max(parallel_vec_f))
            if max_index == 0:
                top_base = "A"
            elif max_index == 1:
                top_base = "C"
            elif max_index == 2:
                top_base = "G"
            elif max_index == 3:
                top_base = "T"

            if calling_base == "A":
                calling_base_seq_num = parallel_vec_f[0]
            elif calling_base == "C":
                calling_base_seq_num = parallel_vec_f[1]
            elif calling_base == "G":
                calling_base_seq_num = parallel_vec_f[2]
            elif calling_base == "T":
                calling_base_seq_num = parallel_vec_f[3]
            # if calling base is not equal to top base, error corrected? need reference to check
            if calling_base != top_base:
                print("well that was a waste of time")
    return

def index_file(read_path, write_path):
    # array to save offsets
    read_line_offset = [0]
    # indices to output
    write_index = 1
    read_index = 1
    # open files read and write files
    with open(read_path) as fr:
        with open(write_path, 'a') as fw:
            # go line by line in read
            read_line = fr.readline()
            while read_line:
                # get the current offset and add to array
                read_line_offset.append(fr.tell())
                # every million save array to file and clear array
                if read_index % 1_000_000 == 0:
                    for offset in read_line_offset:
                        write_line = "{:021d}\n".format(offset)
                        fw.write(write_line)
                        write_index += 1
                    read_line_offset.clear()
                    print("indexed {} records".format(read_index))
                read_index += 1
                read_line = fr.readline()
    return

def print_topology_cut_scores():
    # arrays to save the result
    error_counts = [0] * 93
    all_counts = [0] * 93
    path = "./data/train_file.txt"
    file = open(path, "r")
    for index, line in enumerate(file):
        split_txt = line.split(" ")
        if len(split_txt) != 11:
            continue
        result = split_txt[0]
        position = int(split_txt[6])
        all_counts[position] += 1
        if result == "true":
            error_counts[position] += 1
        if index % 100000 == 0:
            print("Running line {}".format(index))
    print(error_counts)
    print(all_counts)
    return

def calculate_binomial(n, k, prob):
    binomial_coeff = scipy.special.binom(n, k)
    success = np.power(prob, k)
    failure = np.power(1.00 - prob, n - k)
    return (binomial_coeff * success * failure)

def old_data_loader(path, start, length, get_random):
    file = open(path, "r")
    label_tensor = torch.empty((length, 1), dtype = torch.float32)
    input_tensor = torch.empty((length, 17), dtype = torch.float32)
    index = 0
    tensor_pos = 0
    for line in file:
        
        # only get the specified section in file
        if index < start:
            if line != "\n":
                index += 1
            continue
        # if random is required break when tensor is full
        elif get_random:
            if tensor_pos >= length:
                break
        # if random is not required then break when after the specified length
        elif index > start and index >= start + length:
            break
        # continue if random choice is true
        if get_random and random.choice([True, True, True, True, False]):
            continue
        if line != "\n":
            split_txt = line.split(" ")
            if len(split_txt) != 11:
                index += 1
                continue
            # get three base context in one hot encoded
            encoded_bases = one_hot_encoding_bases(split_txt[1][0]) + one_hot_encoding_bases(split_txt[1][1]) + one_hot_encoding_bases(split_txt[1][2])
            # get quality in float
            quality = float(split_txt[3]) / 100
            # get the num of parallel bases in float
            parallel_vec_s = [split_txt[7], split_txt[8], split_txt[9], split_txt[10]]
            char_remov = ["]", "[", ",", "\n"]
            for char in char_remov:
                for index_s in range(len(parallel_vec_s)):
                    temp = parallel_vec_s[index_s].replace(char, "")
                    parallel_vec_s[index_s] = temp
            parallel_vec_f = []
            for parallel in parallel_vec_s:
                parallel_vec_f.append(float(parallel))
            # rearrange so that the calling base num first and rest in decending order
            sorted_vec = rearrange_sort_parallel_bases(parallel_vec_f, split_txt[1][1])
            # make and append to the input tensor,
            input_tensor[tensor_pos] = torch.tensor([encoded_bases + [quality] + sorted_vec])
            # append to result tensor,
            result = split_txt[0]
            if result == "false":
                label_tensor[tensor_pos] = torch.tensor([[1.0]])
                # continue if we want errors and its not a error
            else:
                label_tensor[tensor_pos] = torch.tensor([[0.00]])
            
            if index % 10000 == 0 and get_random == False:
                print("Going through line: {} getting data point: {}/{}".format(index, index - start, length))
            elif tensor_pos % 100 == 0 and get_random == True:
                print("Going through line: {} getting data point: {}/{}".format(index, tensor_pos, length))
            index += 1
            tensor_pos += 1
    file.close()
    return (input_tensor, label_tensor)

def rearrange_sort_parallel_bases(parallel_vec, base):
    selected_base = parallel_vec[0]
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
    parallel_vec.sort(reverse = True)
    parallel_vec = [selected_base] + parallel_vec
    return parallel_vec

def one_hot_encoding_bases(base):
    one_hot_base = [0.0] * 4
    if base == "A":
        one_hot_base[0] = 1.0
    elif base == "C":
        one_hot_base[1] = 1.0
    elif base == "G":
        one_hot_base[2] = 1.0
    elif base == "T":
        one_hot_base[3] = 1.0
    return one_hot_base