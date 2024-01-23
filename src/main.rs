extern crate bio;
mod poabanded;
use bio::alignment::pairwise::banded::Aligner as BandedDP;
use crate::poabanded::Aligner;
use petgraph::{Graph, Directed, graph::NodeIndex};
use petgraph::visit::Topo;
use petgraph::Direction::Outgoing;
use std::cmp;
use petgraph::{Incoming, Direction};
use std::thread;
use std::fs::create_dir_all;
use std::{fs::OpenOptions, io::{prelude::*}};
use rust_htslib::{bam, bam::Read};

const PRINT_ALL: bool = false;
const USEPACBIODATA: bool = true;
const NUM_OF_ITER_FOR_PARALLEL: usize = 10;
const GAP_OPEN: i32 = -2;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const BAND_SIZE: i32 = 100;
const MAX_NODES_IN_POA: usize = 75_000;
const SKIP_SCORE: i32 = 6_000;

use std::io::{self};
use std::env;
use std::str;

fn main() {
    // get the arguments 1. num of threads 2. read bam
    let args: Vec<String> = env::args().collect();
    let num_of_threads = args[1].clone().parse::<usize>().unwrap();
    let read_file_dir = args[2].clone();
    let mut test = false;
    // check whether it is test file or not
    if read_file_dir == "./sample_files/test.ccs.bam" {
        test = true;
    }
    // read bam and run threads to poa, parallel output
    thread_runner(read_file_dir, num_of_threads, test);
}

fn thread_runner (read_file_dir: String, num_of_threads: usize, test: bool) {
    // change the locations if test
    let subread_loc;
    let ip_loc;
    let pw_loc;
    let sn_loc;
    if test {
        subread_loc = 9;
        ip_loc = 13;
        pw_loc = 15;
        sn_loc = 19;
    }
    else {
        subread_loc = 9;
        ip_loc = 12;
        pw_loc = 14;
        sn_loc = 18;
    }
    let mut threads_used = 0;
    let mut children = vec![];
    let mut bam = bam::Reader::from_path(&read_file_dir).unwrap();
    let mut read_count = 0;
    let mut subreads_vec: Vec<String> = vec![];
    let mut ip_vec: Vec<Vec<usize>> = vec![];
    let mut pw_vec: Vec<Vec<usize>> = vec![];
    let mut sn_vec: Vec<f32> = vec![];
    let mut last_sub_read_name: String = "".to_string();
    // bam iterator
    let mut bam_iter = bam.records().into_iter();
    let mut read_name_set = ("".to_string(), "".to_string(), "".to_string());
    let mut get_next_read = true;
    loop {
        if get_next_read {
            let record_option = bam_iter.next();
            match record_option {
                Some(x) => {
                    let record = x.unwrap();
                    read_name_set = ((String::from_utf8(record.seq().as_bytes()).unwrap()), (String::from_utf8(record.qual().to_vec()).unwrap()), (String::from_utf8(record.qname().to_vec()).unwrap()));
                },
                None => {
                    break;
                }
            }
            read_count += 1;
            get_next_read = false;
        }
        let mut input = String::new();
        match io::stdin().read_line(&mut input) {
            Ok(len) => if len == 0 {
                return;
            } else {
                // get all the required info
                let parts = input.split("\t").collect::<Vec<&str>>();
                let name_full = parts[0].to_string();
                let sub_read = parts[subread_loc].to_string();
                let mut temp_subread_ip_vec = vec![];
                let mut temp_subread_pw_vec = vec![];
                let mut temp_sn_vec = vec![];
                let mut ip_collection: Vec<&str> = parts[ip_loc].split(",").collect();
                ip_collection.remove(0);
                for ip in ip_collection {
                    temp_subread_ip_vec.push(ip.parse::<usize>().unwrap());
                }
                // process pw, parse to usize
                let mut pw_collection: Vec<&str> = parts[pw_loc].split(",").collect();
                pw_collection.remove(0);
                for pw in pw_collection {
                    temp_subread_pw_vec.push(pw.parse::<usize>().unwrap());
                }
                let mut sn_collection: Vec<&str> = parts[sn_loc].split(",").collect();
                sn_collection.remove(0);
                for sn in sn_collection {
                    temp_sn_vec.push(sn.parse::<f32>().unwrap());
                }
                // get the read name
                let subread_read_name = name_full.split("/").collect::<Vec<&str>>()[1];
                println!("{}", name_full);
                let current_read_name = read_name_set.2.split("/").collect::<Vec<&str>>()[1];
                println!("current read name {}  subread_name {}", current_read_name, subread_read_name);
                println!("previous subread name {} subread name {}", last_sub_read_name, subread_read_name);
                if current_read_name == subread_read_name {
                    println!("CURRENT MATCH ADDING");
                    // add data to vector
                    subreads_vec.push(sub_read);
                    pw_vec.push(temp_subread_pw_vec);
                    ip_vec.push(temp_subread_ip_vec);
                    sn_vec = temp_sn_vec;
                }
                else if current_read_name == last_sub_read_name {
                    println!("PREV MATCH PROCESSING");
                    println!("Main Thread: processing {} sub reads {}: read_number {}", current_read_name, subreads_vec.len(), read_count);
                    let name = read_name_set.2.clone();
                    let read = read_name_set.0.clone();
                    let quality = read_name_set.1.clone();
                    let subreads = subreads_vec.clone();
                    let ip_stuff = ip_vec.clone();
                    let pw_stuff = pw_vec.clone();
                    let sn_stuff = sn_vec.clone();
                    children.push(thread::spawn(move || {
                        one_function(name, read, quality, subreads, ip_stuff, pw_stuff, sn_stuff, threads_used);
                    }));
                    threads_used += 1;
                    if threads_used == num_of_threads {
                        println!("Main Thread: Waiting for threads to finish......");
                        for child_index in 0..children.len() {
                            let _ = children.pop().unwrap().join();
                            println!("Main Thread: Threads done {}/{}", child_index, num_of_threads);
                        }
                        threads_used = 0;
                    }
                    // clear the vector, add the data
                    subreads_vec.clear();
                    pw_vec.clear();
                    ip_vec.clear();
                    subreads_vec.push(sub_read.clone());
                    pw_vec.push(temp_subread_pw_vec.clone());
                    ip_vec.push(temp_subread_ip_vec.clone());
                    sn_vec = temp_sn_vec.clone();
                }
                else {
                    println!("NO MATCH CLEARING");
                    subreads_vec.clear();
                    pw_vec.clear();
                    ip_vec.clear();
                    subreads_vec.push(sub_read.clone());
                    pw_vec.push(temp_subread_pw_vec.clone());
                    ip_vec.push(temp_subread_ip_vec.clone());
                    sn_vec = temp_sn_vec.clone();
                }
                last_sub_read_name = subread_read_name.to_string();
            }
            Err(error) => {
                eprintln!("error: {}", error);
                return;
            }
        }
    }
}

fn one_function (file_name: String, read: String, quality: String, mut sub_reads: Vec<String>, mut ip_vec: Vec<Vec<usize>>, mut pw_vec: Vec<Vec<usize>>, sn_vec: Vec<f32>, thread_id: usize) {
    // graph!!
    // filter out the long reads and rearrange the reads
    (sub_reads, pw_vec, ip_vec) = reverse_complement_subreads_ip_pw(&sub_reads, pw_vec, ip_vec);
    // reverse if score is too low
    (sub_reads, pw_vec, ip_vec) = check_the_scores_and_change_alignment_subreads_pw_ip(sub_reads, pw_vec, ip_vec, &read);
    if sub_reads.len() == 0 {
        return
    }
    // put the read in first pos
    sub_reads.insert(0, read.clone());
    // do poa with the read and subreads, get the poa and consensus
    let mut sequence_number: usize = 0;
    let mut aligner = Aligner::new(MATCH, MISMATCH, GAP_OPEN, &sub_reads[0].as_bytes().to_vec(), BAND_SIZE);
    for sub_read in &sub_reads {
        if sequence_number != 0 {
            aligner.global(&sub_read.as_bytes().to_vec()).add_to_graph();
        }
        let node_num = aligner.graph().node_count();
        if node_num > MAX_NODES_IN_POA {
            println!("Thread {} : NUM OF NODES {} TOO BIG, SKIPPING", thread_id, node_num);
            return
        }
        sequence_number += 1;
        println!("Thread {} : Sequence {} processed",  thread_id, sequence_number);
    }
    let calculated_graph: &Graph<u8, i32, Directed, usize> = aligner.graph();
    // parallel bases!!
    let (calculated_consensus, calculated_topology) = get_consensus_from_graph(&calculated_graph); //just poa
    let parallel_bases_vec = get_consensus_parallel_bases(sub_reads.len(), &calculated_consensus, &calculated_topology, &calculated_graph, thread_id);
    // align all subreads to ccs
    println!("Thread {} :aligning stuff...", thread_id);
    sub_reads.remove(0);
    let (ip_vec, pw_vec) = align_subreads_to_ccs_read_calculate_avg_ip_pw(&read, sub_reads, ip_vec, pw_vec);
    // match the calculated consensus to the original consensus and get the required indices
    let calc_cons_id = get_redone_consensus_matched_positions(&read, &calculated_consensus);
    let quality_vec_chr = quality.as_bytes().to_vec();
    for (index, pacbio_base) in read.as_bytes().to_vec().iter().enumerate() {
        let char_sequence: Vec<char> = read.chars().collect::<Vec<_>>();
        let mut char_7base_context: Vec<char> = vec![];
        // when 7 base context above 0 this is very dumb re write
        if index >= 3 {
            char_7base_context.push(char_sequence[index - 3]);
        }
        else {
            char_7base_context.push('X');
        }
        if index >= 2 {
            char_7base_context.push(char_sequence[index - 2]);
        }
        else {
            char_7base_context.push('X');
        }
        if index >= 1 {
            char_7base_context.push(char_sequence[index - 1]);
        }
        else {
            char_7base_context.push('X');
        }
        char_7base_context.push(char_sequence[index]);
        // when len is greater than 7 base context
        if read.len() > (1 + index) {
            char_7base_context.push(char_sequence[index + 1]);
        }
        // when len is less than 7 base context
        else {
            char_7base_context.push('X');
        }
        if read.len() > (2 + index) {
            char_7base_context.push(char_sequence[index + 2]);
        }
        else{
            char_7base_context.push('X');
        }
        if read.len() > (3 + index) {
            char_7base_context.push(char_sequence[index + 3]);
        }
        else{
            char_7base_context.push('X');
        }
        let read_sevenbase_context = char_7base_context.iter().collect::<String>();
        let mut pacbio_str = format!("OK({})", *pacbio_base as char);
        let parallel_bases;
        if calc_cons_id[index].1 == 0 {
            parallel_bases = vec![1, 1, 1, 1]; //deletion
            pacbio_str = "DEL()".to_string();
        }
        else if calc_cons_id[index].1 == 1 {
            parallel_bases = parallel_bases_vec[calc_cons_id[index].0].clone(); //normal
        }
        else {
            parallel_bases = parallel_bases_vec[calc_cons_id[index].0].clone(); //subsitution the value corrospond to the sub
            pacbio_str = format!{"SB({})", calc_cons_id[index].1 as char};
        }
        // sn vec padding
        let sn_vec_str = format!("[{:0>7} {:0>7} {:0>7} {:0>7}]", sn_vec[0], sn_vec[1], sn_vec[2], sn_vec[3]);
        // parallel bases padding
        let parallel_bases_str = format!("[{:0>2} {:0>2} {:0>2} {:0>2}]", parallel_bases[0], parallel_bases[1], parallel_bases[2], parallel_bases[3]);
        let write_string = format!("{:0>2} : {:0>5} {:0>5} {} {} {} {:0>3} {:0>3} {}\n", (quality_vec_chr[index] - 33), index, read.len(), read_sevenbase_context, pacbio_str, sn_vec_str, ip_vec[index], pw_vec[index], parallel_bases_str);
        let write_path = format!("{}{}", "./intermediate/", file_name.replace("/", "."));
        write_string_to_file(&write_path, &write_string);
        //print!("Thread_ID {}: {}", thread_id, write_string);
    }
    return
}

pub fn write_string_to_file (file_name: &String, input_string: &String) {
    let path = std::path::Path::new(&file_name);
    let prefix = path.parent().unwrap();
    create_dir_all(prefix).unwrap();
    let mut file = OpenOptions::new().create(true).append(true).open(file_name).unwrap();
    write!(file, "{}", input_string).expect("result file cannot be written");
}

fn reverse_complement_subreads_ip_pw (original_subreads: &Vec<String>, mut pw_vec: Vec<Vec<usize>>, mut ip_vec: Vec<Vec<usize>>) -> (Vec<String>, Vec<Vec<usize>>, Vec<Vec<usize>>) {
    let mut seqvec: Vec<String> = vec![];
    //reverse complement every other line
    let mut index = 0;
    for seq in original_subreads {
        if index % 2 != 0 {
            pw_vec[index].reverse();
            ip_vec[index].reverse();
            let mut tempseq: Vec<char> = vec![];
            let iterator = seq.chars().rev().into_iter();
            for char in iterator{
                tempseq.push(match char {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    _ => ' ',
                });
            }
            seqvec.push(tempseq.iter().cloned().collect::<String>());
        }
        else {
            seqvec.push((*seq.clone()).to_string());
        }
        index += 1;
    }
    // original seqvec
    let seqvec_ori = seqvec.clone();
    //sort the vector by size
    seqvec.sort_by_key(|seq| seq.len());
    //drop the sequences which are > 1.8x median size
    let median_size: f32 = seqvec[(seqvec.len() / 2) - 1].len() as f32;
    let mut drop_index = seqvec.len();
    for index in (seqvec.len() / 2)..(seqvec.len() - 1) {
        if seqvec[index].len() as f32 > (median_size * 1.5) {
            drop_index = index;
            break;
        }
    }
    for _ in drop_index..seqvec.len() {
        seqvec.pop();
    }
    // rearrange the seq vector median first and rest according mediand size difference
    seqvec.sort_by(|a, b| ((a.len() as f32 - median_size).abs()).partial_cmp(&(b.len() as f32 - median_size).abs()).unwrap());
    let mut new_pw_vec: Vec<Vec<usize>> = vec![vec![]; seqvec.len()];
    let mut new_ip_vec: Vec<Vec<usize>> = vec![vec![]; seqvec.len()];
    for i in 0..seqvec_ori.len() {
        for j in 0..seqvec.len() {
            if seqvec_ori[i] == seqvec[j] {
                new_pw_vec[j] = pw_vec[i].clone();
                new_ip_vec[j] = ip_vec[i].clone();
            }
        } 
    }
    (seqvec, new_pw_vec, new_ip_vec)
}

fn check_the_scores_and_change_alignment_subreads_pw_ip (seqvec: Vec<String>, mut pw_vec: Vec<Vec<usize>>, mut ip_vec: Vec<Vec<usize>>, pacbio_consensus: &String) -> (Vec<String>, Vec<Vec<usize>>, Vec<Vec<usize>>) {
    let mut forward_score = 0;
    let mut backward_score = 0;
    // make the pacbio orientation files
    let pacbio_forward = pacbio_consensus.as_bytes().to_vec();
    let pacbio_backward;
    let mut tempseq: Vec<char> = vec![];
    let iterator = pacbio_consensus.chars().rev().into_iter();
    for char in iterator{
        tempseq.push(match char {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => ' ',
        });
    }
    pacbio_backward = tempseq.iter().cloned().collect::<String>().as_bytes().to_vec();
    // check the forward scores for 2 sequences
    for seq in &seqvec {
        let score_func = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let k = 8; // kmer match length
        let w = 20; // Window size for creating the band
        let mut aligner = BandedDP::new(-5, -1, score_func, k, w);
        let alignment = aligner.local(&pacbio_forward, &seq.as_bytes().to_vec());
        let score = alignment.score;
        //println!("forward score: {}", score);
        forward_score += score;
        break;
    }
    // check the backward scores for 2 sequences
    for seq in &seqvec {
        let score_func = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let k = 8; // kmer match length
        let w = 20; // Window size for creating the band
        let mut aligner = BandedDP::new(-5, -1, score_func, k, w);
        let alignment = aligner.local(&pacbio_backward, &seq.as_bytes().to_vec());
        let score = alignment.score;
        //println!("backward score: {}", score);
        backward_score += score;
        break;
    }
    if forward_score < SKIP_SCORE && backward_score < SKIP_SCORE {
        return (vec![], vec![], vec![]);
    }
    else if backward_score > forward_score {
        //println!("Scores are too low, inverting sequences.");
        let mut seqvec2 = vec![];
        //reverse complement every line
        for seq in &seqvec {
            let mut tempseq: Vec<char> = vec![];
            let iterator = seq.chars().rev().into_iter();
            for char in iterator{
                tempseq.push(match char {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    _ => ' ',
                });
            }
            seqvec2.push(tempseq.iter().cloned().collect::<String>());
        }
        for index in 0..pw_vec.len() {
            pw_vec[index].reverse();
        }
        for index in 0..ip_vec.len() {
            ip_vec[index].reverse();
        }
        return (seqvec2, pw_vec, ip_vec);
    }
    else {
        return (seqvec, pw_vec, ip_vec);
    }
}

fn get_consensus_from_graph(graph: &Graph<u8, i32, Directed, usize>) -> (Vec<u8>, Vec<usize>) {
    let mut output: Vec<u8> = vec![];
    let mut topopos: Vec<usize> = vec![];
    let mut topo = Topo::new(graph);
    let mut topo_indices = Vec::new();
    let mut max_index = 0;
    let mut max_score = 0.0;

    while let Some(node) = topo.next(graph) {
        topo_indices.push(node);
        if max_index < node.index(){
            max_index = node.index();
        }
    }
    topo_indices.reverse();
    //define score and nextinpath vectors with capacity of num nodes.
    let mut weight_scores: Vec<i32> = vec![0; max_index + 1];
    let mut scores: Vec<f64> = vec![0.0; max_index + 1];
    let mut next_in_path: Vec<usize> = vec![0; max_index + 1];
    //iterate thorugh the nodes in reverse
    for node in topo_indices{
        let mut best_weight_score_edge: (i32, f64, usize) = (-1 , -1.0, 123456789);
        let mut neighbour_nodes = graph.neighbors_directed(node, Outgoing);
        while let Some(neighbour_node) = neighbour_nodes.next() {
            let mut edges = graph.edges_connecting(node, neighbour_node);
            let mut weight: i32 = 0;
            while let Some(edge) = edges.next() {
                weight += edge.weight().clone();
            }
            let weight_score_edge = (weight, scores[neighbour_node.index()], neighbour_node.index());
            if weight_score_edge > best_weight_score_edge{
                best_weight_score_edge = weight_score_edge;
            }
        }
        //save score and traceback
        if best_weight_score_edge.0 as f64 + best_weight_score_edge.1 > max_score{
            max_score = best_weight_score_edge.0 as f64 + best_weight_score_edge.1;
        }
        scores[node.index()] = best_weight_score_edge.0 as f64 + best_weight_score_edge.1;
        next_in_path[node.index()] = best_weight_score_edge.2;
        weight_scores[node.index()] = best_weight_score_edge.0;
    }
    let mut pos = scores.iter().position(|&r| r == max_score).unwrap();
    //calculate the start weight score
    let mut consensus_started: bool = false;
    let weight_average = scores[pos] / scores.len() as f64;
    let weight_threshold = weight_average as i32 / 2; 
    while pos != 123456789 {
        //continue if starting weight score is too low
        if consensus_started == false && weight_scores[pos] < weight_threshold {
            pos = next_in_path[pos];
            continue;
        }
        //println!("current {} {}", pos, next_in_path[pos]);
        consensus_started = true;
        topopos.push(pos as usize);
        output.push(graph.raw_nodes()[pos].weight);
        pos = next_in_path[pos];
    }
    (output, topopos)
}

pub fn get_consensus_parallel_bases (seq_num: usize, consensus: &Vec<u8>, topology: &Vec<usize>, graph: &Graph<u8, i32, Directed, usize>, thread_id: usize) -> Vec<Vec<usize>> {
    let mut base_count_vec: Vec<Vec<usize>> = vec![];
    //run all the consensus through get indices
    for i in 0..consensus.len() {
        if i % 2000 == 0 {
            println!("Thread {}: progress {}%({} / {})", thread_id,  (i * 100) / consensus.len() , i, consensus.len());
        }
        // skip the indices which are in the passed consensus
        let skip_nodes: Vec<usize> = topology[0 .. i + 1].to_vec();
        // new method using topology cut
        let mut target_node_parent = None;
        let mut target_node_child = None;
        if i != 0{
            target_node_parent = Some(topology[i - 1]);
        }
        if i != consensus.len() - 1 {
            target_node_child = Some(topology[i + 1]);
        }
        let (parallel_nodes, parallel_num_incoming_seq, _) = get_parallel_nodes_with_topology_cut (skip_nodes, seq_num,  topology[i], target_node_parent, target_node_child, graph);
        let base_counts = get_base_counts (parallel_nodes, parallel_num_incoming_seq, graph);
        base_count_vec.push(base_counts);
    }
    base_count_vec
}

pub fn get_parallel_nodes_with_topology_cut (skip_nodes: Vec<usize>, total_seq: usize, target_node: usize, target_node_parent: Option<usize>, target_node_child: Option<usize>, graph: &Graph<u8, i32, Directed, usize>) -> (Vec<usize>, Vec<usize>, Vec<String>) {
    // vector initialization
    let mut debug_strings: Vec<String> = vec![];
    let mut topology = Topo::new(graph);
    let mut topologically_ordered_nodes = Vec::new();
    let mut parallel_nodes: Vec<usize> = vec![];
    let mut parallel_node_parents: Vec<usize> = vec![];
    let mut parallel_num_incoming_seq: Vec<usize> = vec![];
    let mut direction: Option<Direction> = None;
    let temp_string;
    //print stuff
    if PRINT_ALL {
        println!("NODE CHECKING FOR PARALLEL: {}, base {}", target_node, graph.raw_nodes()[target_node].weight as char);
    }

    // make a topologically ordered list
    while let Some(node) = topology.next(graph) {
        topologically_ordered_nodes.push(node.index());
    }
    // find the position of the target node, its child and parent in topology list
    let target_node_topological_position = topologically_ordered_nodes.iter().position(|&r| r == target_node).unwrap();
    let target_child_topological_position = match target_node_child { 
        Some(child_node) => {topologically_ordered_nodes.iter().position(|&r| r == child_node).unwrap()},
        None => {direction = Some(Incoming); target_node_topological_position}
    };
    let target_parent_topological_position = match target_node_parent { 
        Some(parent_node) => {topologically_ordered_nodes.iter().position(|&r| r == parent_node).unwrap()},
        None => {direction = Some(Outgoing); target_node_topological_position}
    };
    // choose a direction with the least amount of intermediate nodes
    if (direction == None) && (topologically_ordered_nodes[target_parent_topological_position..target_node_topological_position].len() > topologically_ordered_nodes[target_node_topological_position..target_child_topological_position].len()) {
        direction = Some(Outgoing);
    }
    else if direction == None {
        direction = Some(Incoming);
    }
    match direction {
        Some(x) => {
            if x == Incoming {
                temp_string = format!("Going backwards");
                if PRINT_ALL {
                    println!("{}", temp_string);
                }
                debug_strings.push(temp_string.clone());
            }
            else {
                temp_string = format!("Going forward");
                if PRINT_ALL {
                    println!("{}", temp_string);
                }
                debug_strings.push(temp_string.clone());
            }
        }
        None => {}
    }
    // check if the target node corrosponds with all the sequences
    let num_seq_through_target_base = find_the_seq_passing_through (target_node, graph);

    if num_seq_through_target_base == total_seq {
        parallel_nodes.push(target_node);
        if USEPACBIODATA {
            parallel_num_incoming_seq.push(num_seq_through_target_base - 1);
        }
        parallel_num_incoming_seq.push(num_seq_through_target_base);
        return (parallel_nodes, parallel_num_incoming_seq, debug_strings);
    }
    // go back skip_count and go forward skip_count + 3 and check if parent and child are before and after target_node_position,
    // iterate skip_count until all sequences are found, break on 5
    let mut seq_found_so_far = num_seq_through_target_base;
    let mut bubble_size = 1;
    while (seq_found_so_far < total_seq)  && (bubble_size < NUM_OF_ITER_FOR_PARALLEL) {
        let temp_debug_strings;
        (parallel_nodes, parallel_node_parents, parallel_num_incoming_seq, seq_found_so_far, temp_debug_strings) = move_in_direction_and_find_crossing_nodes (&skip_nodes, total_seq, direction.unwrap(), parallel_nodes, parallel_node_parents, parallel_num_incoming_seq, seq_found_so_far, target_node, bubble_size, &topologically_ordered_nodes, target_node_topological_position, graph);
        debug_strings = [debug_strings, temp_debug_strings].concat();
        bubble_size += 1;
    }
    if USEPACBIODATA {
        parallel_num_incoming_seq.push(num_seq_through_target_base - 1);
    }
    parallel_num_incoming_seq.push(num_seq_through_target_base);
    parallel_nodes.push(target_node);
    (parallel_nodes, parallel_num_incoming_seq, debug_strings)
}

fn find_the_seq_passing_through (target_node: usize, graph: &Graph<u8, i32, Directed, usize>) -> usize {
    let node_index = NodeIndex::new(target_node);
    //edges directed toward the base
    let incoming_nodes: Vec<NodeIndex<usize>> = graph.neighbors_directed(node_index, Incoming).collect();
    let mut incoming_weight = 0;
    for incoming_node in incoming_nodes {
        let mut edges = graph.edges_connecting(incoming_node, node_index);
        while let Some(edge) = edges.next() {
            incoming_weight += edge.weight().clone();
        }
    }
    //edges directed from the base
    let outgoing_nodes: Vec<NodeIndex<usize>> = graph.neighbors_directed(node_index, Outgoing).collect();
    let mut outgoing_weight = 0;
    for outgoing_node in outgoing_nodes {
        let mut edges = graph.edges_connecting(node_index, outgoing_node);
        while let Some(edge) = edges.next() {
            outgoing_weight += edge.weight().clone();
        }
    }
    cmp::max(outgoing_weight, incoming_weight) as usize
}

fn move_in_direction_and_find_crossing_nodes (skip_nodes: &Vec<usize>, total_seq: usize, direction: Direction, mut parallel_nodes: Vec<usize>, mut parallel_node_parents: Vec<usize>, mut parallel_num_incoming_seq: Vec<usize>, mut seq_found_so_far: usize, focus_node: usize, bubble_size: usize, topologically_ordered_nodes: &Vec<usize>, target_node_position: usize, graph: &Graph<u8, i32, Directed, usize>) -> (Vec<usize>, Vec<usize>, Vec<usize>, usize, Vec<String>) {
    let mut debug_strings: Vec<String> = vec![];
    let mut temp_string: String;
    // get a list of x back_iterations back nodes
    let back_nodes_list = get_xiterations_direction_nodes(direction, bubble_size, vec![], focus_node, graph);
    // get a list of all forward nodes 0..(back_iterations + 3) for all the back_nodes
    let mut edge_nodes_list: Vec<usize> = vec![];
    for back_node in &back_nodes_list {
        let temp_forward_list = get_direction_nodes (direction.opposite(), bubble_size + 3, vec![], *back_node, graph);
        for temp_forward_node in &temp_forward_list {
            if !edge_nodes_list.contains(temp_forward_node) {
                edge_nodes_list.push(*temp_forward_node);
            }
        }
    }
    temp_string = format!("Iteration: {} BackNodes: {:?} CheckNodes: {:?}", bubble_size, back_nodes_list, edge_nodes_list);
    if PRINT_ALL {
        println!("{}", temp_string);
    }
    debug_strings.push(temp_string.clone());
    
    // get the two slices of topologically_ordered_list back front
    let mut slice: Vec<Vec<usize>> = [topologically_ordered_nodes[0..target_node_position].to_vec(), topologically_ordered_nodes[target_node_position + 1..topologically_ordered_nodes.len()].to_vec()].to_vec();
    // for debugging
    if slice[0].len() > 10 {
        temp_string = format!("Back slice {:?}", slice[0][(slice[0].len() - 10)..slice[0].len()].to_vec());
        if PRINT_ALL {
            println!("{}", temp_string);
        }
        debug_strings.push(temp_string.clone());
    }
    else {
        temp_string = format!("Back slice {:?}", slice[0][0..slice[0].len()].to_vec());
        if PRINT_ALL {
            println!("{}", temp_string);
        }
        debug_strings.push(temp_string.clone());
    }
    if slice[1].len() > 10 {
        temp_string = format!("Front slice {:?}", slice[1][0..10].to_vec());
        if PRINT_ALL {
            println!("{}", temp_string);
        }
        debug_strings.push(temp_string.clone());
    }
    else {
        temp_string = format!("Front slice {:?}", slice[1][0..slice[1].len()].to_vec());
        if PRINT_ALL {
            println!("{}", temp_string);
        }
        debug_strings.push(temp_string.clone());
    }

    if direction == Outgoing {
        slice.reverse();
    }
    //iterate through edge nodes obtained
    for edge_node in &edge_nodes_list {
        // get the parents of the edge node
        let edge_node_parents = get_direction_nodes (direction, 1, vec![], *edge_node, graph);
        'parent_loop: for edge_node_parent in &edge_node_parents {
            // if the parent is in back section and node is in front section add to parallel nodes or if both parent and target is in intermediate add to parallel loop
            if slice[0].contains(edge_node_parent) && slice[1].contains(edge_node) && (*edge_node_parent != focus_node) {
                // edge node parent check
                if parallel_nodes.contains(edge_node) && parallel_node_parents.contains(edge_node_parent) {
                    // go through the parallel nodes and if there is a match check if the same parent and continue if so
                    for index in 0..parallel_nodes.len() {
                        if (parallel_nodes[index] == *edge_node) && (parallel_node_parents[index] == *edge_node_parent) {
                            continue 'parent_loop;
                        }
                    }
                }
                // target node front of parallel node check
                if direction == Incoming {
                    if get_direction_nodes(Outgoing, 4, vec![],  *edge_node, graph).contains(&focus_node) {
                        continue;
                    }
                    if skip_nodes.contains(edge_node) {
                        continue;
                    }
                }
                else {
                    if get_direction_nodes(Outgoing, 4, vec![], focus_node, graph).contains(&edge_node) {
                        continue;
                    }
                    if skip_nodes.contains(edge_node_parent) {
                        continue;
                    }
                }
                // all found 
                if seq_found_so_far >= total_seq {
                    break;
                }
                parallel_nodes.push(*edge_node);
                parallel_node_parents.push(*edge_node_parent);
                temp_string = format!("success node {} parent/child {}\n", *edge_node, *edge_node_parent);
                if PRINT_ALL {
                    println!("{}", temp_string);
                }
                debug_strings.push(temp_string.clone());
                // get the edge weight and add to seq_found_so_far
                let mut incoming_weight = 0;
                if direction == Incoming {
                    let mut edges = graph.edges_connecting(NodeIndex::new(*edge_node_parent), NodeIndex::new(*edge_node));
                    while let Some(edge) = edges.next() {
                        incoming_weight += edge.weight().clone();
                    }
                }
                else {
                    let mut edges = graph.edges_connecting(NodeIndex::new(*edge_node), NodeIndex::new(*edge_node_parent));
                    while let Some(edge) = edges.next() {
                        incoming_weight += edge.weight().clone();
                    }
                }
                parallel_num_incoming_seq.push(incoming_weight as usize);
                seq_found_so_far += incoming_weight as usize;
            }
            
        }
    }
    (parallel_nodes, parallel_node_parents, parallel_num_incoming_seq, seq_found_so_far, debug_strings)
}

fn get_direction_nodes (direction: Direction, iteration: usize, mut direction_node_list: Vec<usize>, focus_node: usize, graph: &Graph<u8, i32, Directed, usize>) -> Vec<usize> {
    //forward outgoing
    //backward incoming
    if iteration <= 0 {
        return direction_node_list;
    }
    //get the back nodes of the target
    let mut direction_neighbours = graph.neighbors_directed(NodeIndex::new(focus_node), direction);
    //iterate through the neighbours
    while let Some(direction_neighbour) = direction_neighbours.next() {
        if !direction_node_list.contains(&direction_neighbour.index()){
            direction_node_list.push(direction_neighbour.index());
            direction_node_list = get_direction_nodes (direction, iteration - 1, direction_node_list, direction_neighbour.index(), graph);
        }
    }
    direction_node_list
}

fn get_xiterations_direction_nodes (direction: Direction ,iteration: usize, mut direction_node_list: Vec<usize>, focus_node: usize, graph: &Graph<u8, i32, Directed, usize>) -> Vec<usize> {
    if iteration <= 0 {
        return direction_node_list;
    }
    //get the back nodes of the target
    let mut direction_neighbours = graph.neighbors_directed(NodeIndex::new(focus_node), direction);
    //iterate through the neighbours
    while let Some(direction_neighbour) = direction_neighbours.next() {
        if iteration == 1 {
            if !direction_node_list.contains(&direction_neighbour.index()){
                direction_node_list.push(direction_neighbour.index());
            }
        }
        direction_node_list = get_xiterations_direction_nodes (direction, iteration - 1, direction_node_list, direction_neighbour.index(), graph);
    }
    direction_node_list
}

fn get_base_counts(indices_of_parallel_nodes: Vec<usize>, seq_through_parallel_nodes: Vec<usize>, graph: &Graph<u8, i32, Directed, usize>) -> Vec<usize> {
    let base_counts: Vec<usize>;
    let mut base_a_count = 0;
    let mut base_c_count = 0;
    let mut base_g_count = 0;
    let mut base_t_count = 0;
    //find out how many sequences run through each base
    //match the indices to the base and ++
    for index in 0..indices_of_parallel_nodes.len() {
        match graph.raw_nodes()[indices_of_parallel_nodes[index]].weight {
            65 => {
                base_a_count += seq_through_parallel_nodes[index];
            },
            67 => {
                base_c_count += seq_through_parallel_nodes[index];
            },
            71 => {
                base_g_count += seq_through_parallel_nodes[index];
            },
            84 => {
                base_t_count += seq_through_parallel_nodes[index];
            },
            _ => {
                //nothing
                },
        }
    }
    // save the base counts for debug
    base_counts = [base_a_count, base_c_count, base_g_count, base_t_count].to_vec();
    base_counts
}

fn align_subreads_to_ccs_read_calculate_avg_ip_pw(pacbio_ccs_str: &String, subread_vec: Vec<String>, ip_vec: Vec<Vec<usize>>, pw_vec: Vec<Vec<usize>>) -> (Vec<usize>, Vec<usize>) {
    let mut pacbio_ccs_ip_vec: Vec<usize> = vec![0; pacbio_ccs_str.len()];
    let mut pacbio_ccs_pw_vec: Vec<usize> = vec![0; pacbio_ccs_str.len()];

    let pacbio_ccs: Vec<u8> = pacbio_ccs_str.bytes().collect();
    let mut current_sub_read = 0;
    for subread_str in subread_vec {
        let subread: Vec<u8> = subread_str.bytes().collect();
        let score_func = |a: u8, b: u8| if a == b { 2i32 } else { -2i32 };
        let k = 8; // kmer match length
        let w = 20; // Window size for creating the band
        let mut aligner = BandedDP::new(-2, -2, score_func, k, w);
        let alignment = aligner.global(&pacbio_ccs, &subread);
        let mut pacbio_index = alignment.xstart;
        let mut subread_index = alignment.ystart;
        for op in alignment.operations {
            match op {
                bio::alignment::AlignmentOperation::Match => {
                    pacbio_ccs_ip_vec[pacbio_index] += ip_vec[current_sub_read][subread_index];
                    pacbio_ccs_pw_vec[pacbio_index] += pw_vec[current_sub_read][subread_index];
                    pacbio_index += 1;
                    subread_index += 1;
                },
                bio::alignment::AlignmentOperation::Subst => {
                    pacbio_ccs_ip_vec[pacbio_index] += ip_vec[current_sub_read][subread_index];
                    pacbio_ccs_pw_vec[pacbio_index] += pw_vec[current_sub_read][subread_index];
                    pacbio_index += 1;
                    subread_index += 1;
                },
                bio::alignment::AlignmentOperation::Del => {
                    subread_index += 1;
                },
                bio::alignment::AlignmentOperation::Ins => {
                    pacbio_index += 1;
                },
                _ => {},
            }
        }
        current_sub_read += 1;
    }
    // average the sum vectors
    for index in 0..pacbio_ccs_str.len() {
        pacbio_ccs_ip_vec[index] = pacbio_ccs_ip_vec[index] / current_sub_read;
        pacbio_ccs_pw_vec[index] = pacbio_ccs_pw_vec[index] / current_sub_read;
    }
    (pacbio_ccs_ip_vec, pacbio_ccs_pw_vec)
}

fn get_redone_consensus_matched_positions (pacbio_consensus: &String, calculated_consensus: &Vec<u8>) -> Vec<(usize, u8)> {
    let mut consensus_matched_indices: Vec<(usize, u8)> = vec![];
    let pacbio_consensus_vec: Vec<u8> = pacbio_consensus.bytes().collect();
    let score_func = |a: u8, b: u8| if a == b { 2i32 } else { -2i32 };
    let k = 8; // kmer match length
    let w = 20; // Window size for creating the band
    let mut aligner = BandedDP::new(-2, -2, score_func, k, w);
    let alignment = aligner.global(&calculated_consensus, &pacbio_consensus_vec);
    let temp_score = alignment.score;
    //let (alignment, temp_score) = pairwise(&calculated_consensus, &pacbio_consensus_vec, MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND, 0);
    let mut calc_index = 0;
    for op in alignment.operations {
        match op {
            bio::alignment::AlignmentOperation::Match => {
                consensus_matched_indices.push((calc_index, 1));
                calc_index += 1;
            },
            bio::alignment::AlignmentOperation::Subst => {
                let base = calculated_consensus[calc_index];
                consensus_matched_indices.push((calc_index, base));
                calc_index += 1;
            },
            bio::alignment::AlignmentOperation::Del => {
                consensus_matched_indices.push((calc_index, 0));
            },
            bio::alignment::AlignmentOperation::Ins => {
                calc_index += 1;
            },
            _ => {},
        }
    }
    println!("score {}", temp_score);
    consensus_matched_indices
}

pub fn get_required_start_end_positions_from_read (section_length: usize, current_ref_pos: usize, current_read_pos: usize, required_pos: usize, required_len: usize) -> (usize, usize) {
    let p1;
    let p2;
    // case when both end partial match is required
    if (current_ref_pos < required_pos) && (current_ref_pos + section_length > required_pos + required_len) {
        p1 = required_pos - current_ref_pos + current_read_pos;
        p2 = current_read_pos + required_pos + required_len - current_ref_pos;
    }
    // case when start partial matches are required
    else if (current_ref_pos < required_pos) && (current_ref_pos + section_length >= required_pos) {
        p1 = required_pos - current_ref_pos + current_read_pos;
        p2 = current_read_pos + section_length;
    }
    // case when end partial matches are required
    else if (current_ref_pos <= required_pos + required_len) && (current_ref_pos + section_length > required_pos + required_len) {
        p1 = current_read_pos;
        p2 = current_read_pos + required_pos + required_len - current_ref_pos;
    }
    // normal case
    else {
        p1 = current_read_pos;
        p2 = current_read_pos + section_length;
    }
    (p1, p2)
}