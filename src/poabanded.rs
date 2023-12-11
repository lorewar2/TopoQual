
//obtained from rustbio library and modified//
use std::cmp::{max, Ordering};

use petgraph::Direction::Outgoing;
use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;

use petgraph::{Directed, Graph, Incoming};

pub const MIN_SCORE: i32 = -858_993_459; // negative infinity; see alignment/pairwise/mod.rs
pub type POAGraph = Graph<u8, i32, Directed, usize>;

// Unlike with a total order we may have arbitrary successors in the
// traceback matrix. I have not yet figured out what the best level of
// detail to store is, so Match and Del operations remember In and Out
// nodes on the reference graph.
#[derive(Debug, Clone)]
pub enum AlignmentOperation {
    Match(Option<(usize, usize)>),
    Del(Option<(usize, usize)>),
    Ins(Option<usize>),
}


pub struct Alignment {
    pub score: i32,
    //    xstart: Edge,
    pub operations: Vec<AlignmentOperation>,
}

#[derive(Debug, Clone)]
pub struct TracebackCell {
    score: i32,
    op: AlignmentOperation,
}

impl Ord for TracebackCell {
    fn cmp(&self, other: &TracebackCell) -> Ordering {
        self.score.cmp(&other.score)
    }
}

impl PartialOrd for TracebackCell {
    fn partial_cmp(&self, other: &TracebackCell) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for TracebackCell {
    fn eq(&self, other: &TracebackCell) -> bool {
        self.score == other.score
    }
}

//impl Default for TracebackCell { }

impl Eq for TracebackCell {}

pub struct Traceback {
    cols: usize,

    // store the last visited node in topological order so that
    // we can index into the end of the alignment when we backtrack
    last: NodeIndex<usize>,
    matrix: Vec<(Vec<TracebackCell>, usize, usize)>,
}

impl Traceback {
    /// Create a Traceback matrix with given maximum sizes
    ///
    /// # Arguments
    ///
    /// * `m` - the number of nodes in the DAG
    /// * `n` - the length of the query sequence
    fn with_capacity(m: usize, n: usize, gap_open: i32) -> Self {
        // each row of matrix contain start end position and vec of traceback cells
        let mut matrix: Vec<(Vec<TracebackCell>, usize, usize)> = vec![(vec![], 0, n); m + 1];
        for j in 0..=n {
            matrix[0].0.push(TracebackCell {
                score: (j as i32) * gap_open,
                op: AlignmentOperation::Ins(None),
            });
        }
        matrix[0].0[0] = TracebackCell {
            score: 0,
            op: AlignmentOperation::Match(None),
        };
        Traceback {
            cols: n,
            last: NodeIndex::new(0),
            matrix,
        }
    }
    fn new_row(&mut self, row: usize, size: usize, gap_open: i32, start: usize, end: usize) {
        //println!("row {} start {} end {}", row, start, end);
        self.matrix[row].1 = start;
        self.matrix[row].2 = end;
        if start == 0 {
            self.matrix[row].0.push(TracebackCell {
                score: (row as i32) * gap_open,
                op: AlignmentOperation::Del(None),
            });
        }
        else {
            self.matrix[row].0.push(TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Match(None),
            });
        }
        for _ in 1..=size {
            self.matrix[row].0.push(TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Match(None),
            });
        }
    }
    fn new() -> Self {
        Traceback {
            cols: 0,
            last: NodeIndex::new(0),
            matrix: Vec::new(),
        }
    }

    fn set(&mut self, i: usize, j: usize, cell: TracebackCell) {
        // if j is less than start, if j is greater than end, do nothing
        if !(self.matrix[i].1 > j || self.matrix[i].2 <= j) {
            let real_position = j - self.matrix[i].1;
            self.matrix[i].0[real_position] = cell;
        }
    }

    fn get(&self, i: usize, j: usize) -> &TracebackCell {
        if !(self.matrix[i].1 > j || self.matrix[i].2 <= j) {
            let real_position = j - self.matrix[i].1;
            return &self.matrix[i].0[real_position];
        }
        else if j == 0 {
            return &TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Del(None),
            };
        }
        else if j > self.matrix[i].2 {
            return &TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Ins(None),
            };
        }
        else {
            return &TracebackCell {
                score: MIN_SCORE,
                op: AlignmentOperation::Match(None),
            };
        }
    }
    pub fn alignment(&self) -> Alignment {
        // optimal AlignmentOperation path
        let mut ops: Vec<AlignmentOperation> = vec![];
        
        // Now backtrack through the matrix to construct an optimal path
        let mut i = self.last.index() + 1;
        let mut j = self.cols;
        let mut test_vector: Vec<usize> = vec![0; 6];
        while i > 0 || j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            ops.push(self.get(i,j).op.clone());
            match self.get(i,j).op {
                AlignmentOperation::Match(Some((p, _))) => {
                    i = p + 1;
                    j -= 1;
                    test_vector[0] += 1;
                    //println!("MATCHSOME");
                }
                AlignmentOperation::Del(Some((p, _))) => {
                    i = p + 1;
                    test_vector[1] += 1;
                    //println!("DELSOM");
                }
                AlignmentOperation::Ins(Some(p)) => {
                    i = p + 1;
                    j -= 1;
                    test_vector[2] += 1;
                    //println!("INSSOM");
                }
                AlignmentOperation::Match(None) => {
                    i -= 1; // break;
                    j -= 1;
                    test_vector[3] += 1;
                    //println!("MATCHNON");
                }
                AlignmentOperation::Del(None) => {
                    i -= 1; // j -= 1;
                    test_vector[4] += 1;
                    //println!("DELNON");
                }
                AlignmentOperation::Ins(None) => {
                    j -= 1; // i -= 1;
                    test_vector[5] += 1;
                    //println!("INSNON");
                }
            }
        }
        //println!("MS DS IS MN DN IN {:?}", test_vector);
        ops.reverse();
        Alignment {
            score: self.get(self.last.index() + 1, self.cols).score,
            operations: ops,
        }
    }
}

/// A partially ordered aligner builder
///
/// Uses consuming builder pattern for constructing partial order alignments with method chaining
pub struct Aligner {
    traceback: Traceback,
    query: Vec<u8>,
    pub poa: BandedPoa,
}

impl Aligner {
    pub fn new(match_score: i32, mismatch_score: i32, gap_open_score: i32, reference: &Vec<u8>, band_size: i32) -> Self {
        Aligner {
            traceback: Traceback::new(),
            query: reference.to_vec(),
            poa: BandedPoa::from_string(match_score, mismatch_score, gap_open_score, reference, band_size),
        }
    }

    /// Add the alignment of the last query to the graph.
    pub fn add_to_graph(&mut self) -> &mut Self { 
        let alignment = self.traceback.alignment();
        self.poa.add_alignment(&alignment, &self.query);
        self
    }
    /// Globally align a given query against the graph.
    pub fn global(&mut self, query: &Vec<u8>) -> &mut Self {
        self.query = query.to_vec();
        self.traceback = self.poa.global(query);
        let mut topo =  Topo::new(&self.poa.graph);
        let mut topo_indices = Vec::new();
        while let Some(node) = topo.next(&self.poa.graph) {
            topo_indices.push(node);
        }
        self
    }
    pub fn graph(&self) -> &POAGraph {
        &self.poa.graph
    } 
}

/// A partially ordered alignment graph
///
/// A directed acyclic graph datastructure that represents the topology of a
/// traceback matrix.
pub struct BandedPoa {
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    band_size: i32,
    pub graph: POAGraph,
}

impl BandedPoa {
    pub fn from_string(match_score: i32, mismatch_score: i32, gap_open_score: i32, seq: &Vec<u8>, band_size: i32) -> Self {
        let mut graph: Graph<u8, i32, Directed, usize> =
            Graph::with_capacity(seq.len(), seq.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(seq[0]);
        let mut node: NodeIndex<usize>;
        for base in seq.iter().skip(1) {
            node = graph.add_node(*base);
            graph.add_edge(prev, node, 1);
            prev = node;
        }
        BandedPoa { match_score: match_score, mismatch_score: mismatch_score, gap_open_score: gap_open_score, graph, band_size}
    }

    /// A global Needleman-Wunsch aligner on partially ordered graphs.
    ///
    /// # Arguments
    /// * `query` - the query TextSlice to align against the internal graph member
    pub fn global(&self, query: &Vec<u8>) -> Traceback {
        assert!(self.graph.node_count() != 0);
        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        let mut traceback = Traceback::with_capacity(m, n, self.gap_open_score);
        // construct the score matrix (O(n^2) space)
        let mut topo = Topo::new(&self.graph);
        let mut topo_index = 0;
        // the band required nodes
        let mut band_required_node: Vec<(usize, Vec<usize>)> = vec![];
        while let Some(node) = topo.next(&self.graph) {
            let mut start = 0;
            let mut end = query.len();
            if (topo_index != 0) && (self.band_size != 0) {
                let position_option = band_required_node.iter().position(|r| r.0 == node.index());
                match position_option {
                    // find the node in the band_required_node
                    Some(x) => {
                        let mut min = usize::MAX;
                        let mut max = usize::MIN;
                        // get the min and max required query position
                        for query_pos in &band_required_node[x].1 {
                            if query_pos < &min {
                                min = *query_pos;
                            }
                            if query_pos > &max {
                                max = *query_pos;
                            }
                        }
                        // calculate start and end values
                        if min > self.band_size as usize {
                            start = min - self.band_size as usize;
                        }
                        else {
                            start = 0;
                        }
                        end = max + self.band_size as usize;
                        if topo_index > m - (m / 10) {
                            if end < topo_index + 1 {
                                end = topo_index + 1;
                            }
                        }
                        band_required_node.remove(x);
                    },
                    // if not found do nothing
                    None => {}
                }
            }
            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1;
            //println!("query len {}", query.len());
            traceback.new_row(i, (end - start) + 1, self.gap_open_score, start, end);
            traceback.last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> = self.graph.neighbors_directed(node, Incoming).collect();
            let mut score_position: Vec<(i32, usize)> = vec![];
            let mut max_score = MIN_SCORE;
            for (j_p, q) in query.iter().enumerate().skip(start) {
                let j = j_p + 1;
                if topo_index != 0 {
                    if j_p > end {
                        break;
                    }
                }
                // match and deletion scores for the first reference base
                let max_cell = if prevs.is_empty() {
                    let temp_score;
                    if r == *q {
                        temp_score = self.match_score;
                    }
                    else {
                        temp_score = self.mismatch_score;
                    }
                    TracebackCell {
                        score: traceback.get(0, j - 1).score + temp_score,
                        op: AlignmentOperation::Match(None),
                    }
                } else {
                    let mut max_cell = TracebackCell {
                        score: MIN_SCORE,
                        op: AlignmentOperation::Match(None),
                    };
                    for prev_node in &prevs {
                        let i_p: usize = prev_node.index() + 1; // index of previous node
                        let temp_score;
                        if r == *q {
                            temp_score = self.match_score;
                        }
                        else {
                            temp_score = self.mismatch_score;
                        }
                        max_cell = max(
                            max_cell,
                            max(
                                TracebackCell {
                                    score: traceback.get(i_p, j - 1).score
                                        + temp_score,
                                    op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                },
                                TracebackCell {
                                    score: traceback.get(i_p, j).score + self.gap_open_score,
                                    op: AlignmentOperation::Del(Some((i_p - 1, i))),
                                },
                            ),
                        );
                    }
                    max_cell
                };
    
                let score = max(
                    max_cell,
                    TracebackCell {
                        score: traceback.get(i, j - 1).score + self.gap_open_score,
                        op: AlignmentOperation::Ins(Some(i - 1)),
                    },
                );
                if score.score > max_score {
                    max_score = score.score;
                    score_position.push((score.score, j_p));
                }
                traceback.set(i, j, score);
            }
            //println!("");
            topo_index += 1;
            // find connected nodes from this node, outgoing and save in band required node vec
            let mut neighbour_nodes = self.graph.neighbors_directed(node, Outgoing);
            if self.band_size != 0 {
                score_position.sort_by(|a, b| b.0.cmp(&a.0));
                let mut index = 0;
                let mut temp_vec = vec![];
                for temp in score_position {
                    temp_vec.push(temp.1);
                    index += 1;
                    if index > 20 {
                        break;
                    }
                }
                while let Some(neighbour_node) = neighbour_nodes.next() {
                    let position_option = band_required_node.iter().position(|r| r.0 == neighbour_node.index());
                    match position_option {
                        Some(x) => {
                            band_required_node[x].1 = [band_required_node[x].1.clone(), temp_vec.clone()].concat();
                        },
                        None => {
                            band_required_node.push((neighbour_node.index(), temp_vec.clone()));
                        }
                    }
                }
            }
        }
        traceback
    }
    
    /// Incorporate a new sequence into a graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment of the new sequence to the graph
    /// * `seq` - The sequence being incorporated
    pub fn add_alignment(&mut self, aln: &Alignment, seq: &Vec<u8>) {
        let head = Topo::new(&self.graph).next(&self.graph).unwrap();
        let mut prev: NodeIndex<usize> = NodeIndex::new(head.index());
        let mut i: usize = 0;
        let mut edge_not_connected: bool = false;
        for op in aln.operations.iter() {
            match op {
                AlignmentOperation::Match(None) => {
                    let node: NodeIndex<usize> = NodeIndex::new(0);
                    if (seq[i] != self.graph.raw_nodes()[head.index()].weight) && (seq[i] != b'X') {
                        let node = self.graph.add_node(seq[i]);
                        prev = node;
                    }
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                        edge_not_connected = false;
                    }
                    i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(*p);
                    if (seq[i] != self.graph.raw_nodes()[*p].weight) && (seq[i] != b'X') {
                        let node = self.graph.add_node(seq[i]);
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                    } else {
                        // increment node weight
                        match self.graph.find_edge(prev, node) {
                            Some(edge) => {
                                *self.graph.edge_weight_mut(edge).unwrap() += 1;
                            }
                            None => {
                                if prev.index() != head.index() {
                                    self.graph.add_edge(prev, node, 1);
                                }
                            }
                        }
                        prev = NodeIndex::new(*p);
                    }
                    i += 1;
                }
                AlignmentOperation::Ins(None) => {
                    let node = self.graph.add_node(seq[i]);
                    if edge_not_connected {
                        self.graph.add_edge(prev, node, 1);    
                    }
                    prev = node;
                    edge_not_connected = true;
                    i += 1;
                }
                AlignmentOperation::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    self.graph.add_edge(prev, node, 1);
                    prev = node;
                    i += 1;
                }
                AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes
            }
        }
    }
}
