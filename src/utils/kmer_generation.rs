use bstr::ByteVec;
use handlegraph::handle::{Direction, Handle};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use itertools::Itertools;
use serde::Deserialize;
use serde::Serialize;
use std::cmp::min;
use substring::Substring;

#[derive(Clone, Copy, Eq, Debug, Serialize, Deserialize, Ord, PartialOrd, PartialEq)]
/// Represents the orientation of a SeqPos. (Forward is 0 in dozyg, Reverse is 1)
pub enum SeqOrient {
    Forward,
    Reverse,
}

/// Represents an oriented position, mostly used to represent positions
/// on the forward/reverse linearization of the kmers.
#[derive(Clone, Copy, Eq, Debug, Serialize, Deserialize, Ord, PartialOrd, PartialEq)]
pub struct SeqPos {
    pub orient: SeqOrient,
    pub position: u64,
}

impl SeqPos {
    pub fn new(orient: SeqOrient, position: u64) -> Self {
        SeqPos { orient, position }
    }

    pub fn new_from_bool(is_reverse: bool, position: u64) -> Self {
        let orient: SeqOrient = match is_reverse {
            false => SeqOrient::Forward,
            true => SeqOrient::Reverse,
        };
        SeqPos::new(orient, position)
    }
}

/// Represents a kmer in the graph
#[derive(Debug, Clone, PartialEq)]
pub struct GraphKmer {
    /// The sequence of the kmer
    pub(crate) seq: String,
    /// The start position relative to the first_handle
    pub(crate) begin_offset: SeqPos,
    /// The end position relative to the last_handle
    pub(crate) end_offset: SeqPos,
    /// All the handles that were used to create the kmer
    pub(crate) all_handles: Vec<Handle>,
    /// The first handle of the kmer
    pub(crate) first_handle: Handle,
    /// The last handle of the kmer
    pub(crate) last_handle: Handle,
    /// The orientation of the handles
    pub(crate) handle_orient: bool,
    /// The number of forks the kmer has been involved in
    pub(crate) forks: u64,
}

/*
impl PartialEq for Kmer {
    fn eq(&self, other: &Self) -> bool {
        self.seq == other.seq && self.begin == other.begin && self.end == other.end
        && self.first == other.first && self.last == other.last &&
        self.handle_orient == other.handle_orient
    }
}
 */

impl GraphKmer {
    /// This function allows for a kmer to be extended, used when the
    /// kmer is on more than one handle
    fn extend_kmer(&mut self, new_seq: String, new_handle: Handle) {
        self.seq.push_str(&new_seq);
        self.end_offset = SeqPos::new_from_bool(new_handle.is_reverse(), new_seq.len() as u64);
        self.last_handle = new_handle;
        self.all_handles.push(new_handle);
    }

    fn add_handle_to_complete(&mut self, new_handle: Handle) {
        self.last_handle = new_handle;
    }
}

/// Generate kmers having size k from a given HashGraph. Both the number of overall visited edges
/// and maximum (outgoing) degree for a node can be limited, which should help with bigger graphs.
pub fn generate_kmers_from_graph(
    graph: &HashGraph,
    k: u64,
    edge_max: Option<u64>,
    degree_max: Option<u64>,
) -> Vec<GraphKmer> {
    let mut complete_kmers: Vec<GraphKmer> = Vec::new();

    let sorted_graph_handles: Vec<Handle> = graph.handles_iter().sorted().collect();

    for forward_handle in sorted_graph_handles {
        let mut handle: Handle;
        let mut orient: bool;

        // Kmers will be generated from both forward and reverse, as handlegraph is capable
        // of storing both
        for handle_orient in &[true, false] {
            match handle_orient {
                true => {
                    handle = forward_handle;
                    orient = true;
                }
                false => {
                    handle = forward_handle.flip();
                    orient = false;
                }
            }

            // Check if the handle has more outgoing edges than the maximum allowed,
            // and if that's the case, skip the current handle
            if let Some(degree_max) = degree_max {
                let mut curr_count: u64 = 0;
                graph
                    .handle_edges_iter(handle, Direction::Right)
                    .for_each(|_| {
                        curr_count += 1;
                    });
                if curr_count > degree_max {
                    // Skip current orientation
                    continue;
                }
            }

            // Get current node/handle
            let mut handle_seq = graph.sequence(handle).into_string_lossy();
            let mut handle_length = handle_seq.len() as u64;

            // This will store kmers that have not yet reached size k, which
            // will have to be completed from the neighbours of the handle
            let mut incomplete_kmers: Vec<GraphKmer> = Vec::new();

            // Try generating the "internal" kmers from the given handle
            for i in 0..handle_length {
                let begin = i;
                let end = min(i + k, handle_length);
                let kmer = GraphKmer {
                    seq: handle_seq
                        .substring(begin as usize, end as usize)
                        .to_string(),
                    begin_offset: SeqPos::new_from_bool(handle.is_reverse(), begin),
                    end_offset: SeqPos::new_from_bool(handle.is_reverse(), end),
                    all_handles: vec![handle],
                    first_handle: handle,
                    last_handle: handle,
                    handle_orient: orient,
                    forks: 0,
                };

                // Ignore Ns in kmer generation
                if kmer.seq.contains('N') {
                    continue;
                }

                // If the kmer has already reached size k...
                // NOTE: this implies that the sequence encoded by the current handle
                // has size >= k
                if (kmer.seq.len() as u64) == k {
                    //if !complete_kmers.contains(&kmer) {
                    complete_kmers.push(kmer);
                    //}
                } else {
                    // The kmer is incomplete, thus will have to be completed to reach size k

                    // Check that eventual limits have not been reach yet, otherwise ignore
                    // this kmer

                    let mut next_count: u64 = 0;
                    if edge_max.is_some() || degree_max.is_some() {
                        graph
                            .handle_edges_iter(handle, Direction::Right)
                            .for_each(|_| {
                                next_count += 1;
                            });
                    }

                    if (degree_max.is_none() && edge_max.is_none())
                        || (degree_max.is_some() && next_count < degree_max.unwrap())
                        || (edge_max.is_some() && kmer.forks < edge_max.unwrap())
                    {
                        // Create a copy of the incomplete kmer for each neighbour handle,
                        // so that they can be completed

                        for neighbor in graph.handle_edges_iter(handle, Direction::Right) {
                            let mut inc_kmer = kmer.clone();
                            inc_kmer.last_handle = neighbor;

                            if next_count > 1 {
                                inc_kmer.forks += 1;
                            }

                            incomplete_kmers.push(inc_kmer);
                        }
                    }
                }
            }

            // Then complete all incomplete kmers
            while let Some(mut incomplete_kmer) = incomplete_kmers.pop() {
                handle = incomplete_kmer.last_handle;
                handle_seq = graph.sequence(handle).into_string_lossy();
                handle_length = handle_seq.len() as u64;

                let end = min(k - (incomplete_kmer.seq.len() as u64), handle_length);
                let str_to_add = handle_seq.substring(0, end as usize).to_string();
                incomplete_kmer.extend_kmer(str_to_add, handle);

                // Ignore Ns during kmer generation
                if incomplete_kmer.seq.contains('N') {
                    continue;
                }

                if (incomplete_kmer.seq.len() as u64) == k {
                    //if !complete_kmers.contains(&incomplete_kmer) {
                    complete_kmers.push(incomplete_kmer);
                    //}
                } else {
                    // NOTE: if there is no neighbor, the kmer does not get re-added
                    // to the incomplete ones, so that the external loop can end
                    for neighbor in graph.handle_edges_iter(handle, Direction::Right) {
                        let mut next_count: u64 = 0;
                        if edge_max.is_some() || degree_max.is_some() {
                            graph
                                .handle_edges_iter(handle, Direction::Right)
                                .for_each(|_| {
                                    next_count += 1;
                                });
                        }

                        if (degree_max.is_none() && edge_max.is_none())
                            || (degree_max.is_some() && next_count < degree_max.unwrap())
                            || (edge_max.is_some() && incomplete_kmer.forks < edge_max.unwrap())
                        {
                            let mut inc_kmer = incomplete_kmer.clone();
                            inc_kmer.add_handle_to_complete(neighbor);

                            if next_count > 1 {
                                inc_kmer.forks += 1;
                            }

                            incomplete_kmers.push(inc_kmer);
                        }
                    }
                }
            }

            // IMPORTANT NOTE: a single iteration (of the most external for loop) completes
            // ALL possible kmers that start in the given handle. If the graph "ends" but the
            // kmer cannot reach size k (i.e. it is incomplete) it gets discarded.
            assert!(incomplete_kmers.is_empty());
        }
    }

    // Sort the kmers so that equal kmers (= having the same sequence) are close to each other
    // Note that the same kmer can appear in different places
    // (e.g. CACTTCAC -> CAC and CAC must be consecutive in the ordering)
    complete_kmers.sort_by(|x, y| x.seq.cmp(&y.seq));
    // Also dedup the vec as exact duplicates only waste space. Also note that dedup only works
    // on consecutive duplicates, so only by sorting beforehand it works correctly.
    complete_kmers.dedup();

    complete_kmers
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::read_generation::{
        reads_from_multiple_kmers, reads_to_fasta_file, GeneratedRead,
    };
    use handlegraph::handle::Edge;
    use handlegraph::mutablehandlegraph::MutableHandleGraph;
    use handlegraph::pathgraph::PathHandleGraph;

    /// This function creates a simple graph, used for debugging
    ///            | 2: CT \
    /// FWD  1: AAA         4: GCAC
    ///            \ 3: GA |
    fn create_simple_graph() -> HashGraph {
        let mut graph: HashGraph = HashGraph::new();

        let h1 = graph.create_handle("AAA".as_bytes(), 1);
        let h2 = graph.create_handle("CT".as_bytes(), 2);
        let h3 = graph.create_handle("GA".as_bytes(), 3);
        let h4 = graph.create_handle("GCAC".as_bytes(), 4);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));

        let p1 = graph.create_path_handle("P1".as_bytes(), false);
        graph.append_step(&p1, h1);
        graph.append_step(&p1, h2);
        graph.append_step(&p1, h4);

        let p2 = graph.create_path_handle("P2".as_bytes(), false);
        graph.append_step(&p2, h1);
        graph.append_step(&p2, h3);
        graph.append_step(&p2, h4);

        graph
    }

    #[test]
    fn test_generate_reads() {
        let graph = create_simple_graph();
        let fwd_rev_kmers = generate_kmers_from_graph(&graph, 3, Some(10), Some(10));
        let fwd_kmers: Vec<GraphKmer> = fwd_rev_kmers
            .into_iter()
            .filter(|k| k.handle_orient == true)
            .collect();
        //println!("Kmers: {:#?}", fwd_kmers);
        let reads: Vec<GeneratedRead> = reads_from_multiple_kmers(&fwd_kmers);
        println!("Reads: {:#?}", reads);
        reads_to_fasta_file(&reads, "reads.fasta").ok();
    }
}
