use crate::utils::kmer_generation::GraphKmer;
use bstr::{ByteSlice, ByteVec};
use handlegraph::handle::{Direction, Handle};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use itertools::Itertools;
use rand::prelude::*;
use rand::{thread_rng, Rng};
//use simple_sds::bit_vector::rank_support::RankSupport;
//use simple_sds::bit_vector::*;
//use simple_sds::raw_vector::{AccessRaw, PushRaw, RawVector};
use std::cmp::min;
use std::fs::{read, File};
use std::io::Write;
use crate::utils::kmer_generation::SeqPos;
use substring::Substring;

#[derive(Debug, Clone)]
pub struct GeneratedRead {
    name: String,
    seq: String,
    path_id: Option<i64>,
    handles: Vec<Handle>,
    length: usize,
    first_node_offset: usize,
    last_node_offset: usize,
    errors: Vec<usize>,
}

impl GeneratedRead {
    pub fn new() -> Self {
        GeneratedRead {
            name: String::new(),
            seq: String::new(),
            path_id: Some(0),
            length: 0,
            handles: vec![],
            first_node_offset: 0,
            last_node_offset: 0,
            errors: vec![],
        }
    }

    pub fn new_from_kmer(kmer: &GraphKmer, name: Option<&str>, err_rate: Option<f64>) -> Self {
        let mut rng = thread_rng();
        let alphabet = ['A', 'C', 'G', 'T'];

        let read_name = match name {
            Some(input_name) => String::from(input_name),
            _ => {
                let mut rng = rand::thread_rng();
                format!("read-{}", rng.gen::<u8>())
            }
        };

        let mut read_seq = kmer.seq.clone();
        let mut err_positions: Vec<usize> = Vec::new();

        if let Some(error_rate) = err_rate {
            let seq_with_errors = kmer
                .seq
                .chars()
                .enumerate()
                .map(|(i, og_base)| match rng.gen_bool(error_rate) {
                    true => {
                        let index = rng.gen_range(0..4);
                        let new_base: char = alphabet.get(index).unwrap().clone();
                        println!(
                            "[{}]: Replacing {} with {} (pos: {})",
                            read_name, og_base, new_base, i
                        );
                        err_positions.push(i);
                        new_base
                    }
                    false => og_base,
                })
                .collect();

            err_positions.sort();
            read_seq = seq_with_errors;

            //println!("Kmer: {}", kmer.seq);
            //println!("Read: {}", read_seq);
            //println!("Errors: {:?}", err_positions);
        }

        GeneratedRead {
            name: read_name,
            seq: read_seq,
            path_id: None,
            length: kmer.seq.len(),
            handles: kmer.all_handles.clone(),
            first_node_offset: kmer.begin_offset.position as usize,
            last_node_offset: kmer.end_offset.position as usize,
            errors: err_positions,
        }
    }

    fn extend_read(&mut self, new_seq: &str) {
        self.seq = format!("{}{}", self.seq, new_seq);
        self.length += new_seq.len();
    }

    pub fn to_fasta(&self) -> String {
        let header = format!(
            ">{} {{\"nodes\":{:?},\"start_offset\":{},\"end_offset\":{},\"errors\":{:?}}}",
            self.name,
            self.handles
                .iter()
                .map(|x| x.unpack_number())
                .collect::<Vec<u64>>(),
            self.first_node_offset,
            self.last_node_offset,
            self.errors
        );
        format!("{}\n{}", header, self.seq)
    }

    pub fn fasta_from_kmer(kmer: &GraphKmer, name: Option<&str>, err_rate: Option<f64>) -> String {
        GeneratedRead::new_from_kmer(kmer, name, err_rate).to_fasta()
    }
}

pub fn exact_reads_from_path(graph: &HashGraph, path_id: i64, size: usize) -> Vec<GeneratedRead> {
    let mut reads: Vec<GeneratedRead> = Vec::new();

    let path = graph.paths.get(&path_id).unwrap();
    let path_nodes = path.nodes.clone();

    let mut incomplete_read: Option<GeneratedRead> = None;
    let mut read_count = 0;
    for curr_handle in path_nodes {
        let curr_handle_seq = String::from_utf8(graph.sequence(curr_handle)).unwrap();

        let mut start = 0;
        let mut end = 0;

        // First try to complete incomplete reads
        if let Some(mut read) = incomplete_read.clone() {
            end = min(size - read.length, curr_handle_seq.len());

            read.extend_read(&curr_handle_seq[start..end]);
            read.handles.push(curr_handle.clone());
            read.last_node_offset = end - 1;

            if read.length == size {
                reads.push(read.clone());
                incomplete_read = None;
            }

            start = end;
        }

        // Then try creating reads from the current node
        while start < curr_handle_seq.len() {
            let end = min(start + size, curr_handle_seq.len());

            //println!("Read seq is: {:#?}", read_seq);
            let read = GeneratedRead {
                name: format!("read-{}", read_count),
                seq: curr_handle_seq[start..end].to_string(),
                path_id: Some(path_id),
                handles: vec![curr_handle],
                length: end - start,
                first_node_offset: start,
                last_node_offset: end - 1,
                errors: vec![],
            };

            if read.length == size {
                reads.push(read);
            } else {
                incomplete_read = Some(read);
            }

            read_count += 1;
            start += end;
        }
    }

    reads
}

pub fn split_sequence_into_kmers(sequence: &str, k: usize) -> Vec<String> {
    (0..sequence.len() - k + 1)
        .map(|i| String::from(&sequence[i..i + k]))
        .collect()
}

pub fn split_sequence_into_disjoint_substrings(sequence: &str, k: usize) -> Vec<String> {
    (0..sequence.len() - k + 1)
        .step_by(k)
        .map(|i| String::from(&sequence[i..i + k]))
        .collect()
}

pub fn reads_to_fasta_file(reads: &Vec<GeneratedRead>, file_name: &str) -> std::io::Result<()> {
    let read_strings: Vec<String> = reads.iter().map(|read| read.to_fasta()).collect();
    println!("Creating fasta file...");
    let mut file = File::create(&file_name)?;
    println!("Writing to fasta file...");
    file.write_all(&read_strings.join("\n").as_bytes())?;
    Ok(())
}

pub fn reads_from_multiple_kmers(
    kmers: &Vec<GraphKmer>,
    errors: bool,
    err_rate: Option<f64>,
) -> Vec<GeneratedRead> {
    kmers
        .iter()
        .enumerate()
        .map(|(i, kmer)| GeneratedRead::new_from_kmer(kmer, Some(&format!("Read-{}", i)), err_rate))
        .collect()
}

fn write_kmer_to_file(kmer: GraphKmer, read_count: &mut u64, err_rate: &Option<f64>, file: &mut Option<File>) -> std::io::Result<()> {
    let read_name = format!("read-{}", read_count);
    let fasta = format!("{}\n",GeneratedRead::fasta_from_kmer(&kmer, Some(&read_name), *err_rate));

    if let Some(file_handle) = file {
        file_handle.write(fasta.as_bytes());
    } else {
        print!("{}", fasta);
    }

    *read_count += 1;

    Ok(())
}

pub fn extract_read_direct(graph: &HashGraph, k: u64, n_reads: Option<u64>, only_forward: bool, err_rate: Option<f64>, fasta_file: Option<&str>) -> std::io::Result<()> {
    println!("Creating fasta file...");
    let mut file: Option<File> = match fasta_file {
        Some(filename) => File::create(filename).ok(),
        None => None
    };

    let mut read_count: u64 = 0;
    println!("Generating reads...");
    
    let sorted_graph_handles: Vec<Handle> = graph.handles_iter().sorted().collect();

    'outer: for forward_handle in sorted_graph_handles {
        let mut handle: Handle;
        let mut orient: bool;

        let possible_orients = match only_forward {
            true => vec![true],
            false => vec![true, false],
        };

        for handle_orient in possible_orients {
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

            let starting_handle = handle.clone();

            // Get current node/handle
            let mut handle_seq = graph.sequence(handle).into_string_lossy();
            let mut handle_length = handle_seq.len() as u64;

            // This will store kmers that have not yet reached size k, which
            // will have to be completed from the neighbours of the handle
            let mut incomplete_kmers: Vec<GraphKmer> = Vec::new();

            /*
            println!("Node length: {}", handle_seq.len());
            println!("A - Complete kmers: {}", complete_kmers.len());
            println!("A - Incomplete kmers: {}", incomplete_kmers.len());
            */

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
                    write_kmer_to_file(kmer, &mut read_count, &err_rate, &mut file).ok();

                    if let Some(max_n_reads) = n_reads {
                        if read_count >= max_n_reads {
                            break 'outer;
                        }
                    }
                } else {
                    // The kmer is incomplete, thus will have to be completed to reach size k

                    // Create a copy of the incomplete kmer for each neighbour handle,
                    // so that they can be completed

                    //let close : Vec<Handle> = graph.handle_edges_iter(handle, Direction::Right).collect();
                    //println!("A - Creating incomplete kmers for {} handles", close.len());
                    for neighbor in graph.handle_edges_iter(handle, Direction::Right) {

                        if orient == true && starting_handle.unpack_number() > neighbor.unpack_number()
                        || orient == false && starting_handle.unpack_number() < neighbor.unpack_number() {
                            //println!("Dropping kmer: {:?}", kmer);
                            continue
                        }

                        let mut inc_kmer = kmer.clone();
                        inc_kmer.last_handle = neighbor;

                        incomplete_kmers.push(inc_kmer);
                    }
                
                }
            }
            
            /*
            println!("B - Complete kmers: {}", complete_kmers.len());
            println!("B - Incomplete kmers: {}", incomplete_kmers.len());
            */

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
                    
                    write_kmer_to_file(incomplete_kmer, &mut read_count, &err_rate, &mut file).ok();

                    if let Some(max_n_reads) = n_reads {
                        if read_count >= max_n_reads {
                            break 'outer;
                        }
                    }

                    if let Some(max_n_reads) = n_reads {
                        if read_count >= max_n_reads {
                            break 'outer;
                        }
                    }

                } else {
                    // NOTE: if there is no neighbor, the kmer does not get re-added
                    // to the incomplete ones, so that the external loop can end
                    //let close : Vec<Handle> = graph.handle_edges_iter(handle, Direction::Right).collect();
                    //println!("B - Creating incomplete kmers for {} handles", close.len());
                    for neighbor in graph.handle_edges_iter(handle, Direction::Right) {

                        if orient == true && starting_handle.unpack_number() > neighbor.unpack_number()
                            || orient == false && starting_handle.unpack_number() < neighbor.unpack_number() {
                                //println!("Dropping kmer: {:?}", incomplete_kmer);
                                continue
                        }

                        let mut inc_kmer = incomplete_kmer.clone();
                        inc_kmer.add_handle_to_complete(neighbor);

                        incomplete_kmers.push(inc_kmer);
                    }
                }
            }

            /*
            println!("C - Complete kmers: {}", complete_kmers.len());
            println!("C - Incomplete kmers: {}\n", incomplete_kmers.len());
            */

            // IMPORTANT NOTE: a single iteration (of the most external for loop) completes
            // ALL possible kmers that start in the given handle. If the graph "ends" but the
            // kmer cannot reach size k (i.e. it is incomplete) it gets discarded.
            assert!(incomplete_kmers.is_empty());
        }
    }



    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use handlegraph::handle::Edge;
    use handlegraph::mutablehandlegraph::MutableHandleGraph;
    use handlegraph::pathgraph::PathHandleGraph;

    #[test]
    fn test_kmers_split() {
        let sequence = String::from("AAACGT");
        let kmers = split_sequence_into_kmers(&sequence, 3);
        //println!("Kmers: {:#?}", kmers);
        assert_eq!(kmers.len(), 4);
        assert_eq!(kmers, vec!["AAA", "AAC", "ACG", "CGT"]);
    }

    #[test]
    fn test_kmers_disjoint() {
        let sequence = String::from("AAACGT");
        let kmers = split_sequence_into_disjoint_substrings(&sequence, 3);
        //println!("Kmers: {:#?}", kmers);
        assert_eq!(kmers.len(), 2);
        assert_eq!(kmers, vec!["AAA", "CGT"]);
    }

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
        let reads = exact_reads_from_path(&graph, 0, 3);
        println!("Reads: {:#?}", reads);
        assert_eq!(reads.len(), 3);
        assert_eq!(reads.get(0).unwrap().seq, "AAA");
        assert_eq!(reads.get(1).unwrap().seq, "CTG");
        assert_eq!(reads.get(2).unwrap().seq, "CAC");

        println!(
            "fasta reads: {:#?}",
            reads.iter().map(|r| r.to_fasta()).collect::<Vec<String>>()
        );
        reads_to_fasta_file(&reads, "reads.fasta").ok();
    }
}
