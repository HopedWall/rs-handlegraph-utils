use crate::utils::kmer_generation::GraphKmer;
use bstr::{ByteSlice, ByteVec};
use handlegraph::handle::{Direction, Handle};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use itertools::Itertools;
use rand::prelude::*;
use simple_sds::bit_vector::rank_support::RankSupport;
use simple_sds::bit_vector::*;
use simple_sds::raw_vector::{AccessRaw, PushRaw, RawVector};
use std::cmp::min;
use std::fs::File;
use std::io::Write;

#[derive(Debug, Clone)]
pub struct GeneratedRead {
    name: String,
    seq: String,
    path_id: Option<i64>,
    handles: Vec<Handle>,
    length: usize,
    first_node_offset: usize,
    last_node_offset: usize,
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
        }
    }

    pub fn new_from_kmer(kmer: &GraphKmer, name: Option<&str>) -> Self {
        let read_name = match name {
            Some(input_name) => String::from(input_name),
            _ => {
                let mut rng = rand::thread_rng();
                format!("read-{}", rng.gen::<u8>())
            }
        };

        GeneratedRead {
            name: read_name,
            seq: kmer.seq.clone(),
            path_id: None,
            length: kmer.seq.len(),
            handles: kmer.all_handles.clone(), //TODO: fix
            first_node_offset: kmer.begin_offset.position as usize,
            last_node_offset: kmer.end_offset.position as usize,
        }
    }

    fn extend_read(&mut self, new_seq: &str) {
        self.seq = format!("{}{}", self.seq, new_seq);
        self.length += new_seq.len();
    }

    fn to_fasta(&self) -> String {
        let header = format!(
            ">{} {{nodes:{:?},start_offset:{},end_offset:{}}}",
            self.name,
            self.handles
                .iter()
                .map(|x| x.unpack_number())
                .collect::<Vec<u64>>(),
            self.first_node_offset,
            self.last_node_offset
        );
        format!("{}\n{}", header, self.seq)
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
    let mut file =
        File::create(&file_name).unwrap_or_else(|_| panic!("Couldn't create file {}", &file_name));
    file.write_all(&read_strings.join("\n").as_bytes())
        .unwrap_or_else(|_| panic!("Couldn't write to file {}", &file_name));
    Ok(())
}

pub fn reads_from_multiple_kmers(kmers: &Vec<GraphKmer>) -> Vec<GeneratedRead> {
    kmers
        .iter()
        .enumerate()
        .map(|(i, kmer)| GeneratedRead::new_from_kmer(kmer, Some(&format!("Read-{}", i))))
        .collect()
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
