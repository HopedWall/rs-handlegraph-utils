use clap::{App, Arg};
use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use handlegraph::hashgraph::HashGraph;
use handlegraph_utils::utils::kmer_generation::{generate_kmers_from_graph, GraphKmer};
use handlegraph_utils::utils::read_generation::{
    reads_from_multiple_kmers, reads_to_fasta_file, GeneratedRead,
};
use std::path::PathBuf;

fn main() {
    let matches = App::new("handlegraph-utils")
        .version("1.0")
        .author("Francesco Porto <francesco.porto97@gmail.com>")
        .about("Generates reads in FASTA format from a GFA")
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("in")
                .value_name("FILE.GFA")
                .help("Sets the graph to use")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("out")
                .value_name("FILE.FA")
                .help("Sets the file that will contain the reads")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("length")
                .short("l")
                .long("len")
                .value_name("INT")
                .help("Sets the length of the reads")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("n-reads")
                .short("n")
                .long("n-reads")
                .value_name("INT")
                .help("Sets the number of reads to be generated")
                .required(true)
                .takes_value(true),
        )
        .get_matches();

    let in_gfa_path: &str = matches
        .value_of("input")
        .expect("Could not parse argument --input");

    let out_fasta_path: &str = matches
        .value_of("output")
        .expect("Could not parse argument --output");

    let read_length: u64 = matches
        .value_of("length")
        .expect("Could not parse argument --length")
        .parse::<u64>()
        .unwrap();

    let n_reads: u64 = matches
        .value_of("n-reads")
        .expect("Could not parse argument --n-reads")
        .parse::<u64>()
        .unwrap();

    match generate_reads_from_graph(&in_gfa_path, read_length, n_reads, out_fasta_path) {
        Err(e) => panic!("{}", e),
        _ => println!("Reads stored correctly in {}!", out_fasta_path),
    }
}

pub fn generate_reads_from_graph(
    graph_file: &str,
    read_length: u64,
    n_reads: u64,
    fasta_file: &str,
) -> std::io::Result<()> {
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(&PathBuf::from(graph_file)).unwrap();
    let graph = HashGraph::from_gfa(&gfa);

    let fwd_rev_kmers = generate_kmers_from_graph(&graph, read_length, Some(100), Some(100));
    let fwd_kmers: Vec<GraphKmer> = fwd_rev_kmers
        .into_iter()
        .filter(|k| k.handle_orient == true)
        .collect();
    let n_fwd_kmers: Vec<GraphKmer> = fwd_kmers[0..n_reads as usize].to_vec();

    let reads: Vec<GeneratedRead> = reads_from_multiple_kmers(&n_fwd_kmers);
    reads_to_fasta_file(&reads, fasta_file)
}
