use clap::{App, Arg};
use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use handlegraph::hashgraph::HashGraph;
use handlegraph_utils::utils::kmer_generation::{generate_kmers_from_graph, GraphKmer};
use handlegraph_utils::utils::read_generation::{
    reads_from_multiple_kmers, reads_to_fasta_file, GeneratedRead, extract_read_direct
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
                .required(false)
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
                .required(false)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("errors")
                .short("e")
                .long("errors")
                .help("Add errors to the generated reads")
                .required(false)
                .takes_value(false),
        )
        .arg(
            Arg::with_name("error-rate")
                .short("r")
                .long("error-rate")
                .value_name("FLOAT")
                .help("Sets the ratio of errors of the generated reads")
                .required(false)
                .takes_value(true),
        )
        .get_matches();

    let in_gfa_path: &str = matches
        .value_of("input")
        .expect("Could not parse argument --input");

    let out_fasta_path: Option<&str> = matches.value_of("output");

    let read_length: u64 = matches
        .value_of("length")
        .expect("Could not parse argument --length")
        .parse::<u64>()
        .unwrap();

    let n_reads: Option<u64> = match matches.is_present("n-reads") {
        true => matches
            .value_of("n-reads")
            .expect("Could not parse argument --n-reads")
            .parse::<u64>()
            .ok(),
        false => None,
    };

    let errors = matches.is_present("errors");
    let error_rate: Option<f64> = match errors {
        true => matches
            .value_of("error-rate")
            .expect("Could not parse argument --error-rate")
            .parse::<f64>()
            .ok(),
        false => None,
    };

    match generate_reads_from_graph(
        &in_gfa_path,
        read_length,
        n_reads,
        out_fasta_path,
        errors,
        error_rate,
    ) {
        Err(e) => panic!("{}", e),
        _ => println!(
            "Reads stored correctly in {}!",
            out_fasta_path.unwrap_or("standard output")
        ),
    }
}

pub fn generate_reads_from_graph(
    graph_file: &str,
    read_length: u64,
    n_reads: Option<u64>,
    fasta_file: Option<&str>,
    errors: bool,
    err_rate: Option<f64>,
) -> std::io::Result<()> {
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(&PathBuf::from(graph_file)).unwrap();
    let graph = HashGraph::from_gfa(&gfa);

    /*
    println!("Generating kmers...");
    let fwd_kmers =
        generate_kmers_from_graph(&graph, read_length, Some(100), Some(100), n_reads, true);

    /*
    let n_fwd_kmers: Vec<GraphKmer> = match n_reads {
        Some(amount) => fwd_kmers[0..amount as usize].to_vec(),
        None => fwd_kmers
    };
    let reads: Vec<GeneratedRead> = reads_from_multiple_kmers(&n_fwd_kmers, errors, err_rate);
     */

    println!("Generating reads...");
    let reads: Vec<GeneratedRead> = reads_from_multiple_kmers(&fwd_kmers, errors, err_rate);
    match fasta_file {
        Some(output_file) => reads_to_fasta_file(&reads, output_file),
        _ => Ok(reads.iter().for_each(|r| println!("{}", r.to_fasta()))),
    }
    */

    extract_read_direct(&graph, read_length, n_reads, true, err_rate, fasta_file)

}