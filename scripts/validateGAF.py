'''
Script for the validation of GAFs produced from the reads (generated by this crate, or by other tools).
USAGE:
First, generate some reads, i.e.:
    handlegraph-utils -i graph.gfa -n 1000 --len 100 --out reads.fa
Then, align them with your favorite aligner and obtain a GAF. For example, we can use GraphAligner:
    GraphAligner -g graph.gfa -f reads.fa -a alignments.gaf -x vg
Finally, call this script on the resulting GAF. One of the required parameters is the treshold to be used in order
for the alignment to be considered "correct".
    python3 validateGAF.py alignments.gaf 0.5 graph.og reads.fa yes
Q&A:
Q - Can I use this script with reads generated by other tools other than handlegraph-utils?
A - YES, if the reads were generated with another tool (i.e. pbsim2) then you will need to pass "no"
    as the final argument of the validation script. However, when using another tool, the resulting
    validation will only include the edit similarity and not the node comparison.
'''

import argparse
import pandas as pd
import re
import json
import sys
# Normalized because we want a similarity between 0 and 1.
# Non-normalized would be the actual number of single-character edits.
from strsimpy.normalized_levenshtein import NormalizedLevenshtein
normalized_levenshtein = NormalizedLevenshtein()

# Parse GAF to a dataframe
def parse_GAF(path_to_file):
    gaf_fields = ["name", "qlen", "qstart", "qend",
                  "strand", "path", "plen", "pstart", "pend",
                  "residue", "alblock", "quality", "extra",
                  "extra1", "extra2", "extra3", "extra4"]
    my_gaf = pd.read_csv(path_to_file, sep='\t', names=gaf_fields)
    return my_gaf

# Parse FASTA
class FASTARecord:

    def __init__(self, name, seq, has_truth=False):
        self.name = name
        self.seq = seq

        if has_truth:
            name_split = self.name.split(" {")
            assert(len(name_split) == 2)
            self.actual_name = name_split[0]
            self.mappings = json.loads(' {'+name_split[1])

def parse_FASTA(path_to_file, has_truth):
    records = []
    with open(path_to_file) as f:
        for ext_name, seq in zip(f, f):
            name = ext_name.strip()[1:]
            #name = ext_name.split(" {")[0][1:]
            records.append(FASTARecord(name, seq, has_truth))
    return records

# Parse GFA
def parse_GFA_nodes(path_to_file):
    nodes = dict()
    with open(path_to_file) as f:
        for line in f:
            if line.lstrip().startswith("S"):
                _, nodeid, seq = line.strip().split("\t")
                nodes[int(nodeid)] = seq
    return nodes

# Convert string to boolean
def str2bool(s):
    return s.lower() in ("yes", "true", "t", "1", "y")

# CLI arguments parsing
parser = argparse.ArgumentParser(description='Returns the accuracy of a gaf GAF')
parser.add_argument('GAF', help='Path to the GAF file')
parser.add_argument("Threshold", help="Threshold that should be used to consider a read mapped correctly")
#parser.add_argument('Graph', help='Path to graph in .odgi/.og format')
parser.add_argument('Graph', help='Path to graph in .gfa format')
parser.add_argument('Reads', help='Path to reads in .fasta/.fa format')
parser.add_argument("Has-truth", help="Were the reads generated by handlegraph-utils? (true/false)")
args = vars(parser.parse_args())

print("Command used: {}".format(" ".join(sys.argv)))
#print("GAF file: {}".format(args["GAF"]))

# Read GAF file
my_gaf = parse_GAF(args["GAF"])

# Print out some debug information
threshold = float(args["Threshold"])
print("Threshold is: {}".format(threshold))

# Is the truth there?
has_truth = str2bool(args["Has-truth"])
print("Has truth?", "Yes" if has_truth else "No")

# Read FA file
my_reads = parse_FASTA(args["Reads"], has_truth)

# Load ODGI graph
#g = odgi.graph()
#g.load(args["Graph"])
GFA_nodes = parse_GFA_nodes(args["Graph"])

count_correct = 0
count_correct_nodes = 0

similarity_distribution = []

# For each GAF record
for i in range(len(my_gaf)):

    ##### Obtain truth values (node ids, offsets, errors)
    my_gaf_row = my_gaf.iloc[i]
    my_gaf_name = my_gaf_row["name"]
    name = my_gaf_name.strip()

    '''
    # STEP 0 (Optional): If truth is present (= reads were generated by handlegraph-utils) then obtain it
    if has_truth:
            #TODO: improve this, since jsonarrays can contain " ", split is not what I should do (maybe try REs?)
            name = name.split(" {")[0]
            metadata = ' {'+my_gaf_row["name"].split(" {")[1]
            read_gen_data = json.loads(metadata)

            truth_nodes = [int(node) for node in read_gen_data["nodes"]]
            start_offset = read_gen_data["start_offset"]
            end_offset = read_gen_data["end_offset"]
            if "errors" in read_gen_data:
                    errors = read_gen_data["errors"]

            ##### Get sequence from truth
            truth_as_string = ""
            truth_as_list = []
            for pos,nodeid in enumerate(truth_nodes):
                    handle = g.get_handle(nodeid)
                    handle_seq = g.get_sequence(handle)
                    
                    # Find which substring of the current handle should be considered
                    start_pos = 0
                    end_pos = len(handle_seq)
                    if pos == 0:
                            start_pos = start_offset
                    if pos == len(truth_nodes) - 1:
                            end_pos = end_offset

                    handle_seq = handle_seq[start_pos:end_pos]
                    truth_as_list.append(handle_seq)


            truth_as_string = "".join(truth_as_list)
    '''

    ####### STEP 1: Get sequence from read #######
    ### Here a sequence is obtained by extracting the substring of the query
    ### (aka the READ) that was actually aligned according to the GAF. At first,
    ### it extracts the whole read from the fasta, and then it finds the
    ### correct substring.

    #print("FASTA name is", my_reads[0].name)
    #print("GAF name is", name)
    reads_same_name = list(filter(lambda read: read.name == name, my_reads))

    # Try checking if the aligner has cut part of the FASTA name
    if len(reads_same_name) == 0 and has_truth:
        # Convert the list of FASTA records s.t. the name does not contain the JSON with the truth
        # i.e. "Read-0 {...}" -> "Read-0"
        reads_same_name = list(filter(lambda read: read.actual_name == name, my_reads))

    if len(reads_same_name) == 0:
        print("Couldn't find read {} in FASTA file, skipped".format(name))
        continue
    elif len(reads_same_name) > 1:
        print("Found multiple matching reads for {} in FASTA file, skipped".format(name))
        continue

    # Obtain whole read seq
    curr_read = reads_same_name[0]
    curr_read_seq = curr_read.seq

    if my_gaf_row["qstart"] == '*' or my_gaf_row["qend"] == '*':
        print("Could not find qstart/qend for {}, skipped".format(name))
        continue

    qstart = int(my_gaf_row["qstart"])
    qend = int(my_gaf_row["qend"])
    read_as_string = curr_read_seq[qstart:qend]

    ####### STEP 2: Get sequence from aln #######
    ### Here a sequence is obtained by extracting the node sequence
    ### of each node in the "path matching" of the GAF. This
    ### extracts the sequence from the reference aka the GRAPH.

    my_gaf_nodes_str = my_gaf_row["path"]
    aln_start = int(my_gaf_row["pstart"])
    aln_end = int(my_gaf_row["pend"])

    if my_gaf_nodes_str == '*':
        continue

    my_gaf_tuples = re.findall("(>|<)([0-9]+)", my_gaf_nodes_str)
    my_gaf_int = list(map(lambda x: int(x[1]), my_gaf_tuples))

    aln_as_string = ""
    aln_as_list = []
    for pos,nodeid in enumerate(my_gaf_int):

        #handle = g.get_handle(nodeid)
        #handle_seq = g.get_sequence(handle)

        if nodeid not in GFA_nodes.keys():
            print("Could not find node sequence for id {} in GAF file, skipped".format(nodeid))
            continue

        handle_seq = GFA_nodes[nodeid]

        # Find which substring of the current handle should be considered
        start_pos = 0
        end_pos = len(handle_seq)

        # Note: double if in case the alignment is only 1-node long
        if pos == 0:
            start_pos = aln_start
        if pos == len(my_gaf_int) - 1:
            end_pos = aln_end

        handle_seq = handle_seq[start_pos:end_pos+1]
        aln_as_list.append(handle_seq)

    aln_as_string = "".join(aln_as_list)

    ####### STEP 3: Compute edit distance between read_seq and aln_seq #######
    ### Compute the edit distance between the sequence found in step 1 (READ) and
    ### the one found in step 2 (GRAPH). If this distance is over a certain threshold,
    ### then it is considered to be aligned correctly.

    edit_similarity = normalized_levenshtein.similarity(read_as_string, aln_as_string)

    similarity_distribution.append({"name": curr_read.name if not has_truth else curr_read.actual_name, "similarity": edit_similarity})
    #print("{} - Edit similarity: {}".format(name, edit_similarity))
    if edit_similarity >= threshold:
        count_correct += 1
        #print("{} - OK, is: {}".format(name, edit_similarity))
        #print("ALN seqs are: {}".format(aln_as_list))
        #print("READ seqs are: {}".format(read_as_string))
    else:
        print("{} - Under threshold, is: {} (threshold: {})".format(name, edit_similarity, threshold))
        #print("ALN seqs are: {}".format(aln_as_list))
        print("ALN seq is: {}".format(aln_as_string))
        print("READ seq is: {}".format(read_as_string))
        #print("TRUTH nodes are: {}, ALN nodes are: {}".format(truth_nodes, my_gaf_int))
        #print("TRUTH seqs are: {}".format(truth_as_list))
        #print("ALN seqs are: {}".format(aln_as_list))

    if has_truth:
        truth_nodes = [int(node) for node in curr_read.mappings["nodes"]]

        if my_gaf_int == truth_nodes:
            count_correct_nodes += 1
        else:
            print("{} - Node mismatch, aln is: {} (truth: {})".format(name, my_gaf_int, truth_nodes))

print("There were {}/{} correct reads (edit similarity)".format(count_correct, len(my_reads)))

if has_truth:
    print("There were {}/{} correct reads (nodes matching)".format(count_correct_nodes, len(my_reads)))

print("Similarity distribution: {}".format(similarity_distribution))