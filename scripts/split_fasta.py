import argparse
import sys
import threading
import math
import subprocess
import mmap
from pathlib import Path
import os
import gzip

def parse_nreads(x):
	# upper() just in case 'K' or 'M' is lowercase, doesn't do anything on actual numbers
	if x.upper().endswith('K') and len(x) > 1:
		return int(x.upper().replace('K', '')) * 1000
	elif x.upper().endswith('M') and len(x) > 1:
		return int(x.upper().replace('M', '')) * 1000000
	else:
		return int(x)

def write_to_file(contents, prefix, i, gzipped=False):
	if not gzipped:
		filename = "{}-{}.fa".format(prefix, i)
		with open(filename, "w") as output_fp:
				output_fp.write(''.join(curr_lines))
	else:
		filename = "{}-{}.fa.gz".format(prefix, i)
		with gzip.open(filename, "wb") as output_fp:
				output_fp.write(''.join(curr_lines).encode())


parser = argparse.ArgumentParser(description='Splits a multifasta into multiple multifasta')
parser.add_argument('-i', '--input', help='Path to the input multifasta file', required=True)
parser.add_argument('-o', '--out-prefix', help="Prefix to be used to store output multifasta")
parser.add_argument('-n', '--n-reads', help="Number of reads that each output file should contain", required=True)
parser.add_argument('-g', '--gzip-output', action="store_true", help="Gzip output multifasta")	# This is a yes/no flag
#parser.add_argument('-t', '--threads', type=int, default=1, help="Number of threads to be used")
args = vars(parser.parse_args())

input_file = Path(args["input"])
if not input_file.exists():
	raise Exception("Input file does not exist")

ext = "".join(input_file.suffixes)
if args["out_prefix"]:
	out_prefix = args["out_prefix"]
else:
	# Remove all suffixes from input_file
	out_prefix = input_file
	while out_prefix.suffix:
		out_prefix = out_prefix.with_suffix('')

#n_input_reads = sum(1 for i in open(input_file, 'rb')) // 2 #/2 since header+seq, // = floor division
n_input_reads = 0
if ext == ".fa" or ext == "fasta":
	n_input_reads = int(subprocess.run(["wc", "-l", input_file], check=True, capture_output=True).stdout.split()[0]) // 2
elif ext == ".fa.gz" or ext == ".fasta.gz":
	zcat = subprocess.run(['zcat', input_file], check=True, capture_output=True)
	n_input_reads = int(subprocess.run(["wc", "-l"], input=zcat.stdout, check=True, capture_output=True).stdout.split()[0]) // 2
else:
	raise Exception("Unknown format {} for input file! (Note: FASTQ not currently supported)".format(ext))

n_reads_per_file = parse_nreads(args["n_reads"])
n_output_files = math.ceil(n_input_reads / n_reads_per_file)

print("Found {} reads, splitting into {} files, each with {} reads (or less)".format(n_input_reads, n_output_files, n_reads_per_file))

#intervals = [range(i*n_reads_per_file*2, min(i*n_reads_per_file*2 + n_reads_per_file*2 - 1, n_input_reads*2) + 1) for i in range(n_output_files)]

#in_prefix, ext = os.path.splitext(input_file)
#out_prefix = args["out_prefix"] if args["out_prefix"] else in_prefix

i = 0
curr_lines = []

if ext == ".fa" or ext == ".fasta":
	with open(input_file, "r") as input_fp:
		for name, seq in zip(input_fp, input_fp):
			if len(curr_lines)//2 == n_reads_per_file:
				write_to_file(curr_lines, out_prefix, i, args["gzip_output"])
				i+=1
				curr_lines = []
			curr_lines.append(name)
			curr_lines.append(seq)

elif ext == ".fa.gz" or ext == ".fasta.gz":
	with gzip.open(input_file, "rb") as input_fp:
		for name, seq in zip(input_fp, input_fp):
			if len(curr_lines)//2 == n_reads_per_file:
				write_to_file(curr_lines, out_prefix, i, args["gzip_output"])
				i+=1
				curr_lines = []
			curr_lines.append(name.decode())
			curr_lines.append(seq.decode())	
else:
	raise Exception("Unknown format {} for input file! (Note: FASTQ not currently supported)".format(ext))

if curr_lines:
	write_to_file(curr_lines, out_prefix, i, args["gzip_output"])