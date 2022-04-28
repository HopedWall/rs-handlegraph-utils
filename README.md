# rs-handlegraph-utils
A collection of algorithms and other utilities useful when working with Variation Graphs.

## Algorithms
- graph visits
- finding all the kmers in the graph
- finding all paths between two given nodes
- more to come...

## Utils

```
## Clone this repo
git clone https://github.com/HopedWall/rs-handlegraph-utils
cd rs-handlegraph-utils

## Install it
cargo install --path .

## Generate 1000 exact reads of length 100 from the graph, and put them in reads.fa
handlegraph-utils -i graph.gfa -n 1000 --len 100 --out reads.fa

## Generate 1000 reads of length 100 from the graph, with an error rate of 5%, and put them in reads-error.fa
handlegraph-utils -i graph.gfa -n 1000 --len 100 --out reads-errors.fa --errors -r 0.05

## For each read, the header contains:
## - nodes used to generate the read
## - start offset on the first node
## - end offset on the last node
## - positions on the read that contain errors (i.e. mutations)
cat reads-errors.fa
>Read-0 {"nodes":[13],"start_offset":6,"end_offset":106,"errors":[1, 82]}
...
```