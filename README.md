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
```