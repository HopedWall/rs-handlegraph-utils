use handlegraph::handle::{Direction, Handle};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::*;
use itertools::Itertools;

// Run a BFS visit of the graph
// TODO: add parameter for reverse handles (max_handle -> max_handle - 1 -> ...)
pub fn bfs(graph: &HashGraph) -> Vec<Handle> {
    let mut bfs: Vec<Handle> = Vec::new();

    // First get all the handles in the order of NodeIDs
    // TODO: maybe only get the first handle
    let handles: Vec<Handle> = graph.handles_iter().sorted().collect();

    if let Some(first_handle) = handles.first() {
        let mut bfs_front = vec![first_handle.clone()];
        bfs.extend(bfs_front.clone().into_iter());

        while !bfs_front.is_empty() {
            bfs_front = advance_bfs(graph, bfs_front);
            bfs.extend(bfs_front.clone().into_iter());
        }
    }

    bfs
}

// Only keep forward edges
fn find_neighbors_to_add(graph: &HashGraph, curr_handle: Handle) -> Vec<Handle> {
    graph
        .handle_edges_iter(curr_handle, Direction::Right)
        .filter_map(|x| match x.as_integer() > curr_handle.as_integer() {
            true => Some(x),
            _ => None,
        })
        .collect()
}

// Advance the current front of the BFS
fn advance_bfs(graph: &HashGraph, curr_front: Vec<Handle>) -> Vec<Handle> {
    curr_front
        .iter()
        .flat_map(|x| find_neighbors_to_add(graph, *x))
        .sorted()
        .dedup()
        .collect()
}

/*
// Advance the current front of the BFS in parallel
fn advance_bfs_parallel(graph: &HashGraph, curr_front: Vec<Handle>) -> Vec<Handle> {
    let mut new_front: Vec<Handle> =
        curr_front
            .iter()
            .flat_map(|x| find_neighbors_to_add(graph, *x))
            .collect();

    new_front.par_sort_unstable();
    new_front.dedup();

    new_front
}
 */
