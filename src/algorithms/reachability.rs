use gfa::gfa::Orientation;
use handlegraph::handle::{Direction, Handle, NodeId};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use std::collections::{HashSet, VecDeque};

/// Return a list of paths between two given nodes, where each path is represented
/// as a list of NodeIds
pub fn find_all_paths_between(
    graph: &HashGraph,
    start_node_id: u64,
    end_node_id: u64,
    max_edges: u64,
) -> Vec<Vec<u64>> {
    let mut target_paths_list: Vec<Vec<u64>> = Vec::new();

    // Put a limit on the maximum amount of edges that can be traversed
    // this should prevent excessive memory usage
    let mut curr_edges = 0;
    let mut edges_limit_reached = false;

    // Keep a set of visited nodes so that loops are avoided
    let mut visited_node_id_set: HashSet<u64> = HashSet::new();

    // Create queue
    // NOTE: a queue is used to avoid stack overflows
    let mut queue: VecDeque<u64> = VecDeque::new();
    queue.push_back(start_node_id);
    target_paths_list.push(vec![start_node_id]);

    while !queue.is_empty() {
        let curr_node_id = queue.pop_front().unwrap();

        if curr_node_id == end_node_id {
            continue;
        }

        visited_node_id_set.insert(curr_node_id);
        let current_handle = Handle::new(curr_node_id, Orientation::Forward);

        // Get all paths that end in curr_node
        let mut curr_paths_list: Vec<_> = target_paths_list.clone();
        curr_paths_list.retain(|x| x.ends_with(&[curr_node_id]));

        // Only keep those that don't
        target_paths_list.retain(|x| !x.ends_with(&[curr_node_id]));

        for neighbor in graph.handle_edges_iter(current_handle, Direction::Right) {
            // Append, for each current_path, this neighbor
            let mut temp = curr_paths_list.clone();
            temp.iter_mut()
                .for_each(|x| x.push(neighbor.unpack_number()));
            target_paths_list.append(&mut temp);

            // Add new node to queue
            if !visited_node_id_set.contains(&neighbor.unpack_number())
                && !queue.contains(&neighbor.unpack_number())
            {
                queue.push_back(neighbor.unpack_number());
            }

            // Break if too many edges have been visited
            curr_edges += 1;
            if curr_edges > max_edges {
                edges_limit_reached = true;
                break;
            }
        }

        if edges_limit_reached {
            break;
        }
    }

    // Only keep paths that end in end_node_id
    // start_node_id does not have to be checked
    // TODO: maybe not needed?
    target_paths_list.retain(|x| x.ends_with(&[end_node_id]));

    target_paths_list
}

#[cfg(test)]
mod tests {
    use super::*;
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

        graph
    }

    /// This function creates a more complex graph
    fn create_complex_graph() -> HashGraph {
        let mut graph: HashGraph = HashGraph::new();

        let h1 = graph.create_handle("AAA".as_bytes(), 1);

        let h2 = graph.create_handle("CT".as_bytes(), 2);
        let h3 = graph.create_handle("GA".as_bytes(), 3);

        let h4 = graph.create_handle("GCAC".as_bytes(), 4);
        let h5 = graph.create_handle("GTGC".as_bytes(), 5);

        let h6 = graph.create_handle("GCAC".as_bytes(), 6);
        let h7 = graph.create_handle("GTGC".as_bytes(), 7);

        let h8 = graph.create_handle("TTTT".as_bytes(), 8);
        let h9 = graph.create_handle("TCCC".as_bytes(), 9);

        let h10 = graph.create_handle("AAA".as_bytes(), 10);

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));

        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h2, h5));

        graph.create_edge(&Edge(h3, h6));
        graph.create_edge(&Edge(h3, h7));

        graph.create_edge(&Edge(h4, h8));
        graph.create_edge(&Edge(h5, h8));

        graph.create_edge(&Edge(h6, h9));
        graph.create_edge(&Edge(h7, h9));

        graph.create_edge(&Edge(h8, h10));
        graph.create_edge(&Edge(h9, h10));

        // Add some back edges to increase complexity
        graph.create_edge(&Edge(h6, h2));
        graph.create_edge(&Edge(h10, h3));

        graph
    }

    #[test]
    fn test_simple_paths() {
        let graph = create_simple_graph();
        let paths_between_1_and_4 = find_all_paths_between(&graph, 1, 4, u64::MAX);
        println!("paths: {:?}", paths_between_1_and_4);
        assert_eq!(paths_between_1_and_4.len(), 2);
    }

    #[test]
    fn test_complex_paths() {
        let graph = create_complex_graph();
        let paths_between_1_and_10 = find_all_paths_between(&graph, 1, 10, u64::MAX);
        println!("paths: {:?}", paths_between_1_and_10);
        assert_eq!(paths_between_1_and_10.len(), 4);
    }
}
