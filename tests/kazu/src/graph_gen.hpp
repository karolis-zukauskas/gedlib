#ifndef SRC_KAZU_GRAPH_GEN_HPP_
#define SRC_KAZU_GRAPH_GEN_HPP_

#include "util.hpp"
#include <numeric>

namespace ged {

GEDGraph::GraphID graph_gen_random(GxlGEDEnv& env, std::mt19937& rng, size_t num_nodes, double edge_density,
  std::vector<std::string> node_labels, std::vector<std::string> edge_labels) {

  constexpr char const* NODE_LABEL_KEY = "label";
  constexpr char const* EDGE_LABEL_KEY = "label";

  size_t const num_node_labels = node_labels.size();
  size_t const num_edge_labels = edge_labels.size();

  GEDGraph::GraphID const graph_id = env.add_graph();
  GXLLabel label;

  // Add nodes
  for (size_t i = 0; i < num_nodes; ++i) {
    label.clear();
    label[NODE_LABEL_KEY] = node_labels[rng() % num_node_labels];

    env.add_node(graph_id, std::to_string(i), label);
  }

  // Add edges
  std::vector<bool> edge_matrix(num_nodes * num_nodes, false);
  size_t const num_edges = static_cast<size_t>(edge_density * static_cast<double>(num_nodes * num_nodes - num_nodes) / 2.0);

  for (size_t i = 0; i < num_edges; ++i) {
    while (true) {
      size_t const u_id = rng() % num_nodes;
      size_t const v_id = rng() % num_nodes;
      size_t const uv_edge_index = u_id * num_nodes + v_id;
      size_t const vu_edge_index = v_id * num_nodes + u_id;

      if (u_id != v_id && edge_matrix[uv_edge_index] == false && edge_matrix[vu_edge_index] == false) {
        edge_matrix[uv_edge_index] = true;
        edge_matrix[vu_edge_index] = true;

        label.clear();
        label[EDGE_LABEL_KEY] = edge_labels[rng() % num_edge_labels];

        env.add_edge(graph_id, std::to_string(u_id), std::to_string(v_id), label, false);
        break;
      }
    }
  }

  return graph_id;
}

}

#endif
