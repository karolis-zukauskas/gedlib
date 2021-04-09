#ifndef SRC_KAZU_GRAPH_GEN_HPP_
#define SRC_KAZU_GRAPH_GEN_HPP_

#include "util.hpp"
#include <numeric>

namespace ged {

GEDGraph::GraphID graph_gen_random(GxlGEDEnv& env, std::mt19937& rng, size_t num_nodes, double edge_density,
  std::vector<std::string> const& node_labels, std::vector<std::string> const& edge_labels) {

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

size_t _find_node_index (size_t r, std::vector<size_t> const& prefix_sums) {
  size_t low = 0;
  size_t high = prefix_sums.size() - 1;
  while (low < high)
    {
    size_t mid = (low + high) / 2;
    if (r <= prefix_sums[low])
      return low;
    else if (low == mid)
      return high;
    else if (r <= prefix_sums[mid])
      high = mid;
    else
      low = mid;
    }

  return low;
}

GEDGraph::GraphID graph_gen_power(GxlGEDEnv& env, std::mt19937& rng, size_t num_nodes, size_t edges_per_node,
  std::vector<std::string> const& node_labels, std::vector<std::string> const& edge_labels) {

  constexpr char const* NODE_LABEL_KEY = "label";
  constexpr char const* EDGE_LABEL_KEY = "label";

  size_t const num_node_labels = node_labels.size();
  size_t const num_edge_labels = edge_labels.size();

  GEDGraph::GraphID const graph_id = env.add_graph();
  GXLLabel label;

  size_t const initial_node_count = edges_per_node + 2;
  for (size_t i = 0; i < initial_node_count; ++i) {
    label.clear();
    label[NODE_LABEL_KEY] = node_labels[rng() % num_node_labels];
    env.add_node(graph_id, std::to_string(i), label);
  }

  for (size_t i = 0; i < initial_node_count - 1; ++i) {
    label.clear();
    label[EDGE_LABEL_KEY] = edge_labels[rng() % num_edge_labels];
    env.add_edge(graph_id, std::to_string(i), std::to_string(i + 1), label, false);
  }

  label.clear();
  label[EDGE_LABEL_KEY] = edge_labels[rng() % num_edge_labels];
  env.add_edge(graph_id, std::to_string(0), std::to_string(initial_node_count - 1), label, false);

  std::vector<size_t> prefix_sums;
  for (size_t i = 0; i < initial_node_count; ++i) {
    size_t const degree = env.node_degree(graph_id, std::to_string(i));

    if (i == 0) {
      prefix_sums.push_back(degree);
    } else {
      prefix_sums.push_back(prefix_sums[i - 1]);
      prefix_sums[i] += degree;
    }
  }

  for (size_t i = initial_node_count; i < num_nodes; ++i) {
    label.clear();
    label[NODE_LABEL_KEY] = node_labels[rng() % num_node_labels];
    env.add_node(graph_id, std::to_string(i), label);
    prefix_sums.push_back(prefix_sums[i - 1]);

    size_t edges_added = 0;
    while (edges_added < edges_per_node) {
      size_t r = (rng() % (prefix_sums[i - 1] - 1)) + 1;
      size_t index = _find_node_index(r, prefix_sums);

      bool const edge_exists = env.is_edge(graph_id, std::to_string(i), std::to_string(index));
      if (edge_exists)
        continue;

      label.clear();
      label[EDGE_LABEL_KEY] = edge_labels[rng() % num_edge_labels];

      assert(i != index && "self loops should not happen");
      env.add_edge(graph_id, std::to_string(i), std::to_string(index), label, false);

      edges_added++;

      for (size_t j = index; j <= i; ++j)
        prefix_sums[j]++;
      prefix_sums[i]++;
    }
  }

  return graph_id;
}

GEDGraph::GraphID graph_gen_cluster(GxlGEDEnv& env, std::mt19937& rng, size_t num_nodes, size_t num_clusters,
  double p_inner_edge, double p_outer_edge, std::vector<std::string> const& node_labels, std::vector<std::string> const& edge_labels) {

  assert(num_nodes > 1);
  assert(num_clusters > 0);
  assert(num_nodes % num_clusters == 0);

  constexpr char const* NODE_LABEL_KEY = "label";
  constexpr char const* EDGE_LABEL_KEY = "label";

  size_t const num_node_labels = node_labels.size();
  size_t const num_edge_labels = edge_labels.size();

  GEDGraph::GraphID const graph_id = env.add_graph();
  GXLLabel label;

  for (size_t i = 0; i < num_nodes; ++i) {
    label.clear();
    label[NODE_LABEL_KEY] = node_labels[rng() % num_node_labels];
    env.add_node(graph_id, std::to_string(i), label);
  }

  for (size_t i = 0; i < num_nodes - 1; ++i) {
    label.clear();
    label[EDGE_LABEL_KEY] = edge_labels[rng() % num_edge_labels];
    env.add_edge(graph_id, std::to_string(i), std::to_string(i + 1), label, false);
  }

  label.clear();
  label[EDGE_LABEL_KEY] = edge_labels[rng() % num_edge_labels];
  env.add_edge(graph_id, std::to_string(0), std::to_string(num_nodes - 1), label, false);


  std::uniform_real_distribution<> p_distribution(0.0, 1.0);
  size_t const cluster_size = num_nodes / num_clusters;
  double p_edge = 0.0;

  for (size_t i = 0; i < num_nodes; ++i) {
    for (size_t j = 0; j < num_nodes; ++j) {
      if (i == j)
        continue;

      if (i / cluster_size == j / cluster_size) {
        if ((i + 1) % cluster_size == j % cluster_size)
          p_edge = 1.0;
        else
          p_edge = p_inner_edge;
      } else {
        p_edge = p_outer_edge;
      }

      if (p_distribution(rng) <= p_edge) {
        label.clear();
        label[EDGE_LABEL_KEY] = edge_labels[rng() % num_edge_labels];

        env.add_edge(graph_id, std::to_string(i), std::to_string(j), label, true);
      }
    }
  }

  return graph_id;
}

}

#endif
