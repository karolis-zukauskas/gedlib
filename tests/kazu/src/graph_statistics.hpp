#ifndef SRC_KAZU_GRAPH_STATISTICS_HPP_
#define SRC_KAZU_GRAPH_STATISTICS_HPP_

#include "util.hpp"
#include <numeric>
#include <cassert>

namespace ged {

constexpr size_t infinite_distance = std::numeric_limits<size_t>::max();

std::vector<size_t> _apsp_floyd(GxlExchangeGraph const& g) {
  std::vector<size_t> distances(g.num_nodes * g.num_nodes, infinite_distance);

  assert(g.num_nodes == g.adj_matrix.size());
  for (size_t i = 0; i < g.num_nodes; ++i) {
    for (size_t j = 0; j < g.num_nodes; ++j) {
      assert(g.num_nodes == g.adj_matrix[i].size());

      if (g.adj_matrix[i][j] != 0)
        distances[i * g.num_nodes + j] = 1;
    }
  }

  for (size_t i = 0; i < g.num_nodes; ++i)
    distances[i * g.num_nodes + i] = 0;

  for (size_t k = 0; k < g.num_nodes; ++k) {
    for (size_t i = 0; i < g.num_nodes; ++i) {
      for (size_t j = 0; j < g.num_nodes; ++j) {
        size_t const ij = i * g.num_nodes + j;
        size_t const ik = i * g.num_nodes + k;
        size_t const kj = k * g.num_nodes + j;

        if (distances[ik] != infinite_distance && distances[kj] != infinite_distance && distances[ij] >= distances[ik] + distances[kj]) {
          distances[ij] = distances[ik] + distances[kj];
        }
      }
    }
  }

  return distances;
}

size_t _graph_stat_isJoint(std::vector<size_t> const& distances) {
  for (size_t i = 0; i < distances.size(); ++i) {
    if (distances[i] == infinite_distance)
      return false;
  }

  return true;
}

size_t _graph_stat_diameter(std::vector<size_t> const& distances) {
  size_t max_distance = 0;
  for (size_t i = 0; i < distances.size(); ++i) {
    if (distances[i] > max_distance && distances[i] != infinite_distance)
      max_distance = distances[i];
  }

  return max_distance;
}

size_t _graph_stat_radius(GxlExchangeGraph const& g, std::vector<size_t> const& distances) {
  size_t min_eccentricity = infinite_distance;
  for (size_t i = 0; i < g.num_nodes; ++i) {
    size_t eccentricity = 0;
    for (size_t j = 0; j < g.num_nodes; ++j) {
      if (distances[i * g.num_nodes + j] > eccentricity && distances[i * g.num_nodes + j] != infinite_distance) {
        eccentricity = distances[i * g.num_nodes + j];
      }
    }

    if (min_eccentricity > eccentricity)
      min_eccentricity = eccentricity;
  }

  return min_eccentricity;
}

void _graph_stat_node_degrees(GxlExchangeGraph const& g, size_t& min_degree, size_t& max_degree, double& avg_degree) {
  min_degree = infinite_distance;
  max_degree = 0;
  avg_degree = 0;

  for (size_t i = 0; i < g.adj_matrix.size(); ++i) {
    size_t degree = std::accumulate(g.adj_matrix[i].begin(), g.adj_matrix[i].end(), 0);
    avg_degree += static_cast<double>(degree);

    if (min_degree > degree)
      min_degree = degree;
    if (max_degree < degree)
      max_degree = degree;
  }

  avg_degree /= static_cast<double>(g.num_nodes);
}

struct GraphStatistics {
  GEDGraph::GraphID id;

  size_t num_nodes;
  size_t num_edges;
  double avg_degree;
  bool is_joint;
  size_t diameter;
  size_t radius;
  size_t min_node_degree;
  size_t max_node_degree;
};

GraphStatistics graph_stat_compute(GxlGEDEnv& env, GEDGraph::GraphID graph_id) {
  auto g = env.get_graph(graph_id);
  std::vector<size_t> distances = _apsp_floyd(g);

  GraphStatistics stats;
  stats.id = graph_id;
  stats.num_nodes = g.num_nodes;
  stats.num_edges = g.num_edges;
  stats.is_joint = _graph_stat_isJoint(distances);
  stats.diameter = _graph_stat_diameter(distances);
  stats.radius = _graph_stat_radius(g, distances);
  _graph_stat_node_degrees(g, stats.min_node_degree, stats.max_node_degree, stats.avg_degree);

  return stats;
}

}

#endif