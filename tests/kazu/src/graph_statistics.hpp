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
  size_t min_node_index = infinite_distance;

  for (size_t i = 0; i < g.num_nodes; ++i) {
    size_t eccentricity = 0;

    for (size_t j = 0; j < g.num_nodes; ++j) {
      if (distances[i * g.num_nodes + j] > eccentricity && distances[i * g.num_nodes + j] != infinite_distance) {
        eccentricity = distances[i * g.num_nodes + j];
      }
    }

    if (min_eccentricity == infinite_distance) {
      min_eccentricity = eccentricity;
      min_node_index = i;
      continue;
    }

    // Find radius of the largest component ("disconnected" sub-graph)
    if (distances[min_node_index * g.num_nodes + i] == infinite_distance && min_eccentricity < eccentricity) {
      min_eccentricity = eccentricity;
      min_node_index = i;
    } else if (distances[min_node_index * g.num_nodes + i] != infinite_distance && min_eccentricity > eccentricity) {
      min_eccentricity = eccentricity;
      min_node_index = i;
    }
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
  double edge_density;
  bool is_joint;
  size_t diameter;
  size_t radius;
  size_t min_node_degree;
  size_t max_node_degree;
};

using GraphStatsMap = std::map<GEDGraph::GraphID, GraphStatistics>;


GraphStatistics graph_stat_compute_single(GxlGEDEnv& env, GEDGraph::GraphID graph_id) {
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

  double const max_edges = static_cast<double>(g.num_nodes * (g.num_nodes - 1)) / 2.0;
  stats.edge_density = static_cast<double>(g.num_edges) / max_edges;

  return stats;
}

void graph_stat_compute_all(GxlGEDEnv& env, GraphStatsMap& stats) {
  stats.clear();
  for (GEDGraph::GraphID i = env.graph_ids().first; i != env.graph_ids().second; ++i)
    stats[i] = graph_stat_compute_single(env, i);
}

// Structure for holding comparison information of two given graphs
// Most values are callculated as: min_X(g1, g2) / max_X(g1, g2), where X refers
// to a specific attribute of the graph, e.g. node_count.
struct GraphDiff {
  double node_count;
  double edge_count;
  double avg_node_degree;
  double edge_density;

  // Ratio of matching node label count.
  double node_label_count;
  // Ratio of matching node label count.
  double edge_label_count;

  double diameter;
  double radius;
  bool both_joint;
};


// constexpr double finger_tolerance = M_PI / 100.0;

// struct LabelKey_Finger {
//   double orientation;
// };

// bool operator< (ged::GXLLabel const& lhs, LabelKey_Finger const& rhs) {
//   double const _lhs = std::stod(lhs.at("orient"));
//   double const _rhs = rhs.orientation;

//   return _lhs - _rhs < finger_tolerance;
// }
// bool operator< (LabelKey_Finger const& lhs, ged::GXLLabel const& rhs) {
//   double const _lhs = lhs.orientation;
//   double const _rhs = std::stod(rhs.at("orient"));

//   return _lhs - _rhs < finger_tolerance;
// }

double graph_diff_compute_node_labels(GxlExchangeGraph const& g, GxlExchangeGraph const& h) {
  auto count_labels = [](GxlExchangeGraph const& g) -> std::map<ged::GXLLabel, size_t> {
    std::map<ged::GXLLabel, size_t> label_map;

    for (auto const& label : g.node_labels) {
      auto iterator = label_map.find(label);
      if (iterator != label_map.end())
        (*iterator).second += 1;
      else
        label_map[label] = 1;
    }

    return label_map;
  };

  std::map<ged::GXLLabel, size_t> g_label_map = count_labels(g);
  std::map<ged::GXLLabel, size_t> h_label_map = count_labels(h);

  size_t matching_labels = 0;
  for (auto const& g_label : g_label_map) {
    auto h_it = h_label_map.find(g_label.first);
    if (h_it != h_label_map.end())
      matching_labels += std::min(g_label.second, (*h_it).second);
  }

  size_t const max_nodes = std::max(g.num_nodes, h.num_nodes);
  if (max_nodes == 0)
    return 0;

  return static_cast<double>(matching_labels) / static_cast<double>(max_nodes);
}

double graph_diff_compute_edge_labels(GxlExchangeGraph const& g, GxlExchangeGraph const& h) {
  auto count_labels = [](GxlExchangeGraph const& g) -> std::map<ged::GXLLabel, size_t> {
    std::map<ged::GXLLabel, size_t> label_map;

    for (auto const& label : g.edge_labels) {
      auto iterator = label_map.find(label.second);
      if (iterator != label_map.end())
        (*iterator).second += 1;
      else
        label_map[label.second] = 1;
    }

    return label_map;
  };

  std::map<ged::GXLLabel, size_t> g_label_map = count_labels(g);
  std::map<ged::GXLLabel, size_t> h_label_map = count_labels(h);

  size_t matching_labels = 0;
  for (auto const& g_label : g_label_map) {
    auto h_it = h_label_map.find(g_label.first);
    if (h_it != h_label_map.end())
      matching_labels += std::min(g_label.second, (*h_it).second);
  }

  // We are working with undirected graphs
  assert(matching_labels % 2 == 0);
  matching_labels /= 2;

  size_t const max_edges = std::max(g.num_edges, h.num_edges);
  if (max_edges == 0)
    return 0;

  return static_cast<double>(matching_labels) / static_cast<double>(max_edges);
}

template<typename T>
double _calculate_ratio(T a, T b) {
  double const min = static_cast<double>(std::min(a, b));
  double const max = static_cast<double>(std::max(a, b));

  if (max <= 0)
    return 0;
  else
    return min / max;
}

GraphDiff graph_diff_compute(GxlGEDEnv& env, GraphStatsMap const& stats, GEDGraph::GraphID g_id, GEDGraph::GraphID h_id) {
  assert(stats.find(g_id) != stats.end());
  assert(stats.find(h_id) != stats.end());

  auto const g = env.get_graph(g_id);
  auto const h = env.get_graph(h_id);
  GraphStatistics const& g_stats = stats.at(g_id);
  GraphStatistics const& h_stats = stats.at(h_id);

  GraphDiff diff;
  diff.node_count = _calculate_ratio(g.num_nodes, h.num_nodes);
  diff.edge_count = _calculate_ratio(g.num_edges, h.num_edges);
  diff.avg_node_degree = _calculate_ratio(g_stats.avg_degree, h_stats.avg_degree);
  diff.edge_density = _calculate_ratio(g_stats.edge_density, h_stats.edge_density);
  diff.node_label_count = graph_diff_compute_node_labels(g, h);
  diff.edge_label_count = graph_diff_compute_edge_labels(g, h);
  diff.diameter = _calculate_ratio(g_stats.diameter, h_stats.diameter);
  diff.radius = _calculate_ratio(g_stats.radius, h_stats.radius);
  diff.both_joint = g_stats.is_joint && h_stats.is_joint;

  return diff;
}

}

#endif