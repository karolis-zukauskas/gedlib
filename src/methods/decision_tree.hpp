#ifndef SRC_METHODS_DECISION_TREE_HPP_
#define SRC_METHODS_DECISION_TREE_HPP_

#include <map>
#include <memory>
#include <cassert>

namespace ged {

struct FastGraphDiff {
  // O(1) metrics
  double node_count;
  double edge_count;
  double avg_node_degree;
  double edge_density;

  // O(N) metrics
  double node_label_count;
  double edge_label_count;
};

template<class UserNodeLabel, class UserEdgeLabel>
class DecisionTreeMethod : public GEDMethod<UserNodeLabel, UserEdgeLabel> {
public:
  using GEDMethodPtr = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>;
  using GEDMethodMap = std::map<Options::GEDMethod, GEDMethodPtr>;

public:
  DecisionTreeMethod(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data, GEDEnv<GXLNodeID, UserNodeLabel, UserEdgeLabel>* pEnv);

  FastGraphDiff compute_graph_diff(GEDGraph const& g, GEDGraph const& h);

private:
  virtual Options::GEDMethod pick_method(GEDGraph const& g, GEDGraph const& h) = 0;

  // Member functions inherited from GEDMethod.
  virtual void ged_init_() final;
  virtual void ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) final;
  virtual bool ged_parse_option_(const std::string & option, const std::string & arg) final;
  virtual std::string ged_valid_options_string_() const final;
  virtual void ged_set_default_options_() final;

protected:
  GEDMethodMap m_methods;
  GEDData<UserNodeLabel, UserEdgeLabel> const& m_ged_data;
  GEDEnv<GXLNodeID, UserNodeLabel, UserEdgeLabel>* m_env;
};



template<class UserNodeLabel, class UserEdgeLabel> 
DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>::
DecisionTreeMethod(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data, GEDEnv<GXLNodeID, UserNodeLabel, UserEdgeLabel>* env)
  : GEDMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
    m_ged_data(ged_data),
    m_env(env) {
}

template<class UserNodeLabel, class UserEdgeLabel>
void
DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>::
ged_init_() {
  for (auto const& method : m_methods) {
    method.second->init();
  }
}

template<class UserNodeLabel, class UserEdgeLabel>
void
DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>::
ged_run_(const GEDGraph & g, const GEDGraph & h, Result & result) {
  Options::GEDMethod method_index = this->pick_method(g, h);
  assert(m_methods.find(method_index) != m_methods.end() && "Method is not initialized");

  m_methods[method_index]->run_as_util(g, h, result);
}

template<class UserNodeLabel, class UserEdgeLabel>
bool
DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>::
ged_parse_option_(const std::string & option, const std::string & arg) {
  return true;
}

template<class UserNodeLabel, class UserEdgeLabel>
std::string
DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>::
ged_valid_options_string_() const {
  return "";
}

template<class UserNodeLabel, class UserEdgeLabel>
void
DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>::
ged_set_default_options_() {
}


template<typename T>
double calculate_ratio(T a, T b) {
  double const min = static_cast<double>(std::min(a, b));
  double const max = static_cast<double>(std::max(a, b));

  if (max <= 0)
    return 0;
  else
    return min / max;
}

template<class UserNodeLabel, class UserEdgeLabel>
double graph_diff_compute_node_labels_temp(GEDEnv<GXLNodeID, UserNodeLabel, UserEdgeLabel>& env, GEDGraph const& g, GEDGraph const& h) {
  auto count_labels = [&env](ged::GEDGraph const& g) -> std::map<ged::GXLLabel, size_t> {
    std::map<ged::GXLLabel, size_t> label_map;

    for (auto it = g.nodes().first; it != g.nodes().second; it++) {
      LabelID label_id = g.get_node_label(*it);
      GXLLabel const& label = env.get_node_label_ref(label_id);

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

  size_t const max_nodes = std::max(g.num_nodes(), h.num_nodes());
  if (max_nodes == 0)
    return 0;

  return static_cast<double>(matching_labels) / static_cast<double>(max_nodes);
}

template<class UserNodeLabel, class UserEdgeLabel>
double graph_diff_compute_edge_labels_temp(GEDEnv<GXLNodeID, UserNodeLabel, UserEdgeLabel>& env, GEDGraph const& g, GEDGraph const& h) {
  auto count_labels = [&env](ged::GEDGraph const& g) -> std::map<ged::GXLLabel, size_t> {
    std::map<ged::GXLLabel, size_t> label_map;

    for (auto it = g.edges().first; it != g.edges().second; it++) {
      LabelID label_id = g.get_edge_label(*it);
      GXLLabel const& label = env.get_edge_label_ref(label_id);

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

  size_t const max_edges = std::max(g.num_edges(), h.num_edges());
  if (max_edges == 0)
    return 0;

  return static_cast<double>(matching_labels) / static_cast<double>(max_edges);
}

template<class UserNodeLabel, class UserEdgeLabel>
FastGraphDiff
DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>::
compute_graph_diff(GEDGraph const& g, GEDGraph const& h) {

  FastGraphDiff diff;
  diff.node_count = calculate_ratio(g.num_nodes(), h.num_nodes());
  diff.edge_count = calculate_ratio(g.num_edges(), h.num_edges());

  auto compute_avg_degree = [](GEDGraph const& g) -> double {
    double avg_degree = 0;
    for (auto it = g.nodes().first; it != g.nodes().second; it++) {
      avg_degree += static_cast<double>(g.degree(*it));
    }

    return avg_degree / static_cast<double>(g.num_nodes());
  };

  diff.avg_node_degree = calculate_ratio(compute_avg_degree(g), compute_avg_degree(h));

  auto compute_edge_density = [](GEDGraph const& g) -> double {
    double const max_edges = static_cast<double>(g.num_nodes() * (g.num_nodes() - 1)) / 2.0;
    if (max_edges == 0)
      return 0;

    return static_cast<double>(g.num_edges()) / max_edges;
  };

  diff.edge_density = calculate_ratio(compute_edge_density(g), compute_edge_density(h));

  diff.node_label_count = graph_diff_compute_node_labels_temp(*(this->m_env), g, h);
  diff.edge_label_count = graph_diff_compute_edge_labels_temp(*(this->m_env), g, h);

  return diff;
}

}

#endif