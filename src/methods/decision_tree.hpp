#ifndef SRC_METHODS_DECISION_TREE_HPP_
#define SRC_METHODS_DECISION_TREE_HPP_

#include <map>
#include <memory>
#include <cassert>

namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
class DecisionTreeMethod : public GEDMethod<UserNodeLabel, UserEdgeLabel> {
public:
  using GEDMethodPtr = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>;
  using GEDMethodMap = std::map<Options::GEDMethod, GEDMethodPtr>;

public:
  DecisionTreeMethod(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data);

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
};



template<class UserNodeLabel, class UserEdgeLabel> 
DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>::
DecisionTreeMethod(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data)
  : GEDMethod<UserNodeLabel, UserEdgeLabel>(ged_data),
    m_ged_data(ged_data) {
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




// struct GraphDiff {
//   double node_count;
//   double edge_count;
//   double avg_node_degree;
//   double edge_density;
//   double diameter;
//   double radius;
//   bool both_joint;
// };


}

#endif