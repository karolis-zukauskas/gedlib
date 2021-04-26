#ifndef SRC_METHODS_DECISION_TREE_HPP_
#define SRC_METHODS_DECISION_TREE_HPP_

#include <map>
#include <memory>
#include <cassert>

namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
class DecisionTree {
public:
  using GEDMethodPtr = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>;
  using GEDMethodMap = std::map<Options::GEDMethod, GEDMethodPtr>;

public:
  DecisionTree(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data);
  void generate_node_map(GEDGraph const& g, GEDGraph const& h, Result& result);

private:
  virtual Options::GEDMethod pick_method(GEDGraph const& g, GEDGraph const& h) = 0;

protected:
  GEDMethodMap m_methods;
  GEDData<UserNodeLabel, UserEdgeLabel> const& m_ged_data;
};



template<class UserNodeLabel, class UserEdgeLabel> 
DecisionTree<UserNodeLabel, UserEdgeLabel>::
DecisionTree(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data)
  : m_ged_data(ged_data),
    m_methods(GEDMethodMap()) {
}

template<class UserNodeLabel, class UserEdgeLabel>
void DecisionTree<UserNodeLabel, UserEdgeLabel>::
generate_node_map(GEDGraph const& g, GEDGraph const& h, Result& result) {
  Options::GEDMethod method_index = this->pick_method(g, h);
  assert(m_methods.find(method_index) != m_methods.end() && "Method is not initialized");

  m_methods[method_index]->run_as_util(g, h, result);
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


template<class UserNodeLabel, class UserEdgeLabel>
class DecisionTree_REPTree : public DecisionTree<UserNodeLabel, UserEdgeLabel> {
public:
  DecisionTree_REPTree(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data)
    : DecisionTree<UserNodeLabel, UserEdgeLabel>(ged_data) {

    this->m_methods[Options::GEDMethod::NODE] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new Node<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[Options::GEDMethod::BRANCH_FAST] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new BranchFast<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[Options::GEDMethod::STAR] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new Star<UserNodeLabel, UserEdgeLabel>(ged_data));
  }

private:
  virtual Options::GEDMethod pick_method(GEDGraph const& g, GEDGraph const& h) {

    return Options::GEDMethod::NODE;
  }
};

}

#endif