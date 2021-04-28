#ifndef SRC_METHODS_DECISION_TREE_REP_HPP_
#define SRC_METHODS_DECISION_TREE_REP_HPP_


namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
class DecisionTree_REPTree : public DecisionTreeMethod<UserNodeLabel, UserEdgeLabel> {
public:
  DecisionTree_REPTree(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data, GEDEnv<GXLNodeID, UserNodeLabel, UserEdgeLabel>* env)
    : DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>(ged_data, env) {

    this->m_methods[Options::GEDMethod::NODE] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new Node<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[Options::GEDMethod::BRANCH_FAST] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new BranchFast<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[Options::GEDMethod::BRANCH_UNIFORM] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new BranchUniform<UserNodeLabel, UserEdgeLabel>(ged_data));
  }

private:
  virtual Options::GEDMethod pick_method(GEDGraph const& g, GEDGraph const& h) {
    //  edge_label_count < 0.1
    // |    avg_node_degree < 0.7 : IPFP (BRANCH_FAST) (805/252) [406/155]
    // |    avg_node_degree >= 0.7
    // |   |    edge_density < 0.99 : IPFP (BRANCH_FAST) (1662/1288) [817/613]
    // |   |    edge_density >= 0.99 : IPFP (NODE) (134/75) [65/35]
    //  edge_label_count >= 0.1 : IPFP (BRANCH_UNIFORM) (36479/30711) [18252/15375]

    FastGraphDiff diff = this->compute_graph_diff(g, h);
    if (diff.edge_label_count < 0.1) {
      if (diff.avg_node_degree < 0.7) {
        return Options::GEDMethod::BRANCH_FAST;
      } else {
        if (diff.edge_density < 0.99) {
          return Options::GEDMethod::BRANCH_FAST;
        } else {
          return Options::GEDMethod::NODE;
        }
      }
    } else {
      return Options::GEDMethod::BRANCH_UNIFORM;
    }

    assert(false);
  }
};

}

#endif