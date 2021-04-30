#ifndef SRC_METHODS_DECISION_TREE_REP_HPP_
#define SRC_METHODS_DECISION_TREE_REP_HPP_


namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
class DecisionTree_REPTree : public DecisionTreeMethod<UserNodeLabel, UserEdgeLabel> {
private:
  enum {
    NODE = 0,
    BRANCH_FAST = 1,
    BRANCH_UNIFORM = 2,
  };

public:
  DecisionTree_REPTree(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data, GEDEnv<GXLNodeID, UserNodeLabel, UserEdgeLabel>* env)
    : DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>(ged_data, env) {

    this->m_methods[NODE] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new Node<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[BRANCH_FAST] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new BranchFast<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[BRANCH_UNIFORM] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new BranchUniform<UserNodeLabel, UserEdgeLabel>(ged_data));
  }

private:
  virtual size_t pick_method(GEDGraph const& g, GEDGraph const& h) {
    //  edge_label_count < 0.1
    // |    avg_node_degree < 0.7 : IPFP (BRANCH_FAST) (805/252) [406/155]
    // |    avg_node_degree >= 0.7
    // |   |    edge_density < 0.99 : IPFP (BRANCH_FAST) (1662/1288) [817/613]
    // |   |    edge_density >= 0.99 : IPFP (NODE) (134/75) [65/35]
    //  edge_label_count >= 0.1 : IPFP (BRANCH_UNIFORM) (36479/30711) [18252/15375]

    FastGraphDiff diff = this->compute_graph_diff(g, h);
    if (diff.edge_label_count < 0.1) {
      if (diff.avg_node_degree < 0.7) {
        return BRANCH_FAST;
      } else {
        if (diff.edge_density < 0.99) {
          return BRANCH_FAST;
        } else {
          return NODE;
        }
      }
    } else {
      return BRANCH_UNIFORM;
    }

    assert(false);
  }
};

template<class UserNodeLabel, class UserEdgeLabel>
class DecisionTree_REPTree2 : public DecisionTreeMethod<UserNodeLabel, UserEdgeLabel> {
private:
  enum {
    STAR = 0,
    STAR6 = 1,
    BRANCH_UNIFORM = 2,
  };

public:
  DecisionTree_REPTree2(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data, GEDEnv<GXLNodeID, UserNodeLabel, UserEdgeLabel>* env)
    : DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>(ged_data, env) {

    this->m_methods[STAR] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new Star<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[STAR6] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new Star6<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[BRANCH_UNIFORM] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new BranchUniform<UserNodeLabel, UserEdgeLabel>(ged_data));
  }

private:
  virtual size_t pick_method(GEDGraph const& g, GEDGraph const& h) {
    // Created from generated graph data only

    //  edge_label_count < 0.5 : IPFP (BRANCH_UNIFORM) (2302/1925) [1128/931]
    //  edge_label_count >= 0.5
    // |    node_count < 0.85
    // |   |    edge_density < 0.91 : IPFP (BRANCH_UNIFORM) (1461/1281) [732/633]
    // |   |    edge_density >= 0.91 : IPFP (STAR) (1467/1299) [700/630]
    // |    node_count >= 0.85
    // |   |    node_label_count < 0.64 : IPFP (BRANCH_UNIFORM) (375/321) [196/165]
    // |   |    node_label_count >= 0.64 : IPFP (STAR6) (5011/4462) [2553/2272]

    FastGraphDiff diff = this->compute_graph_diff(g, h);
    if (diff.edge_label_count < 0.5) {
      return BRANCH_UNIFORM;
    } else {
      if (diff.node_count < 0.85) {
        if (diff.edge_density < 0.91) {
          return BRANCH_UNIFORM;
        } else {
          return STAR;
        }
      } else {
        if (diff.node_label_count < 0.64) {
          return BRANCH_UNIFORM;
        } else {
          return STAR6;
        }
      }
    }

    assert(false);
  }
};

}

#endif