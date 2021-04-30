#ifndef SRC_METHODS_DECISION_TREE_J48_HPP_
#define SRC_METHODS_DECISION_TREE_J48_HPP_


namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
class DecisionTree_J48 : public DecisionTreeMethod<UserNodeLabel, UserEdgeLabel> {
private:
  enum {
    BRANCH_FAST = 0,
    REFINE_BRANCH_UNIFORM = 1,
    BP_BEAM_BRANCH_UNIFORM = 2,
    BRANCH_UNIFORM = 3,
    STAR = 4,
  };

public:
  DecisionTree_J48(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data, GEDEnv<GXLNodeID, UserNodeLabel, UserEdgeLabel>* env)
    : DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>(ged_data, env) {

    this->m_methods[BRANCH_FAST] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new BranchFast<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[REFINE_BRANCH_UNIFORM] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new Refine<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[BP_BEAM_BRANCH_UNIFORM] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new BPBeam<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[BRANCH_UNIFORM] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new BranchUniform<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[STAR] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new Star<UserNodeLabel, UserEdgeLabel>(ged_data));

    this->m_methods[REFINE_BRANCH_UNIFORM]->set_options("--threads 1 --initialization-method BRANCH_UNIFORM");
    this->m_methods[BP_BEAM_BRANCH_UNIFORM]->set_options("--threads 1 --initialization-method BRANCH_UNIFORM");
  }

private:
  virtual size_t pick_method(GEDGraph const& g, GEDGraph const& h) {
    //  edge_label_count <= 0.096774: IPFP (BRANCH_FAST) (3889.0/2499.0)
    //  edge_label_count > 0.096774
    // |    node_label_count <= 0.111111
    // |   |    edge_label_count <= 0.518519: IPFP (REFINE) ( --ls-initialization-method BRANCH_UNIFORM) (2002.0/1794.0)
    // |   |    edge_label_count > 0.518519: IPFP (BP_BEAM) ( --ls-initialization-method BRANCH_UNIFORM) (2534.0/2314.0)
    // |    node_label_count > 0.111111
    // |   |    edge_label_count <= 0.574074
    // |   |   |    avg_node_degree <= 0.857963: IPFP (BRANCH_UNIFORM) (7285.0/5767.0)
    // |   |   |    avg_node_degree > 0.857963: IPFP (STAR) (3642.0/3089.0)
    // |   |    edge_label_count > 0.574074: IPFP (BRANCH_UNIFORM) (39268.0/33033.0)

    FastGraphDiff diff = this->compute_graph_diff(g, h);
    if (diff.edge_label_count <= 0.096774) {
      return BRANCH_FAST;
    } else {
      if (diff.node_label_count <= 0.111111) {
        if (diff.edge_label_count <= 0.518519) {
          return REFINE_BRANCH_UNIFORM;
        } else {
          return BP_BEAM_BRANCH_UNIFORM;
        }
      } else {
        if (diff.edge_label_count <= 0.574074) {
          if (diff.avg_node_degree <= 0.857963) {
            return BRANCH_UNIFORM;
          } else {
            return STAR;
          }
        } else {
          return BRANCH_UNIFORM;
        }
      }
    }

    assert(false);
  }
};


template<class UserNodeLabel, class UserEdgeLabel>
class DecisionTree_J48_2 : public DecisionTreeMethod<UserNodeLabel, UserEdgeLabel> {
private:
  enum {
    NODE = 0,
    BRANCH_FAST = 1,
    STAR = 2,
  };

public:
  DecisionTree_J48_2(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data, GEDEnv<GXLNodeID, UserNodeLabel, UserEdgeLabel>* env)
    : DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>(ged_data, env) {

    this->m_methods[BRANCH_FAST] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new BranchFast<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[NODE] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new Node<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[STAR] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new Star<UserNodeLabel, UserEdgeLabel>(ged_data));
  }

private:
  virtual size_t pick_method(GEDGraph const& g, GEDGraph const& h) {
    //  edge_label_count <= 0.090909: IPFP (BRANCH_FAST) (4169.0/2780.0)
    //  edge_label_count > 0.090909
    // |    node_label_count <= 0.068965: IPFP (NODE) (12027.0/8166.0)
    // |    node_label_count > 0.068965: IPFP (STAR) (10657.0/8709.0)

    FastGraphDiff diff = this->compute_graph_diff(g, h);
    if (diff.edge_label_count <= 0.090909) {
      return BRANCH_FAST;
    } else {
      if (diff.node_label_count <= 0.068965) {
        return NODE;
      } else {
        return STAR;
      }
    }

    assert(false);
  }
};



template<class UserNodeLabel, class UserEdgeLabel>
class DecisionTree_J48_3 : public DecisionTreeMethod<UserNodeLabel, UserEdgeLabel> {
private:
  enum {
    NODE = 0,
    STAR = 1,
  };

public:
  DecisionTree_J48_3(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data, GEDEnv<GXLNodeID, UserNodeLabel, UserEdgeLabel>* env)
    : DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>(ged_data, env) {

    this->m_methods[NODE] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new Node<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[STAR] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new Star<UserNodeLabel, UserEdgeLabel>(ged_data));
  }

private:
  virtual size_t pick_method(GEDGraph const& g, GEDGraph const& h) {
    //  avg_node_degree <= 0.833977: IPFP (NODE) (8052.0/5571.0)
    //  avg_node_degree > 0.833977
    // |    edge_count <= 0.5: IPFP (STAR) (6261.0/5251.0)
    // |    edge_count > 0.5: IPFP (NODE) (12540.0/10485.0)

    FastGraphDiff diff = this->compute_graph_diff(g, h);
    if (diff.avg_node_degree <= 0.833977) {
      return NODE;
    } else {
      if (diff.edge_count <= 0.5) {
        return STAR;
      } else {
        return STAR;
      }
    }

    assert(false);
  }
};



template<class UserNodeLabel, class UserEdgeLabel>
class DecisionTree_J48_4 : public DecisionTreeMethod<UserNodeLabel, UserEdgeLabel> {
private:
  enum {
    BRANCH_UNIFORM = 0,
    STAR6 = 1,
  };

public:
  DecisionTree_J48_4(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data, GEDEnv<GXLNodeID, UserNodeLabel, UserEdgeLabel>* env)
    : DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>(ged_data, env) {

    this->m_methods[BRANCH_UNIFORM] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new BranchUniform<UserNodeLabel, UserEdgeLabel>(ged_data));
    this->m_methods[STAR6] = std::unique_ptr<GEDMethod<UserNodeLabel, UserEdgeLabel>>(new Star6<UserNodeLabel, UserEdgeLabel>(ged_data));
  }

private:
  virtual size_t pick_method(GEDGraph const& g, GEDGraph const& h) {
    // Created from generated graph data only

    //  edge_label_count <= 0.50237: IPFP (BRANCH_UNIFORM) (3431.0/2857.0)
    //  edge_label_count > 0.50237
    // |    node_label_count <= 0.642857: IPFP (BRANCH_UNIFORM) (1510.0/1301.0)
    // |    node_label_count > 0.642857
    // |   |    node_count <= 0.883721
    // |   |   |    edge_density <= 0.919192: IPFP (BRANCH_UNIFORM) (2401.0/2139.0)
    // |   |   |    edge_density > 0.919192: IPFP (STAR6) (2600.0/2317.0)
    // |   |    node_count > 0.883721: IPFP (STAR6) (5983.0/5319.0)

    FastGraphDiff diff = this->compute_graph_diff(g, h);
    if (diff.edge_label_count <= 0.50237) {
      return BRANCH_UNIFORM;
    } else {
      if (diff.node_label_count <= 0.642857) {
        return BRANCH_UNIFORM;
      } else {
        if (diff.node_count <= 0.883721) {
          if (diff.edge_density <= 0.919192) {
            return BRANCH_UNIFORM;
          } else {
            return STAR6;
          }
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
