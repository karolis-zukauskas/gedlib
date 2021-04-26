#ifndef SRC_METHODS_DECISION_TREE_REP_HPP_
#define SRC_METHODS_DECISION_TREE_REP_HPP_


namespace ged {

template<class UserNodeLabel, class UserEdgeLabel>
class DecisionTree_REPTree : public DecisionTreeMethod<UserNodeLabel, UserEdgeLabel> {
public:
  DecisionTree_REPTree(GEDData<UserNodeLabel, UserEdgeLabel> const& ged_data)
    : DecisionTreeMethod<UserNodeLabel, UserEdgeLabel>(ged_data) {

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