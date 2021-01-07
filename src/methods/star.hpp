/***************************************************************************
 *                                                                          *
 *   Copyright (C) 2018 by David B. Blumenthal                              *
 *                                                                          *
 *   This file is part of GEDLIB.                                           *
 *                                                                          *
 *   GEDLIB is free software: you can redistribute it and/or modify it      *
 *   under the terms of the GNU Lesser General Public License as published  *
 *   by the Free Software Foundation, either version 3 of the License, or   *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   GEDLIB is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           *
 *   GNU Lesser General Public License for more details.                    *
 *                                                                          *
 *   You should have received a copy of the GNU Lesser General Public       *
 *   License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                          *
 ***************************************************************************/

/*!
 * @file  star.hpp
 * @brief ged::Star class declaration.
 */

#ifndef SRC_METHODS_STAR_HPP_
#define SRC_METHODS_STAR_HPP_

namespace ged {

/*!
 * @brief Computes lower and upper bounds for uniform edit costs.
 * @details Implements the method Star suggested in:
 * - Z. Zeng, A. K. H. Tung, J. Wang, J. Feng, and L. Zhou:
 *   &ldquo;Comparing stars: On approximating graph edit distance&rdquo;,
 *   https://doi.org/10.1016/10.14778/1687627.1687631
 *
 * Supports the following option in addition to the ones supported by ged::LSAPEBasedMethod:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--sort-method STD\|COUNTING</tt> | the employed sorting algorithm | @p COUNTING | @ref ged::util::counting_sort() <br> use counting sort if the number of different edge labels is constant |
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Star : public LSAPEBasedMethod<UserNodeLabel, UserEdgeLabel> {

public:

	virtual ~Star();

	Star(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

public:

	enum SortMethod_ {STD, COUNTING};

	class SortedNodeLabels_ {
	public:
		SortedNodeLabels_(const GEDGraph & g, SortMethod_ sort_method);

		SortedNodeLabels_();

		void operator=(const SortedNodeLabels_ & sorted_edge_labels);

		const std::vector<LabelID> & get_incident_labels(GEDGraph::NodeID) const;

	private:
		std::map<GEDGraph::NodeID, std::vector<LabelID>> sorted_node_labels_;
	};

	SortMethod_ sort_method_;

	std::map<GEDGraph::GraphID, SortedNodeLabels_> sorted_node_labels_;

	// Inherited member functions from LSAPEBasedMethod.

	virtual void lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & lsape_instance);

	virtual bool lsape_parse_option_(const std::string & option, const std::string & arg) final;

	virtual std::string lsape_valid_options_string_() const final;

	virtual void lsape_set_default_options_() final;

	virtual void lsape_init_graph_(const GEDGraph & graph) final;

	virtual double lsape_lower_bound_scaling_factor_(const GEDGraph & g, const GEDGraph & h) final;

	// Helper member functions.

	virtual double compute_substitution_cost_(const GEDGraph & g, const GEDGraph & h, GEDGraph::NodeID i, GEDGraph::NodeID k,
			const SortedNodeLabels_ & sorted_node_labels_g, const SortedNodeLabels_ & sorted_node_labels_h) const;

	virtual double compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const;

	virtual double compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const;
};

/**
 * Modifications:
 * Insertion and deletion cost functions are taken from BIPARTINE/BRANCH.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Star2 : public Star<UserNodeLabel, UserEdgeLabel> {
public:
  Star2(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

  double compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const;
  double compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const;
};

/**
 * Modifications:
 * Insertion and deletion cost functions are taken from BIPARTINE/BRANCH and they are not multiplied by minimum edit costs
 * in lsape_populate_instance_().
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Star3 : public Star<UserNodeLabel, UserEdgeLabel> {
public:
  Star3(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

  virtual void lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & lsape_instance);
  double compute_deletion_cost_(const GEDGraph & g, GEDGraph::NodeID i) const;
  double compute_insertion_cost_(const GEDGraph & h, GEDGraph::NodeID k) const;
};

/**
 * Modifications:
 * Mean node edit consts are used instead of minimum costs in lsape_populate_instance_
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Star4 : public Star<UserNodeLabel, UserEdgeLabel> {
public:
  Star4(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

  virtual void lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & lsape_instance);
};

/**
 * Modifications:
 * 1) Mean node edit consts are used instead of minimum costs in lsape_populate_instance_.
 * 2) If the matrix is square, node insert/delete consts are set to infinity.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Star5 : public Star<UserNodeLabel, UserEdgeLabel> {
public:
  Star5(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

  virtual void lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & lsape_instance);
};

/**
 * Modifications:
 * 1) Mean node edit consts are used instead of minimum costs in lsape_populate_instance_.
 * 2) If the matrix is "almost square", node insert/delete consts are set to infinity. "Almost equal" is calculuted as
 * `abs(n - m) <= threshold` where threshold is calculated as `round(min(n, m) / 10)`.
 */
template<class UserNodeLabel, class UserEdgeLabel>
class Star6 : public Star<UserNodeLabel, UserEdgeLabel> {
public:
  Star6(const GEDData<UserNodeLabel, UserEdgeLabel> & ged_data);

  virtual void lsape_populate_instance_(const GEDGraph & g, const GEDGraph & h, DMatrix & lsape_instance);
};

}

#endif /* SRC_METHODS_STAR_HPP_ */