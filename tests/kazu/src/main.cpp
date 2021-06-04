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

// #define WRITE_STATS_FILE
constexpr int TEST_THREADS = 8;
constexpr bool TEST_ONLY_UNIQUE_PAIRS = false;

#include "util.hpp"
#include "graph_gen.hpp"
#include "graph_statistics.hpp"
#include "tests/ls.hpp"
#include "tests/lsape.hpp"
#include "tests/statistics.hpp"

int main(int argc, char* argv[]) {
#ifdef WRITE_STATS_FILE
  s_stats_filename = ::util::create_stats_file("all");
#endif

  std::vector<std::string> const all_datasets = {
    "Letter_HIGH", "Mutagenicity", "AIDS", "Protein", "GREC", "Fingerprint",
  };
  std::vector<std::string> const sized_datasets = {
    "Mutagenicity", "AIDS", "Protein",
  };

  std::vector<Method> const methods {
    Method (Options::GEDMethod::IPFP, "RANDOM", true),

    Method (Options::GEDMethod::IPFP, "BRANCH_UNIFORM", true),
    Method (Options::GEDMethod::IPFP, "BRANCH_FAST", true),
    Method (Options::GEDMethod::IPFP, "BRANCH", true),
    Method (Options::GEDMethod::IPFP, "NODE", true),
    Method (Options::GEDMethod::IPFP, "STAR", true),

    // NOTE: These two require training ... see Method options
    Method (Options::GEDMethod::IPFP, "WALKS", true),
    Method (Options::GEDMethod::IPFP, "SUBGRAPH", true),

    Method (Options::GEDMethod::IPFP, "BP_BEAM", true, " --ls-initialization-method BRANCH_FAST"),
    Method (Options::GEDMethod::IPFP, "BP_BEAM", true, " --ls-initialization-method BRANCH_UNIFORM"),
    Method (Options::GEDMethod::IPFP, "BP_BEAM", true, " --ls-initialization-method BIPARTITE"),
    Method (Options::GEDMethod::IPFP, "BP_BEAM", true, " --ls-initialization-method NODE"),
    Method (Options::GEDMethod::IPFP, "BP_BEAM", true, " --ls-initialization-method STAR"),

    Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method BRANCH_FAST"),
    Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method BRANCH_UNIFORM"),
    Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method BIPARTITE"),
    Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method NODE"),
    Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method STAR"),

    Method (Options::GEDMethod::IPFP, "BRANCH_UNIFORM2", true),
    Method (Options::GEDMethod::IPFP, "BRANCH_UNIFORM3", true),
    Method (Options::GEDMethod::IPFP, "BRANCH_UNIFORM4", true),

    Method (Options::GEDMethod::IPFP, "STAR2", true),
    Method (Options::GEDMethod::IPFP, "STAR3", true),
    Method (Options::GEDMethod::IPFP, "STAR4", true),
    Method (Options::GEDMethod::IPFP, "STAR6", true),
    Method (Options::GEDMethod::IPFP, "STAR7", true),
    Method (Options::GEDMethod::IPFP, "STAR8", true),

    Method (Options::GEDMethod::IPFP, "REP_TREE", true),
    Method (Options::GEDMethod::IPFP, "REP_TREE_2", true),
    Method (Options::GEDMethod::IPFP, "REP_TREE_3", true),
    Method (Options::GEDMethod::IPFP, "J48", true),
    Method (Options::GEDMethod::IPFP, "J48_2", true),
    Method (Options::GEDMethod::IPFP, "J48_3", true),
    Method (Options::GEDMethod::IPFP, "J48_4", true),
    Method (Options::GEDMethod::IPFP, "J48_5", true),
  };

  // compute_statistics()
  // test_lsape_all(lsape_methods, all_datasets);

  test_ls_all_datasets(methods, all_datasets);
  test_ls_graph_sizes(methods, all_datasets);

  // ================================================================
  // GENERATED - CLUSTER
  // ================================================================
  {
    size_t const num_graphs = 50;
    size_t const node_variance = 5;
    double const pInner = 0.8;
    double const pOuter = 0.02;

    test_ls_cluster_graphs(methods, num_graphs, node_variance, { 20, 30, 40, 50, 60, 70 }, { 2 }, pInner, pOuter);
    test_ls_cluster_graphs(methods, num_graphs, node_variance, { 20, 30, 40, 50, 60, 70 }, { 3 }, pInner, pOuter);
    test_ls_cluster_graphs(methods, num_graphs, node_variance, { 20, 30, 40, 50, 60, 70 }, { 5 }, pInner, pOuter);

    return 0;
  }

  // ================================================================
  // GENERATED - POWER
  // ================================================================
  {
    size_t const num_graphs = 50;
    size_t const node_variance = 5;
    std::vector<size_t> const graph_sizes { 10, 20, 30, 40, 50, 60 };

    test_ls_power_graphs(methods, num_graphs, node_variance, graph_sizes, { 3 });
    test_ls_power_graphs(methods, num_graphs, node_variance, graph_sizes, { 5 });
    test_ls_power_graphs(methods, num_graphs, node_variance, { 20, 30, 40, 50, 60 }, { 10 });
  }

  // ================================================================
  // GENERATED - RANDOM
  // ================================================================
  {
    size_t const num_graphs = 100;
    size_t const node_variance = 5;
    std::vector<double> const edge_densities {
      0.1, 0.2, 0.3, 0.4, 0.5, 0.6
    };

    test_ls_rand_graphs(methods, num_graphs, node_variance, { 10 }, edge_densities);
    test_ls_rand_graphs(methods, num_graphs, node_variance, { 30 }, edge_densities);
    test_ls_rand_graphs(methods, num_graphs, node_variance, { 50 }, edge_densities);
  }

  return 0;
}
