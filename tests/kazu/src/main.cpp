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
 * @file tests/vldbj2019/src/test_best_methods.cpp
 * @brief Runs tests for VLDB J. submission.
 * @details The binary built from this file was used for the experiments in the following paper:
 * - D. B. Blumenthal, N. Boria, J. Gamper, S. Bougleux L. Brun:
 *   &ldquo;Comparing heuristics for graph edit distance computation&rdquo;,
 *   VLDB J. 2019
 */

//#define WRITE_STATS_FILE
constexpr int TEST_THREADS = 2;
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
    Method (Options::GEDMethod::IPFP, "BRANCH_UNIFORM", true),
    Method (Options::GEDMethod::IPFP, "BRANCH_FAST", true),
    Method (Options::GEDMethod::IPFP, "BRANCH", true),
    Method (Options::GEDMethod::IPFP, "NODE", true),
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
    Method (Options::GEDMethod::IPFP, "STAR", true),
    Method (Options::GEDMethod::IPFP, "STAR4", true),
    Method (Options::GEDMethod::IPFP, "STAR6", true),
  };

  // compute_statistics();

  // {
  //   GxlGEDEnv env;
  //   std::mt19937 rng;
  //   std::vector<std::string> const node_labels { "A", "B", "C" };
  //   std::vector<std::string> const edge_labels { "0", "1" };

  //   auto id1 = graph_gen_cluster(env, rng, 19, 4, 1.0, 0.00, node_labels, edge_labels);
  //   env.save_as_gxl_graph(id1, "../results/_test_graph1.gxl");

  //   auto id2 = graph_gen_cluster(env, rng, 43, 4, 0.8, 0.00, node_labels, edge_labels);
  //   env.save_as_gxl_graph(id2, "../results/_test_graph2.gxl");
  // }

  {
    // GxlGEDEnv env;
    // ::util::setup_environment("Fingerprint", false, env);

    // GraphStatsMap stats;
    // graph_stat_compute_all(env, stats);

    // auto const g = env.get_graph(0);
    // auto const h = env.get_graph(1);

    // for (auto it = env.graph_ids().first; it != env.graph_ids().second - 1; it++) {
    //   GraphDiff diff_1 = graph_diff_compute(env, stats, it, it + 1);
    //   GraphDiff diff_2 = graph_diff_compute_temp(env, it, it + 1);

    //   assert(diff_1.node_count == diff_2.node_count);
    //   assert(diff_1.edge_count == diff_2.edge_count);
    //   assert(diff_1.avg_node_degree == diff_2.avg_node_degree);
    //   assert(diff_1.edge_density == diff_2.edge_density);
    //   assert(diff_1.node_label_count == diff_2.node_label_count);
    //   assert(diff_1.edge_label_count == diff_2.edge_label_count);
    // }
  }

  {
    std::vector<std::string> const all_datasets = {
      "Mutagenicity", "AIDS"
    };
    std::vector<Method> const methods {
      Method (Options::GEDMethod::IPFP, "J48_4", true),
      Method (Options::GEDMethod::IPFP, "REP_TREE_2", true),
    };

    test_ls_all_datasets(methods, all_datasets);
    return 0;
  }

  // ================================================================
  // GENERATED - CLUSTER
  // ================================================================

  double p_inner_edge = 0.75;
  double p_outer_edge = 0.02;
  {
    size_t const num_graphs = 50;
    size_t const node_variance = 5;
    std::vector<size_t> const graph_sizes {
      20, 30,
    };
    std::vector<size_t> const cluster_count {
      3, 5,
    };
    test_ls_cluster_graphs(methods, num_graphs, node_variance, graph_sizes, cluster_count, p_inner_edge, p_outer_edge);
  }

  {
    size_t const num_graphs = 50;
    size_t const node_variance = 20;
    std::vector<size_t> const graph_sizes {
      30,
    };
    std::vector<size_t> const cluster_count {
      5
    };
    test_ls_cluster_graphs(methods, num_graphs, node_variance, graph_sizes, cluster_count, p_inner_edge, p_outer_edge);
  }

  // ================================================================
  // GENERATED - POWER
  // ================================================================

  {
    size_t const num_graphs = 50;
    size_t const node_variance = 5;
    std::vector<size_t> const graph_sizes {
      // 10, 20, 30, 40
      10, 30,
    };
    std::vector<size_t> const edges_per_node {
      // 2, 3, 4, 5,
      2, 3,
    };
    test_ls_power_graphs(methods, num_graphs, node_variance, graph_sizes, edges_per_node);
  }

  // ================================================================
  // GENERATED - RANDOM
  // ================================================================

  {
    size_t const num_graphs = 50;
    size_t const node_variance = 5;
    std::vector<size_t> const graph_sizes {
      10, 30,
    };
    std::vector<double> const edge_densities {
      0.1, 0.2,
    };
    test_ls_rand_graphs(methods, num_graphs, node_variance, graph_sizes, edge_densities);
  }

  // test_ls_graph_sizes();
  // test_ls_rand_graphs();
  // test_lsape_all();

  return 0;
}
