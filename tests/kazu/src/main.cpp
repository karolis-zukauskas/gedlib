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

#define WRITE_STATS_FILE
constexpr int TEST_THREADS = 8;
constexpr bool TEST_ONLY_UNIQUE_PAIRS = true;

#include "util.hpp"
#include "graph_gen.hpp"
#include "graph_statistics.hpp"
#include "tests/ls_all.hpp"
#include "tests/ls_rand.hpp"
#include "tests/ls_sizes.hpp"
#include "tests/lsape.hpp"
#include "tests/statistics.hpp"

int main(int argc, char* argv[]) {
  test_ls_all_datasets();

  // compute_statistics();

  // GxlGEDEnv env;

  // std::mt19937 rng;
  // std::vector<std::string> const node_labels { "A", "B", "C" };
  // std::vector<std::string> const edge_labels { "0", "1" };

  // for (size_t n = 10; n < 25; n += 5) {
  //   for (size_t i = 3; i < 6; ++i) {
  //     std::stringstream ss;
  //     ss << "../results/" << "power_" << std::to_string(n) << "_" << std::to_string(i) << ".gxl";

  //     std::string filename = ss.str();
  //     auto graph_id = graph_gen_power(env, rng, n, i, node_labels, edge_labels);
  //     env.save_as_gxl_graph(graph_id, filename);
  //   }
  // }

  // auto create_cluster = [&](size_t num_nodes, size_t num_clusters, double p_inner_edge, double p_outer_edge) -> void {
  //   std::stringstream ss;
  //   ss << "../results/" << "cluster_" << std::to_string(num_nodes) << "_" << std::to_string(num_clusters) << 
  //     "_" << std::to_string(p_inner_edge) << "_" << std::to_string(p_outer_edge) << ".gxl";

  //   std::string filename = ss.str();
  //   auto graph_id = graph_gen_cluster(env, rng, num_nodes, num_clusters, p_inner_edge, p_outer_edge, node_labels, edge_labels);
  //   env.save_as_gxl_graph(graph_id, filename);
  // };

  // create_cluster(20, 4, 0.8, 0);
  // create_cluster(20, 4, 0.8, 0.05);
  // create_cluster(40, 4, 0.8, 0);
  // create_cluster(40, 4, 0.5, 0);

  // create_cluster(40, 4, 0.75, 0.01);
  // create_cluster(40, 4, 0.75, 0.02);
  // create_cluster(40, 4, 0.75, 0.03);
  // create_cluster(40, 4, 0.75, 0.04);
  // create_cluster(40, 4, 0.75, 0.05);

  // auto graph_ids = ::util::setup_environment("AIDS", false, env);
  // for (auto& id : graph_ids) {
  //   GraphStatistics stats = graph_stat_compute_single(env, id);
  //   std::cout << env.get_graph_name(id) << " D: " << stats.diameter << " R: " << stats.radius << " joint: " << stats.is_joint <<
  //     " min_pow: " << stats.min_node_degree << " max_pow: " << stats.max_node_degree << std::endl;
  // }

  // test_ls_graph_sizes();
  // test_ls_all_datasets();
  // test_ls_rand_graphs();
  // test_lsape_all();

  return 0;
}
