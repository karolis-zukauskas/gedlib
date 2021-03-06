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
 * @file tests/vldbj2019/src/util.hpp
 * @brief Provides utility functions for tests of VLDB J. submission.
 * @details The binary built from this file was used for the experiments in the following paper:
 * - D. B. Blumenthal, N. Boria, J. Gamper, S. Bougleux L. Brun:
 *   &ldquo;Comparing heuristics for graph edit distance computation&rdquo;,
 *   VLDB J. 2019
 */

#ifndef SRC_TESTS_VLDBJ2018_UTIL_HPP_
#define SRC_TESTS_VLDBJ2018_UTIL_HPP_

#define GXL_GEDLIB_SHARED
#include "../../../src/env/ged_env.hpp"

using namespace ged;

#include <functional>
#include <vector>
#include <numeric>
#include <cassert>

using GxlExchangeGraph = ged::ExchangeGraph<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>;
using GxlGEDEnv = ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>;

static std::string s_stats_filename;

namespace util {

bool is_chemical_dataset(const std::string & dataset) {
  return ((dataset == "AIDS") or (dataset == "Mutagenicity") or (dataset == "acyclic") or (dataset == "alkane") or (dataset == "mao") or (dataset == "pah"));
}

bool is_letter_dataset(const std::string & dataset) {
  return ((dataset == "Letter_HIGH") or (dataset == "Letter_LOW") or (dataset == "Letter_MED"));
}

void check_dataset(const std::string & dataset) {
  if (not (is_chemical_dataset(dataset) or is_letter_dataset(dataset) or (dataset == "CMU-GED") or (dataset == "Fingerprint") or (dataset == "GREC") or (dataset == "Protein"))) {
    throw ged::Error(std::string("Dataset \"") + dataset + "\" does not exists.");
  }
}

std::string graph_dir(const std::string & dataset) {
  std::string root_dir("../../../data/datasets/");
  if ((dataset == "AIDS") or (dataset == "Fingerprint") or (dataset == "GREC") or (dataset == "Protein") or (dataset == "Mutagenicity")) {
    return (root_dir + dataset + "/data/");
  }
  else if ((dataset == "Letter_HIGH")) {
    return (root_dir + "Letter/HIGH/");
  }
  else if ((dataset == "Letter_LOW")) {
    return (root_dir + "Letter/LOW/");
  }
  else if ((dataset == "Letter_MED")) {
    return (root_dir + "Letter/MED/");
  }
  else if (dataset == "CMU-GED") {
    return (root_dir + dataset + "/CMU/");
  }
  else if ((dataset == "acyclic") or (dataset == "alkane") or (dataset == "mao") or (dataset == "pah")) {
    return (root_dir + dataset + "/");
  }
  else {
    throw ged::Error(std::string("Dataset \"") + dataset + "\" does not exists.");
  }
  return "";
}

std::string train_collection(const std::string & dataset) {
  std::string root_dir("../collections/");
  check_dataset(dataset);
  if (is_letter_dataset(dataset)) {
    return (root_dir + "Letter_train.xml");
  }
  return root_dir + dataset + "_train.xml";
}

std::string test_collection(const std::string & dataset) {
  std::string root_dir("../collections/");
  check_dataset(dataset);
  if (is_letter_dataset(dataset)) {
    return (root_dir + "Letter_test.xml");
  }
  return root_dir + dataset + "_test.xml";
}

std::string size_constrained_collection(const std::string & dataset, std::size_t max_size_div_10) {
  std::string suffix(std::to_string((max_size_div_10 * 10) - 9) + "-" + std::to_string(max_size_div_10 * 10));
  std::string root_dir("../collections/");
  check_dataset(dataset);
  return root_dir + dataset + suffix + ".xml";
}

std::string config_prefix(const std::string & dataset) {
  check_dataset(dataset);
  return std::string("../ini/" + dataset + "_");
}

std::string init_options(const std::string & dataset, const std::string & config_suffix, const std::string & data_suffix = "", bool save_train = false, bool load_train = false, std::size_t threads = 8) {
  check_dataset(dataset);
  std::string options("--threads ");
  options += std::to_string(threads) + " --save ../ini/";
  options += dataset + "_" + config_suffix + ".ini";
  if (save_train) {
    if (load_train) {
      throw ged::Error("Training data cannot be both saved and loaded.");
    }
    options += " --save-train ../ini/" + dataset + "_" + data_suffix + ".data";
  }
  if (load_train) {
    options += " --load-train ../ini/" + dataset + "_" + data_suffix + ".data";
  }
  return options;
}

std::string ground_truth_option(const std::string & dataset) {
  check_dataset(dataset);
  //if (is_letter_dataset(dataset)) {
  //  return std::string(" --ground-truth-method EXACT");
  //}
  return std::string(" --ground-truth-method IPFP");
}

ged::Options::EditCosts edit_costs(const std::string & dataset) {
  if (is_chemical_dataset(dataset)) {
    return ged::Options::EditCosts::CHEM_2;
  }
  else if (is_letter_dataset(dataset)) {
    return ged::Options::EditCosts::LETTER;
  }
  else if (dataset == "CMU-GED") {
    return ged::Options::EditCosts::CMU;
  }
  else if (dataset == "Fingerprint") {
    return ged::Options::EditCosts::FINGERPRINT;
  }
  else if (dataset == "GREC") {
    return ged::Options::EditCosts::GREC_2;
  }
  else if (dataset == "Protein") {
    return ged::Options::EditCosts::PROTEIN;
  }
  else {
    throw ged::Error(std::string("Dataset \"") + dataset + "\" does not exists.");
  }
  return ged::Options::EditCosts::CONSTANT;
}

ged::Options::GXLNodeEdgeType node_type(const std::string & dataset) {
  check_dataset(dataset);
  if ((dataset == "Fingerprint") or (dataset == "CMU-GED")) {
    return ged::Options::GXLNodeEdgeType::UNLABELED;
  }
  return ged::Options::GXLNodeEdgeType::LABELED;
}

ged::Options::GXLNodeEdgeType edge_type(const std::string & dataset) {
  check_dataset(dataset);
  if (is_letter_dataset(dataset)) {
    return ged::Options::GXLNodeEdgeType::UNLABELED;
  }
  return ged::Options::GXLNodeEdgeType::LABELED;
}

std::unordered_set<std::string> irrelevant_node_attributes(const std::string & dataset) {
  check_dataset(dataset);
  std::unordered_set<std::string> irrelevant_attributes;
  if ((dataset == "AIDS")) {
    irrelevant_attributes.insert({"x", "y", "symbol"});
  }
  else if (dataset == "Protein") {
    irrelevant_attributes.insert("aaLength");
  }
  return irrelevant_attributes;
}

std::unordered_set<std::string> irrelevant_edge_attributes(const std::string & dataset) {
  check_dataset(dataset);
  std::unordered_set<std::string> irrelevant_attributes;
  if ((dataset == "GREC")) {
    irrelevant_attributes.insert({"angle0", "angle1"});
  }
  else if (dataset == "Protein") {
    irrelevant_attributes.insert({"distance0", "distance1"});
  }
  else if (dataset == "Fingerprint") {
    irrelevant_attributes.insert("angle");
  }
  return irrelevant_attributes;
}

ged::Options::InitType init_type(const std::string & dataset) {
  if (is_chemical_dataset(dataset) or (dataset == "Protein")) {
    return ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES;
  }
  return ged::Options::InitType::LAZY_WITHOUT_SHUFFLED_COPIES;
}

std::vector<ged::GEDGraph::GraphID> setup_environment(const std::string & dataset, bool train, ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & env) {
  std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(graph_dir(dataset), (train ? train_collection(dataset) : test_collection(dataset)), node_type(dataset), edge_type(dataset), irrelevant_node_attributes(dataset), irrelevant_edge_attributes(dataset)));
  env.set_edit_costs(edit_costs(dataset));
  env.init(init_type(dataset));
  return graph_ids;
}

std::vector<ged::GEDGraph::GraphID> setup_environment(const std::string & dataset, std::size_t max_size_div_10, ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & env) {
  std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs(graph_dir(dataset), size_constrained_collection(dataset, max_size_div_10), node_type(dataset), edge_type(dataset), irrelevant_node_attributes(dataset), irrelevant_edge_attributes(dataset)));
  env.set_edit_costs(edit_costs(dataset));
  env.init(init_type(dataset));
  return graph_ids;
}

void setup_generated_environment(ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & env) {
  // Costs: [node_ins_cost, node_del_cost, node_rel_cost, edge_ins_cost, edge_del_cost, edge_rel_cost]
  std::vector<double> edit_costs { 2.5, 2.5, 5, 1.0, 1.0, 2.0 };

  env.set_edit_costs(ged::Options::EditCosts::CONSTANT, edit_costs);
  env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
}

void setup_datasets(std::vector<std::string> & datasets) {
  datasets = {"Letter_HIGH", "Mutagenicity", "AIDS", "Protein", "GREC", "Fingerprint"};
}

void setup_size_test_datasets(std::vector<std::string> & datasets) {
  datasets = {"Mutagenicity", "AIDS", "Protein"};
}

std::string create_stats_file(std::string const& base_name) {
  std::time_t t = std::time(0);
  std::tm* now = std::localtime(&t);
  std::stringstream ss;
  ss << "../results/stats" << "_" << base_name << "_" << now->tm_mon << "-" << now->tm_mday << "_"
    << now->tm_hour << "-" << now->tm_min << "-" << now->tm_sec << ".csv";

  std::string filename = ss.str();
  std::ofstream stats_file(filename.c_str(), std::ios_base::trunc);

  stats_file << "method, ub, runtime, dataset, g_id, h_id, node_count, edge_count, avg_node_degree, edge_density, node_label_count, edge_label_count, diameter, radius, both_joint" << std::endl;
  stats_file.close();

  return filename;
}

std::string create_result_file(std::string const& base_name, std::string const& dataset) {
  std::time_t t = std::time(0);
  std::tm* now = std::localtime(&t);
  std::stringstream ss;
  ss << "../results/" << base_name << "_" << dataset << "_" << now->tm_mon << "-" << now->tm_mday << "_"
    << now->tm_hour << "-" << now->tm_min << "-" << now->tm_sec << ".csv";

  std::string filename = ss.str();
  std::ofstream result_file(filename.c_str(), std::ios_base::trunc);

  result_file << "lb, ub, runtime, init_solutions_t, ls_iterations_t, ls_iterations, num_nodes, method" << std::endl;
  result_file.close();

  return filename;
}

}

#endif /* SRC_TESTS_VLDBJ2018_UTIL_HPP_ */
