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

#include "util.hpp"

using GxlGEDEnv = ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>;

class Method {
private:
  ged::Options::GEDMethod ged_method_;
  std::string initialization_method_;
  bool use_optimal_lsape_initialization_;
  std::string extra_options_;

  std::string options_() const {
    // TODO: Do we want to test multithreaded?
    std::string options("--threads 1");

    options += " --optimal_initialization " + use_optimal_lsape_initialization_ ? "TRUE" : "FALSE";
    options += " --initialization-method " + initialization_method_;

    return options + extra_options_;
  }

public:
  Method(ged::Options::GEDMethod ged_method, std::string const& initialization_method, bool use_optimal_lsape_initialization, std::string const& extra_options = "")
    : ged_method_{ged_method},
      initialization_method_{initialization_method},
      use_optimal_lsape_initialization_{use_optimal_lsape_initialization},
      extra_options_(extra_options) {}

    std::string name() const {
      std::stringstream name;

      if (ged_method_ == ged::Options::GEDMethod::IPFP) {
        name << "IPFP";
      } else if (ged_method_ == ged::Options::GEDMethod::REFINE) {
        name << "REFINE";
      } else {
        assert(false);
      }

      name << " (" << initialization_method_ << ")";

      return name.str();
    }

    void run_on_dataset(const std::string & dataset, ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> & env, bool ensure_n_greater_m,
      double & avg_lb, double & avg_ub, double & avg_runtime, double& avg_init_time, double& avg_ls_iterations_time, double& avg_num_ls_iterations) const {

      env.set_method(ged_method_, options_());
      if (dataset != "Protein" or ged_method_ != ged::Options::GEDMethod::PARTITION) {
        env.init_method();
      }

      std::size_t num_runs{env.graph_ids().second * env.graph_ids().second};
      ged::ProgressBar progress_bar(num_runs);

      avg_runtime = 0;
      avg_ub = 0;
      avg_lb = 0;
      avg_init_time = 0;
      avg_ls_iterations_time = 0;
      avg_num_ls_iterations = 0;

      std::cout << "\n\tDataset: " << dataset << "   Method: " << name() << std::endl;

      for (ged::GEDGraph::GraphID g_id = env.graph_ids().first; g_id != env.graph_ids().second; g_id++) {
        for (ged::GEDGraph::GraphID h_id = env.graph_ids().first; h_id != env.graph_ids().second; h_id++) {
          ged::GEDGraph::GraphID _g_id = g_id;
          ged::GEDGraph::GraphID _h_id = h_id;

          // Some greedy LSAPE solvers require graph nodes to be n >= m
          if (ensure_n_greater_m && env.get_num_nodes(_g_id) < env.get_num_nodes(_h_id))
            std::swap(_g_id, _h_id);

          env.run_method(_g_id, _h_id);

          avg_lb += env.get_lower_bound(_g_id, _h_id);
          avg_ub += env.get_upper_bound(_g_id, _h_id);
          avg_runtime += env.get_runtime(_g_id, _h_id);

          auto pLSMethod = dynamic_cast<ged::LSBasedMethod<ged::GXLLabel, ged::GXLLabel>*>(env.get_method());
          assert(pLSMethod != nullptr);

          avg_init_time += pLSMethod->initial_sulutions_time.count();
          avg_ls_iterations_time += pLSMethod->ls_iterations_time.count();
          avg_num_ls_iterations += static_cast<double>(pLSMethod->num_ls_iterations);

          progress_bar.increment();
          std::cout << "\r\t" << name() << ": " << progress_bar << std::flush;
        }
      }

      avg_lb /= static_cast<double>(num_runs);
      avg_ub /= static_cast<double>(num_runs);
      avg_runtime /= static_cast<double>(num_runs);
      avg_init_time /= static_cast<double>(num_runs);
      avg_ls_iterations_time /= static_cast<double>(num_runs);
      avg_num_ls_iterations /= static_cast<double>(num_runs);

      std::cout << "\n";
    }
};

std::string create_result_file(std::string const& base_name, std::string const& dataset) {
  std::time_t t = std::time(0);
  std::tm* now = std::localtime(&t);
  std::stringstream ss;
  ss << "../results/" << base_name << "_" << dataset << "_" << now->tm_mon << "-" << now->tm_mday << "_"
    << now->tm_hour << "-" << now->tm_min << "-" << now->tm_sec << ".csv";

  std::string filename = ss.str();
  std::ofstream result_file(filename.c_str());

  result_file << "lb, ub, runtime, init_solutions_t, ls_iterations_t, ls_iterations, num_nodes, method" << std::endl;
  result_file.close();

  return filename;
}

void run_methods(std::vector<Method> const& methods, GxlGEDEnv& env, std::string const& dataset, bool ensure_n_greater_m, std::string const& result_filename) {
  double avg_ub{0};
  double avg_lb{0};
  double avg_runtime{0};
  double avg_init_time{0};
  double avg_ls_iterations_time{0};
  double avg_num_ls_iterations = 0;

  for (auto & method : methods) {
    method.run_on_dataset(dataset, env, ensure_n_greater_m, avg_lb, avg_ub, avg_runtime, avg_init_time, avg_ls_iterations_time, avg_num_ls_iterations);

    std::ofstream result_file(result_filename.c_str(), std::ios_base::app);

    result_file << avg_lb << "," << avg_ub << "," << avg_runtime << "," << avg_init_time << "," << avg_ls_iterations_time << ",";
    result_file << avg_num_ls_iterations << "," << env.get_avg_num_nodes() << "," << method.name() << "\n";
    result_file.close();
  }
}

void run_on_sized_dataset(std::vector<Method> const& methods, std::string const& dataset, bool ensure_n_greater_m, std::string const& result_filename) {
  std::cout << "\n\n=== " << dataset << " ===\n\n";

  std::size_t max_max_size_div_10 {0};
  if (dataset == "AIDS") {
    max_max_size_div_10 = 8;
  }
  else if (dataset == "Protein") {
    max_max_size_div_10 = 6;
  }
  else if (dataset == "Mutagenicity") {
    max_max_size_div_10 = 10;
  }

  for (std::size_t max_size_dev_10{1}; max_size_dev_10 <= max_max_size_div_10; max_size_dev_10++) {
    std::cout << "\n@@@@@\n@@@@@ DATASET PARTITION: " << max_size_dev_10 << "\n@@@@@\n";

    ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
    util::setup_environment(dataset, max_size_dev_10, env);

    run_methods(methods, env, dataset, ensure_n_greater_m, result_filename);
  }
}

void run_on_test_dataset(std::vector<Method> const& methods, std::string const& dataset, bool ensure_n_greater_m, std::string const& result_filename) {
  std::cout << "\n\n=== " << dataset << " ===\n\n";

  ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
  util::setup_environment(dataset, false, env);

  run_methods(methods, env, dataset, ensure_n_greater_m, result_filename);
}

void test_ls_all_datasets() {
  std::vector<std::string> datasets = {
    "Letter_HIGH", "Mutagenicity", "AIDS", "Protein", "GREC", "Fingerprint",
  };

  std::vector<Method> methods {
    Method (ged::Options::GEDMethod::IPFP, "STAR", true),
    Method (ged::Options::GEDMethod::IPFP, "BRANCH_UNIFORM", true),
    Method (ged::Options::GEDMethod::IPFP, "BRANCH_FAST", true),
    Method (ged::Options::GEDMethod::IPFP, "BRANCH", true),
    Method (ged::Options::GEDMethod::IPFP, "NODE", true),
    Method (ged::Options::GEDMethod::IPFP, "BIPARTITE", true),
    Method (ged::Options::GEDMethod::IPFP, "RANDOM", true),

    Method (ged::Options::GEDMethod::REFINE, "STAR", true),
    Method (ged::Options::GEDMethod::REFINE, "BRANCH_UNIFORM", true),
    Method (ged::Options::GEDMethod::REFINE, "BRANCH_FAST", true),
    Method (ged::Options::GEDMethod::REFINE, "BRANCH", true),
    Method (ged::Options::GEDMethod::REFINE, "NODE", true),
    Method (ged::Options::GEDMethod::REFINE, "BIPARTITE", true),
    Method (ged::Options::GEDMethod::REFINE, "RANDOM", true),
  };

  for (auto dataset : datasets) {
    try {
      std::string results_filename = create_result_file("ls_all_datasets", dataset);

      run_on_test_dataset(methods, dataset, false, results_filename);
    }
    catch (const std::exception & error) {
      std::cerr << error.what() << ". " << "Error on test_ls_all_datasets: " << dataset << ".\n";
    }
  }
}

void test_ls_graph_sizes() {
  std::vector<std::string> datasets = {
    "Mutagenicity", "AIDS", "Protein",
  };

  std::vector<Method> methods {
    Method (ged::Options::GEDMethod::IPFP, "STAR", true),
    Method (ged::Options::GEDMethod::IPFP, "BIPARTITE", true),
    Method (ged::Options::GEDMethod::IPFP, "BRANCH_UNIFORM", true),
    Method (ged::Options::GEDMethod::IPFP, "BRANCH_FAST", true),
    Method (ged::Options::GEDMethod::IPFP, "BRANCH", true),
    Method (ged::Options::GEDMethod::IPFP, "NODE", true),
    Method (ged::Options::GEDMethod::IPFP, "RANDOM", true),
  };

  for (auto dataset : datasets) {
    try {
      std::string results_filename = create_result_file("ls_graph_sizes", dataset);

      run_on_sized_dataset(methods, dataset, false, results_filename);
    }
    catch (const std::exception & error) {
      std::cerr << error.what() << ". " << "Error on test_ls_graph_sizes: " << dataset << ".\n";
    }
  }
}

void test_greedy_lsape() {

}

void test_iterations() {

}

int main(int argc, char* argv[]) {
  test_ls_all_datasets();
  test_ls_graph_sizes();

  return 0;
}

    //{ ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --iterations 100 --initialization-method BRANCH_UNIFORM" },
    //{ ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --iterations 50 --initialization-method BRANCH_UNIFORM" },
    //{ ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --iterations 33 --initialization-method BRANCH_UNIFORM" },
    //{ ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --iterations 20 --initialization-method BRANCH_UNIFORM" },
    //{ ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --iterations 10 --initialization-method BRANCH_UNIFORM" },

    //{ ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1" },
    //{ ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --initialization-method BRANCH_UNIFORM" },
    //{ ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --initialization-method BRANCH_FAST" },
    //{ ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --initialization-method BRANCH" },
    //{ ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --initialization-method NODE" },

    //{ ged::Options::GEDMethod::REFINE, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --max-swap-size 2" },
    //{ ged::Options::GEDMethod::REFINE, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --max-swap-size 2 --initialization-method BRANCH_UNIFORM" },
    //{ ged::Options::GEDMethod::REFINE, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --max-swap-size 2 --initialization-method BRANCH_FAST" },
    //{ ged::Options::GEDMethod::REFINE, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --max-swap-size 2 --initialization-method BRANCH" },
    //{ ged::Options::GEDMethod::REFINE, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --max-swap-size 2 --initialization-method NODE" },

    //{ ged::Options::GEDMethod::REFINE, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --max-swap-size 3" },
    //{ ged::Options::GEDMethod::REFINE, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --max-swap-size 3 --initialization-method BRANCH_UNIFORM" },
    //{ ged::Options::GEDMethod::REFINE, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --max-swap-size 3 --initialization-method BRANCH_FAST" },
    //{ ged::Options::GEDMethod::REFINE, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --max-swap-size 3 --initialization-method BRANCH" },
    //{ ged::Options::GEDMethod::REFINE, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --max-swap-size 3 --initialization-method NODE" },

    // { ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --initialization-method SUBGRAPH" },
    // { ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --initialization-method WALKS" },

    //{ ged::Options::GEDMethod::IPFP, " --initialization-method BIPARTITE_ML" },
    //{ ged::Options::GEDMethod::IPFP, " --initialization-method BIPARTITE" },
    //{ ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --initialization-method RING_ML" },
    //{ ged::Options::GEDMethod::IPFP, " --initial-solutions 1 --ratio-runs-from-initial-solutions 1 --initialization-method RING" },

    // { ged::Options::GEDMethod::NODE, "" },
    // { ged::Options::GEDMethod::BRANCH_FAST, "" },
    // { ged::Options::GEDMethod::BRANCH_UNIFORM, "" },
    // { ged::Options::GEDMethod::RING, "" }, // AI based - really slow training ..
    // { ged::Options::GEDMethod::STAR, "" },
    // { ged::Options::GEDMethod::BLP_NO_EDGE_LABELS, "" },
