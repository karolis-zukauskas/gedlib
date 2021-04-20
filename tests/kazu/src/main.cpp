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
#include "graph_gen.hpp"
#include "graph_statistics.hpp"
#include <functional>

using namespace ged;
using GxlExchangeGraph = ExchangeGraph<GXLNodeID, GXLLabel, GXLLabel>;
using GxlGEDEnv = GEDEnv<GXLNodeID, GXLLabel, GXLLabel>;

constexpr int TEST_THREADS = 1;

class Method {
private:
  Options::GEDMethod ged_method_;
  std::string initialization_method_;
  bool use_optimal_lsape_initialization_;
  std::string extra_options_;

  std::string options_(const std::string & dataset) const {
    // TODO: Do we want to run algorithms multithreaded?
    std::string options("--threads 1");

    if (ged_method_ == Options::GEDMethod::REFINE || ged_method_ == Options::GEDMethod::IPFP) {
      //options += std::string(" --optimal_initialization ") + (use_optimal_lsape_initialization_ ? "TRUE" : "FALSE");
      options += " --initialization-method " + initialization_method_;

      if (initialization_method_ == "WALKS") {
        options += " --load ../ini/" + dataset + "_walks.ini";
      } else if (initialization_method_ == "SUBGRAPH") {
        options += " --load ../ini/" + dataset + "_walks.ini";
      }
    }

    if (ged_method_ == Options::GEDMethod::WALKS) {
      options += " --load ../ini/" + dataset + "_walks.ini";
    }
    if (ged_method_ == Options::GEDMethod::SUBGRAPH) {
      options += " --load ../ini/" + dataset + "_subgraph.ini";
    }

    return options + extra_options_;
  }

public:
  Method(Options::GEDMethod ged_method, std::string const& initialization_method, bool use_optimal_lsape_initialization, std::string const& extra_options = "")
    : ged_method_{ged_method},
      initialization_method_{initialization_method},
      use_optimal_lsape_initialization_{use_optimal_lsape_initialization},
      extra_options_(extra_options) {}

    std::string name() const {
      std::stringstream name;

      if (ged_method_ == Options::GEDMethod::IPFP) {
        name << "IPFP" << " (" << initialization_method_ << ")";

        // Could format this better ..
        if (!extra_options_.empty())
          name << " (" << extra_options_ << ")";
      } else if (ged_method_ == Options::GEDMethod::REFINE) {
        name << "REFINE" << " (" << initialization_method_ << ")";
      } else if (ged_method_ == Options::GEDMethod::BIPARTITE) {
        name << "BP";
      } else if (ged_method_ == Options::GEDMethod::BP_BEAM) {
        name << "BP_BEAM";
      } else if (ged_method_ == Options::GEDMethod::BRANCH_FAST) {
        name << "BRANCHFAST";
      } else if (ged_method_ == Options::GEDMethod::BRANCH) {
        name << "BRANCH";
      } else if (ged_method_ == Options::GEDMethod::BRANCH_UNIFORM) {
        name << "BRANCH_UNIFORM";
      } else if (ged_method_ == Options::GEDMethod::BRANCH_UNIFORM2) {
        name << "BRANCH_UNIFORM2";
      } else if (ged_method_ == Options::GEDMethod::STAR) {
        name << "STAR";
      } else if (ged_method_ == Options::GEDMethod::STAR2) {
        name << "STAR2";
      } else if (ged_method_ == Options::GEDMethod::STAR3) {
        name << "STAR3";
      } else if (ged_method_ == Options::GEDMethod::STAR4) {
        name << "STAR4";
      } else if (ged_method_ == Options::GEDMethod::STAR5) {
        name << "STAR5";
      } else if (ged_method_ == Options::GEDMethod::STAR6) {
        name << "STAR6";
      } else if (ged_method_ == Options::GEDMethod::SUBGRAPH) {
        name << "SUBGRAPH";
      } else if (ged_method_ == Options::GEDMethod::WALKS) {
        name << "WALKS";
      } else if (ged_method_ == Options::GEDMethod::NODE) {
        name << "NODE";
      }

      return name.str();
    }

    void run_on_dataset(const std::string & dataset, GxlGEDEnv & env, bool ensure_n_greater_m,
      double & avg_lb, double & avg_ub, double & avg_runtime, double& avg_init_time, double& avg_ls_iterations_time, double& avg_num_ls_iterations) const {

      env.set_method(ged_method_, options_(dataset));
      env.init_method();

      std::size_t num_runs{env.graph_ids().second * env.graph_ids().second};
      ProgressBar progress_bar(num_runs);

      avg_runtime = 0;
      avg_ub = 0;
      avg_lb = 0;
      avg_init_time = 0;
      avg_ls_iterations_time = 0;
      avg_num_ls_iterations = 0;

      // std::cout << "\n\tDataset: " << dataset << "   Method: " << name() << std::endl;

      for (GEDGraph::GraphID g_id = env.graph_ids().first; g_id != env.graph_ids().second; g_id++) {
        for (GEDGraph::GraphID h_id = env.graph_ids().first; h_id != env.graph_ids().second; h_id++) {
          GEDGraph::GraphID _g_id = g_id;
          GEDGraph::GraphID _h_id = h_id;

          // Some greedy LSAPE solvers require graph nodes to be n >= m
          if (ensure_n_greater_m && env.get_num_nodes(_g_id) < env.get_num_nodes(_h_id))
            std::swap(_g_id, _h_id);

          env.run_method(_g_id, _h_id);

          avg_lb += env.get_lower_bound(_g_id, _h_id);
          avg_ub += env.get_upper_bound(_g_id, _h_id);
          avg_runtime += env.get_runtime(_g_id, _h_id);

          auto pLSMethod = dynamic_cast<LSBasedMethod<GXLLabel, GXLLabel>*>(env.get_method());
          if (pLSMethod != nullptr) {
            avg_init_time += pLSMethod->initial_sulutions_time.count();
            avg_ls_iterations_time += pLSMethod->ls_iterations_time.count();
            avg_num_ls_iterations += static_cast<double>(pLSMethod->num_ls_iterations);
          }

          progress_bar.increment();
          // std::cout << "\r\t" << name() << ": " << progress_bar << std::flush;
        }
      }

      avg_lb /= static_cast<double>(num_runs);
      avg_ub /= static_cast<double>(num_runs);
      avg_runtime /= static_cast<double>(num_runs);
      avg_init_time /= static_cast<double>(num_runs);
      avg_ls_iterations_time /= static_cast<double>(num_runs);
      avg_num_ls_iterations /= static_cast<double>(num_runs);

      // std::cout << "\n";
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

  result_file << "method, num_nodes, edge_density, lb, ub, runtime, ls_iterations, init_solutions_t, ls_iterations_t" << std::endl;
  result_file.close();

  return filename;
}

void run_methods(std::vector<Method> const& methods, std::function<GxlGEDEnv()> const& setup_ged_env, std::string const& dataset,
  double edge_density, bool ensure_n_greater_m, std::string const& result_filename) {

#ifdef _OPENMP
    omp_set_num_threads(TEST_THREADS);
#pragma omp parallel for if(TEST_THREADS > 1) schedule(dynamic)
#endif
  for (std::size_t i = 0; i < methods.size(); ++i) {
    double avg_ub{0};
    double avg_lb{0};
    double avg_runtime{0};
    double avg_init_time{0};
    double avg_ls_iterations_time{0};
    double avg_num_ls_iterations = 0;

    // GEDEnv is tightly coupled with it's graph_data and edit_costs data structures - can't reuse that memory without making a lot of
    // adjustments to the whole gedlib. Initialize a copy of GEDEnv (includes graph and edit cost data) for each run of a method (needed for parallelization).
    GxlGEDEnv env = setup_ged_env();

    methods[i].run_on_dataset(dataset, env, ensure_n_greater_m, avg_lb, avg_ub, avg_runtime, avg_init_time, avg_ls_iterations_time, avg_num_ls_iterations);

#pragma omp critical
    {
      std::ofstream result_file(result_filename.c_str(), std::ios_base::app);

      result_file << methods[i].name() << "," << env.get_avg_num_nodes() << "," << edge_density << ",";
      result_file << avg_lb << "," << avg_ub << "," << avg_runtime << "," << avg_num_ls_iterations << ",";
      result_file << avg_init_time << "," << avg_ls_iterations_time << "\n";
      result_file.close();
    }
  }
}

void run_on_sized_dataset(std::vector<Method> const& methods, std::string const& dataset, bool ensure_n_greater_m, std::string const& result_filename) {
  std::cout << "\n=== " << dataset << " ===\n";

  std::size_t max_max_size_div_10 {0};
  if (dataset == "AIDS") {
    max_max_size_div_10 = 8;
  }
  else if (dataset == "Protein") {
    // NOTE: originally Protein is divided into 6 datasets
    max_max_size_div_10 = 7;
  }
  else if (dataset == "Mutagenicity") {
    max_max_size_div_10 = 10;
  }

  for (std::size_t max_size_dev_10{1}; max_size_dev_10 <= max_max_size_div_10; max_size_dev_10++) {
    std::cout << "\n@@@ Dataset Partition: " << max_size_dev_10 << "\n";

    auto setup_ged_env = [&dataset, &max_size_dev_10]() -> GxlGEDEnv {
      GxlGEDEnv env;
      ::util::setup_environment(dataset, max_size_dev_10, env);

      return env;
    };

    run_methods(methods, setup_ged_env, dataset, 0.0, ensure_n_greater_m, result_filename);
  }
}

void run_on_test_dataset(std::vector<Method> const& methods, std::string const& dataset, bool ensure_n_greater_m, std::string const& result_filename) {
  std::cout << "\n=== " << dataset << " ===\n";

  auto setup_ged_env = [&dataset]() -> GxlGEDEnv {
    GxlGEDEnv env;
    ::util::setup_environment(dataset, false, env);

    return env;
  };

  run_methods(methods, setup_ged_env, dataset, 0.0, ensure_n_greater_m, result_filename);
}

void run_on_generated_dataset(std::vector<Method> const& methods, std::function<void(GxlGEDEnv&)> const& generate_graphs, std::string const& dataset,
  double edge_density, bool ensure_n_greater_m, std::string const& result_filename) {
  std::cout << "\n=== " << dataset << " ===\n";

  auto setup_ged_env = [&generate_graphs]() -> GxlGEDEnv {
    GxlGEDEnv env;
    generate_graphs(env);

    return env;
  };

  run_methods(methods, setup_ged_env, dataset, edge_density, ensure_n_greater_m, result_filename);
}

void test_ls_all_datasets() {
  std::vector<std::string> datasets = {
    "AIDS", //"Letter_HIGH", "Mutagenicity", "Protein", "GREC", "Fingerprint",
  };

  std::vector<Method> const methods {
    // Method (Options::GEDMethod::IPFP, "BRANCH_UNIFORM", true),
    // Method (Options::GEDMethod::IPFP, "BRANCH_FAST", true),
    // Method (Options::GEDMethod::IPFP, "BRANCH", true),
    Method (Options::GEDMethod::IPFP, "BP_BEAM", true, " --ls-initialization-method NODE"),
    // Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method BRANCH_FAST"),
    // Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method BIPARTITE"),
    // Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method NODE"),
    // Method (Options::GEDMethod::IPFP, "WALKS", true),
    // Method (Options::GEDMethod::IPFP, "SUBGRAPH", true),
    // Method (Options::GEDMethod::IPFP, "STAR4", true),
    // Method (Options::GEDMethod::IPFP, "STAR5", true),
    // Method (Options::GEDMethod::IPFP, "STAR6", true),
    // Method (Options::GEDMethod::IPFP, "RANDOM", true),
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
    "AIDS", "Protein", "Mutagenicity",
  };

  std::vector<Method> const methods {
    Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method BRANCH_FAST"),
    Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method BIPARTITE"),
    Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method NODE"),
    // Method (Options::GEDMethod::IPFP, "WALKS", true),
    // Method (Options::GEDMethod::IPFP, "SUBGRAPH", true),
    // Method (Options::GEDMethod::IPFP, "STAR4", true),
    // Method (Options::GEDMethod::IPFP, "STAR4", true),
    // Method (Options::GEDMethod::IPFP, "STAR5", true),
    // Method (Options::GEDMethod::IPFP, "STAR6", true),
    // Method (Options::GEDMethod::IPFP, "BIPARTITE", true),
    // Method (Options::GEDMethod::IPFP, "BRANCH_FAST", true),
    // Method (Options::GEDMethod::IPFP, "BRANCH", true),
    // Method (Options::GEDMethod::IPFP, "NODE", true),
    // Method (Options::GEDMethod::IPFP, "RANDOM", true),
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

void test_ls_rand_graphs() {
  std::vector<Method> const methods {
    Method (Options::GEDMethod::IPFP, "BRANCH_UNIFORM", true),
    Method (Options::GEDMethod::IPFP, "BRANCH_FAST", true),
    Method (Options::GEDMethod::IPFP, "BRANCH", true),
    Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method BRANCH_FAST"),
    //Method (Options::GEDMethod::IPFP, "STAR5", true),
    //Method (Options::GEDMethod::IPFP, "STAR6", true),
  };

  size_t const num_graphs = 100;
  size_t const node_variance = 5;
  std::vector<size_t> const graph_sizes {
    10, 20, 30, 40, 50,
  };
  std::vector<double> const edge_densities {
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
  };
  std::vector<std::string> const node_labels { "A", "B", "C" };
  std::vector<std::string> const edge_labels { "0", "1" };

  std::string results_filename = create_result_file("ls_generated", "rand");
  for (size_t graph_size : graph_sizes) {
    for (double edge_density : edge_densities) {
      std::stringstream ss;
      ss << "random, nodes: " << graph_size << " edge density: " << edge_density;
      std::string dataset = ss.str();

      try {
        auto generate_graphs = [num_graphs, graph_size, edge_density, &node_labels, &edge_labels, node_variance](GxlGEDEnv& env) -> void {
          std::mt19937 size_variance_rng(0);          // For graph_size variance
          std::mt19937 graph_gen_rng(42);             // For generating graphs

          for (size_t i = 0; i < num_graphs; ++i) {
            double multiplier = static_cast<double>((static_cast<long>(size_variance_rng() % 2000)) - 1000) / 1000.0;
            size_t num_nodes = graph_size + static_cast<size_t>(std::round(multiplier * static_cast<double>(node_variance)));

            graph_gen_random(env, graph_gen_rng, num_nodes, edge_density, node_labels, edge_labels);
          }

          ::util::setup_rand_environment(env);
        };

        run_on_generated_dataset(methods, generate_graphs, dataset, edge_density, false, results_filename);
      } catch (const std::exception & error) {
        std::cerr << error.what() << ". " << "Error on test_ls_rand_graphs: " << dataset << ".\n";
      }
    }
  }
}

void test_lsape_all() {
  std::vector<std::string> all_datasets = {
    "Letter_HIGH", "Mutagenicity", "AIDS", "Protein", "GREC", "Fingerprint",
  };
  std::vector<std::string> sized_datasets = {
    "Mutagenicity", "AIDS", "Protein",
  };

  std::vector<Method> const methods {
    Method (Options::GEDMethod::WALKS, "", true),
    Method (Options::GEDMethod::SUBGRAPH, "", true),
    Method (Options::GEDMethod::BIPARTITE, "", true),
    Method (Options::GEDMethod::STAR, "", true),
    // Method (Options::GEDMethod::STAR2, "", true),
    // Method (Options::GEDMethod::STAR3, "", true),
    Method (Options::GEDMethod::STAR4, "", true),
    // Method (Options::GEDMethod::STAR5, "", true),
    Method (Options::GEDMethod::STAR6, "", true),
    // Method (Options::GEDMethod::BRANCH_UNIFORM2, "", true),
    Method (Options::GEDMethod::NODE, "", true),
    Method (Options::GEDMethod::BRANCH, "", true),
    Method (Options::GEDMethod::BRANCH_FAST, "", true),
    Method (Options::GEDMethod::BRANCH_UNIFORM, "", true),
  };

  for (auto dataset : all_datasets) {
    try {
      std::string results_filename = create_result_file("lsape_all_datasets", dataset);

      run_on_test_dataset(methods, dataset, false, results_filename);
    }
    catch (const std::exception & error) {
      std::cerr << error.what() << ". " << "Error on test_lsape_all: " << dataset << ".\n";
    }
  }

  // TEMP
  return;

  for (auto dataset : sized_datasets) {
    try {
      std::string results_filename = create_result_file("lsape_sized_datasets", dataset);

      run_on_sized_dataset(methods, dataset, false, results_filename);
    }
    catch (const std::exception & error) {
      std::cerr << error.what() << ". " << "Error on test_lsape_all: " << dataset << ".\n";
    }
  }

  size_t const num_graphs = 100;
  size_t const node_variance = 5;
  std::vector<size_t> const graph_sizes {
    10, 20, 30, 40, 50,
  };
  std::vector<double> const edge_densities {
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
  };
  std::vector<std::string> const node_labels { "A", "B", "C" };
  std::vector<std::string> const edge_labels { "0", "1" };

  std::string results_filename = create_result_file("lsape_generated", "rand");
  for (size_t graph_size : graph_sizes) {
    for (double edge_density : edge_densities) {
      std::stringstream ss;
      ss << "random, nodes: " << graph_size << " edge density: " << edge_density;
      std::string dataset = ss.str();

      try {
        auto generate_graphs = [num_graphs, graph_size, edge_density, &node_labels, &edge_labels, node_variance](GxlGEDEnv& env) -> void {
          std::mt19937 size_variance_rng(graph_size);         // For graph_size variance
          std::mt19937 graph_gen_rng(graph_size);             // For generating graphs

          for (size_t i = 0; i < num_graphs; ++i) {
            double multiplier = static_cast<double>((static_cast<long>(size_variance_rng() % 2000)) - 1000) / 1000.0;
            size_t num_nodes = graph_size + static_cast<size_t>(std::round(multiplier * static_cast<double>(node_variance)));

            graph_gen_random(env, graph_gen_rng, num_nodes, edge_density, node_labels, edge_labels);
          }

          ::util::setup_rand_environment(env);
        };

        run_on_generated_dataset(methods, generate_graphs, dataset, edge_density, false, results_filename);
      } catch (const std::exception & error) {
        std::cerr << error.what() << ". " << "Error on test_lsape_all: " << dataset << ".\n";
      }
    }
  }
}

void test_greedy_lsape() {

}

void test_iterations() {

}

void compute_statistics() {
  std::vector<std::string> all_datasets = {
    "Letter_HIGH", "Mutagenicity", "AIDS", "Protein", "GREC", "Fingerprint",
  };
  std::vector<std::string> sized_datasets = {
    "Mutagenicity", "AIDS", "Protein",
  };
  std::vector<std::string> const node_labels { "A", "B", "C" };
  std::vector<std::string> const edge_labels { "0", "1" };

  auto print_statistics = [](GxlGEDEnv& env, std::string const& dataset) -> void {
    double is_joint = 0.0;
    double num_nodes = 0.0;
    double num_edges = 0.0;
    double avg_degree = 0.0;
    double diameter = 0.0;
    double radius = 0.0;
    double min_node_degree = 0.0;
    double max_node_degree = 0.0;

    for (GEDGraph::GraphID i = env.graph_ids().first; i != env.graph_ids().second; ++i) {
      GraphStatistics stats = graph_stat_compute(env, i);

      is_joint += static_cast<double>(stats.is_joint);
      num_nodes += static_cast<double>(stats.num_nodes);
      num_edges += static_cast<double>(stats.num_edges);
      avg_degree += stats.avg_degree;
      diameter += static_cast<double>(stats.diameter);
      radius += static_cast<double>(stats.radius);
      min_node_degree += static_cast<double>(stats.min_node_degree);
      max_node_degree += static_cast<double>(stats.max_node_degree);
    }

    is_joint /= static_cast<double>(env.num_graphs());
    num_nodes /= static_cast<double>(env.num_graphs());
    num_edges /= static_cast<double>(env.num_graphs());
    avg_degree /= static_cast<double>(env.num_graphs());
    diameter /= static_cast<double>(env.num_graphs());
    radius /= static_cast<double>(env.num_graphs());
    min_node_degree /= static_cast<double>(env.num_graphs());
    max_node_degree /= static_cast<double>(env.num_graphs());


    std::cout << std::setprecision(3) << dataset << "\tis_joint: " << is_joint << "\tdiameter: " << diameter
              << "\tradius: " << radius << "\tmin_node_degree: " << min_node_degree << "\tmax_node_degree: " << max_node_degree
              << "\tnum_nodes: " << num_nodes << "\tnum_edges: " << num_edges << "\tavg_degree: " << avg_degree << std::endl;
  };

  for (auto const& dataset : all_datasets) {
    GxlGEDEnv env;
    ::util::setup_environment(dataset, false, env);
    print_statistics(env, dataset);
  }

  std::cout << std::endl;

  for (auto const& dataset : sized_datasets) {
    std::size_t max_max_size_div_10 = 0;
    if (dataset == "AIDS") 
      max_max_size_div_10 = 8;
    else if (dataset == "Protein")
      max_max_size_div_10 = 7;
    else if (dataset == "Mutagenicity")
      max_max_size_div_10 = 10;

    for (std::size_t max_size_dev_10 = 1; max_size_dev_10 <= max_max_size_div_10; max_size_dev_10++) {
      GxlGEDEnv env;
      ::util::setup_environment(dataset, max_size_dev_10, env);
      print_statistics(env, dataset);
    }
  }

  std::cout << std::endl;

  for (size_t num_nodes = 5; num_nodes <= 40; num_nodes += 5) {
    GxlGEDEnv env;
    std::mt19937 rng;

    for (size_t i = 0; i < 100; ++i)
      graph_gen_power(env, rng, num_nodes, 2, node_labels, edge_labels);

    print_statistics(env, "power" + std::to_string(num_nodes));
  }

  std::cout << std::endl;

  for (size_t num_nodes = 20; num_nodes <= 40; num_nodes += 20) {
    GxlGEDEnv env;
    std::mt19937 rng;

    for (size_t i = 0; i < 100; ++i)
      graph_gen_cluster(env, rng, num_nodes, 4, 0.8, 0.02, node_labels, edge_labels);

    print_statistics(env, "cluster" + std::to_string(num_nodes));
  }
}

int main(int argc, char* argv[]) {
  test_ls_all_datasets();

  //compute_statistics();


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
  //   GraphStatistics stats = graph_stat_compute(env, id);
  //   std::cout << env.get_graph_name(id) << " D: " << stats.diameter << " R: " << stats.radius << " joint: " << stats.is_joint <<
  //     " min_pow: " << stats.min_node_degree << " max_pow: " << stats.max_node_degree << std::endl;
  // }

  // test_ls_graph_sizes();
  // test_ls_all_datasets();
  // test_ls_rand_graphs();
  // test_lsape_all();

  return 0;
}
