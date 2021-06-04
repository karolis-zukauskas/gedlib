#ifndef SRC_KAZU_METHOD_HPP_
#define SRC_KAZU_METHOD_HPP_

#include "util.hpp"

struct MethodRunResult {
  GEDGraph::GraphID g_id;
  GEDGraph::GraphID h_id;
  double ub = 0.0;
  double runtime = 0.0;

  MethodRunResult() {}

  MethodRunResult(GEDGraph::GraphID g_id, GEDGraph::GraphID h_id, double ub, double runtime)
    : g_id(g_id),
      h_id(h_id),
      ub(ub),
      runtime(runtime) {}
};

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
      } else if (ged_method_ == Options::GEDMethod::BP_BEAM) {
        name << "BP_BEAM" << " (" << initialization_method_ << ")";
      } else if (ged_method_ == Options::GEDMethod::BIPARTITE) {
        name << "BP";
      }  else if (ged_method_ == Options::GEDMethod::BRANCH_FAST) {
        name << "BRANCHFAST";
      } else if (ged_method_ == Options::GEDMethod::BRANCH) {
        name << "BRANCH";
      } else if (ged_method_ == Options::GEDMethod::BRANCH_UNIFORM) {
        name << "BRANCH_UNIFORM";
      } else if (ged_method_ == Options::GEDMethod::BRANCH_UNIFORM2) {
        name << "BRANCH_UNIFORM2";
      } else if (ged_method_ == Options::GEDMethod::BRANCH_UNIFORM3) {
        name << "BRANCH_UNIFORM3";
      } else if (ged_method_ == Options::GEDMethod::BRANCH_UNIFORM4) {
        name << "BRANCH_UNIFORM4";
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
      } else if (ged_method_ == Options::GEDMethod::STAR7) {
        name << "STAR7";
      } else if (ged_method_ == Options::GEDMethod::STAR8) {
        name << "STAR8";
      } else if (ged_method_ == Options::GEDMethod::SUBGRAPH) {
        name << "SUBGRAPH";
      } else if (ged_method_ == Options::GEDMethod::WALKS) {
        name << "WALKS";
      } else if (ged_method_ == Options::GEDMethod::NODE) {
        name << "NODE";
      }
      else {
        assert(false);
      }

      return name.str();
    }

    void run_on_dataset(const std::string & dataset, GxlGEDEnv & env, std::vector<MethodRunResult>& results, bool ensure_n_greater_m, bool test_only_unique_pairs,
      double & avg_lb, double & avg_ub, double & avg_runtime, double& avg_init_time, double& avg_ls_iterations_time, double& avg_num_ls_iterations) const {

      env.set_method(ged_method_, options_(dataset));
      env.init_method();

      size_t num_runs;
      if (test_only_unique_pairs)
        num_runs = static_cast<size_t>((env.graph_ids().second * (env.graph_ids().second - 1)) / 2.0);
      else
        num_runs = env.graph_ids().second * env.graph_ids().second;
      ProgressBar progress_bar(num_runs);

      results.clear();
      results.reserve(num_runs);

      avg_runtime = 0;
      avg_ub = 0;
      avg_lb = 0;
      avg_init_time = 0;
      avg_ls_iterations_time = 0;
      avg_num_ls_iterations = 0;

      auto run_ged = [&](GEDGraph::GraphID g_id, GEDGraph::GraphID h_id) -> void {
        GEDGraph::GraphID _g_id = g_id;
        GEDGraph::GraphID _h_id = h_id;

        if (ensure_n_greater_m && env.get_num_nodes(_g_id) < env.get_num_nodes(_h_id))
            std::swap(_g_id, _h_id);

        env.run_method(_g_id, _h_id);
        results.emplace_back(_g_id, _h_id, env.get_upper_bound(_g_id, _h_id), env.get_runtime(_g_id, _h_id));

        avg_lb += env.get_lower_bound(_g_id, _h_id);
        avg_ub += env.get_upper_bound(_g_id, _h_id);
        avg_runtime += env.get_runtime(_g_id, _h_id);

        auto pLSMethod = dynamic_cast<LSBasedMethod<GXLLabel, GXLLabel>*>(env.get_method());
        if (pLSMethod != nullptr) {
          avg_init_time += pLSMethod->initial_sulutions_time.count();
          avg_ls_iterations_time += pLSMethod->ls_iterations_time.count();
          avg_num_ls_iterations += static_cast<double>(pLSMethod->num_ls_iterations);
        }

        // progress_bar.increment();
        // std::cout << "\r\t" << name() << ": " << progress_bar << std::flush;
      };

      // std::cout << "\n\tDataset: " << dataset << "   Method: " << name() << std::endl;
      if (test_only_unique_pairs) {
        for (GEDGraph::GraphID g_id = env.graph_ids().first; g_id != env.graph_ids().second; g_id++)
          for (GEDGraph::GraphID h_id = g_id + 1; h_id != env.graph_ids().second; h_id++)
            run_ged(g_id, h_id);
      } else {
        for (GEDGraph::GraphID g_id = env.graph_ids().first; g_id != env.graph_ids().second; g_id++)
          for (GEDGraph::GraphID h_id = env.graph_ids().first; h_id != env.graph_ids().second; h_id++)
            run_ged(g_id, h_id);
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

void run_methods(std::vector<Method> const& methods, std::function<GxlGEDEnv()> const& setup_ged_env, std::string const& dataset,
  double edge_density, bool ensure_n_greater_m, bool test_only_unique_pairs, std::string const& result_filename,
  std::map<std::string, std::vector<double>>* pMethodUbMap = nullptr) {

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

    GraphStatsMap stats;
    graph_stat_compute_all(env, stats);

    std::vector<MethodRunResult> results;

    methods[i].run_on_dataset(dataset, env, results, ensure_n_greater_m, test_only_unique_pairs,
      avg_lb, avg_ub, avg_runtime, avg_init_time, avg_ls_iterations_time, avg_num_ls_iterations);

#pragma omp critical
    {
    if (pMethodUbMap != nullptr) {
      std::string methodName = methods[i].name();
      assert(pMethodUbMap->find(methodName) != pMethodUbMap->end());
      pMethodUbMap->at(methodName).push_back(avg_ub);
    }

#ifdef WRITE_STATS_FILE
      std::ofstream stats_file(s_stats_filename.c_str(), std::ios_base::app);
      for (auto const& result : results) {
        GraphDiff diff = graph_diff_compute(env, stats, result.g_id, result.h_id);

        stats_file << methods[i].name() << "," << result.ub << "," << result.runtime << "," << dataset << "," << result.g_id << "," << result.h_id << ","
                   << diff.node_count << ","
                   << diff.edge_count << ","
                   << diff.avg_node_degree << ","
                   << diff.edge_density << ","
                   << diff.node_label_count << ","
                   << diff.edge_label_count << ","
                   << diff.diameter << ","
                   << diff.radius << ","
                   << (diff.both_joint ? "joint" : "disjoint") << "\n";
      }
      stats_file.close();
#endif

      // lb, ub, runtime, init_solutions_t, ls_iterations_t, ls_iterations, num_nodes, method
      std::ofstream result_file(result_filename.c_str(), std::ios_base::app);

      result_file << avg_lb << "," << avg_ub << "," << avg_runtime << "," << avg_init_time << ",";
      result_file << avg_ls_iterations_time << "," << avg_num_ls_iterations << ",";
      result_file << env.get_avg_num_nodes() << "," << methods[i].name() << std::endl;

      result_file.close();
    }
  }
}

void run_on_sized_dataset(std::vector<Method> const& methods, std::string const& dataset, bool ensure_n_greater_m, bool test_only_unique_pairs,
  std::string const& result_filename) {

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

    run_methods(methods, setup_ged_env, dataset, 0.0, ensure_n_greater_m, test_only_unique_pairs, result_filename);
  }
}

void run_on_test_dataset(std::vector<Method> const& methods, std::string const& dataset, bool ensure_n_greater_m, bool test_only_unique_pairs,
  std::string const& result_filename) {

  std::cout << "\n=== " << dataset << " ===\n";

  auto setup_ged_env = [&dataset]() -> GxlGEDEnv {
    GxlGEDEnv env;
    ::util::setup_environment(dataset, false, env);

    return env;
  };

  run_methods(methods, setup_ged_env, dataset, 0.0, ensure_n_greater_m, test_only_unique_pairs, result_filename);
}

void run_on_generated_dataset(std::vector<Method> const& methods, std::function<void(GxlGEDEnv&)> const& generate_graphs, std::string const& dataset,
  double edge_density, bool ensure_n_greater_m, bool test_only_unique_pairs, std::string const& result_filename) {

  std::cout << "\n=== " << dataset << " ===\n";

  auto setup_ged_env = [&generate_graphs]() -> GxlGEDEnv {
    GxlGEDEnv env;
    generate_graphs(env);

    return env;
  };

  run_methods(methods, setup_ged_env, dataset, edge_density, ensure_n_greater_m, test_only_unique_pairs, result_filename);
}

#endif
