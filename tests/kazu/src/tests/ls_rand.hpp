#ifndef SRC_KAZU_TEST_LS_RAND_HPP_
#define SRC_KAZU_TEST_LS_RAND_HPP_

#include "../util.hpp"
#include "../method.hpp"
#include "../graph_gen.hpp"

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

  std::string stats_filename = ::util::create_stats_file("ls_generated");
  std::string results_filename = ::util::create_result_file("ls_generated", "rand");
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

        run_on_generated_dataset(methods, generate_graphs, dataset, edge_density, false, TEST_ONLY_UNIQUE_PAIRS, results_filename, stats_filename);
      } catch (const std::exception & error) {
        std::cerr << error.what() << ". " << "Error on test_ls_rand_graphs: " << dataset << ".\n";
      }
    }
  }
}

#endif
