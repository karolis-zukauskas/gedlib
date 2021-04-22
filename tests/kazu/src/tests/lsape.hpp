#ifndef SRC_KAZU_TEST_LSAPE_HPP_
#define SRC_KAZU_TEST_LSAPE_HPP_

#include "../util.hpp"
#include "../method.hpp"

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

  std::string stats_filename;

  stats_filename = ::util::create_stats_file("lsape_all_datasets");
  for (auto dataset : all_datasets) {
    try {
      std::string results_filename = ::util::create_result_file("lsape_all_datasets", dataset);

      run_on_test_dataset(methods, dataset, false, TEST_ONLY_UNIQUE_PAIRS, results_filename, stats_filename);
    }
    catch (const std::exception & error) {
      std::cerr << error.what() << ". " << "Error on test_lsape_all: " << dataset << ".\n";
    }
  }

  stats_filename = ::util::create_stats_file("lsape_sized_datasets");
  for (auto dataset : sized_datasets) {
    try {
      std::string results_filename = ::util::create_result_file("lsape_sized_datasets", dataset);

      run_on_sized_dataset(methods, dataset, false, TEST_ONLY_UNIQUE_PAIRS, results_filename, stats_filename);
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

  stats_filename = ::util::create_stats_file("lsape_generated");
  std::string results_filename = ::util::create_result_file("lsape_generated", "rand");
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

        run_on_generated_dataset(methods, generate_graphs, dataset, edge_density, false, TEST_ONLY_UNIQUE_PAIRS, results_filename, stats_filename);
      } catch (const std::exception & error) {
        std::cerr << error.what() << ". " << "Error on test_lsape_all: " << dataset << ".\n";
      }
    }
  }
}

#endif
