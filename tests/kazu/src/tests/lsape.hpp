#ifndef SRC_KAZU_TEST_LSAPE_HPP_
#define SRC_KAZU_TEST_LSAPE_HPP_

#include "../util.hpp"
#include "../method.hpp"

void test_lsape_all(std::vector<Method> const& methods, std::vector<std::string> const& datasets) {
  for (auto dataset : datasets) {
    try {
      std::string results_filename = ::util::create_result_file("lsape_all_datasets", dataset);

      run_on_test_dataset(methods, dataset, false, TEST_ONLY_UNIQUE_PAIRS, results_filename);
    }
    catch (const std::exception & error) {
      std::cerr << error.what() << ". " << "Error on test_lsape_all: " << dataset << ".\n";
    }
  }
}

void test_lsape_sized(std::vector<Method> const& methods, std::vector<std::string> const& datasets) {
  for (auto dataset : datasets) {
    try {
      std::string results_filename = ::util::create_result_file("lsape_sized_datasets", dataset);

      run_on_sized_dataset(methods, dataset, false, TEST_ONLY_UNIQUE_PAIRS, results_filename);
    }
    catch (const std::exception & error) {
      std::cerr << error.what() << ". " << "Error on test_lsape_all: " << dataset << ".\n";
    }
  }
}

void test_lsape_rand(std::vector<Method> const& methods, size_t num_graphs, size_t node_variance,
  std::vector<size_t> const& graph_sizes, std::vector<double> const& edge_densities) {

  // size_t const num_graphs = 100;
  // size_t const node_variance = 5;
  // std::vector<size_t> const graph_sizes {
  //   10, 20, 30, 40, 50,
  // };
  // std::vector<double> const edge_densities {
  //   0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
  // };
  std::vector<std::string> const node_labels { "A", "B", "C" };
  std::vector<std::string> const edge_labels { "0", "1" };

  std::string results_filename = ::util::create_result_file("lsape_generated", "rand");
  for (size_t graph_size : graph_sizes) {
    for (double edge_density : edge_densities) {
      std::stringstream ss;
      ss << "random n: " << graph_size << " ed: " << edge_density;
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

          ::util::setup_generated_environment(env);
        };

        run_on_generated_dataset(methods, generate_graphs, dataset, edge_density, false, TEST_ONLY_UNIQUE_PAIRS, results_filename);
      } catch (const std::exception & error) {
        std::cerr << error.what() << ". " << "Error on test_lsape_all: " << dataset << ".\n";
      }
    }
  }
}

#endif
