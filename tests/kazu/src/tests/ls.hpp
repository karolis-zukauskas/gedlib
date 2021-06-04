#ifndef SRC_KAZU_TEST_LS_ALL_HPP_
#define SRC_KAZU_TEST_LS_ALL_HPP_

#include "../util.hpp"
#include "../method.hpp"

void test_ls_all_datasets(std::vector<Method> const& methods, std::vector<std::string> const& datasets) {
  for (auto dataset : datasets) {
    try {
      std::string results_filename = ::util::create_result_file("ls_all_datasets", dataset);

      run_on_test_dataset(methods, dataset, false, TEST_ONLY_UNIQUE_PAIRS, results_filename);
    }
    catch (const std::exception & error) {
      std::cerr << error.what() << ". " << "Error on test_ls_all_datasets: " << dataset << ".\n";
    }
  }
}

void test_ls_graph_sizes(std::vector<Method> const& methods, std::vector<std::string> const& datasets) {
  for (auto dataset : datasets) {
    try {
      std::string results_filename = ::util::create_result_file("ls_graph_sizes", dataset);

      run_on_sized_dataset(methods, dataset, false, TEST_ONLY_UNIQUE_PAIRS, results_filename);
    }
    catch (const std::exception & error) {
      std::cerr << error.what() << ". " << "Error on test_ls_graph_sizes: " << dataset << ".\n";
    }
  }
}

void test_ls_rand_graphs(std::vector<Method> const& methods, size_t num_graphs, size_t node_variance,
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

  std::string results_filename = ::util::create_result_file("ls_generated", "rand");
  for (size_t graph_size : graph_sizes) {
    for (double edge_density : edge_densities) {
      std::stringstream ss;
      ss << "random n: " << graph_size << " ed: " << edge_density;
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

          ::util::setup_generated_environment(env);
        };

        run_on_generated_dataset(methods, generate_graphs, dataset, edge_density, false, TEST_ONLY_UNIQUE_PAIRS, results_filename);
      } catch (const std::exception & error) {
        std::cerr << error.what() << ". " << "Error on test_ls_rand_graphs: " << dataset << ".\n";
      }
    }
  }
}

auto print_statistics = [](GxlGEDEnv& env, std::string const& dataset) -> void {
    double is_joint = 0.0;
    double num_nodes = 0.0;
    double num_edges = 0.0;
    double avg_degree = 0.0;
    double edge_density = 0.0;
    double diameter = 0.0;
    double radius = 0.0;
    double min_node_degree = 0.0;
    double max_node_degree = 0.0;

    for (GEDGraph::GraphID i = env.graph_ids().first; i != env.graph_ids().second; ++i) {
      GraphStatistics stats = graph_stat_compute_single(env, i);

      is_joint += static_cast<double>(stats.is_joint);
      num_nodes += static_cast<double>(stats.num_nodes);
      num_edges += static_cast<double>(stats.num_edges);
      avg_degree += stats.avg_degree;
      edge_density += stats.edge_density;
      diameter += static_cast<double>(stats.diameter);
      radius += static_cast<double>(stats.radius);
      min_node_degree += static_cast<double>(stats.min_node_degree);
      max_node_degree += static_cast<double>(stats.max_node_degree);
    }

    is_joint /= static_cast<double>(env.num_graphs());
    num_nodes /= static_cast<double>(env.num_graphs());
    num_edges /= static_cast<double>(env.num_graphs());
    avg_degree /= static_cast<double>(env.num_graphs());
    edge_density /= static_cast<double>(env.num_graphs());
    diameter /= static_cast<double>(env.num_graphs());
    radius /= static_cast<double>(env.num_graphs());
    min_node_degree /= static_cast<double>(env.num_graphs());
    max_node_degree /= static_cast<double>(env.num_graphs());


    std::cout  /*<< std::setprecision(3) */ << dataset << "\tis_joint: " << is_joint << "\tdiameter: " << diameter
              << "\tradius: " << radius << "\tmin_node_degree: " << min_node_degree << "\tmax_node_degree: " << max_node_degree
              << "\tnum_nodes: " << num_nodes << "\tnum_edges: " << num_edges << "\tavg_degree: " << avg_degree
              << "\tedge_density: " << edge_density << std::endl << std::endl;
  };

void test_ls_power_graphs(std::vector<Method> const& methods, size_t num_graphs, size_t node_variance,
  std::vector<size_t> const& graph_sizes, std::vector<size_t> const& edges_per_node) {

  std::vector<std::string> const node_labels { "A", "B", "C" };
  std::vector<std::string> const edge_labels { "0", "1" };

  std::string results_filename = ::util::create_result_file("ls_generated", "power");
  for (size_t graph_size : graph_sizes) {
    for (size_t edges : edges_per_node) {
      std::stringstream ss;
      ss << "power n: " << graph_size << " epn: " << edges;
      std::string dataset = ss.str();

      try {
        auto generate_graphs = [num_graphs, graph_size, edges, &node_labels, &edge_labels, node_variance](GxlGEDEnv& env) -> void {
          std::mt19937 size_variance_rng(graph_size * edges);          // For graph_size variance
          std::mt19937 graph_gen_rng(graph_size * edges);              // For generating graphs

          for (size_t i = 0; i < num_graphs; ++i) {
            double multiplier = static_cast<double>((static_cast<long>(size_variance_rng() % 2000)) - 1000) / 1000.0;
            size_t num_nodes = graph_size + static_cast<size_t>(std::round(multiplier * static_cast<double>(node_variance)));

            graph_gen_power(env, graph_gen_rng, num_nodes, edges, node_labels, edge_labels);
          }

          ::util::setup_generated_environment(env);
        };

        // GxlGEDEnv env;
        // generate_graphs(env);
        // print_statistics(env, dataset);

        run_on_generated_dataset(methods, generate_graphs, dataset, static_cast<double>(edges), false, TEST_ONLY_UNIQUE_PAIRS, results_filename);
      } catch (const std::exception & error) {
        std::cerr << error.what() << ". " << "Error on test_ls_power_graphs: " << dataset << ".\n";
      }
    }
  }
}

void test_ls_cluster_graphs(std::vector<Method> const& methods, size_t num_graphs, size_t node_variance,
  std::vector<size_t> const& graph_sizes, std::vector<size_t> const& num_clusters, double p_inner_edge, double p_outer_edge) {

  std::vector<std::string> const node_labels { "A", "B", "C" };
  std::vector<std::string> const edge_labels { "0", "1" };

  std::string results_filename = ::util::create_result_file("ls_generated", "cluster");
  for (size_t graph_size : graph_sizes) {
    for (size_t clusters : num_clusters) {
      std::stringstream ss;
      ss << "cluster n: " << graph_size << " c: " << clusters;
      std::string dataset = ss.str();

      try {
        auto generate_graphs = [num_graphs, graph_size, clusters, &node_labels, &edge_labels, node_variance, p_inner_edge, p_outer_edge](GxlGEDEnv& env) -> void {
          std::mt19937 size_variance_rng(graph_size * clusters);          // For graph_size variance
          std::mt19937 graph_gen_rng(graph_size * clusters);              // For generating graphs

          for (size_t i = 0; i < num_graphs; ++i) {
            double multiplier = static_cast<double>((static_cast<long>(size_variance_rng() % 2000)) - 1000) / 1000.0;
            size_t num_nodes = graph_size + static_cast<size_t>(std::round(multiplier * static_cast<double>(node_variance)));

            graph_gen_cluster(env, graph_gen_rng, num_nodes, clusters, p_inner_edge, p_outer_edge, node_labels, edge_labels);
          }

          ::util::setup_generated_environment(env);
        };

        // GxlGEDEnv env;
        // generate_graphs(env);
        // print_statistics(env, dataset);

        run_on_generated_dataset(methods, generate_graphs, dataset, static_cast<double>(0.0), false, TEST_ONLY_UNIQUE_PAIRS, results_filename);
      } catch (const std::exception & error) {
        std::cerr << error.what() << ". " << "Error on test_ls_power_graphs: " << dataset << ".\n";
      }
    }
  }
}

#endif
