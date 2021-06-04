#ifndef SRC_KAZU_TEST_STATISTICS_HPP_
#define SRC_KAZU_TEST_STATISTICS_HPP_

#include "../util.hpp"
#include "../graph_statistics.hpp"

void min_max_statistics() {
  std::vector<std::string> all_datasets = {
    "AIDS", "Mutagenicity", "Protein", "GREC", "Letter_HIGH", "Fingerprint"
  };

  auto print_statistics = [](GxlGEDEnv& env, std::string const& dataset) -> void {
    size_t min_nodes = std::numeric_limits<size_t>::max();
    size_t max_nodes = 0;
    double avg_nodes = 0;

    size_t min_edges = std::numeric_limits<size_t>::max();
    size_t max_edges = 0;
    double avg_edges = 0;

    for (GEDGraph::GraphID i = env.graph_ids().first; i != env.graph_ids().second; ++i) {
      auto g = env.get_graph(i);

      if (max_nodes < g.num_nodes)
        max_nodes = g.num_nodes;
      if (min_nodes > g.num_nodes)
        min_nodes = g.num_nodes;

      if (max_edges < g.num_edges)
        max_edges = g.num_edges;
      if (min_edges > g.num_edges)
        min_edges = g.num_edges;

      avg_nodes += static_cast<double>(g.num_nodes);
      avg_edges += static_cast<double>(g.num_edges);
    }

    avg_nodes /= static_cast<double>(env.num_graphs());
    avg_edges /= static_cast<double>(env.num_graphs());

    std::cout << dataset << ", " << std::setprecision(3) << min_nodes << ", " << max_nodes << ", " << avg_nodes
              << ", " << min_edges << ", " << max_edges << ", " << avg_edges << std::endl;
  };

  std::cout << "Dataset, min_nodes, max_nodes, avg_nodes, min_edges, max_edges, avg_edges" << std::endl;

  for (auto const& dataset : all_datasets) {
    GxlGEDEnv env;
    ::util::setup_environment(dataset, false, env);
    print_statistics(env, dataset);
  }
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


    std::cout << std::setprecision(3) << dataset << "\tis_joint: " << is_joint << "\tdiameter: " << diameter
              << "\tradius: " << radius << "\tmin_node_degree: " << min_node_degree << "\tmax_node_degree: " << max_node_degree
              << "\tnum_nodes: " << num_nodes << "\tnum_edges: " << num_edges << "\tavg_degree: " << avg_degree
              << "\tedge_density: " << edge_density << std::endl;
  };

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

  for (auto const& dataset : all_datasets) {
    GxlGEDEnv env;
    ::util::setup_environment(dataset, false, env);
    print_statistics(env, dataset);
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

#endif
