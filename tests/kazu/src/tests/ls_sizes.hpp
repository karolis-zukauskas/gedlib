#ifndef SRC_KAZU_TEST_LS_SIZES_HPP_
#define SRC_KAZU_TEST_LS_SIZES_HPP_

#include "../util.hpp"
#include "../method.hpp"

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

  std::string stats_filename = ::util::create_stats_file("ls_graph_sizes");
  for (auto dataset : datasets) {
    try {
      std::string results_filename = ::util::create_result_file("ls_graph_sizes", dataset);

      run_on_sized_dataset(methods, dataset, false, TEST_ONLY_UNIQUE_PAIRS, results_filename, stats_filename);
    }
    catch (const std::exception & error) {
      std::cerr << error.what() << ". " << "Error on test_ls_graph_sizes: " << dataset << ".\n";
    }
  }
}

#endif
