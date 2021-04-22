#ifndef SRC_KAZU_TEST_LS_ALL_HPP_
#define SRC_KAZU_TEST_LS_ALL_HPP_

#include "../util.hpp"
#include "../method.hpp"

void test_ls_all_datasets() {
  std::vector<std::string> datasets = {
    "AIDS", "Letter_HIGH", "Mutagenicity", "Protein", "GREC", "Fingerprint",
    // "Fingerprint", "AIDS", "Letter_HIGH", "GREC"
  };

  std::vector<Method> const methods {
    Method (Options::GEDMethod::IPFP, "BRANCH_UNIFORM", true),
    Method (Options::GEDMethod::IPFP, "BRANCH_FAST", true),
    Method (Options::GEDMethod::IPFP, "BRANCH", true),
    Method (Options::GEDMethod::IPFP, "NODE", true),
    Method (Options::GEDMethod::IPFP, "BP_BEAM", true, " --ls-initialization-method BRANCH_FAST"),
    Method (Options::GEDMethod::IPFP, "BP_BEAM", true, " --ls-initialization-method BIPARTITE"),
    Method (Options::GEDMethod::IPFP, "BP_BEAM", true, " --ls-initialization-method NODE"),
    Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method BRANCH_FAST"),
    Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method BIPARTITE"),
    Method (Options::GEDMethod::IPFP, "REFINE", true, " --ls-initialization-method NODE"),
    // Method (Options::GEDMethod::IPFP, "WALKS", true),
    // Method (Options::GEDMethod::IPFP, "SUBGRAPH", true),
    Method (Options::GEDMethod::IPFP, "STAR", true),
    Method (Options::GEDMethod::IPFP, "STAR4", true),
    // Method (Options::GEDMethod::IPFP, "STAR5", true),
    Method (Options::GEDMethod::IPFP, "STAR6", true),
    Method (Options::GEDMethod::IPFP, "RANDOM", true),
  };

  std::string stats_filename = ::util::create_stats_file("ls_all_datasets");
  for (auto dataset : datasets) {
    try {
      std::string results_filename = ::util::create_result_file("ls_all_datasets", dataset);

      run_on_test_dataset(methods, dataset, false, TEST_ONLY_UNIQUE_PAIRS, results_filename, stats_filename);
    }
    catch (const std::exception & error) {
      std::cerr << error.what() << ". " << "Error on test_ls_all_datasets: " << dataset << ".\n";
    }
  }
}

#endif
