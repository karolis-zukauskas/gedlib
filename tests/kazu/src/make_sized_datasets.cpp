#include <dirent.h>
#include <cstdio>
#include <iostream>
#include <string>
#include "util.hpp"

using GxlGEDEnv = ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel>;

void create_dataset_file(std::string const& filename, GxlGEDEnv const& ged_env, std::vector<std::pair<std::string, size_t>> const& graphs) {
  std::ofstream file(filename.c_str());

  file << "<?xml version=\"1.0\"?>" << std::endl;
  file << "<GraphCollection>" << std::endl;

  size_t count = 0;
  for (auto& pair : graphs) {
    if (count == 100)
      break;

    file << "  <graph file=\"" << std::get<0>(pair) << "\" class=\"" << ged_env.get_graph_class(std::get<1>(pair)) << "\" />" << std::endl;
    count++;
  }

  file << "</GraphCollection>" << std::endl;
  file.close();
}

int main() {
  std::vector<std::string> datasets = {
    "AIDS", "Mutagenicity", "Protein",
  };

  for (std::string const& dataset : datasets) {
    std::vector<std::vector<std::pair<std::string, size_t>>> graph_bins(10);
    std::string graph_dir = "../../../data/datasets/" + dataset + "/data/";
    GxlGEDEnv ged_env;

    DIR *pDir;
    dirent* pEnt;

    if ((pDir = opendir (graph_dir.c_str())) != nullptr) {
      while ((pEnt = readdir (pDir)) != nullptr) {
        std::string filename(pEnt->d_name);
        if (filename.find(".gxl") == -1)
          continue;

        size_t g_id = ged_env.load_gxl_graph(graph_dir, filename, util::node_type(dataset), util::edge_type(dataset),
          util::irrelevant_node_attributes(dataset), util::irrelevant_edge_attributes(dataset), ged::undefined(), "any");

        auto num_nodes = ged_env.get_num_nodes(g_id);
        auto bin_id = static_cast<uint>(num_nodes) / 10;
        if (bin_id >= 10)
          continue;

        graph_bins[bin_id].push_back(std::make_pair(filename, g_id));
      }

      std::cout << "\n" << dataset << "\n";

      size_t bin_id = 1;
      for (auto& bin : graph_bins) {
        std::mt19937 urng;
        urng.seed(time(0));
        std::shuffle(bin.begin(), bin.end(), urng);

        std::cout << "bin_size: " << bin.size();
        double avg_num_nodes = 0;
        for (auto& pair : bin)
          avg_num_nodes += static_cast<double>(ged_env.get_num_nodes(std::get<1>(pair)));

        avg_num_nodes /= static_cast<double>(bin.size());
        std::cout << "\t avg_nodes: " << avg_num_nodes << std::endl;

      std::stringstream ss;
      ss << "../results/" << dataset << bin_id << "-" << bin_id + 9 << ".xml";
      create_dataset_file(ss.str(), ged_env, bin);
      bin_id += 10;
      }

      closedir (pDir);
    } else {
      perror ("ERROR reading directory");
      return -1;
    }
  }

  return 0;
}
