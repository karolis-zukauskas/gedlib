#//////////////////////////////////////////////////////////////////////////#
#                                                                          #
#   Copyright (C) 2018 by David B. Blumenthal                              #
#                                                                          #
#   This file is part of GEDLIB.                                           #
#                                                                          #
#   GEDLIB is free software: you can redistribute it and/or modify it      #
#   under the terms of the GNU Lesser General Public License as published  #
#   by the Free Software Foundation, either version 3 of the License, or   #
#   (at your option) any later version.                                    #
#                                                                          #
#   GEDLIB is distributed in the hope that it will be useful,              #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           #
#   GNU Lesser General Public License for more details.                    #
#                                                                          #
#   You should have received a copy of the GNU Lesser General Public       #
#   License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>. #
#                                                                          #
#//////////////////////////////////////////////////////////////////////////#

##
# @file sample.py
# @brief Python script that generates a random sample of given size from a given dataset.
#
# @details 
# Usage: 
# ```sh
# $ python sample.py \<dataset\> \<sample\> [-h] [--help] [--exclude \<graph collection\>] [--balanced] \<--size \<size\> | --size_ratio \<size ratio\>\>
# ```
# 
# Arguments:
#
# <table>
# <tr><th colspan="2"> positional arguments
# <tr><td> <tt>\<dataset\></tt> <td> path to existing graph collection XML file from which the sample should be drawn; must respect GraphCollection.dtd
# <tr><td> <tt>\<sample\></tt> <td> path to sample file to be generated by the script
# <tr><th colspan="2"> optional arguments
# <tr><td> <tt>-h</tt> <td> show help
# <tr><td> <tt>--balanced</tt> <td> generate sample with equal number of graphs per class
# <tr><td> <tt>--exclude \<graph collection\></tt> <td> path to existing graph collection XML file whose graphs should be excluded from the sample; must respect GraphCollection.dtd
# <tr><td> <tt>--size \<size\></tt> <td> number of graphs listed in the sample; must be between 0 and the number of graphs listed in <tt>\<dataset\></tt>
# <tr><td> <tt>--size_ratio \<size ratio\></tt> <td> number of graphs listed in the sample divided by number of graphs listed in <tt>\<dataset\></tt>; must be between 0 and 1
# </table>
'''
Python script that generates a random sample of given size from a given dataset.
'''

import xml.etree.ElementTree as ET
import argparse
import random

# Parse the input arguments.
parser = argparse.ArgumentParser(description="Generates a random sample of given size from a given dataset.")
parser.add_argument("dataset", help="path to existing dataset file")
parser.add_argument("sample", help="path to sample file to be generated by the script")
parser.add_argument("--exclude", help="path to existing file that list the graphs contained in the dataset which should not appear in the sample")
parser.add_argument("--balanced", help="generate sample with equal number of graphs per class", action="store_true")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--size", help="size of sample; must be greater that 0; if larger than size of and the size of the dataset", type=int)
group.add_argument("--size_ratio", help="size of sample divided by size of dataset; must be between 0 and 1", type=float)
args = parser.parse_args()
if args.dataset == args.sample:
	raise Exception("dataset file equals sample file")

# Collect excluded graphs.
excluded_graphs = set()
if args.exclude:
	tree = ET.parse(args.exclude)
	excluded_dataset = tree.getroot()
	for graph in excluded_dataset:
		excluded_graphs.add(graph.attrib["file"])

# Collect the classes.
dataset = ET.parse(args.dataset).getroot()
classes = set()
graph_classes = {graph.attrib["file"] : graph.attrib["class"] for graph in dataset}
for graph in dataset:
	classes.add(graph.attrib["class"])
num_classes = len(classes)

# Collect the candidate graphs and group them w.r.t. their classes.
candidate_graphs = {cl : [] for cl in classes}
for graph in dataset:
	if not graph.attrib["file"] in excluded_graphs:
		candidate_graphs[graph.attrib["class"]].append(graph.attrib["file"])
candidate_sizes = {cl : len(candidate_graphs[cl]) for cl in classes}
total_candidate_size = sum([candidate_sizes[cl] for cl in candidate_sizes])
min_candidate_sizes = min([candidate_sizes[cl] for cl in candidate_sizes])

# Determine the number of sampled graphs per class.
if args.size_ratio:
	if args.size_ratio < 0 or args.size_ratio > 1:
		raise Exception("SIZE_RATIO must be between 0 and 1")
	if args.balanced:
		sample_sizes = {cl : min(min_candidate_sizes, int((total_candidate_size * args.size_ratio) / num_classes)) for cl in classes}
	else:
		sample_sizes = {cl : int(candidate_sizes[cl] * args.size_ratio) for cl in classes}
else:
	if args.size < 0:
		raise Exception("SIZE must be greater than 0")
	if args.size > total_candidate_size:
		args.size = total_candidate_size
		
# Sample the graphs.
sampled_graphs = []
if args.balanced:
    sample_sizes = {cl : min(min_candidate_sizes, int(args.size / num_classes)) for cl in classes}
    sampled_graphs = [graph for cl in classes for graph in random.sample(candidate_graphs[cl], sample_sizes[cl])]
else:
    sampled_graphs = random.sample([graph for cl in classes for graph in candidate_graphs[cl]], args.size)

# Write sampled graphs to XML file.
file = open(args.sample, "w")
file.write("<?xml version=\"1.0\"?>")
file.write("\n<!DOCTYPE GraphCollection SYSTEM \"http://www.inf.unibz.it/~blumenthal/dtd/GraphCollection.dtd\">")
file.write("\n<GraphCollection>")
for graph in sampled_graphs:
	file.write("\n\t<graph file=\"" + graph + "\" class=\"" + graph_classes[graph] + "\"/>")
file.write("\n</GraphCollection>")
file.close()