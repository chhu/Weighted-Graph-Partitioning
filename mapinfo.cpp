#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <valarray>
#include <cassert>
#include "wpart.hpp"

using namespace std;
using namespace wpart;

int main(int argc, char **argv) {
	int entries, i = 0;
	int n_part = 0;
	int n_iter = 50000;
	if (argc < 3) {
		cerr << "Call: gpart <map_file> <grf_file> \n";
		exit(888);
	}
	ifstream map(argv[1], ifstream::in);
	map >> entries;
	valarray<int32_t> part_from_map(entries);
	for (int i = 0; i < entries; i++) {
		int idx, part;
		map >> idx; map >> part;
		assert(idx < entries);
		part_from_map[idx] = part;
	}

	CAPart dc(entries, part_from_map.max() + 1, argv[2]);

	cerr << "Input graph has " << dc.n_nodes << " nodes, map contains " << dc.n_partitions << " partitions." << endl;
	dc.partitioning = part_from_map;
	dc.init(true);
	dc.no_alterations = true;
	dc.step();
	cerr << "Max imbalance: " << 100. * abs(dc.pressure).max() / dc.desired[0] << "% "<< " Halo: " << 100. * dc.border_count.sum() / (float)dc.n_nodes << "%" << endl;
	cerr << "Per domain:\n";
	for (int i = 0; i < dc.border_count.size(); i++)
		cerr << "Domain " << i << ": Imbalance: " << dc.pressure[i] << " nodes / " << 100. * dc.pressure[i] / dc.desired[i] << "%, halo-count: " << dc.border_count[i] << " nodes / " <<  100. * dc.border_count[i] / (float)dc.desired[i] << "%" << endl;

	return 0;
}


