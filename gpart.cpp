#include <stdio.h>
#include <iostream>
#include <vector>
#include <valarray>
#include <cassert>
#include "wpart.hpp"

using namespace std;
using namespace wpart;

int main(int argc, char **argv) {
	int dummy, entries, i = 0;
	int n_part = 2;
	int n_iter = 50000;
	if (argc < 4) {
		cerr << "Call: gpart <n_partitions> <iterations> <grf_file> >map \n";
		exit(888);
	}
	n_part = atoi(argv[1]);
	n_iter = atoi(argv[2]);
	CAPart dc(dummy, n_part, argv[3]);
	cerr << "Input graph has " << dc.n_nodes << " nodes." << endl;
	dc.init(false);
	unsigned alts = 0; i = 0;
//	while ((abs(dc.pressure).max() > 5) || (dc.pressure.sum() != 0)) {
	while (i < n_iter) {
		i++;
		alts += dc.step();
		if (i%1000 == 0) {
			cerr << alts << " PMAX: " << 100. * abs(dc.pressure).max() / dc.desired[0] << "% "<< " Halo: " << 100. * dc.border_count.sum() / (float)dc.n_nodes << "% PotMAX " << dc.potential.max()<< endl;
			if (alts < 100)
				break;
			alts = 0;
		}
	}
	// Write to std out
	cout << dc.n_nodes << endl;
	for (int i = 0; i < dc.n_nodes; i++)
		cout << i << " " << dc.partitioning[i] << endl;

	cerr << "Final worst imbalance: " << 100. * abs(dc.pressure).max()/ dc.desired[0] << "% Halo: " << 100. * dc.border_count.sum() / (float)dc.n_nodes << "% Iterations: " << i << endl;

	return 0;
}


