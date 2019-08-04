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
	int n_iter_coarse = 2000;
	const int n_iter_max = 10000;
	float max_imbalance = 0.01;

	if (argc < 4) {
		cerr << "Call: gpart <n_partitions> <max-final-imbalance> <grf_file> >map \n";
		exit(888);
	}
	n_part = atoi(argv[1]);
	max_imbalance = atof(argv[2]);
	CAPart dc(dummy, n_part, argv[3]);
	cerr << "Input graph has " << dc.n_nodes << " nodes." << endl;
	dc.init(false);
	unsigned alts = 0; i = 0;

	if (dc.n_nodes / n_part > 10000 /*&& dc.n_nodes > 1000000*/) {

		dc.createSplit(50, 0.1);
		cerr << "Split took " << dc.split->iteration << " iterations\n";
	//	while ((abs(dc.pressure).max() > 5) || (dc.pressure.sum() != 0)) {
		while (i < n_iter_coarse) {
			i++;
			alts += dc.coarse->step();
			if (i%100 == 0) {
				cerr << "COARSE " << alts << " PMAX: " << 100. * abs(dc.coarse->pressure).max() / dc.coarse->desired[0] << "% "<< " Halo: " << 100. * dc.coarse->border_count.sum() / (float)dc.coarse->n_nodes << "% PotMAX " << dc.coarse->potential.max()<< endl;
				if (alts < 100)
					break;
				alts = 0;
			}
		}
		dc.coarse->project(dc.split, &dc);
	}
	alts = 0; i = 0;
//	while ((abs(dc.pressure).max() > 5) || (dc.pressure.sum() != 0)) {
	while (i < n_iter_max) {
		i++;
		alts += dc.step();
		if (i%100 == 0) {
			float max_imb = 100. * abs(dc.pressure).max() / dc.desired[0];
			cerr << "FINAL " << alts << " PMAX: " << max_imb << "% "<< " Halo: " << 100. * dc.border_count.sum() / (float)dc.n_nodes << "% PotMAX " << dc.potential.max()<< endl;
			if (alts < 100 || max_imb < max_imbalance)
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


