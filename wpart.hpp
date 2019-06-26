/*
 * wpart.hpp
 *
 *	Weighted graph partitioning using cellular automaton rules
 *
 *  Created on: Jun 18, 2019
 *      Author: C Huettig (christian.huettig@dlr.de)
 *
 */

#ifndef WPART_HPP_
#define WPART_HPP_

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include <vector>
#include <set>
#include <cassert>
#include <list>
#include <algorithm>
#include <valarray>
#include <string>
#include <map>
#include <functional>
#include <numeric>
#include <cstring>
#include <limits>
#include <math.h>
#include <time.h>
#include <stdint.h>

using namespace std;

namespace wgcapart {

struct CAPart {

	// INPUT -- These two vectors must be set prior to init() call
	vector<vector<int32_t> > graph;	// Connections between nodes, first vector's size determines N nodes.
	vector<double> weights;			// Size of that vector determines how many partitions will be created. n_partitions = weights.size()

	// OUTPUT
	valarray<uint32_t> partitioning;

	// INTERNAL
	valarray<double>  potential;		// Size: n_nodes, holds pot field
	valarray<int32_t> balance;		// Size: n_partitions + 1, holds deviation from target N
	valarray<int32_t> target;		// Size: n_partitions + 1, holds target domain size.
	valarray<bool>    endangered;   // Size: n_partitions + 1,

	size_t	n_nodes, n_partitions, iteration;


	void init(bool pre_seeded) {
		n_nodes = graph.size();
		n_partitions = weights.size();

		prev_potential.resize(n_nodes, 0.);
		potential.resize(n_nodes, 0.);
		prev_partitioning.resize(n_nodes, 0);
		partitioning.resize(n_nodes, 0);

		// Translate weights to node counts
		valarray<double> weights_va(weights.data(), weights.size());
		if (abs(1 - weights_va.sum()) > 1e-6) {	// Equal distribution, ignore contents, just use length
			cout << "Using equal distribution for " << n_partitions << " domains. Nodes per domain: " << n_nodes / n_partitions << endl;
			prev_pressure.resize(n_partitions + 1, n_nodes / n_partitions);
		} else {
			prev_pressure.resize(n_partitions + 1);
			for (int i = 1; i < n_partitions + 1; i++) {
				prev_pressure[i] = n_nodes * weights_va[i-1];
			}
		}
		prev_pressure[0] = -double(n_nodes);
		pressure = prev_pressure;

		if (pre_seeded) {
			prev_partitioning = partitioning;
			return;
		}
		uint32_t step = n_nodes/ n_partitions - 1;
		for (int i = 0; i < n_partitions; i++)	// Equally distribute seed within nodes
			prev_partitioning[step * i] = i + 1;
		prev_potential[prev_partitioning > 0U] = 10;
		iteration = 0;
	}

	uint32_t step() {
		uint32_t alterations = 0;
		for (size_t i = 0; i < n_nodes; i++) {
			vector<int32_t> &neighbors = graph[i];
			unsigned own_count = 0;
			double own_sum = 0, own_max = 0;
			unsigned other_count = 0, other_max_domain = 0;
			double other_sum = 0, other_max = 0;

			unsigned current_domain = prev_partitioning[i];
			double pot_current = potential[i];

			for (unsigned i_n = 0; i_n < neighbors.size(); i_n++) {
				uint32_t ne = neighbors[i_n];
				uint32_t neighbor_domain = prev_partitioning[ne];
				double neighbor_potential = prev_potential[ne];

				if (neighbor_domain == current_domain) {
					own_sum += neighbor_potential;
					own_max = max(own_max, neighbor_potential);
					own_count++;
				} else {
					other_sum += neighbor_potential;
					if (neighbor_potential > other_max) {
						other_max = neighbor_potential;
						other_max_domain = neighbor_domain;
					}
					other_count++;
				}
			}

			double new_potential = current_domain == 0 ? 0 :
					(max(prev_pressure[current_domain], 0.) + 0.5 * pot_current + 0.5 * (own_sum - other_sum) / neighbors.size());// - 0.2 * foreign_sum / 4/*neighbors.length*/);// - (foreign_count ? foreign_sum / foreign_count : 0);

			potential[i] = new_potential;

			if (other_count == 0){
				continue;
			}
			uint32_t new_domain = own_max >= other_max ? current_domain : other_max_domain;
			if (new_domain == 0) // Never let vacuum grow
				new_domain = current_domain;
			if (new_domain != current_domain) {
				pressure[new_domain]--;
				pressure[current_domain]++;
				partitioning[i] = new_domain;
				alterations++;
			}
		}
		prev_potential = potential;
		prev_partitioning = partitioning;
		prev_pressure = pressure;
		iteration++;
		return alterations;
	}

	double sdev() {
		valarray<double> pot_peak(n_partitions);
		pot_peak = 0;
		for (uint32_t i = 0; i < n_nodes; i++) {
			if (partitioning[i] == 0)
				continue;
			uint32_t domain = partitioning[i];
			double pot = potential[i];
			if (pot_peak[domain-1] < pot)
				pot_peak[domain-1] = pot;
		}
		double avg = pot_peak.sum() / pot_peak.size();
		return sqrt(1./(pot_peak.size() - 1) * pow(pot_peak - avg, 2).sum()) / (avg == 0 ? 1 : avg);
	}

};


};




#endif /* WPART_HPP_ */
