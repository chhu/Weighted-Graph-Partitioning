/*
 * wpart.hpp
 *
 *	Weighted graph partitioning using a cellular automaton rule.
 *
 *
 *  Created on: Jun 18, 2019
 *      Author: C Huettig (christian.huettig@dlr.de)
 *
 */

#ifndef WPART_HPP_
#define WPART_HPP_

#include <vector>
#include <set>
#include <cassert>
#include <valarray>
#include <math.h>
#include <stdint.h>
#include <iostream>
#include <fstream>

using namespace std;

namespace wpart {

#define likely(x)	__builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

struct CAPart {
	valarray<double> pressure;		// Contains imbalance of each domain. VAL at index X means: Domain X need to grow / shrink VAL nodes
	valarray<double> desired;		// Contains calculated node counts from weights

	valarray<uint32_t> border_count; // Contains number of cells for each partition that have a foreign neighbor in graph
	size_t	n_nodes, n_partitions;

	valarray<double> potential;

	// INPUT
	vector<vector<int32_t> > graph;	// Connections between nodes, first vector's size determines N nodes.
	valarray<double> weights;		// Size of that vector determines how many partitions will be created. Zero value weight means auto-calc weight from remainder

	// OUTPUT ( / INPUT if pre-seeded)
	valarray<int32_t> partitioning; // Domain index for each node

	uint32_t iteration;
	bool	no_alterations;
	double	bias;

	CAPart(uint32_t n_nodes_, uint32_t n_partitions_, const char* grf_name = NULL) {
		ifstream grf;
		if (grf_name) { // Open grf file and read n_nodes
			grf.open(grf_name, ifstream::in);
			grf >> n_nodes_;grf >> n_nodes_;
		}
		n_nodes = n_nodes_;
		n_partitions = n_partitions_;
		graph.resize(n_nodes);
		weights.resize(n_partitions_, 0.);
		iteration = 0;
		pressure.resize(n_partitions, 0);
		desired.resize(n_partitions, 0);
		border_count.resize(n_partitions, 0);
		potential.resize(n_nodes, -1.);
		partitioning.resize(n_nodes, -1);
		no_alterations = false;
		bias = 0;
		split = coarse = NULL;
		if (grf.is_open()) { // Read graph from file
			int dummy, entries, i = 0;
			grf >> dummy;grf >> dummy;grf >> dummy; // Values on graph info should probably not be ignored
			assert(dummy == 0);
			while (grf >> entries) {
				assert( i < n_nodes);
				for (int j = 0; j < entries; j++) {
					grf >> dummy;
					graph[i].push_back(dummy);
					assert(dummy < n_nodes);
				}
				i++;
			}
			grf.close();
		}
	};
/*
	~CAPart() {
		if (split != NULL) delete split;
		if (coarse != NULL) delete coarse;
	}
*/
	// preseeded (or pre-partitioned) expects either a domain decomposition present
	// in the partitioning valarray; it will adjust / update that according to weights.
	// Another possibility is that only seeds are present (single index), the rest
	// left to void value (-1).
	void init(bool preseeded) {
		assert(n_nodes > 10);
		assert(n_partitions >= 2);
		assert((n_partitions * 10) < n_nodes);

		if (preseeded) {
			assert(partitioning.min() == 0);	// There must be a partitioning present
			assert(partitioning.max() == (int)n_partitions - 1);
		} else {
			partitioning = -1;
			// Place seeds
			uint32_t step = n_nodes/ n_partitions - 1;
			for (int i = 0; i < n_partitions; i++)	// distribute with random offset within linear (1D) partition
				partitioning[step * i + rand() % step] = i;
		}
		potential[partitioning >= 0] = 1 / ((partitioning == -1).sum() + 1);

		// Translate weights to node counts
		double total_weights = weights.sum();
		assert(total_weights <= 1 && total_weights >= 0);

		double remainder = (1 - total_weights);

		if (remainder > 1e-4) { // Auto balance, fill zero weights with what remains
			double unset_count = 0;
			for (int i = 0; i < n_partitions; i++) unset_count += weights[i] == 0. ? 1 : 0;
			remainder /= unset_count;
			weights[weights == 0.] = remainder;
		}
		bias = 1/(weights.min() * n_nodes);
		for (int i = 0; i < n_partitions; i++) {
			pressure[i] = round(n_nodes * weights[i]);
		}
		pressure[0] += n_nodes - pressure.sum(); // Total must be equal to n_nodes, rounding errors get added to first domain
		desired = pressure;
		// Adjust pressure according to actual partition or seeds.
		for (int i = 0; i < n_nodes; i++)
			if (partitioning[i] >= 0)
				pressure[partitioning[i]]--;
	}

	uint32_t step() {
		uint32_t alterations = 0;
		border_count = 0;

		for (size_t i = 0; i < n_nodes; i++) {
			const vector<int32_t> &neighbors = graph[i];
			double own_sum = 0;
			int other_count = 0, other_max_domain = 0;
			double other_sum = 0, other_max = -1e8;

			int current_domain = partitioning[i];

			// Gather sum of neighbor potential, split in "own" and "other"
			for (unsigned i_n = 0; i_n < neighbors.size(); i_n++) {
				const uint32_t &ne = neighbors[i_n];
				const int32_t &neighbor_domain = partitioning[ne];
				const double &neighbor_potential = potential[ne];

				if (likely(neighbor_domain == current_domain)) {
					own_sum += neighbor_potential;
				} else {
					other_sum += neighbor_potential;
					if (neighbor_potential > other_max) {
						other_max = neighbor_potential; // Remember domain of strongest foreign potential
						other_max_domain = neighbor_domain;
					}
					other_count++;
				}
			}

			if (likely(current_domain >= 0)) {	// Update potential field
				double power = (pressure[current_domain] / n_nodes) / (weights[current_domain]);
				if (power < 0)
					power = 0;
				else
					power += bias;

				potential[i] = power + (own_sum - other_sum) / neighbors.size();

				// Endangered domain?
				if (unlikely((desired[current_domain] - pressure[current_domain]) < (0.5 * desired[current_domain])))
					continue;
			}

	        if (likely(other_count == 0)) // Shortcut, surrounded by own domain, nothing more to do.
	        	continue;

	        if (likely(current_domain >= 0))	// Mark as halo cell
	        	border_count[current_domain]++; // Just because 'others' are neighboring

	        if (potential[i] < 0 && !no_alterations) { // We alter domain for node i, likely-hood of branch 50/50.
		        if (likely(current_domain >= 0))
		        	pressure[current_domain]++;
		        pressure[other_max_domain]--;

		        partitioning[i] = other_max_domain;
		        alterations++;
	        }
		}
		iteration++;
		return alterations;
	}

	// For recursive MG operations
	struct CAPart* split;
	struct CAPart* coarse;

	void createSplit(unsigned nodes_per_partition) {
		assert(nodes_per_partition > 10);
		unsigned n_partitions_split = n_nodes / nodes_per_partition;
		split = new CAPart(n_nodes, n_partitions_split);
		split->graph = graph;
		split->init(false);

		while (split->pressure.max() > 0.1 * nodes_per_partition)
			split->step();
		coarse = split->compact(n_partitions);
		coarse->weights = weights;
		coarse->init(false);
	}

	// This method assumes there is a fully partitioned graph present (no voids)
	CAPart* compact(uint32_t n_partitions_) {
		assert(n_partitions_ < n_partitions / 10);
		CAPart* result = new CAPart(n_partitions, n_partitions_);
		// Build a graph that builds upon current partitions and their neighbors
		vector<set<int32_t> > compact(n_partitions);
		for (int i = 0; i < n_nodes; i++) {
			const vector<int32_t> &neighbors = graph[i];
			int current_domain = partitioning[i];
			for (unsigned i_n = 0; i_n < neighbors.size(); i_n++) {
				int32_t neighbor_domain = partitioning[neighbors[i_n]];
				if (neighbor_domain != current_domain)
					compact[current_domain].insert(neighbor_domain);
			}
		}
		// Convert vector<set> to vector<vector>
		for (int i = 0; i < n_partitions; i++)
			for (set<int32_t>::iterator it = compact[i].begin(); it != compact[i].end(); ++it)
				result->graph[i].push_back(*it);

		return result;
	}

	// Projects mapping from this graph onto target using split as map
	void project(CAPart* split, CAPart* target) {
		assert(split->n_nodes == target->n_nodes);
		assert(target->n_partitions == n_partitions);
		for (int i = 0; i < target->n_nodes; i++)
			target->partitioning[i] = partitioning[split->partitioning[i]];
		target->init(true);	// Forces recalc of pressure
		for (int i = 0; i < target->n_nodes; i++) // Project potential as well
			target->potential[i] = potential[split->partitioning[i]]; // TODO scaling is wrong

	}


	// Helpers

	// Return max potential (centroid) positions for domains
	valarray<size_t> peaks() {
		valarray<size_t> result(n_partitions);result = 0;
		valarray<double> pot_peak(0., n_partitions);
		for (uint32_t i = 0; i < n_nodes; i++) {
			if (partitioning[i] == 0)
				continue;
			uint32_t domain = partitioning[i];
			double pot = potential[i];
			if (pot > pot_peak[domain-1]) {
				pot_peak[domain-1] = pot;
				result[domain - 1] = i;
			}
		}
		return result;
	}

	// Takes current peak-pot locations as seed and resets
	void restart() {
		valarray<size_t> pks = peaks();
		partitioning = -1;
		for (int i = 0; i < pks.size(); i++) {
			partitioning[pks[i]] = i;
		}
		init(true);
	}
};

struct CAPartMG {

};

};

#endif /* WPART_HPP_ */
