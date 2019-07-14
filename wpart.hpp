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
#include <cassert>
#include <valarray>
#include <math.h>
#include <stdint.h>

using namespace std;

namespace wpart {

#define likely(x)	__builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

struct CAPart {
	valarray<double> pressure;		// Contains imbalances of each domain. Positive VAL at index X means: Domain X need to grow VAL nodes
	valarray<double> desired;		// Contains calculated n nodes from weights

	valarray<uint32_t> border_count; // Contains number of cells for each partition that have a foreign neighbor in graph
	size_t	n_nodes, n_partitions;

	valarray<double> potential;

	// INPUT
	vector<vector<int32_t> > graph;	// Connections between nodes, first vector's size determines N nodes.
	valarray<double> weights;		// Size of that vector determines how many partitions will be created. Zero weight means auto-calc weight from remainder

	// OUTPUT ( / INPUT if pre-seeded)
	valarray<int32_t> partitioning; // Domain index for each node

	uint32_t iteration;

	CAPart(uint32_t n_nodes_, uint32_t n_partitions_) {
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
	};

	// preseeded (or pre-partitioned) expects either a domain decomposition present
	// in the partitioning valarray; it will adjust / update that according to weights.
	// Another possibility is that only seeds are present (single index), the rest
	// left to void (-1).
	void init(bool preseeded) {
		assert(n_nodes > 10);
		assert(n_partitions >= 2);
		assert((n_partitions * 10) < n_nodes);

		if (preseeded) {
			assert((partitioning >= 0).sum());	// There must be a
			assert(partitioning.max() < (int)n_partitions);
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
			vector<int32_t> &neighbors = graph[i];
			int own_count = 0;
			double own_sum = 0, own_max = 0;
			int other_count = 0, other_max_domain = 0;
			double other_sum = 0, other_max = 0;

			int current_domain = partitioning[i];
			double current_potential = potential[i];

			for (unsigned i_n = 0; i_n < neighbors.size(); i_n++) {
				uint32_t ne = neighbors[i_n];
				int32_t neighbor_domain = partitioning[ne];
				double neighbor_potential = potential[ne];

				if (likely(neighbor_domain == current_domain)) {
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

			if (likely(current_domain >= 0)) {	// Update potential field
				double power = max((pressure[current_domain] / n_nodes) / (weights[current_domain]), 0.);
				potential[i] = max(power + 0.5 * current_potential + 0.5 * (own_sum - other_sum) / neighbors.size(), 0.);

				// Endangered domain?
				if (unlikely((desired[current_domain] - pressure[current_domain]) < (0.5 * desired[current_domain])))
					continue;
			}

	        if (likely(other_count == 0)) // Shortcut
	          continue;

	        if (likely(current_domain >= 0))
	        	border_count[current_domain]++; // Just because 'others' are neighboring

	        if (own_max < other_max) { // We alter domain for node i, likely-hood of branch 50/50.
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
};

#endif /* WPART_HPP_ */
