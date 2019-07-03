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

struct CAPart {
	valarray<double> pressure;		// Contains imbalances of each domain. Positive VAL at index X means: Domain X need to grow VAL nodes
	valarray<double> desired;		// Contains calculated n nodes from weights

	valarray<uint32_t> border_count; // Contains number of cells for each partition that have a foreign neighbor in graph
	size_t	n_nodes, n_partitions;

	valarray<double> potential;

	// INPUT
	vector<vector<int32_t> > graph;	// Connections between nodes, first vector's size determines N nodes.
	valarray<double> weights;		// Size of that vector determines how many partitions will be created. Zero weight means auto-calc weight from remainder
	valarray<size_t> seeds;			// Node indices of seed points. Must be same size as weights or seeds will be auto-distributed.

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

	// preseeded expects a domain decomposition present in partitioning. It will adjust that according to weights.
	void init(bool preseeded) {
		assert(n_nodes > 10);
		assert(n_partitions >= 2);
		assert((n_partitions * 10) < n_nodes);

		if (preseeded) {
			assert((partitioning >= 0).sum());
			assert((partitioning < (int)n_partitions).sum());
		} else {
			partitioning = -1;
			// Place seed.
			if (seeds.size() != n_partitions) {
				// Auto-seed
				seeds.resize(n_partitions);
				uint32_t step = n_nodes/ n_partitions - 1;
				for (int i = 0; i < n_partitions; i++)	// distribute linear within nodes
					seeds[i]  = step * i + rand() % step;
			}
			for (int i = 0; i < n_partitions; i++)
				partitioning[seeds[i]] = i;
		}
		potential[partitioning >= 0] = 1;

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
		if (preseeded) {	// Adjust pressure based on current partitioning
			for (int i = 0; i < n_nodes; i++)
				pressure[partitioning[i]]--;
		} else
			pressure -= 1;  // Because one seed exists per domain
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
			double pot_current = potential[i];

			for (unsigned i_n = 0; i_n < neighbors.size(); i_n++) {
				uint32_t ne = neighbors[i_n];
				int32_t neighbor_domain = partitioning[ne];
				double neighbor_potential = potential[ne];

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
	        bool endangered = false;
			if (current_domain >= 0) {
				double power = (pressure[current_domain] / n_nodes) / (weights[current_domain]);
				if (power < 0)
					power = 0;
				double pot_new = power + 0.5 * pot_current + 0.5 * (own_sum - other_sum) / neighbors.size();
				if (pot_new < 0)
					pot_new = 0;

				potential[i] = pot_new;
				endangered = (desired[current_domain] - pressure[current_domain]) < (0.5 * desired[current_domain]);
			}

	        if (other_count == 0 || endangered) { // Shortcut
	          continue;
	        }

	        if (current_domain >= 0)
	        	border_count[current_domain]++; // Just because 'others' are neighboring

	        if (own_max < other_max) { // We alter domain for node i
		        if (current_domain >= 0) pressure[current_domain]++;
		        pressure[other_max_domain]--;

		        partitioning[i] = other_max_domain;
		        alterations++;
	        }
		}
		iteration++;
		return alterations;
	}

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

	// Takes current peak-pot locations as seeds and resets
	void restart() {
		seeds = peaks();
		init(false);

	}
};
};

#endif /* WPART_HPP_ */
