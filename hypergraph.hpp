#include "graph.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <functional>
#include <fstream>
#include <cstring>
#include <random>
#include <climits>
#include <omp.h>

using namespace std;

typedef pair<int,int> ii;
typedef pair<float,float> ff;
typedef vector<int> vi;
typedef vector<bool> vb;
typedef vector<float> vf;
typedef vector<ii> vii;
typedef vector<vi> vvi;
typedef vector<unordered_set<int>> vus;
typedef vector<vii> vvii;

// compute approximation of mitigation of S_M using 30K monte carlo simulations
double estimateMitigation(Graph &g, int n, vector<int> &seeds, ff &global_reward_data) 
{
	long long int global_running_total = 0;
	global_reward_data.first = 0.0;
	global_reward_data.second = 0.0;
	long long int num = 30000;

	#pragma omp parallel 
	{
		ii reward_data;
		vb visit(n,false);
		vb seed(n,false);
		vi visit_index(n,0);
		vi aw_length(n,-1);
		vi aw_close(n,0);
		vb adoption(n,false);
		vb fr_visit(n,false);
		vi fr_index(n,0);
		vb fake_reachable(n,false);
		vvi visit_neighbours(n);
		vvi parents(n);
		vvi parent_arrivals(n);
		vvi parent_permutation(n);
		priority_queue<ii, vii, greater<ii> > pq;

		for (unsigned int i = 0; i < seeds.size(); i++) {
			seed[seeds[i]] = true;
		}

		long long int running_total = 0;
		long long int ones = 0;
		long long int twos = 0;
		#pragma omp for
		for (int i = 0; i < num; i++) {
	    	running_total += g.generateMitigationSample(seeds, seed, reward_data, visit, visit_index, aw_length, aw_close, adoption, fr_visit, fr_index, fake_reachable, visit_neighbours, parents, parent_arrivals, parent_permutation, pq);
			ones += reward_data.first;
			twos += reward_data.second;
		}

		#pragma omp critical
		{
			global_running_total += running_total;
			global_reward_data.first += ones;
			global_reward_data.second += twos;
		}
	}

	global_reward_data.first /= global_running_total;
	global_reward_data.second *= 2;
	global_reward_data.second /= global_running_total;

	return (float)global_running_total / num;
}

/*
* generate hyperedges in parallel following TCIC model
*/
void addHyperedgeParallel(Graph &g, HyperGraph &hg, long long num, bool sa_upper, bool rs)
{
	int numNodes = g.getNumNodes();
	bool empty = false;

	vvii hyperedges;

	#pragma omp parallel 
	{
		ii root_data;
		ii phase_one_traversal_data;
		ii phase_two_traversal_data;
		vb phase_one_visit(numNodes,false);
		vb phase_two_visit(numNodes,false);
		vb dead_visit(numNodes,false);
		vb phase_three_visit(numNodes,false);
		vi phase_one_visit_index(numNodes,0);
		vi phase_two_visit_index(numNodes,0);
		vi dead_visit_index(numNodes,0);
		vi dist(numNodes,INT_MAX);
		vi delayed_dist(numNodes,INT_MAX);
		vi aw_length(numNodes,-1);
		vb overlap(numNodes,false);
		vi tb_nodes(numNodes,0);
		vi tb_index(numNodes,0);
		vb adopt_fake(numNodes,false);
		vb adopt_true(numNodes,false);
		vvi visit_neighbours(numNodes);
		vvi visit_neighbours_meet_len(numNodes);
		vvi parent_permutation(numNodes);
		vvi phase_three_parents(numNodes);
		vvi phase_three_parent_arrivals(numNodes);
		vvi dead_parents(numNodes);
		vii hyperedge;
		vvii hyperedges_priv;
		priority_queue<ii, vii, greater<ii> > pq;

		#pragma omp for
		for (int i = 0; i < num; i++) {
			if (rs) {
				empty = hg.phaseOneRS(g, root_data, phase_one_traversal_data, pq, phase_one_visit, phase_one_visit_index, dead_visit, dead_visit_index, dist, aw_length, dead_parents);
			} else {
				empty = hg.phaseOne(g, root_data, phase_one_traversal_data, pq, phase_one_visit, phase_one_visit_index, dead_visit, dead_visit_index, dist, aw_length, dead_parents);
			}
	    	if (!empty) {
				hg.phaseTwo(g, root_data, phase_two_traversal_data, pq, phase_one_visit, parent_permutation, phase_two_visit, phase_two_visit_index, 
								delayed_dist, aw_length, overlap, tb_nodes, visit_neighbours, visit_neighbours_meet_len, dead_parents, hyperedge, sa_upper);
				if (phase_two_traversal_data.first > 0) hg.phaseThree(g, root_data.first, phase_two_traversal_data.first, pq, tb_nodes, aw_length, visit_neighbours, visit_neighbours_meet_len, parent_permutation, 
								phase_three_parents, phase_three_parent_arrivals, adopt_fake, adopt_true, phase_three_visit, tb_index, hyperedge, sa_upper);
				hyperedges_priv.push_back(hyperedge);
			}
			hg.reset(g, phase_one_traversal_data, phase_two_traversal_data, phase_one_visit_index, phase_one_visit, dead_visit_index, dead_visit, phase_two_visit_index, visit_neighbours, 
							visit_neighbours_meet_len, parent_permutation, aw_length, dead_parents, hyperedge);
		}

		#pragma omp critical
		hyperedges.insert(hyperedges.end(), hyperedges_priv.begin(), hyperedges_priv.end());
	}

	for (unsigned int i = 0; i < hyperedges.size(); i++) {
		hg.addEdge(hyperedges[i]);
	}
}

/*
* generate hyperedges in serial following TCIC model
*/
void addHyperedgeSerial(Graph &g, HyperGraph &hg, long long num, bool sa_upper, bool rs)
{	
	int numNodes = g.getNumNodes();
	long long iter = 0;
	bool empty = false;

	ii root_data;
	ii phase_one_traversal_data;
	ii phase_two_traversal_data;
	vb phase_one_visit(numNodes,false);
	vb phase_two_visit(numNodes,false);
	vb dead_visit(numNodes,false);
	vb phase_three_visit(numNodes,false);
	vi phase_one_visit_index(numNodes,0);
	vi phase_two_visit_index(numNodes,0);
	vi dead_visit_index(numNodes,0);
	vi dist(numNodes,INT_MAX);
	vi delayed_dist(numNodes,INT_MAX);
	vi aw_length(numNodes,-1);
	vb overlap(numNodes,false);
	vi tb_nodes(numNodes,0);
	vi tb_index(numNodes,0);
	vb adopt_fake(numNodes,false);
	vb adopt_true(numNodes,false);
	vvi visit_neighbours(numNodes);
	vvi visit_neighbours_meet_len(numNodes);
	vvi parent_permutation(numNodes);
	vvi phase_three_parents(numNodes);
	vvi phase_three_parent_arrivals(numNodes);
	vvi dead_parents(numNodes);
	vii hyperedge;
	priority_queue<ii, vii, greater<ii> > pq;

	while (iter < num) {
    	if (rs) {
			empty = hg.phaseOneRS(g, root_data, phase_one_traversal_data, pq, phase_one_visit, phase_one_visit_index, dead_visit, dead_visit_index, dist, aw_length, dead_parents);
		} else {
			empty = hg.phaseOne(g, root_data, phase_one_traversal_data, pq, phase_one_visit, phase_one_visit_index, dead_visit, dead_visit_index, dist, aw_length, dead_parents);
		}
    	if (!empty) {
			hg.phaseTwo(g, root_data, phase_two_traversal_data, pq, phase_one_visit, parent_permutation, phase_two_visit, phase_two_visit_index, 
							delayed_dist, aw_length, overlap, tb_nodes, visit_neighbours, visit_neighbours_meet_len, dead_parents, hyperedge, sa_upper);
			if (phase_two_traversal_data.first > 0) hg.phaseThree(g, root_data.first, phase_two_traversal_data.first, pq, tb_nodes, aw_length, visit_neighbours, visit_neighbours_meet_len, parent_permutation, 
							phase_three_parents, phase_three_parent_arrivals, adopt_fake, adopt_true, phase_three_visit, tb_index, hyperedge, sa_upper);
			hg.addEdge(hyperedge);
		}
		hg.reset(g, phase_one_traversal_data, phase_two_traversal_data, phase_one_visit_index, phase_one_visit, dead_visit_index, dead_visit, phase_two_visit_index, visit_neighbours, 
						visit_neighbours_meet_len, parent_permutation, aw_length, dead_parents, hyperedge);
		iter++;
	}
}

/*
* linear pass over coverages to find node with maximum marginal coverage
* also maintains top k marginals for computation of improved upper bound
*/
int getMaxIndex(int n, vector<int> &node_weight, vector<int> &k_max_mc) {
	int max_ind = -1;
	int max_cov = 0;

	for (int i = 0; i < n; i++) {
		if (node_weight[i] > max_cov) {
			max_ind = i;
			max_cov = node_weight[i];
		}
		if (node_weight[i] > k_max_mc[0]) {
			k_max_mc[0] = node_weight[i];
			sort(k_max_mc.begin(), k_max_mc.end());
		}
	}

	return max_ind;
}

/*
* greedy algorithm for weighted max cover over collection of RDR sets w/ improved UB computation
*/
int buildSeedSet(Graph &g, HyperGraph &hg, unsigned int n, unsigned int k, vector<int> &seeds)
{	
	unsigned int i, j;
	int diff, coverage_weight, cur_cov_ub, max_index;
	bool alive;
	vector<pair<int,int> > edge_list, node_list;
	vector<int> k_max_mc(k,0);

	vector<int> node_weight(n,0);
	for (i = 0; i < n; i++) {
		node_weight[i] = hg.getNodeWeight(i);
	}

	int cur_coverage = 0;
	int improved_cov_ub = INT_MAX;
	long long numEdge = hg.getNumEdge();

	// check if an edge is removed
	vector<bool> edge_removed(numEdge, false);
	
	unsigned int cur_seed = 0;
	// building each seed at a time
	while(cur_seed < k) {
		max_index = getMaxIndex(n, node_weight, k_max_mc);
		if (max_index == -1) break; // all sets have been covered 

		cur_cov_ub = cur_coverage;
		for (i = 0; i < k; i++) {
			cur_cov_ub += k_max_mc[i];
			k_max_mc[i] = 0; // reset for next iteration
		}
		if (cur_cov_ub < improved_cov_ub) improved_cov_ub = cur_cov_ub;

		seeds.push_back(max_index);
		cur_coverage += node_weight[max_index];
		edge_list = hg.getNode(max_index);
		for (i = 0; i < edge_list.size(); i++) {
			if (edge_removed[edge_list[i].first]) continue;
			alive = false;
			coverage_weight = edge_list[i].second;
			node_list = hg.getEdge(edge_list[i].first);
			for (j = 0; j < node_list.size(); j++) {
				if (node_list[j].second > 0) { // only process those nodes with a non-zero marginal for this hyperedge
					diff = node_list[j].second - coverage_weight;
					if (diff <= 0) {
						node_weight[node_list[j].first] -= node_list[j].second;
						node_list[j].second = 0;
					} else {
						node_weight[node_list[j].first] -= diff;
						node_list[j].second = diff;
					}
					if (!alive && node_list[j].second > 0) alive = true; // still marginal gain achievable from this hyperedge
				}
			}
			if (!alive) edge_removed[edge_list[i].first] = true;
		}
		cur_seed++;
	}

	getMaxIndex(n, node_weight, k_max_mc);
	cur_cov_ub = cur_coverage;
	for (i = 0; i < k; i++) {
		cur_cov_ub += k_max_mc[i];
	}
	if (cur_cov_ub < improved_cov_ub) improved_cov_ub = cur_cov_ub;

	return improved_cov_ub;
}

/*
* linear pass over coverages to find node with maximum marginal coverage
*/
int getMaxIndexBaseline(int n, vector<int> &node_weight) {
	int max_ind = -1;
	int max_cov = 0;

	for (int i = 0; i < n; i++) {
		if (node_weight[i] > max_cov) {
			max_ind = i;
			max_cov = node_weight[i];
		}
	}

	return max_ind;
}

/*
* greedy algorithm for weighted max cover over collection of RDR sets w/o improved UB computation
*/
int buildSeedSetBaseline(Graph &g, HyperGraph &hg, unsigned int n, unsigned int k, vector<int> &seeds)
{	
	unsigned int i, j;
	int diff, coverage_weight, max_index;
	bool alive;
	vector<pair<int,int> > edge_list, node_list;

	vector<int> node_weight(n,0);
	for (i = 0; i < n; i++) {
		node_weight[i] = hg.getNodeWeight(i);
	}

	int coverage = 0;
	long long numEdge = hg.getNumEdge();

	// check if an edge is removed
	vector<bool> edge_removed(numEdge, false);
	
	unsigned int cur_seed = 0;
	// building each seed at a time
	while(cur_seed < k) {
		max_index = getMaxIndexBaseline(n, node_weight);
		if (max_index == -1) break; // all sets have been covered 
		seeds.push_back(max_index);
		coverage += node_weight[max_index];
		edge_list = hg.getNode(max_index);
		for (i = 0; i < edge_list.size(); i++) {
			if (edge_removed[edge_list[i].first]) continue;
			alive = false;
			coverage_weight = edge_list[i].second;
			node_list = hg.getEdge(edge_list[i].first);
			for (j = 0; j < node_list.size(); j++) {
				if (node_list[j].second > 0) { // only process those nodes with a non-zero marginal for this hyperedge
					diff = node_list[j].second - coverage_weight;
					if (diff <= 0) {
						node_weight[node_list[j].first] -= node_list[j].second;
						node_list[j].second = 0;
					} else {
						node_weight[node_list[j].first] -= diff;
						node_list[j].second = diff;
					}
					if (!alive && node_list[j].second > 0) alive = true; // still marginal gain achievable from this hyperedge
				}
			}
			if (!alive) edge_removed[edge_list[i].first] = true;
		}
		cur_seed++;
	}

	return coverage;
}
