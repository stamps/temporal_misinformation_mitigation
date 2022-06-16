#ifndef _GRAPH_H
#define _GRAPH_H
#include <vector>
#include <unordered_set>
#include <random>
#include <functional>
#include <queue>
#include "sfmt/SFMT.h"

typedef uint32_t UI;

class Graph
{
	friend class HyperGraph;
	private:
		UI UI_MAX = 4294967295U;
		// number of nodes
		unsigned int num_nodes;
		// number of edges
		unsigned int num_edges;
		// fake seeds
		std::vector<int> fake_seeds_index;
		unsigned int num_fake_seeds;
		std::vector<bool> fake_seed;
		// dynamics parameters
		int aw;
		int ml;
		float rp;
		int br;
		bool ego;
		// adjacency lists
		std::vector<std::vector<int> > adj_list;
		std::vector<std::vector<int> > rev_adj_list;
		std::vector<int> node_deg;
		std::vector<int> rev_node_deg;
		std::vector<std::vector<float> > weights;
		std::vector<std::vector<float> > rev_weights;
		sfmt_t sfmt_seed;
	
	public:
		Graph(int aw, int ml, float rp, int br, bool ego);
		// get a vector of neighbours of node u in G
		const std::vector<int> & getOutNeighbours(int u) const;
		// get a vector of neighbours of node u in G^T
		const std::vector<int> & getInNeighbours(int u) const;
		// return weights of neighbours of node u in G
		const std::vector<float> & getOutWeights(int u) const;
		// return weights of neighbours of node u in G^T
		const std::vector<float> & getInWeights(int u) const;
		// return a vector of the fake seeds of G
		const std::vector<int> & getFakeSeeds() const;
		// get out degree of node u
		unsigned int getOutDegree(int u) const;
		// get in degree of node u
		unsigned int getInDegree(int u) const;
		// get size of the graph
		unsigned int getNumNodes() const;
		// get number of edges
		unsigned int getNumEdges() const;
		// get number of fake seeds
		unsigned int getNumFakeSeeds() const;
		// determine if node v is a fake seed
		bool isFakeSeed(int v) const;
		// read graph from a file
		void readGraph(const char * filename, bool fixed, float w);
		// read fake seeds from a file
		void readFakeSeeds(const char* filename);
		// generate a random activation window length
		int getActivationWindow();
		// generate a random meeting length
		int generateMeetingLength(int u);
		// compute mitigation lower bound via depth 1 MIA
		double computeMitigationLowerBound(unsigned int n, unsigned int k);
		// generate a single forward monte carlo estimate of the influence of input seed
		int generateInfluenceSample(std::vector<bool> &visit, std::vector<int> &visit_index, int root);
		// generate a single forward monte carlo estimate of the influence of F
		int generateFakeInfluenceSample(std::vector<bool> &visit, std::vector<int> &visit_index);
		// generate a single forward monte carlo estimate of the mitigation of S_M
		int generateMitigationSample(std::vector<int> &seeds, std::vector<bool> &seed, std::pair<int,int> &reward_data, std::vector<bool> &visit, std::vector<int> &visit_index, std::vector<int> &aw_length, std::vector<int> &aw_close, 
			std::vector<bool> &adoption, std::vector<bool> &fr_visit, std::vector<int> &fr_index, std::vector<bool> &fake_reachable, std::vector<std::vector<int>> &visit_neighbours, std::vector<std::vector<int>> &parents, std::vector<std::vector<int>> &parent_arrivals, 
			std::vector<std::vector<int>> &parent_permutation, std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int,int>> > &pq);
};

class HyperGraph
{
	private:
		// store the hyperedges that a node is incident on together with the corresponding coverage weight 
		std::vector<std::vector<std::pair<int,int> > > node_hyperedges;
		// store hyperedges
		std::vector<std::vector<std::pair<int,int> > > hyperedges;
		// store the weight of hyperedges that a node is covered by
		std::vector<int> node_hyperedge_weight;
		//unsigned int cur_hyperedge;
		sfmt_t sfmtSeed;
		
	public:
		HyperGraph(unsigned int n);
		void addEdge(std::vector<std::pair<int,int> > & edge);
		const std::vector<std::pair<int,int> > & getEdge(int e) const;
		const std::vector<std::pair<int,int> > & getNode(int u) const;
		int getNodeWeight(int u) const;
        int getNumEdge() const;
		void clearEdges();

		bool phaseOneRS(Graph &g, std::pair<int,int> &root_data, std::pair<int,int> &traversal_data, 
			std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int,int>> > &pq,
			std::vector<bool> &visit, std::vector<int> &visit_index, std::vector<bool> &dead_visit, std::vector<int> &dead_visit_index, 
			std::vector<int> &dist, std::vector<int> &aw_length, std::vector<std::vector<int>> &dead_parents);
		bool phaseOne(Graph &g, std::pair<int,int> &root_data, std::pair<int,int> &traversal_data, 
			std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int,int>> > &pq,
			std::vector<bool> &visit, std::vector<int> &visit_index, std::vector<bool> &dead_visit, std::vector<int> &dead_visit_index, 
			std::vector<int> &dist, std::vector<int> &aw_length, std::vector<std::vector<int>> &dead_parents);
		void phaseTwo(Graph &g, std::pair<int,int> &root_data, std::pair<int,int> &traversal_data, 
			std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int,int>> > &pq, 
			std::vector<bool> &phase_one_visit, std::vector<std::vector<int>> &parent_permutation, std::vector<bool> &visit, 
			std::vector<int> &visit_index, std::vector<int> &delayed_dist, std::vector<int> &aw_length, std::vector<bool> &overlap, 
			std::vector<int> &tb_nodes, std::vector<std::vector<int>> &visit_neighbours, std::vector<std::vector<int>> &visit_neighbours_meet_len, 
			std::vector<std::vector<int>> &dead_parents, std::vector<std::pair<int,int>> &hyperedge, bool sa_upper);
		void phaseThree(Graph &g, int root, int num_tb, std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int,int>> > &pq, 
			std::vector<int> &tb_nodes, std::vector<int> &aw_length, std::vector<std::vector<int>> &visit_neighbours, std::vector<std::vector<int>> &visit_neighbours_meet_len, 
			std::vector<std::vector<int>> &parent_permutation, std::vector<std::vector<int>> &parents, std::vector<std::vector<int>> &parent_arrivals, 
			std::vector<bool> &adopt_fake, std::vector<bool> &adopt_true, std::vector<bool> &visit, std::vector<int> &tb_index, std::vector<std::pair<int,int>> &hyperedge, bool sa_upper);
		void reset(Graph &g, std::pair<int,int> &phase_one_traversal_data, std::pair<int,int> &phase_two_traversal_data, std::vector<int> &phase_one_visit_index, 
			std::vector<bool> &phase_one_visit, std::vector<int> &dead_visit_index, std::vector<bool> &dead_visit, std::vector<int> &phase_two_visit_index, std::vector<std::vector<int>> &visit_neighbours, 
			std::vector<std::vector<int>> &visit_neighbours_meet_len, std::vector<std::vector<int>> &parent_permutation, std::vector<int> &aw_length, 
			std::vector<std::vector<int>> &dead_parents, std::vector<std::pair<int,int>> &hyperedge);
};

float getCurrentMemoryUsage();

#endif
