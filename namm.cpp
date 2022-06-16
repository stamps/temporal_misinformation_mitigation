#include "option.h"
#include "hypergraph.hpp"
#include "graph.h"
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <fstream>
#include <cstring>
#include <random>
#include <omp.h>

using namespace std;

// equation 5
double computeThetaMax(int n, double epsilon, double epsilon_prime, double delta, double delta_prime, double mit_lower_bound, int k) {
	double Delta = delta - delta_prime;
	double log_n_choose_k = lgammal(n+1)-lgammal(k+1)-lgammal(n-k+1);
	double denom = epsilon * (1+epsilon_prime) - 2 * epsilon_prime * (1-1/exp(1));
	double numerator = 8 * n * (1-1/exp(1)) * (3+epsilon_prime) * (log(27.0 / (4 * Delta)) + log_n_choose_k);
	return numerator / (3 * mit_lower_bound * denom * denom);
}

/*
* read fake seeds
*/
float readFakeInfluence(const char* filename)
{
    float fi;
    ifstream in(filename);
    in >> fi;
    in.close();
    return fi;
}

// compute eps-delta approx of INF_F using generalized SRA
double estimateFakeInfluence(Graph &g, int n, double epsilon_prime, double delta_prime) {
	int b = n - g.getNumFakeSeeds();
	double eps = epsilon_prime * (1 - (epsilon_prime * b)/((2 + 2.0/3.0 * epsilon_prime) * log(2.0 / delta_prime) * b));
	long long int gamma = (1 + epsilon_prime) * (2 + 2.0/3.0 * eps) * log(2.0 / delta_prime) * (1.0 / (eps*eps)) * b;
	long long int counter = 0;
	long long int running_total = 0;

	vector<bool> visit(n,false);
	vector<int> visit_index(n,0);

	while (running_total < gamma) {
		running_total += g.generateFakeInfluenceSample(visit, visit_index);
		counter++;
	}

	return (float)running_total / counter;
}

// equation 6
double computeLowerBound(HyperGraph &hg, unsigned int n, unsigned int k, long long num_edge, double epsilon, double delta, double inf_f, vector<int> &seeds) {
	unsigned int i, j;
	int diff, coverage_weight;
	bool alive;
	vector<pair<int,int> > edge_list, node_list;

	vector<int> node_weight(n,0);
	for (i = 0; i < n; i++){
		node_weight[i] = hg.getNodeWeight(i);
	}

	vector<bool> edge_removed(num_edge, false);	
	unsigned int cur_seed = 0;
	int coverage = 0;
	
	while(cur_seed < seeds.size()) {
		coverage += node_weight[seeds[cur_seed]];
		edge_list = hg.getNode(seeds[cur_seed]);
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

	double weighted_coverage = coverage * inf_f / n;
	double term_two = n / ((1+epsilon) * num_edge);
	double root_one = sqrt(weighted_coverage + ((25.0/36.0) * term_one));
	double root_two = sqrt(term_one);
	double term_three = root_one - root_two;
	return (term_three * term_three - (term_one / 36)) * term_two;
}

// equation 7
double computeUpperBound(int n, int k, long long int num_edge, double epsilon, double delta, double inf_f, int cov_ub) {
	double weighted_coverage = cov_ub * inf_f / n;
	double term_one = log(1.0/delta);
	double term_two = n / ((1-epsilon) * num_edge);
	double root_one = sqrt(weighted_coverage + term_one);
	double root_two = sqrt(term_one);
	double term_three = root_one + root_two;
	return term_three * term_three * term_two;
}

// lemma 12
double computeBetaDenom(HyperGraph &hg, unsigned int n, unsigned int k, long long num_edge, double inf_f, vector<int> &seeds) {
	unsigned int i, j;
	int diff, coverage_weight;
	bool alive;
	vector<pair<int,int> > edge_list, node_list;

	vector<int> node_weight(n,0);
	for (i = 0; i < n; i++){
		node_weight[i] = hg.getNodeWeight(i);
	}

	vector<bool> edge_removed(num_edge, false);	
	unsigned int cur_seed = 0;
	int coverage = 0;
	
	while(cur_seed < seeds.size()) {
		coverage += node_weight[seeds[cur_seed]];
		edge_list = hg.getNode(seeds[cur_seed]);
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

	double mit = coverage * inf_f / num_edge;
	return mit;
}

int main(int argc, char ** argv)
{
	srand(time(NULL));
	
	OptionParser op(argc, argv);
	if (!op.validCheck()){
		printf("Parameters error, please check the readme.txt file for correct format!\n");
		return -1;
	}

	char * inFile = op.getPara("-i");
	if (inFile == NULL){
		inFile = (char*)"network";
	}

	char * outFile = op.getPara("-o");
	if (outFile == NULL){
		outFile = (char*)"results.txt";
	}

	char * fakeSeedsFile = op.getPara("-fakeseeds");
	if (fakeSeedsFile == NULL){
		fakeSeedsFile = (char*)"fake.seeds";
	}

	char * fakeInfFile = op.getPara("-fakeinf");
	if (fakeInfFile == NULL){
		fakeInfFile = (char*)"fake.inf";
	}

	char * tmp = op.getPara("-epsilon");
	float epsilon = 0.1;
	if (tmp != NULL){
		epsilon = atof(tmp);
	}

	int aw = 30;
	tmp = op.getPara("-aw");
	if (tmp != NULL){
		aw = atoi(tmp);
	}

	float rp = 0.6;
	tmp = op.getPara("-rp");
	if (tmp != NULL){
		aw = atof(tmp);
	}

	int ml = 6;
	tmp = op.getPara("-ml");
	if (tmp != NULL){
		ml = atoi(tmp);
	}

	int k = 1;
	tmp = op.getPara("-k");
	if (tmp != NULL){
		k = atoi(tmp);
	}

	int br = 200;
	tmp = op.getPara("-br");
	if (tmp != NULL){
		br = atoi(tmp);
	}

	bool nub = false;
	tmp = op.getPara("-nub");
	if (tmp != NULL){
		nub = atoi(tmp);
	}

	bool ego = false;
	tmp = op.getPara("-ego");
	if (tmp != NULL){
		ego = atoi(tmp);
	}

	float perturb = -1.0;
	tmp = op.getPara("-perturb");
	if (tmp != NULL){
		perturb = atof(tmp);
	}
	bool robust = (perturb < 0.0) ? false : true;

	bool cic = false;
	tmp = op.getPara("-cic");
	if (tmp != NULL){
		cic = atoi(tmp);
	}

	float ew = -1.0;
	tmp = op.getPara("-ew");
	if (tmp != NULL){
		ew = atof(tmp);
	}
	bool fixed = (ew < 0.0) ? false : true;

	cout << "\n*******************" << endl;
	cout << "\tSTART" << endl;
	cout << "*******************\n" << endl;

	Graph g_prime(aw, ml, rp, br, ego);
	g_prime.readGraph(inFile, fixed, ew);
	g_prime.readFakeSeeds(fakeSeedsFile);

	if (cic) {
		cout << "Ignoring temporal paramters." << endl;
		rp = 1.0;
		ml = 1;
	}

	if (robust) {
		cout << "Perturbing ground truth temporal paramters." << endl;
		float sd;
		random_device rd;
		mt19937 gen(rd());
		// perturb AW
		sd = aw * perturb;
		normal_distribution<float> d_aw(0.0, sd);
		aw = (int)(round(aw + d_aw(gen)));
		cout << "perturbed aw: " << aw << endl;
		// perturb RP
		sd = rp * perturb;
		normal_distribution<float> d_rp(0.0, sd);
		rp = rp + d_rp(gen);
		cout << "perturbed rp: " << rp << endl;
		// perturb ML
		sd = ml * perturb;
		normal_distribution<float> d_ml(0.0, sd);
		ml = (int)(round(ml + d_ml(gen)));
		cout << "perturbed ml: " << ml << endl;
		// perturb BR
		sd = br * perturb;
		normal_distribution<float> d_br(0.0, sd);
		br = (int)(round(br + d_br(gen)));
		cout << "perturbed br: " << br << endl;
	}

	Graph g(aw, ml, rp, br, ego);
	g.readGraph(inFile, fixed, ew);
	g.readFakeSeeds(fakeSeedsFile);

	int n = g.getNumNodes();

	float delta = 1.0/n;
	tmp = op.getPara("-delta");
    if (tmp != NULL){
    	delta = atof(tmp);
    }

    float precision = 1-1/exp(1);
    tmp = op.getPara("-precision");
    if (tmp != NULL){
    	precision = atof(tmp);
    }
	
	HyperGraph coll_one(n);
	HyperGraph coll_two(n);

    double epsilon_prime = epsilon / 2.0;
    double delta_prime = delta / 9.0;

    double interval;
    double time_pre = 0.0;
    double time_namm = 0.0;
    double time_est_mit = 0.0;
    double start_total = omp_get_wtime();
    double start = omp_get_wtime();

    cout << "\nComputing mitigation lower bound." << endl;
    double mit_lower_bound = g.computeMitigationLowerBound(n, k);
    interval = omp_get_wtime()-start;
    time_pre += interval;
    cout << "Time to compute mitigation lower bound: " << interval << "s" << endl;
    cout << "k = " << k << "\tLB = " << mit_lower_bound << endl; 

    start = omp_get_wtime();
    double inf_f = readFakeInfluence(fakeInfFile);
    interval = omp_get_wtime()-start;
    time_pre += interval;
    cout << "INF_F = " << inf_f << endl;

    cout << "fake seed set: ";
    const vi &fs = g.getFakeSeeds();
	for (unsigned int s = 0; s < g.getNumFakeSeeds(); s++) {
		cout << fs[s] << " ";
	}
	cout << endl;

    double theta_max = computeThetaMax(n, epsilon, epsilon_prime, delta, delta_prime, mit_lower_bound, k);
    double theta_nought = theta_max * (epsilon*epsilon) * mit_lower_bound / n;
    cout << "theta_max = " << theta_max << endl; 
    cout << "theta_0 = " << theta_nought << endl; 
	
	int i_max = ceil(log2(theta_max / theta_nought));
	cout << "i_max = " << i_max << endl; 

	double delta_iter = (delta - delta_prime) / (3.0*i_max);
	cout << "delta_iter = " << delta_iter << endl; 

	vector<int> seeds;
	ff lower_reward_data;
	ff upper_reward_data;
	double lower;
	double upper;
	double alpha;
	int cov_ub;

	long long int cur_samples = (long long) theta_nought;
	int i = 1;

	// LOWER BOUND FOR SA
	cout << "\n*******************" << endl;
	cout << "LOWER BOUND FOR SA" << endl;
	cout << "*******************" << endl;
	start = omp_get_wtime();

	addHyperedgeParallel(g, coll_one, cur_samples, false, false);
	addHyperedgeParallel(g, coll_two, cur_samples, false, false);
	while (true) {
		cout << "\nIteration " << i << endl;
		seeds.clear();
		if (nub) {
			cov_ub = buildSeedSetBaseline(g, coll_one, n, k, seeds);
			cov_ub /= (1-1/exp(1));
		} else {
			cov_ub = buildSeedSet(g, coll_one, n, k, seeds);
		}
		lower = computeLowerBound(coll_two, n, k, cur_samples, epsilon_prime, delta_iter, inf_f, seeds);
		cout << "lower is " << lower << endl;

		upper = computeUpperBound(n, k, cur_samples, epsilon_prime, delta_iter, inf_f, cov_ub);
		cout << "upper is " << upper << endl;
		alpha = lower / upper;
		cout << "alpha is " << alpha << endl;
		if ((alpha >= (precision - epsilon)) || (i == i_max)) {
			break;
		}
		cout << "Generating up to " << 2*cur_samples << " samples " << endl;
		addHyperedgeParallel(g, coll_one, cur_samples, false, false);
		addHyperedgeParallel(g, coll_two, cur_samples, false, false);
		cur_samples *= 2;
		i++;
	}
	interval = omp_get_wtime()-start;
    time_namm += interval;
	cout << "\nTime to execute NAMM: " << interval << "s" << endl;

	cout << "seed set: ";
	for (unsigned int s = 0; s < seeds.size(); s++) {
		cout << seeds[s] << " ";
	}
	cout << endl;

	start = omp_get_wtime();
	double sa_mit_lower = estimateMitigation(g_prime, n, seeds, lower_reward_data);
	interval = omp_get_wtime()-start;
	time_est_mit += interval;
	cout << "\nTime to estimate lower mitigation via MC: " << interval << "s" << endl;
	cout << "lower mitigation is " << sa_mit_lower << endl;
	cout << "lower reward breakdown: " << lower_reward_data.first << " " << lower_reward_data.second << endl;

	// UPPER BOUND FOR SA
	cout << "\n*******************" << endl;
	cout << "UPPER BOUND FOR SA" << endl;
	cout << "*******************" << endl;
	coll_one = HyperGraph(n);
	coll_two = HyperGraph(n);
	cur_samples = (long long) theta_nought;
	i = 1;
	start = omp_get_wtime();

	addHyperedgeParallel(g, coll_one, cur_samples, true, false);
	addHyperedgeParallel(g, coll_two, cur_samples, true, false);
	while (true) {
		cout << "\nIteration " << i << endl;
		seeds.clear();
		if (nub) {
			cov_ub = buildSeedSetBaseline(g, coll_one, n, k, seeds);
			cov_ub /= (1-1/exp(1));
		} else {
			cov_ub = buildSeedSet(g, coll_one, n, k, seeds);
		}
		lower = computeLowerBound(coll_two, n, k, cur_samples, epsilon_prime, delta_iter, inf_f, seeds);
		cout << "lower is " << lower << endl;

		upper = computeUpperBound(n, k, cur_samples, epsilon_prime, delta_iter, inf_f, cov_ub);
		cout << "upper is " << upper << endl;
		alpha = lower / upper;
		cout << "alpha is " << alpha << endl;
		if ((alpha >= (precision - epsilon)) || (i == i_max)) {
			break;
		}
		cout << "Generating up to " << 2*cur_samples << " samples " << endl;
		addHyperedgeParallel(g, coll_one, cur_samples, true, false);
		addHyperedgeParallel(g, coll_two, cur_samples, true, false);
		cur_samples *= 2;
		i++;
	}
	interval = omp_get_wtime()-start;
    time_namm += interval;
	cout << "\nTime to execute NAMM: " << interval << "s" << endl;

	cout << "seed set: ";
	for (unsigned int s = 0; s < seeds.size(); s++) {
		cout << seeds[s] << " ";
	}
	cout << endl;

	start = omp_get_wtime();
	double sa_mit_upper = estimateMitigation(g_prime, n, seeds, upper_reward_data);
	interval = omp_get_wtime()-start;
	time_est_mit += interval;
	cout << "\nTime to estimate upper mitigation via MC: " << interval << "s" << endl;
	cout << "upper mitigation is " << sa_mit_upper << endl;
	cout << "upper reward breakdown: " << upper_reward_data.first << " " << upper_reward_data.second << endl;

	// generate additional RDR sets before computing upper bound estimate for beta calculation
	addHyperedgeParallel(g, coll_two, cur_samples, true, false);
	cur_samples *= 2;
	double beta_denom = computeBetaDenom(coll_two, n, k, cur_samples, inf_f, seeds);
	cout << "beta_denom is " << beta_denom << endl;
	double beta = sa_mit_upper / beta_denom;
	cout << "beta is " << beta << endl;

	cout << "\n*******************" << endl;
	cout << "\tSA RESULT" << endl;
	cout << "*******************" << endl;
	double time_total = omp_get_wtime()-start_total;
	cout << "\nTotal Time: " << time_total << "s" << endl;
	ofstream out(outFile, ios::app);
	if (sa_mit_lower > sa_mit_upper) {
		cout << "SA-NAMM yields mitigation = " << sa_mit_lower << endl;
		out << "mit " << sa_mit_lower << endl;
		out << "reward breakdown " << lower_reward_data.first << " " << lower_reward_data.second << endl;
	} else {
		cout << "SA-NAMM yields mitigation = " << sa_mit_upper << endl;
		out << "mit " << sa_mit_upper << endl;
		out << "reward breakdown " << upper_reward_data.first << " " << upper_reward_data.second << endl;
	}
	out << "beta " << beta << endl;
	out << "time total " << time_total << endl;
	out << "time pre-processing " << time_pre << endl;
	out << "time namm " << time_namm << endl;
	out << "time estimating mit " << time_est_mit << endl;
	out << endl;
	out.close();

	cout << "\n*******************" << endl;
	cout << "\tALL DONE" << endl;
	cout << "*******************" << endl;

	//
	// ALL DONE
	//
}