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
#include <algorithm>
#include <random>
#include <chrono> 
#include <omp.h>

using namespace std;

void readSeeds(Graph &g, int k, const char* seedFile, vector<int> &seeds) {
	ifstream in(seedFile);
	int i, v;

	i = 0;
	while (i < k) {
		in >> v;
		if ( !g.isFakeSeed(v) ) {
			seeds.push_back(v);
		}
		i++;
	}

	in.close();
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

	int aw = 30;
	char * tmp = op.getPara("-aw");
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

	int br = 200;
	tmp = op.getPara("-br");
	if (tmp != NULL){
		br = atoi(tmp);
	}

	bool ego = false;
	tmp = op.getPara("-ego");
	if (tmp != NULL){
		ego = atoi(tmp);
	}

	float ew = -1.0;
	tmp = op.getPara("-ew");
	if (tmp != NULL){
		ew = atof(tmp);
	}
	bool fixed = (ew < 0.0) ? false : true;

	Graph g(aw, ml, rp, br, ego);
	g.readGraph(inFile, fixed, ew);
	g.readFakeSeeds(fakeSeedsFile);

	int n = g.getNumNodes();

	int k = 1;
	tmp = op.getPara("-k");
	if (tmp != NULL){
		k = atoi(tmp);
	}

	char * seedFile = op.getPara("-s");
	if (seedFile == NULL){
		seedFile = (char*)"input.seeds";
	}

	vector<int> seeds;
	readSeeds(g, k, seedFile, seeds);

	double start = omp_get_wtime();
	unsigned int i;
	cout << "seed set: ";
	for (i = 0; i < seeds.size(); i++) {
		cout << seeds[i] << " ";
	}
	cout << endl;
	
	ff reward_data;
	double mit = estimateMitigation(g, n, seeds, reward_data);
	cout << "mitigation is " << mit << endl;
	cout << "reward breakdown: " << reward_data.first << " " << reward_data.second << endl;

	ofstream out(outFile, ios::app);
	out << "\nmit: " << mit << endl;
	out << "reward breakdown " << reward_data.first << " " << reward_data.second << endl;
	out.close();

	//
	// ALL DONE
	//
}