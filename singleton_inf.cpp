#include "option.h"
#include "graph.h"
#include <algorithm>
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

typedef pair<float,int> fi;
typedef vector<fi> vfi;

// estimate influence of singleton seeds using 10K monte carlo simulations
void estimateInfluence(Graph &g, int n, int r, vfi &inf) {

	#pragma omp parallel 
	{
		long long int running_total;
		float influence;

		vector<bool> visit(n,false);
		vector<int> visit_index(n,0);

		#pragma omp for
		for (int i = 0; i < n; i++) {
			running_total = 0;
			for (int j = 0; j < r; j++) {
				running_total += g.generateInfluenceSample(visit, visit_index, i);
			}
			influence = (float)running_total / r;
			inf[i] = fi(influence, i);
		}
	}

	sort(inf.begin(), inf.end());
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
		outFile = (char*)"singleton.inf";
	}

	float ew = -1.0;
	char * tmp = op.getPara("-ew");
	if (tmp != NULL){
		ew = atof(tmp);
	}
	bool fixed = (ew < 0.0) ? false : true;

	Graph g(0, 0, 0.0, 0, false);
	g.readGraph(inFile, fixed, ew);

	int n = g.getNumNodes();

	int r = 10000;
	tmp = op.getPara("-r");
	if (tmp != NULL){
		r = atoi(tmp);
	}

	vfi inf(n);

    double start = omp_get_wtime();
    estimateInfluence(g, n, r, inf);
    cout << "Time to estimate singleton influences: " << (omp_get_wtime()-start) << "s" << endl;

    ofstream of(outFile);
	for (int i = n-1; i >= 0; i--){
		of << inf[i].second << "\t" << inf[i].first << endl;
	}
	of.close();

	//
	// ALL DONE
	//
}