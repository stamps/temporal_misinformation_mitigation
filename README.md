# Misinformation Mitigation under Differential Propagation Rates and Temporal Penalties
Code for the paper "Misinformation Mitigation under Differential Propagation Rates and Temporal Penalties" in VLDB 2022.

Information:
--------------------------------------------------------
Version 1.0: Implementation of NAMM Algorithm for Misinformation Mitigation under Temporal Competitive Independent (TCIC). For more details about NAMM, please read our paper: "Simpson, M., Hashemi, F., Lakshmanan, L.S., Misinformation Mitigation under Differential Propagation Rates and Temporal Penalties, VLDB 2022"

Contact Authors: Michael Simpson (michaelesimp@gmail.com)


Requirements:
--------------------------------------------------------
In order to compile all the tools, it requires GCC 4.7.2 and later (which also support OpenMP 3.1).


Compile:
--------------------------------------------------------
Use `make' command to compile everything


How to use:
--------------------------------------------------------
This package offers a set of functions to use in order to find a seed set of given size k:

0. (Optional) Computing edge weights (probabilities) as described in the experiments:
		`./format <input file> <output file> 1`

		<input file>: the path to text file in edge list format with no weights: the first line contains the number of nodes n and number of edges m, each of the next m lines describes an edge following the format: <src> <dest>.
		<output file>: the path to text output file with edge probabilities
		The last parameter (1) means the input graph is considered as directed.

1. Convert from a graph file in weighted edge list to two binary files that encode the graph and its transpose
        `./el2bin <input file> <output file> <transpose output file>`

    	<input file>: the path to text file in weighted edge list format: the first line contains the number of nodes n and number of edges m, each of the next m lines describes an edge following the format: <src> <dest> <weight>.
    	<output file>: the path to binary output file for the input graph
    	<transpose output file>: the path to binary output file for the transpose graph

2. (Optional) Compute influence of singleton seed sets
		`./singleton_inf [Options]`

	Options:

		-i <binary graph file>
            specify the path to the binary graph file (default: network.bin)

        -o <seed output file>
            specify the path to the output file containing singleton influence values (default: singleton.inf)

		-r <number of MC simulations>
			number of monte carlo simulations used to estimate the influence of each node in the graph (default = 10,000)

3. Generate fake seeds
		`./fake_seeds [Options]`

	Options:

		-n <number of nodes>
			number of nodes in the graph (mandatory)

        -o <seed output file>
            specify the path to the output file containing selected seeds (default: fake.seeds)

        -k <number of seeds>
            number of fake seed nodes (default = 1)

		-m <method>
			method for generating fake seeds (top or random, default = top)

		-f <fraction>
        	fake seeds will be drawn from the top f-th fraction

		-s <singleton influence file>
			specify the path to the file containing the singleton influence values (default: singleton.inf)

4. Compute influence of fake seeds
		`./fake_inf [Options]`

	Options:

		-i <binary graph file>
            specify the path to the binary graph file (default: network.bin)

        -o <seed output file>
            specify the path to the output file containing fake seeds influence (default: fake.inf)

        -fakeseeds <fake seeds file>
            specify the path to the fake seeds file (default: fake.seeds)

		-epsilon <epsilon value used>
            epsilon value in (epsilon,delta)-approximation (see our paper for more details, default = 0.1)

4. Run NAMM to find the seed sets
        `./namm [Options]`

    Options:

        -i <binary graph file>
            specify the path to the binary graph file (default: network)

        -o <seed output file>
            specify the path to the output file containing mitigation and timing results (default: results.txt)

        -fakeseeds <fake seeds file>
        	specify the path to the input file containing the fake seeds (default: fake.seeds)

       	-fakeinf <fake inf file>
        	specify the path to the input file containing the fake seeds influence (default: fake.inf)

        -k <number of seeds>
            number of selected seed nodes in the seed set (default = 1)

        -aw <activation window>
            parameterization for binomial distribution from which activation window lengths are drawn (see our paper for more details, default = 30)

        -rp <reading probability>
            parameterization for reading probabilities (see our paper for more details, default = 0.6)

        -ml <meeting length>
            parameterization for binomial distribution from which meeting lengths are drawn (see our paper for more details, default = 6)

        -br <base rate>
            parameterization for base propagation rate of misinformation (see our paper for more details, default = 200)

        -epsilon <epsilon value used>
            epsilon value in (epsilon,delta)-approximation (see our paper for more details, default = 0.1)

        -delta <delta value used>
            delta value in (epsilon,delta)-approximation (default: 1/n)

     Output format:
        The outputs are printed on standard output stream in the following order

                Seed Nodes: <list of selected seed nodes>
                Mitigation: <Mitigation reward of the select seed set>
                Time: <running time in seconds>

5. (Optional) Verify mitigation of a seed set - returns a (epsilon, 1/n)-estimate of the mitigation:
        
		`./estimate_mit <binary graph file> <seed file> <epsilon>`

********************************************************************************************************

Examples on a toy network: find seed set having 2 seed nodes on the graph network.txt.

The sample network network.txt, in this case, contains only 4 nodes and 4 edges and is formated as follows:
		
		4 4
		0 1 0.3
		0 2 0.4
		0 3 0.2
		1 2 0.5

1. Convert to binary file:
	
		`./el2bin network.txt network.bin network_rev.bin`

2. (Optional) Compute influence of singleton seed sets using 10K MC simulations:
	
		`./singleton_inf -i network -o singleton.inf -r 100000`

3. Generate a fake seed at random:
	
		`./fake_seeds -n 4 -o fake.seeds -k 1 -m random`
		`./fake_seeds -n 4 -o fake.seeds -k 1 -m top -f 0.25 -s singleton.inf`

4. Estimate influence of fake seed:
	
		`./fake_inf -o fake.inf -fakeseeds fake.seeds`

5. Run NAMM with k=2, epsilon=0.1, delta=0.01:
	
		`./namm -i network -fakeseeds fake.seeds -fakeinf fake.inf -k 2 -epsilon 0.1 -delta 0.01`

6. Verify mitigation reward with epsilon=0.01, assume that the seed nodes are put in network.seeds:
	
		`./estimate_mit network.bin network.seeds 0.01`

********************************************************************************************************

