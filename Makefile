PARA = -std=c++11 -Wall -O0
com: option.cpp graph.cpp namm.cpp el2bin.cpp format.cpp singleton_inf.cpp fake_seeds.cpp fake_inf.cpp estimate_mit.cpp
	g++ -c option.cpp -o option.o $(PARA)
	g++ -c -g graph.cpp -o graph.o $(PARA)
	g++ -Xpreprocessor -fopenmp -lomp -g namm.cpp graph.o option.o -o namm $(PARA) sfmt/SFMT.c
	g++ el2bin.cpp -o el2bin $(PARA)
	g++ format.cpp -o format $(PARA)
	g++ -Xpreprocessor -fopenmp -lomp singleton_inf.cpp graph.o option.o -o singleton_inf $(PARA) sfmt/SFMT.c
	g++ fake_seeds.cpp graph.o option.o -o fake_seeds $(PARA) sfmt/SFMT.c
	g++ -Xpreprocessor -fopenmp -lomp fake_inf.cpp graph.o option.o -o fake_inf $(PARA) sfmt/SFMT.c
	g++ -Xpreprocessor -fopenmp -lomp estimate_mit.cpp graph.o option.o -o estimate_mit $(PARA) sfmt/SFMT.c