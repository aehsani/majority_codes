#include <iostream>
#include "nlcode.h"
#include <fstream>



int main() {

	// first parameter is k
	// second parameter is n
	// third parameter is subset length
	// fourth parameter is epsilon

	
	double epsilon = 0.5;
	int seed = 7;
	
	int k = 500000;
	int n = 1000000;
	int m = 3;

	
	NLCodeEngine nlcode = NLCodeEngine(k, n, m, epsilon, seed);
	nlcode.iterate(50);
	return 0;
	
	
	
	

	//debugging
	
	/*
	NLCodeEngine nlcode_debug = NLCodeEngine(7, 3, 5, 0.5, 9);
	
	nlcode_debug.debug();
	nlcode_debug.iterate(1);
	nlcode_debug.debugSubsetUpdate();
	return 0;
	
	*/






	/*
	int n = 150000;
	int k = 100000;
	int m = 7;
	myfile.open("data_file")
	for (int i = 0; i < 21; i++) {
		eps = i*0.05;
		NLCodeEngine nlcode = NLCodeEngine(k, n, m, eps);
		error = nlcode.iterate(50);
  		myfile << (1-eps)*n/(double) k << " " << error << std::endl;
	}
	myfile.close();
	return 0;
	*/
}