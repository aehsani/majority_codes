#include "nlcode.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cassert>



NLCodeEngine::NLCodeEngine(int k_, int n_, int m_, double epsilon_, int seed) {
	srand(seed);
	epsilon = epsilon_;
	k = k_;
	n = n_;
	m = m_;
	iters = 0;

	data_bits = new int [k];
	code_bits = new int [n];

	// generate random data bits
	for (int i = 0; i < k; i++) {
		data_bits[i] = rand() % 2;
	}

	std::cout << "Generated data bits." << std::endl;

	int x;
	bool same;
	for (int i = 0; i < n; i++) {
		// new subset
		int *new_subset = new int [m];

		// x is the data bit assigned to subset i (there are m such data bits)
		for (int j = 0; j < m; j++) {
			do {
				same = false;
				x = rand() % k;
				for (int q = 0; q < j; q++) {
					if (x == new_subset[q]) {
						same = true;
					}
				}
			} while (same);
			new_subset[j] = x;
			membership[x].push_back(i);
		}
		subsets[i] = new_subset;

		// calculate code bit
		int sum = 0;
		for (int j = 0; j < m; j++) {
			sum += data_bits[subsets[i][j]];
		} 
		code_bits[i] = (int) (sum > m/2.0);

		// erase code bith with a certain probability
		double p = rand() / (RAND_MAX + 1.);
		if (p < epsilon) {
			erased.insert(i);
		}
	}

	std::cout << "Generated code bits and erasures." << std::endl;

	sb_mssg_0 = new double [n];
	sb_mssg_1 = new double [n];

	dt_mssg_0 = new double [k];
	dt_mssg_1 = new double [k];

	for (int i = 0; i < k; i++) {
		dt_mssg_0[i] = 0.5;
		dt_mssg_1[i] = 0.5;
	}

	for (int i = 0; i < n; i++) {
		if (erased.find(i) != erased.end()) {
			// code bit was erased
			sb_mssg_0[i] = 0.5;
			sb_mssg_1[i] = 0.5;
		}
		else if (code_bits[i] == 1) {
			sb_mssg_0[i] = 0;
			sb_mssg_1[i] = 1;
		}
		else {
			sb_mssg_0[i] = 1;
			sb_mssg_1[i] = 0;
		}
	}

	std::cout << "Initialized BP messages memory." << std::endl;

	generateMajoritySequences();

	std::cout << "Generated Majority Sequences" << std::endl;

	// generate probabilities
	// Previously tried to generate proability of f given x
	p_fa_gv_xa = 0.5 + pow(0.5, m)*nChoosek(m-1, (m-1)/2);
	p_fa_gv_xb = 0.5 - pow(0.5, m)*nChoosek(m-1, (m-1)/2);
	

	int count;

	count = (m-1)/2;
	num_f_eq_x = 0;
	while (count <= m-1) {
		num_f_eq_x += nChoosek(m-1, count);
		count++;
	}

	count = (m-1)/2 + 1;
	num_f_neq_x = 0;
	while (count <= m-1) {
		num_f_neq_x += nChoosek(m-1, count);
		count++;
	}


	std::cout << "Finished Initialization" << std::endl;
	std::cout << std::endl;
}

void NLCodeEngine::generateMajoritySequences() {
	maj_seqs = new int* [(int) (pow(2.0, m) + 0.5)];
	n_seqs = 0;

	int *arr = new int [m];
	generateBinaryStrings(0, arr, 0);

}

void NLCodeEngine::generateBinaryStrings(int ones, int *arr, int i) {
	if (i == m) {
        if (ones > m/2.0) {
        	maj_seqs[n_seqs] = new int [m];
        	for (int j = 0; j < m; j++) {
        		maj_seqs[n_seqs][j] = arr[j];
        	}
        	n_seqs++;
        }
        return; 
    } 
  
    // First assign "0" at ith position 
    // and try for all other permutations 
    // for remaining positions 
    arr[i] = 0; 
    generateBinaryStrings(ones, arr, i + 1); 
  
    // And then assign "1" at ith position 
    // and try for all other permutations 
    // for remaining positions 
    arr[i] = 1; 
    generateBinaryStrings(ones+1, arr, i + 1); 
}

double NLCodeEngine::iterate(int steps) {
	int same = 0;
	double prev_error = 1, bit_error = 1;

	std::cout << "Error beforehand: " << bitError() << std::endl;

	for (int i = 0; i < steps; i++) {
		iters++;
		std::cout << "Iteration " << iters << " ";
		updateData();
		updateSubsets();
		prev_error = bit_error;
		bit_error = bitError();
		std::cout << "error: " << bit_error << std::endl;
		if (fabs(bit_error - prev_error) < threshold) {
			same += 1;
			if (same >= 2) {
				std::cout << "Error rate converged to " << bit_error;
				std::cout << " after " << iters << " iterations." << std::endl;
				return bit_error;
			}
		}
		else {
			same = 0;
		}
	}
	std::cout << std::endl;
}

double NLCodeEngine::bitError() {
	double sum = 0.0;
	for (int i = 0; i < k; i++) {
		if (dt_mssg_1[i] > dt_mssg_0[i]) {
			sum += (data_bits[i] == 0);
		}
		else {
			sum += (data_bits[i] == 1);
		}
	}
	return sum/k;
}

void NLCodeEngine::updateSubsets() {
	double p_f_1;
	double p_config;
	int *members, *current_seq;
	int data_index;
	// iterate through the subset bis
	for (int i = 0; i < n; i++) {
		// if observed subset, don't update
		
		
		if (erased.find(i) == erased.end()) {
			// code bit was not erased
			// so don't update probabilites
			continue;
		}
		


		// probability that the subset equals 1
		p_f_1 = 0;
		members = subsets[i]; // data bit indices

		// iterate through each sequence of majority ones.
		for(int j = 0; j < n_seqs; j++) {
			current_seq = maj_seqs[j];
			p_config = 1;
			for (int h = 0; h < m; h++) {
				data_index = members[h];
				if (current_seq[h] == 1) {
					p_config *= dt_mssg_1[data_index];
				}
				else {
					p_config *= dt_mssg_0[data_index];
				}
			}
			p_f_1 += p_config;
		}
		assert(p_f_1 <= 1);
		sb_mssg_1[i] = p_f_1;
		sb_mssg_0[i] = 1 - p_f_1;
	}
}

void NLCodeEngine::updateData() {
	/*
	double p_data_1, p_data_0;
	int subset, code_bit;
	// iterate throught the data bits
	for (int i = 0; i < k; i++) {
		p_data_1 = p_data_0 = 1;
		for(int j=0; j < membership[i].size(); j++) {
			subset = membership[i][j];
			if (erased.find(subset) == erased.end()) {
				code_bit = code_bits[subset];
				if (code_bit == 1) {
					p_data_1 *= p_fa_gv_xa;
					p_data_0 *= p_fa_gv_xb;
				}
				else {
					p_data_1 *= p_fa_gv_xb;
					p_data_0 *= p_fa_gv_xa;
				}
			}
			else {
				p_data_1 *= (p_fa_gv_xa*sb_mssg_1[subset]
						+ p_fa_gv_xb*sb_mssg_0[subset]);
				
				p_data_0 *= (p_fa_gv_xb*sb_mssg_1[subset]
						+ p_fa_gv_xa*sb_mssg_0[subset]);
			}
		}
		dt_mssg_0[i] = p_data_0/(p_data_0 + p_data_1);
		dt_mssg_1[i] = p_data_1/(p_data_0 + p_data_1);
	}
	*/
	
	
	
	double p_data_1, p_data_0;
	int subset, code_bit, likely_maj;
	// iterate throught the data bits
	for (int i = 0; i < k; i++) {
		//std::cout << i << ": ";
		p_data_1 = p_data_0 = 1;
		for(int j=0; j < membership[i].size(); j++) {
			subset = membership[i][j];
			//std::cout << subset << " ";
			/*
			p_data_1 *= num_f_eq_x*sb_mssg_1[subset] + 
					num_f_neq_x*sb_mssg_0[subset];

			p_data_0 *= num_f_eq_x*sb_mssg_0[subset] + 
					num_f_neq_x*sb_mssg_1[subset];
			*/
			
			// what bit is the subset most likely?
			if (sb_mssg_1[subset] > sb_mssg_0[subset]) {
				p_data_1 *= num_f_eq_x;
				p_data_0 *= num_f_neq_x;
			}
			else {
				p_data_1 *= num_f_neq_x;
				p_data_0 *= num_f_eq_x;
			}
			

		}
		//std::cout << std::endl;
		dt_mssg_0[i] = p_data_0/(p_data_0 + p_data_1);
		dt_mssg_1[i] = p_data_1/(p_data_0 + p_data_1);
	}
	


}





int NLCodeEngine::nChoosek(int a, int b) {
    if (b > a) return 0;
    if (b * 2 > a) b = a-b;
    if (b == 0) return 1;

    int result = a;
    for( int i = 2; i <= b; ++i ) {
        result *= (a-i+1);
        result /= i;
    }
    return result;
}

void NLCodeEngine::debug() {
	std::cout << "DEBUG" << std::endl;
	std::cout << "data bits" << std::endl;
	for (int i = 0; i < k; i++) {
		std::cout << data_bits[i];
	}
	std::cout << std::endl;

	std::cout << "code bits" << std::endl;
	for (int i = 0; i < n; i++) {
		std::cout << code_bits[i];
	}
	std::cout << std::endl;

	std::cout << "subsets" << std::endl;
	for (auto const& pair: subsets) {
		std::cout << pair.first << ": ";

		std::cout << "{indices: ";
		for (int i = 0; i < m; i++) {
			std::cout << pair.second[i] << " ";
		}
		std::cout << "} ";

		std::cout << "{data: ";
		for (int i = 0; i < m; i++) {
			std::cout << data_bits[pair.second[i]] << " ";
		}
		std::cout << "} ";

		std::cout << "bit: " << code_bits[pair.first];
		std::cout << std::endl;
	}

	std::cout << "erased code bits" << std::endl;
	for (auto j = erased.begin(); j != erased.end(); j++) {
		std::cout << (*j) << " ";
	}
	std::cout << std::endl;

	
	std::cout << "membership" << std::endl;
	for (auto const& pair: membership) {
		std::cout << pair.first << ": ";

		std::cout << "{ ";
		for(int i=0; i < pair.second.size(); i++) {
   			std::cout << pair.second[i] << " ";
		}
		std::cout << "} \n";
	}
	

	/*
	std::cout << "majority seqs" << std::endl;
	for (int i = 0; i < n_seqs; i++) {
		for (int j = 0; j < m; j++) {
			std::cout << maj_seqs[i][j];
		}
		std::cout << std::endl;
	}
	*/

	std::cout << "num configs f = a if x = a: " << num_f_eq_x << std::endl;
	std::cout << "num configs f = a if x = b: " << num_f_neq_x << std::endl;

	std::cout << std::endl;

}

void NLCodeEngine::debugSubsetUpdate() {
	// print data messages
	std::cout << "data messages 1 and 0" << std::endl;
	std::cout << std::setprecision(3);
	for (int i = 0; i < k; i++) {
		std::cout << dt_mssg_1[i] << " ";
	}
	std::cout << std::endl;
	for (int i = 0; i < k; i++) {
		std::cout << dt_mssg_0[i] << " ";
	}
	std::cout << std::endl;

	updateSubsets();

	std::cout << "subset messages 1 and 0" << std::endl;
	std::cout << std::setprecision(3);
	for (int i = 0; i < n; i++) {
		std::cout << sb_mssg_1[i] << " ";
	}
	std::cout << std::endl;
	for (int i = 0; i < n; i++) {
		std::cout << sb_mssg_0[i] << " ";
	}
	std::cout << std::endl;

}




