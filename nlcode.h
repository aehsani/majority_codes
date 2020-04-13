#ifndef NLCODE_H
#define NLCODE_H

#include <map>
#include <vector>
#include <unordered_map>
#include <unordered_set>

class NLCodeEngine {
	private:
		// iterations
		int iters = 0;

		//convergence threshold
		double threshold = 0.0005;

		// data length
		int k;

		// code length
		int n;

		// subset size
		int m;

		// probability or erasure
		double epsilon;

		// map from subsets to the bits each subset contains
		std::unordered_map<int, std::vector<int>> membership;

		// map from bits to the subsets that contain that bit
		std::unordered_map<int, int*> subsets;

		// bits in f that were erased (later flipped)
		std::unordered_set<int> erased;

		// data bits
		int *data_bits;

		// code bits
		int *code_bits;

		// subset bit messages (similar to probabilities)
		double *sb_mssg_0;
		double *sb_mssg_1;

		// data bit messages (similar to probabilities)
		double *dt_mssg_0;
		double *dt_mssg_1;

		// sequences of length m that have more than m 1s
		int **maj_seqs;
		int n_seqs;

		double p_fa_gv_xa; 
		double p_fa_gv_xb;

		int num_f_eq_x;
		int num_f_neq_x;

	public:
		// new messages (used in bp)
		//std::unordered_map<

		// old messages (used in bp)
		//std::unordered_map

		//Constructor
		// k is data bits
		// n is code length
		// m is the number of bits in a subset
		// epsilon is probability of erasure
		NLCodeEngine(int k, int n, int m, double epsilon, int seed);

		void generateMajoritySequences();
		void generateBinaryStrings(int ones, int *arr, int i);

		double iterate(int steps);
		double bitError();

		// use these for "fake" loopy BP algorithm
		void updateSubsets();
		void updateData();

		int nChoosek(int a, int b);

		void debug();

		void testUpdateSubsets();

		void testUpdateData();

		void debugSubsetUpdate();
};

#endif