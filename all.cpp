#include "EGM96.h"
using namespace std;
int main(void){
	try{
		string input_file = "EGM96.grv";
		Parser parser(input_file);
		parser.calc_coeff();
		
		Solver solver(&parser);
		solver.set_max(25);
		
		int N = 100; //number of points in a row (total N^2 points)
		
		string output_file = "potential.txt";
		
		Analyser analyser(N, &solver);
		analyser.set_bound(MIN_LON, MAX_LON, MIN_LAT, MAX_LAT);
		analyser.make_table();
		
		int ** minimums, ** maximums, num_min, num_max;
		analyser.get_extremums(60, -50, maximums, minimums, num_max, num_min);
		
		Presenter presenter(output_file, &analyser);
		presenter.write_indulation();
		presenter.write_extremums(maximums, minimums, num_max, num_min);
		/*for(int i = 0; i < num_max; i++)
			free(maximums[i]);
		for(int i = 0; i < num_min; i++)
			free(minimums[i]);*/
	}catch(Error err){
		cout << "Error. See log.txt\n";
		err.msg();
	}
	return 0;
}

/*
 * "EGM96.grv" -> parser-> coeff
 *  coeff -> solver -> func get_indul
 *  func get_indul -> analyser -> extremums & table
 *  extremums & table -> presenter -> "potential.txt"
 *  "potential.txt" -> Presenter_flat.m -> plot indulation & plot extremums
 * 
 * ]
 * LON N*N LAT N*N u N*N
 * [i, j]
 */
