#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>
#include <queue>

#define A 6.3781363e6
#define MU 3.986004415e14
#define Delta 1e-8
#define MAX_LAT 90
#define MIN_LAT -90
#define MAX_LON 180
#define MIN_LON -180
#define EPS 1e-6

using namespace std;
class Error{
	string log = "log.txt";
	string Msg;
public:
	enum Err{ERROR_OPEN, ERROR_EMPTY, ERROR_NO_BEGIN, ERROR_NO_END, ERROR_NO_DATA, ERROR_NO_MEMORY,
			 ERROR_BAD_LONGITUDE, ERROR_BAD_LATITUDE, ERROR_BAD_N};
	Error(Err err){
		const char * ERR_MSG[] = {"Error opening file\n", "File is empty\n", "No begin keyword in file\n", 
		"No end keyword in file\n", "No data in opened file\n", "No memory\n", "Bad latitude interval\n", 
		"Bad longitude interval\n", "Bad N (less or equal to 0)\n"};
		Msg = ERR_MSG[err];
	}
	~Error(){}
    void msg(){
		ofstream out(log);
		out << Msg;
	}
};

class Parser{
	double ** coeff;
	int size;
	string filename;
	void split_input(const string line, string & sn, string & sm, string & sc, string & ss) const;
public:
	Parser(string fname = ""){
		size = 1;
		filename = fname;
		coeff = new double * [size];
		for(int i = 0; i < size; i++)
			coeff[i] = new double[2];
		if(!coeff) throw Error(Error::ERROR_NO_MEMORY);
	}
	/*Parser(const Parser & p){
		filename = p.filename;
		size = p.size;
		coeff = new double * [size];
		for(int i = 0; i < size; i++){
			coeff[i] = new double[2];
			coeff[i][0] = p.get_coeff()[i][0];
			coeff[i][1] = p.get_coeff()[i][1];
		}
	}*/
	~Parser(){
		for(int i = 0; i < size; i++) delete [] coeff[i];   //here
	}
	int get_size() const{return size; }
	int get_max() const;
	void calc_coeff();
	double ** get_coeff() const;
	Parser * operator= (Parser * p){
		//allocate memory
		this->coeff = new double * [size];
		for(int i = 0; i < size; i++){
			this->coeff[i] = new double[2];
			this->coeff[i][0] = p->get_coeff()[i][0];
			this->coeff[i][1] = p->get_coeff()[i][1];
		}
		if(!coeff) throw Error(Error::ERROR_NO_MEMORY);
		this->size = p->size;
		this->filename = p->filename;
		return this;
	}
};

class Solver{
	Parser * parser;
	int size;
	int max;
	long long C(int n, int k);
	double assoc_1egendre(unsigned int n, unsigned int m, double x);
	double fact(unsigned int n);
	double legendre(unsigned int n, unsigned int m, double x);
public:
	Solver(){}
	Solver(Parser * p){
		parser = p;
		size = p->get_size();
	}
	/*Solver(const Solver & s){
		Parser tmp(*s.parser);
		parser = &tmp;
		size = s.size;
		max = s.max;
	}*/
	~Solver(){}
	void set_max(int n);
	double get_u(double r, double phi, double theta);
	double get_v(double r, double phi, double theta);
	double get_indul(double r, double phi, double theta);
	Solver *  operator= (Solver * s){
		//allocate memory
		size = s->size;
		Parser tmp(*(s->parser));
		parser = &tmp;
		max = s->max;
		return this;
	}
};

class Analyser{
	int N;
	double min_lon, max_lon, min_lat, max_lat;
	double *** table;
	Solver * solver;
	int* search_next(double ** mat, int ** col);
	int ** find(double ** mat, int type, int &size);
public:
	Analyser(){}
	Analyser(int n, Solver * s){
		solver = s;
		N = n;
		if(n <= 0) throw Error(Error::ERROR_BAD_N);
		table = new double ** [N];
		for(int i = 0; i < N; i++){
			table[i] = new double * [N];
			for(int j = 0; j < N; j++)
				table[i][j] = new double [3];
		}
		if(!table) throw Error(Error::ERROR_NO_MEMORY);
	}
	/*Analyser(const Analyser & a): solver(a.solver){
		table = a.table;
		min_lon = a.min_lon;
		max_lon = a.max_lon;
		min_lat = a.min_lat;
		max_lat = a.max_lat;
		N = a.N;
	}*/
	~Analyser(){
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++)
				delete [] table[i][j];				//and here..
		}
	}
	
	//getters (trivial)
	double *** get_table(){ return table; }
	int get_N(){ return N; }
	
	//set bounders: minimum and maximum for longtude and latitude
	void set_bound(double mlon, double Mlon, double mlat, double Mlat);
	//make table NxN 3 values in a cell: longitude, latitude and indulation
	void make_table();
	//get extremums. answer gets by link
	void get_extremums(double threshold_max, double threshold_min, int ** & maximums,  int ** & minimums, int & num_max, int & num_min);
	Analyser * operator= (Analyser * a){
		//table = a->table;
		min_lon = a->min_lon;
		max_lon = a->max_lon;
		min_lat = a->min_lat;
		max_lat = a->max_lat;
		N = a->N;
		//copy table
		table = new double ** [N];
		if(!table) throw Error(Error::ERROR_NO_MEMORY);
		for(int i = 0; i < N; i++){
			table[i] = new double * [N];
			if(!table[i]) throw Error(Error::ERROR_NO_MEMORY);
			for(int j = 0; j < N; j++){
				table[i][j] = new double [3];
				table[i][j][0] = a->table[i][j][0];
				table[i][j][1] = a->table[i][j][1];
				table[i][j][2] = a->table[i][j][2];
			}
		}
		
		return this;
	}
};


class Presenter{
	string outfile;
	Analyser * analyser;
public:	
	Presenter(){}
	Presenter(string outf, Analyser * a){
		analyser = a;
		outfile = outf;
	}
	~Presenter(){}
	//write table from analyser to the output file
	void write_indulation();
	//and also write extremums 
	void write_extremums(int ** maximums, int ** minimums, int num_max, int num_min);
};



