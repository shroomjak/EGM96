#include "EGM96.h"
using namespace std;
void Parser::split_input(const string line, string & sn, string & sm, string & sc, string & ss) const{
	stringstream str(line);			//make artifical stream
	str >> sn >> sm >> sc >> ss; 	//and read from it
}
int Parser::get_max() const{
	// (n+1)*(n+2)/2 = s n+1-?
	return int(0.5*(-1+sqrt(1+8*size)));
}

long long Solver::C(int n, int k){
	long long a = 1, b = 1;
	for(int i = 1; i <= k; i++){
		a *= (n+1-i);
		b *= i;
	}
	return a/b;
}

double Solver::assoc_1egendre(unsigned int n, unsigned int m, double x){
	double res = 0.0;
	for(unsigned int i = m; i < n; i++)
		res += C(n, i)*C((n+i-1)/2, n)*pow(x, i-m)*fact(i)/fact(i-m);
	return (1<<n)*pow((1-x*x), m/2.0)*((m%2)?(-1):(1))*res;
}

void Parser::calc_coeff(){
		string startkeyword  = "BEGIN Coefficients";
		string endkeyword    =  "END Coefficients";
			ifstream file(filename);
			if(!file) throw Error(Error::ERROR_OPEN);
			file.seekg(0, std::ios::end);
			if(file.tellg() == 0) throw Error(Error::ERROR_EMPTY);
			file.seekg(0, std::ios::beg);
				int n;
				int max_n = 0;
				string line;	
				while (getline(file, line) and line != startkeyword);
				if(line != startkeyword) throw Error(Error::ERROR_NO_BEGIN);
				while (getline(file, line) and line != endkeyword){
					string sn, sm, sc, ss;
					split_input(line, sn, sm, sc, ss);
					n = stoi(sn);
					if(n > max_n) max_n = n;
				}
				if(line != endkeyword) throw Error(Error::ERROR_NO_END);
				if(max_n == 0) throw Error(Error::ERROR_NO_DATA);
				file.close();
				file.open(filename);
				size = ((max_n+2)*(max_n+1))/2;
				cout << size << endl;
				coeff = new double * [size];
				for(int i = 0; i < size; i++) coeff[i] = new double [2];
				if(!coeff) throw Error(Error::ERROR_NO_MEMORY);
				coeff[0][0] = 0, coeff[0][1] = 0;
				coeff[1][0] = 0, coeff[1][1] = 0;
				
				while (getline(file, line) and line != startkeyword);
				while (getline(file, line) and line != endkeyword){
					string sn, sm, sc, ss;
					split_input(line, sn, sm, sc, ss);
					int n = stoi(sn);
					int m = stoi(sm);
					double C = stod(sc);
					double S = stod(ss);
					int index = (n*(n+1))/2-1+m;
					coeff[index][0] = C;
					coeff[index][1] = S; 
				}
				file.close();
				//get = 1;
}
double ** Parser::get_coeff() const{
	return coeff;
}

double Solver::fact(unsigned int n){
	double res = 1;
	if(n == 0) return 1;
	for(unsigned int i = 2; i <= n; i++) res *= (double) i;
	return res;
}
double Solver::legendre(unsigned int n, unsigned int m, double x){
	return sqrt(double((2-(m==0))*(2*n+1)*fact(n-m))/double(fact(n+m)))*assoc_legendre(n, m, x);
}

void Solver::set_max(int n){
	max = n;
}

double Solver::get_u(double r, double phi, double theta){
	double ** coeff = parser->get_coeff();
	double tmp = 0.0;
	double res = 0.0;
	for(int n = 1; n <= max; n++){
		tmp = 0.0;
		for(int m = 0; m <= n; m++){
			double C = coeff[(n*(n+1))/2-1+m][0],
				   S = coeff[(n*(n+1))/2-1+m][1];
			//cout << "C[" << n << ", " << m << "] = " << C << " " << "S[" << n << ", " << m << "] = " << S << endl;
			tmp += legendre(n, m, sin(theta))*(C*cos(m*phi) + S*sin(m*phi));
		}
		res += pow(double(A)/r, n)*tmp;
	}
	return -(1.0+res)*(MU/r);
}
double Solver::get_v(double r, double phi, double theta){
	double ** coeff = parser->get_coeff();
	double res = 0.0;
	for(int n = 2; n <= max; n += 2){
		double J = coeff[(n*(n+1))/2-1][0];
		res += pow(double(A)/r, n)*J*legendre(n, 0, sin(theta));
	}
	return -(1.0+res)*(MU/r);
}
double Solver::get_indul(double r, double phi, double theta){
	return double(A*A)/MU*(-get_u(r, phi, theta)+get_v(r, phi, theta));
}
void Analyser::set_bound(double mlon, double Mlon, double mlat, double Mlat){
	min_lon = mlon;
	max_lon = Mlon;
	min_lat = mlat;
	max_lat = Mlat;
	if(min_lon > max_lon || abs(min_lon - max_lon) < EPS) throw Error(Error::ERROR_BAD_LONGITUDE);
	if(min_lat > max_lat || abs(min_lat - max_lat) < EPS) throw Error(Error::ERROR_BAD_LATITUDE);
}
void Analyser::make_table(){
	double lon_step = (max_lon-min_lon)/(N-1),
		   lat_step = (max_lat-min_lat)/(N-1),
		   r = A;
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			double  lon = (min_lon+lon_step*j)*M_PI/180.0,
					lat = (min_lat+lat_step*i)*M_PI/180.0;
			table[i][j][0] = lon;
			table[i][j][1] = lat; 	
			table[i][j][2] = solver->get_indul(r, lon, lat);
		}
	}
	
}
void Presenter::write_indulation(){
	ofstream out(outfile);
	int N = analyser->get_N();
	double *** table = analyser->get_table();
	out << N << endl;
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			out << table[i][j][0] << " " << table[i][j][1] << " " << table[i][j][2] << " ";
		}
		out << endl;
	}
	out.close();
}

void Presenter::write_extremums(int ** maximums, int ** minimums, int num_max, int num_min){
	ofstream out;
	out.open(outfile, ios::app);
	out << endl;
	out << num_max << endl;
	for(int i = 0; i < num_max; i++) out << maximums[0][i] << " " << maximums[1][i] << endl;
	out << num_min << endl;
	for(int i = 0; i < num_min; i++) out << minimums[0][i] << " " << minimums[1][i] << endl;
	out.close();
}
int * Analyser::search_next(double ** mat, int ** col){
	int* ans = new int[2];
	ans[0] = ans[1] = -1;

	int flag = 1;
	for(int i = 0; i < N && flag; i++)
		for(int j = 0; j < N && flag; j++) 
			if(mat[i][j] != 0 && col[i][j] == 0){
				flag = 0;
				ans[0] = i; ans[1] = j;
				return ans;
			}
	return ans;
}

int ** Analyser::find(double ** mat, int type, int &size){  
	int ** col = new int * [N];
	for(int i = 0; i < N; i++) col[i] = new int[N];
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++) col[i][j] = 0;
			
	queue <int*> q;

	int * start;
	double max;
	int * max_ind, * v, * tmp;
	int ** ans = new int *[2];
	ans[0] = ans[1] = nullptr;
	size = 0;
	tmp = new int[2];
	while ((start = search_next(mat, col))[0] != -1){
		max = mat[start[0]][start[1]] * type;
		max_ind = start;
		q.push(start);
		while(!q.empty()){
			v = q.front(); 
			col[v[0]][v[1]] = 1;
			
			if(mat[v[0]][v[1]]*type - max > EPS){
				max = mat[v[0]][v[1]]*type;
				max_ind = v;
			}
			
			int ngbr[][2] = {{-1, -1}, {-1, 0}, {-1, +1}, {0, -1}, {0, +1}, {+1, -1}, {+1, 0}, {+1, +1}};
			for(int i = 0; i < 8; i++){
				if(v[0]+ngbr[i][0] >= 0 && v[0]+ngbr[i][0] < N && v[1]+ngbr[i][1] >= 0 && v[1]+ngbr[i][1] < N && col[v[0]+ngbr[i][0]][v[1]+ngbr[i][1]] == 0 && fabs(mat[v[0]+ngbr[i][0]][v[1]+ngbr[i][1]]) > EPS){
					tmp = new int[2];
					col[v[0]+ngbr[i][0]][v[1]+ngbr[i][1]] = -1;
					tmp[0] = v[0]+ngbr[i][0];
					tmp[1] = v[1]+ngbr[i][1];
					q.push(tmp);
				}
			}
			//delete[] v;
			q.pop();
		}
		++size;
		ans[0] = (int*)realloc(ans[0], size*sizeof(int));
		if(!ans[0]) throw Error(Error::ERROR_NO_MEMORY);
		ans[1] = (int*)realloc(ans[1], size*sizeof(int));
		if(!ans[1]) throw Error(Error::ERROR_NO_MEMORY);
		ans[0][size-1] = max_ind[0];
		ans[1][size-1] = max_ind[1];    
	}
	for(int i = 0; i < N; i++) delete [] col[i];
	while (!q.empty()) q.pop();
	return ans;
}

void Analyser::get_extremums(double threshold_max, double threshold_min, int ** & maximums,  int ** & minimums, int & num_max, int & num_min){
	double ** mat_max = new double * [N];
    double ** mat_min = new double * [N];
    for(int i = 0; i < N; i++){
		mat_max[i] = new double [N];
		mat_min[i] = new double [N];
    }
    for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++){
			mat_min[i][j] = mat_max[i][j] = 0;
			if(table[i][j][2] > threshold_max) mat_max[i][j] = table[i][j][2]; 
			if(table[i][j][2] < threshold_min) mat_min[i][j] = table[i][j][2]; 
      }
    minimums = find(mat_min, -1, num_min);
    maximums = find(mat_max, +1, num_max);
    for(int i = 0; i < N; i++){
		delete [] mat_max[i];
		delete [] mat_min[i];
	}
}
