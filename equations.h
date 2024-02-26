#pragma once
using namespace std;
//Данные в новом временном слое 
extern double p1;
extern double c1;
extern double err_n,
ezz_n,
e00_n,
erz_n;
extern double delta_err_n,
delta_ezz_n,
delta_e00_n,
delta_erz_n;
extern double Drr_n,
Dzz_n,
D00_n,
Drz_n;
extern double E_n;
extern double sigma_rr_n,
sigma_zz_n,
sigma_00_n,
sigma_rz_n;
extern double rn, zn;
extern double vrn, vzn;
extern double A1;
extern double r01;
extern double z01;
extern double V1;
extern double ro1;
extern double deltaT;
extern double deltaT12;
extern double A12;
extern double V12;
extern double V12tV;
extern double qkomb_n,
qrr_n,
qzz_n,
q00_n,
qrz_n;

const double PI = 3.141592653589793;
void Geometry_nach(vector<Part>& Pr, vector<Mat>& Mt, vector<Yacheika>& Ych, vector<Uzel>& Uz);
void Geometry(int i, double& A1, double& r01, double& z01, double& V1, double& ro1,
	vector<Part>& Pr, vector<Mat>& Mt, vector<Yacheika>& Ych, vector<Uzel>& Uz);
void Kurant(int i, double kr, double& deltaT, vector<Yacheika>& Ych, vector<Uzel>& Uz);
void SkorKoord(double& vrn, double& vzn, double& rn, double& zn, vector<Yacheika>& Ych,
	vector<Uzel>& Uz, vector<Part>& Pr, Pressure& Press, double time);
void Defor(int i, double& A12, double& V12, double& V12tV, double& err_n, double& ezz_n,
	double& e00_n, double& erz_n, double& delta_err_n, double& delta_ezz_n,
	double& delta_e00_n, double& delta_erz_n, vector<Yacheika>& Ych,
	vector<Uzel>& Uz, double t);
void Psevdov(int i, double lamdal, double lamdak, double& qkomb_n, double& qrr_n,
	double& qzz_n, double& q00_n, double& qrz_n, vector<Yacheika>& Ych);
void Naprag(int i, double& Drr_n, double& Dzz_n, double& D00_n, double& Drz_n,
	double& sigma_rr_n, double& sigma_zz_n, double& sigma_00_n, double& sigma_rz_n,
	vector<Yacheika>& Ych, vector<Uzel>& Uz, vector <Mat> Mt, vector <Part> Pr);
void Energy(int i, double& E_n, double Drr_n, double Dzz_n, double D00_n, double Drz_n,
	double err_n, double ezz_n, double e00_n, double erz_n, double V12, double V1,
	double deltaT, double qkomb_n, vector<Yacheika>& Ych);
void Perezap(int i, vector<Yacheika>& Ych);
void Obnul();



