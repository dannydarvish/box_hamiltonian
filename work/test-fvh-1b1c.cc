#include "../src/fvh_solver.h"
#include <armadillo>
#include <iostream>
#include <vector>
#include <complex>
#include <fstream>

double g_sigma_pi_pi = 1.6380;
double c_sigma_pi_pi = 1.0200; // fm
double G_pi_pi_pi_pi = 0.5560;
double d_pi_pi = 0.5140; //fm
double m_pi = 134.976;
vector<double> m_bare = {700.0};
vector<pair<double,double> > m_tp = {{m_pi,m_pi}};
int beta = 0;
double fm_times_MeV = 0.005067731;
double hbar_c_over_fm = 197.327;

complex<double> g_00(double k)
{
    return g_sigma_pi_pi/sqrt(m_pi)/(1+sqr(c_sigma_pi_pi*k*fm_times_MeV));
}
// vector<vector<complex<double> > > g(double k)
// {
//     vector<vector<complex<double > > > ret {{g_00(k)}};
//     return ret;
// }
vector<vector<complex<double> (*)(double)> > g = {{g_00}};

complex<double> v_00(double k, double kp)
{
    return G_pi_pi_pi_pi/sqr(m_pi) * 1.0/sqr(1+sqr(d_pi_pi*k*fm_times_MeV)) *
        1.0/sqr(1+sqr(d_pi_pi*kp*fm_times_MeV));
}

vector<vector<complex<double> (*)(double, double)> > v = {{v_00}};

int main(int argc, char** argv)
{
    int nsq_max;
    double L;
    try
    {
        nsq_max = stoi(argv[1]);
        L = stod(argv[2]);
    }
    catch(logic_error)
    {
        cerr << "Usage: test-1b1c nsq_max L" << endl;
        exit(1);
    }
    ofstream fout("L" + to_string(int(L)) + ".txt");
    FVSpectrumSolver fvss(nsq_max, L / hbar_c_over_fm, g, v, m_bare, m_tp);
    vector<double> spectrum = fvss.get_spectrum();
    for (auto e : spectrum)
        fout << e << endl;
    fout.close();
    return 0;
}