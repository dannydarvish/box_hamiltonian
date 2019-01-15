#include "../src/fvh_solver.h"
#include <armadillo>
#include <iostream>
#include <vector>
#include <complex>
#include <fstream>

using namespace std;

inline double k_from_om(double om, double m1, double m2)
{
    return sqrt(sqr(1.0/(2*om)*(sqr(om)+sqr(m1)-sqr(m2)))-sqr(m2));
}

double g_sigma_pi_pi = 2.0000;
double c_sigma_pi_pi = 0.6722; // fm
double G_pi_pi_pi_pi = 2.4998;
double d_pi_pi = 0.2440; //fm
double m_pi = 134.976;
double m_k = 497.648;
double g_sigma_k_kbar = 0.6451;
double c_sigma_k_kbar = 1.0398;
double G_k_kbar_k_kbar = 0.0200;
double d_k_kbar = 0.1000;
double G_pi_pi_k_kbar = 0.3500;

vector<double> m_bare = {700.0};
vector<pair<double,double> > m_tp = {{m_pi,m_pi},
                                        {m_k, m_k}};
double fm_times_MeV = 0.005067731;
double hbar_c_over_fm = 197.327;

vector<double> E_vec;

complex<double> g_00(double k)
{
    return g_sigma_pi_pi/sqrt(m_pi)/(1+sqr(c_sigma_pi_pi*k*fm_times_MeV));
}
complex<double> g_01(double k)
{
    return g_sigma_k_kbar/sqrt(m_pi)/(1+sqr(c_sigma_k_kbar*k*fm_times_MeV));
}
vector<vector<complex<double> (*)(double)> > g = {{g_00, g_01}};

complex<double> v_00(double k, double kp)
{
    return G_pi_pi_pi_pi/sqr(m_pi) * 1.0/sqr(1+sqr(d_pi_pi*k*fm_times_MeV)) *
        1.0/sqr(1+sqr(d_pi_pi*kp*fm_times_MeV));
}
complex<double> v_01(double k, double kp)
{
    return G_pi_pi_k_kbar/sqr(m_pi) * 1.0/sqr(1+sqr(d_pi_pi*k*fm_times_MeV)) *
        1.0/sqr(1+sqr(d_k_kbar*kp*fm_times_MeV));
}
complex<double> v_10(double k, double kp)
{
    return G_pi_pi_k_kbar/sqr(m_pi) * 1.0/sqr(1+sqr(d_k_kbar*k*fm_times_MeV)) *
        1.0/sqr(1+sqr(d_pi_pi*kp*fm_times_MeV));
}
complex<double> v_11(double k, double kp)
{
    return G_k_kbar_k_kbar/sqr(m_pi) * 1.0/sqr(1+sqr(d_k_kbar*k*fm_times_MeV)) *
        1.0/sqr(1+sqr(d_k_kbar*kp*fm_times_MeV));
}
vector<vector<complex<double> (*)(double, double)> > v = {{v_00, v_01},
                                                          {v_10, v_11}};

int main(int argc, char** argv)
{
    int nsq_max;
    try
    {
        nsq_max = stoi(argv[1]);
    }
    catch(logic_error)
    {
        cerr << "Usage: test-1b2c nsq_max" << endl;
        exit(1);
    }
    // double L;
    // try
    // {
    //     nsq_max = stoi(argv[1]);
    //     L = stod(argv[2]);
    // }
    // catch(logic_error)
    // {
    //     cerr << "Usage: test-1b2c nsq_max L" << endl;
    //     exit(1);
    // }
    ofstream fout("levels.txt");
    double L = 4;
    FVSpectrumSolver fvss(nsq_max, L / hbar_c_over_fm, g, v, m_bare, m_tp);
    vector<double> spectrum;
    for (L = 4; L <= 10; L += .01)
    {
        fvss.set_L(L / hbar_c_over_fm);
        spectrum = fvss.get_spectrum();
        for (auto e : spectrum)
            fout << e << " ";
        fout << endl;
    }
    fout.close();
    return 0;
}
