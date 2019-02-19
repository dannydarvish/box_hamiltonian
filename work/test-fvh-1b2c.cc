#include "../src/fvh_solver.h"
#include <armadillo>
#include <iostream>
#include <vector>
#include <complex>
#include <fstream>

using namespace std;

double m_pi = 134.976;
double m_k = 497.648;
double hbar_c_over_fm = 197.327;
vector<double> m_bare = {700.0};
vector<pair<double,double> > m_tp = {{m_pi,m_pi},
                                        {m_k, m_k}};

class gDerived : public gBase
{
public:
    cdouble operator() (uint i, uint alpha, double k) const override
    {
        if (alpha == 0)
            return g_00(k);
        else if (alpha == 1)
            return g_01(k);
        else
            throw(range_error("i or alpha is out of bounds for g."));
    }
private:
    double g_sigma_pi_pi = 2.0000;
    double c_sigma_pi_pi = 0.6722; // fm
    double G_pi_pi_pi_pi = 2.4998;
    double d_pi_pi = 0.2440; //fm
    double g_sigma_k_kbar = 0.6451;
    double c_sigma_k_kbar = 1.0398;
    double G_k_kbar_k_kbar = 0.0200;
    double d_k_kbar = 0.1000;
    double G_pi_pi_k_kbar = 0.3500;
    double fm_times_MeV = 0.005067731;
    cdouble g_00(double k) const
    {
        return g_sigma_pi_pi / sqrt(m_pi) / (1 + sqr(c_sigma_pi_pi * k * fm_times_MeV));
    }
    cdouble g_01(double k) const
    {
        return g_sigma_k_kbar / sqrt(m_pi) / (1 + sqr(c_sigma_k_kbar * k * fm_times_MeV));
    }
};

class vDerived : public vBase
{
public:
    cdouble operator() (uint alpha, uint beta, double k, double kp) const override
    {
        if (alpha == 0)
        {
            if (beta == 0)
                return v_00(k, kp);
            else if (beta == 1)
                return v_01(k, kp);
            else
                throw(range_error("alpha or beta is out of bounds for v"));
        }
        else if (alpha == 1)
        {
            if (beta == 0)
                return v_10(k, kp);
            else if (beta == 1)
                return v_11(k, kp);
            else
                throw(range_error("alpha or beta is out of bounds for v."));
        }
        else 
            throw(range_error("alpha or beta is out of bounds for v."));
    }
private:
    double g_sigma_pi_pi = 2.0000;
    double c_sigma_pi_pi = 0.6722; // fm
    double G_pi_pi_pi_pi = 2.4998;
    double d_pi_pi = 0.2440; //fm
    double g_sigma_k_kbar = 0.6451;
    double c_sigma_k_kbar = 1.0398;
    double G_k_kbar_k_kbar = 0.0200;
    double d_k_kbar = 0.1000;
    double G_pi_pi_k_kbar = 0.3500;
    double fm_times_MeV = 0.005067731;
    cdouble v_00(double k, double kp) const
    {
        return G_pi_pi_pi_pi/sqr(m_pi) * 1.0/sqr(1+sqr(d_pi_pi*k*fm_times_MeV)) *
            1.0/sqr(1+sqr(d_pi_pi*kp*fm_times_MeV));
    }
    cdouble v_01(double k, double kp) const
    {
        return G_pi_pi_k_kbar/sqr(m_pi) * 1.0/sqr(1+sqr(d_pi_pi*k*fm_times_MeV)) *
            1.0/sqr(1+sqr(d_k_kbar*kp*fm_times_MeV));
    }
    cdouble v_10(double k, double kp) const
    {
        return G_pi_pi_k_kbar/sqr(m_pi) * 1.0/sqr(1+sqr(d_k_kbar*k*fm_times_MeV)) *
            1.0/sqr(1+sqr(d_pi_pi*kp*fm_times_MeV));
    }
    cdouble v_11(double k, double kp) const
    {
        return G_k_kbar_k_kbar/sqr(m_pi) * 1.0/sqr(1+sqr(d_k_kbar*k*fm_times_MeV)) *
            1.0/sqr(1+sqr(d_k_kbar*kp*fm_times_MeV));
    }
};

int main(int argc, char** argv)
{
    gDerived g;
    vDerived v;

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

    ofstream fout("levels-1b2c.txt");
    double L = 4;
    FVSpectrumSolver fvss(nsq_max, L / hbar_c_over_fm, &g, &v, m_bare, m_tp);
    vector<double> spectrum;
    for (L = 4; L <= 10; L += .01)
    {
        fvss.set_L(L / hbar_c_over_fm);
        spectrum = fvss.solve_spectrum();
        for (auto e : spectrum)
            fout << e << " ";
        fout << endl;
    }
    fout.close();
    return 0;
}
