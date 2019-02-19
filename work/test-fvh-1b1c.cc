#include "../src/fvh_solver.h"
#include <armadillo>
#include <iostream>
#include <vector>
#include <complex>
#include <fstream>

const double hbar_c_over_fm = 197.327;
double m_pi = 134.976;
vector<double> m_bare = {700.0};
vector<pair<double,double> > m_tp = {{m_pi,m_pi}};
int beta = 0;

class gDerived : public gBase
{
public:
    cdouble operator() (uint i, uint alpha, double k) const override
    {
        if (i > 0 || alpha > 0)
            throw(range_error("i or alpha is out of bounds for g."));
        return g_00(k);
    }
private:
    double g_sigma_pi_pi = 1.6380;
    double c_sigma_pi_pi = 1.0200; // fm
    double G_pi_pi_pi_pi = 0.5560;
    double d_pi_pi = 0.5140; //fm
    double m_pi = 134.976;
    double fm_times_MeV = 0.005067731;
    double hbar_c_over_fm = 197.327;

    cdouble g_00(double k) const
    {
        return g_sigma_pi_pi/sqrt(m_pi)/(1+sqr(c_sigma_pi_pi*k*fm_times_MeV));
    }
};

class vDerived : public vBase
{
public:
    cdouble operator() (uint alpha, uint beta, double k, double kp) const override
    {
        if (alpha > 0 || beta > 0)
            throw(range_error("alpha or beta is out of bounds for v."));
        return v_00(k, kp);
    }
private:
    double g_sigma_pi_pi = 1.6380;
    double c_sigma_pi_pi = 1.0200; // fm
    double G_pi_pi_pi_pi = 0.5560;
    double d_pi_pi = 0.5140; //fm
    double fm_times_MeV = 0.005067731;

    cdouble v_00(double k, double kp) const
    {
        return G_pi_pi_pi_pi/sqr(m_pi) * 1.0/sqr(1+sqr(d_pi_pi*k*fm_times_MeV)) *
            1.0/sqr(1+sqr(d_pi_pi*kp*fm_times_MeV));
    }
};

int main(int argc, char** argv)
{
    gDerived g;
    vDerived v;

    int nsq_max;
    double L = 4;
    try
    {
        nsq_max = stoi(argv[1]);
    }
    catch(logic_error)
    {
        cerr << "Usage: test-1b1c nsq_max" << endl;
        exit(1);
    }
    ofstream fout("levels-1b1c.txt");
    FVSpectrumSolver fvss(nsq_max, L / hbar_c_over_fm, &g, &v, m_bare, m_tp);
    for (L = 4; L <= 10; L += 0.01)
    {
        fvss.set_L(L / hbar_c_over_fm);
        vector<double> spectrum = fvss.solve_spectrum();
        for (auto e : spectrum)
            fout << e << " ";
        fout << endl;
    }
    fout.close();
    return 0;
}