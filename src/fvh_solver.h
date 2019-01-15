#ifndef FVH_SOLVER
#define FVH_SOLVER
#include <armadillo>
#include <complex>
#include <algorithm>
#include <map>

using namespace std;
using namespace arma;

typedef complex<double> cdouble;
typedef vector<complex<double> > cvec;

inline int sqr(int x){return x*x;}
inline double sqr(double x){return x*x;}

// Notation here follows the Wu et. al paper
// Unless otherwise stated, units are in MeV and fm.
class FVSpectrumSolver
{
public:
    FVSpectrumSolver(int nsq_max, double L, vector<vector<cdouble (*)(double)> > g,
                     vector<vector<cdouble (*)(double, double)> > v,
                     const vector<double> & m_bare, const vector<pair<double, double> > & m_tp) :
                     nsq_max(nsq_max), m_bare(m_bare), m_tp(m_tp), L(L),
                     n_bare(m_bare.size()), n_chan(m_tp.size()), g(g), v(v),
                     finished_calc(false)
    {
        for (int i = 0; i <= sqrt(nsq_max); i++)
            for (int j = 0; j <= i; j++)
                for (int k = 0; k <= j; k++)
                {                    
                    int val = sqr(i)+sqr(j)+sqr(k);
                    if (val <= nsq_max &&
                        find(momentum_ints.begin(), momentum_ints.end(), val) == momentum_ints.end())
                        {
                            momentum_ints.push_back(val);
                        }
                }
        for (int nsq : momentum_ints)
            C3_map.insert(make_pair(nsq, C3(nsq)));
    }
    ~FVSpectrumSolver() {}
    const Mat<cdouble> & get_hamiltonian();
    const vector<double> & get_spectrum();
    void set_L(double L);
private:
    // n_bare: number of bare particles
    // n_chan: number of two-particle channels
    // m_bare: vector of single particle masses
    // m_tp: vector of pairs of two-particle channel masses
    // L: lattice length in fm
    // g: User supplied 1-to-2 coupling (infinite volume)
    // v: User supplied 2-to-2 coupling (infinite volume)
    vector<double> m_bare;
    vector<pair<double, double> > m_tp; // \vb{k}^2 = (2\pi/L)^2 * (n_x^2 + n_y^2 + n _z^2)
    int nsq_max, n_bare, n_chan; // nsq == (n_x^2 + n_y^2 + n _z^2)
    double L;
    bool finished_calc;
    vector<int> momentum_ints;
    map<int, int> C3_map;
    Mat<cdouble> h; // the Hamiltonian
    vector<vector<cdouble (*)(double)> > g;
    vector<vector<cdouble (*)(double, double)> > v;
    void construct_hamiltonian(Mat<cdouble> & M);
    vector<double> spectrum;
    int C3(int nsq);
};
#endif