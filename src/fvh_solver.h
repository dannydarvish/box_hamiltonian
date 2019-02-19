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

class gBase
{
public:
    // The user will implement a function for the 1-to-2 coupling, g_{i, \alpha}(k)
    virtual cdouble operator() (uint i, uint alpha, double k) const = 0;
};

class vBase
{
public:
    // The user will implement a function for the 2-to-2 coupling, v_{\alpha, \beta}(k, k^\prime)
    virtual cdouble operator() (uint alpha, uint beta, double k, double kp) const = 0;
};

// Notation here follows the Wu et. al paper
// Unless otherwise stated, units are in MeV and fm.
class FVSpectrumSolver
{
public:
    FVSpectrumSolver(int nsq_max, double L, const gBase * g_ptr, const vBase * v_ptr,
                    const vector<double> & m_bare, const vector<pair<double, double> > & m_tp) :
                    nsq_max(nsq_max), L(L), g_ptr(g_ptr), v_ptr(v_ptr), m_bare(m_bare), m_tp(m_tp),
                    n_bare(m_bare.size()), n_chan(m_tp.size()), finished_calc(false)
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
    vector<double> solve_spectrum();
    void set_L(double L);
    void set_g_ptr(const gBase * g_ptr);
    void set_v_ptr(const vBase * v_ptr);
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
    const gBase * g_ptr;
    const vBase * v_ptr;
    void construct_hamiltonian(Mat<cdouble> & M);
    int C3(int nsq);
};
#endif