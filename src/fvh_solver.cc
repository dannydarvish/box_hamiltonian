#include "fvh_solver.h"
#include <algorithm>
#include <assert.h>

void FVSpectrumSolver::construct_hamiltonian(Mat<cdouble> & M)
{
    int dim = n_chan * momentum_ints.size() + n_bare;
    M.zeros(dim, dim);
    // Non-interacting part
    int pos = 0;
    for (int i = 0; i < n_bare; i++)
    {
        h(pos, pos) += m_bare[i];
        pos++;
    }
    for (int nsq : momentum_ints)
        for (int alpha = 0; alpha < n_chan; alpha++)
        {
            M(pos, pos) += sqrt(sqr(m_tp[alpha].first)  + sqr(2 * M_PI / L) * nsq) +
                           sqrt(sqr(m_tp[alpha].second) + sqr(2 * M_PI / L) * nsq);
            pos++;
        }
    assert(pos == dim);
    // Interacting part
    int row(0), col;
    for (int i = 0; i < n_bare; i++)
    {
        col = 0;
        for (int j = 0; j < n_bare; j++)
            col++;
        for (int nsq : momentum_ints)
            for (int beta = 0; beta < n_chan; beta++)
            {
                M(row, col) += sqrt(double(C3_map[nsq])/(4.0*M_PI)) * pow(2*M_PI/L,1.5) *
                                (*g_ptr)(i, beta, 2*M_PI*sqrt(nsq)/L);
                col++;
            }
        row++;
    }
    for (int nsq : momentum_ints)
        for (int alpha = 0; alpha < n_chan; alpha++)
        {
            col = 0;
            for (int j = 0; j < n_bare; j++)
            {
                M(row, col) += conj(sqrt(double(C3_map[nsq])/(4.0*M_PI)) * pow(2*M_PI/L,1.5) *
                                (*g_ptr)(j, alpha, 2*M_PI*sqrt(nsq)/L));
                col++;
            }
            for (int msq : momentum_ints)
                for (int beta = 0; beta < n_chan; beta++)
                {
                    M(row, col) += sqrt(C3_map[nsq]/(4*M_PI))*sqrt(C3_map[msq]/(4*M_PI)) *
                                   pow(2*M_PI/L, 3.0) * 
                                   (*v_ptr)(alpha, beta, 2*M_PI/L*sqrt(nsq), 2*M_PI/L*sqrt(msq));
                    col++;
                }
            row++;
        }
    assert(row == dim);
    assert(col == dim);
    finished_calc = true;
}

vector<double> FVSpectrumSolver::solve_spectrum()
{
    construct_hamiltonian(h);
    assert(h.is_symmetric(1e-6*abs(h(0,0))));
    vec eigval;
    eig_sym(eigval, h);
    return conv_to<vector<double> >::from(eigval);
}

int FVSpectrumSolver::C3(int nsq)
{
    int count = 0;
    for (int i = -1*(sqrt(nsq) + 1); i <= sqrt(nsq) + 1; i++)
        for (int j = -1*(sqrt(nsq) + 1); j <= sqrt(nsq) + 1; j++)
            for (int k = -1*(sqrt(nsq) + 1); k <= sqrt(nsq) + 1; k++)
                if (sqr(i)+sqr(j)+sqr(k) == nsq)
                    count++;
    return count;
}

void FVSpectrumSolver::set_L(double L)
{
    finished_calc = false;
    this->L = L;
}

void FVSpectrumSolver::set_g_ptr(const gBase * g_ptr)
{
    this->g_ptr = g_ptr;
}

void FVSpectrumSolver::set_v_ptr(const vBase * v_ptr)
{
    this->v_ptr = v_ptr;
}