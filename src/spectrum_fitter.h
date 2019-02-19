// TO DO: Make use of friend class functionality in Residual
#ifndef SPECTRUM_FITTER
#define SPECTRUM_FITTER
#include <armadillo>
#include <complex>
#include "fvh_solver.h"
#include <vector>

using namespace std;

typedef complex<double> cdouble;
typedef unsigned int uint;

enum Mode {BOOTSTRAP, JACKKNIFE};

class ModelSpectrumBase
{
public:
    virtual vector<double> operator() (const vector<double> & params) = 0;
};

class SpectrumFitter
{
public:
    SpectrumFitter(vector<vector<double> > & spec_data, Mode mode,
                   const vector<double> & initial_guess, ModelSpectrumBase * model_spectrum_ptr) :
        spec_data(spec_data), mode(mode), num_fit_params(initial_guess.size()), initial_guess(initial_guess),
        model_spectrum_ptr(model_spectrum_ptr), calculated_covariance(false)
    {
        if (initial_guess.size() != num_fit_params)
            throw length_error("Initial guess has size " + to_string(initial_guess.size()) +
                                    ", but num_fit_params is " + to_string(num_fit_params) + ".");
        for (vector<double> & spec : spec_data)
            if (spec.size() != spec_data[0].size())
                throw length_error("Number of levels is not the same across all sample mean and all resamplings.");
        cov.resize(spec_data[0].size(), spec_data[0].size());
        inv_cov.resize(spec_data[0].size(), spec_data[0].size());
    }
    ~SpectrumFitter() {};
    vector<double> calc_result_params();
    vector<double> calc_result_covariance();
private:
    bool calculated_covariance;
    // Store the spectrum as a vector of vector<double>s.  Each row is a different sample,
    // each column is a different level.
    // The first vector is the full sample mean of the levels.
    // Every other vector is over the individual resamplings.
    vector<vector<double> > spec_data;
    const Mode mode; // resampling mode
    uint num_fit_params;
    vector<double> initial_guess;
    mat cov;     // Covariance matrix
    mat inv_cov; // Inverse of the covariance matrix
    // This user-implemented functor takes in a vector of parameters
    // and returns the model Hamiltonian spectrum.
    ModelSpectrumBase * model_spectrum_ptr;
    // The first in the pair is the energy, the second is its uncertainty
    bool find_params(uint spec_data_index, vector<double> & result_params);
    void calculate_covariance(); // Caluclate both the covariance and the inverse covariance
};
#endif