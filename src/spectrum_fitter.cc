#include "spectrum_fitter.h"
#include "fvh_solver.h"
#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnPrint.h"

int CHANGE_ME_MAX_ITS = 1000;
double CHANGE_ME_REL_TOL = 1e-2;

ostream & operator<<(ostream & out, const vector<double> & vec)
{
    for (double v : vec)
        out << v << " ";
    return out;
}

class Residual : public ROOT::Minuit2::FCNBase
{
public:
    Residual(ModelSpectrumBase * model_spectrum_ptr, const vector<double> * lat_spec_ptr,
             const mat & inv_cov) :
             model_spectrum_ptr(model_spectrum_ptr), lat_spec_ptr(lat_spec_ptr), inv_cov(inv_cov) {}
    double operator() (const vector<double> & params) const
    {
        cout << params << endl;
        vector<double> model_spec = (*model_spectrum_ptr)(params);
        
        assert(inv_cov.n_rows == lat_spec_ptr->size() && inv_cov.n_rows == inv_cov.n_cols);
        assert(model_spec.size() >= (*lat_spec_ptr).size());
        sort(model_spec.begin(), model_spec.end());
        assert(is_sorted((*lat_spec_ptr).begin(), (*lat_spec_ptr).end()));

        double chisq = 0.0;
        for (uint row = 0; row < lat_spec_ptr->size(); row++)
        {
            for (uint col = 0; col < lat_spec_ptr->size(); col++)
                chisq += (model_spec[row]-(*lat_spec_ptr)[row])*(model_spec[col]-(*lat_spec_ptr)[col]) *
                         inv_cov(row, col);
        }
        return chisq;
    }

    void set_model_spec_func(ModelSpectrumBase * model_spectrum_ptr)
    {this->model_spectrum_ptr = model_spectrum_ptr;}

    void set_lat_spec_ptr(const vector<double> * lat_spec_ptr)
    {this->lat_spec_ptr = lat_spec_ptr;}

    double Up() const {return 1.0;}

private:
    ModelSpectrumBase * model_spectrum_ptr;
    // This needs to be a pointer since we need to be able to reassign it.
    const vector<double> * lat_spec_ptr;
    const mat & inv_cov;
};

bool SpectrumFitter::find_params(uint spec_data_index, vector<double> & result_params)
{
    result_params.resize(num_fit_params);
    if (!calculated_covariance)
        calculate_covariance();

    unsigned int strategylevel=2;  // 0 = low, 1 = med, 2 = high quality
                                   // lower level means faster, higher means
                                   // more reliable minimization
    vector<double> unc(initial_guess.size());
    for (uint i = 0; i < initial_guess.size(); i++)
        unc[i] = 0.1 * initial_guess[i];
    Residual residual(model_spectrum_ptr, &spec_data[spec_data_index], inv_cov);
    ROOT::Minuit2::MnMinimize M(residual, initial_guess, unc, strategylevel);
    ROOT::Minuit2::FunctionMinimum result = M(CHANGE_ME_MAX_ITS, CHANGE_ME_REL_TOL);
    if (result.IsValid())
        for (int i = 0; i < result_params.size(); i++)
            result_params[i] = result.UserParameters().Value(i);
    else
        return false;
    return true;
}

// We are freezing the covariance matrix here. Tests should be done to see if
// this is a valid thing to do.
void SpectrumFitter::calculate_covariance()
{
    int num_levels = spec_data[0].size();
    int num_resamplings = spec_data.size() - 1;
    vector<double> resampled_mean(num_levels, 0.0);
    for (uint i = 0; i < num_resamplings; i++)
        for (uint level = 0; level < num_levels; level++)
            {
                resampled_mean[level] += spec_data[i+1][level] / double(num_resamplings);
            }
    if (mode == BOOTSTRAP)
    {
        for (uint row = 0; row < num_levels; row++)
            for (uint col = 0; col < num_levels; col++)
            {
                for (uint b = 0; b < num_resamplings; b++)
                   cov(row, col) +=  (spec_data[b+1][row] - resampled_mean[row]) *
                                     (spec_data[b+1][col] - resampled_mean[col]);
                cov(row, col) /= double(num_resamplings - 1);
            }
    }
    else if (mode == JACKKNIFE)
    {
        for (uint row = 0; row < num_levels; row++)
            for (uint col = 0; col < num_levels; col++)
            {
                for (uint b = 0; b < num_resamplings; b++)
                   cov(row, col) +=  (spec_data[b+1][row] - resampled_mean[row]) *
                                     (spec_data[b+1][col] - resampled_mean[col]);
                cov(row, col) *=  double(num_resamplings - 1) /
                                  double(num_resamplings);
            }
    }
    inv_cov = cov.i();
    calculated_covariance = true;
}

vector<double> SpectrumFitter::calc_result_params()
{
    vector<double> params;
    bool success = find_params(0, params);
    if(success)
        return params;
    else throw(runtime_error("Could not find a minimum for the full sample mean."));
}