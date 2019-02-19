#include "../src/spectrum_fitter.h"
#include <armadillo>
#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <assert.h>
const double hbar_c_over_fm = 197.327;

class gDerived : public gBase
{
public:
    gDerived(const vector<double> & params) :
        g_sigma_pi_pi(params[0]), c_sigma_pi_pi(params[1]) {}
    cdouble operator() (uint i, uint alpha, double k) const override
    {
        if (i > 0 || alpha > 0)
            throw(range_error("i or alpha is out of bounds for g."));
        return g_00(k);
    }
    void set_params(const vector<double> & params)
    {
        g_sigma_pi_pi = params[0];
        c_sigma_pi_pi = params[1];
    }
private:
    double g_sigma_pi_pi;
    double c_sigma_pi_pi;
    double m_pi = 134.976;
    double fm_times_MeV = 0.005067731;

    cdouble g_00(double k) const
    {
        return g_sigma_pi_pi/sqrt(m_pi)/(1+sqr(c_sigma_pi_pi*k*fm_times_MeV));
    }
};

class vDerived : public vBase
{
public:
    vDerived(const vector<double> & params) :
            G_pi_pi_pi_pi(params[2]), d_pi_pi(params[3]) {}
    cdouble operator() (uint alpha, uint beta, double k, double kp) const override
    {
        if (alpha > 0 || beta > 0)
            throw(range_error("alpha or beta is out of bounds for v."));
        return v_00(k, kp);
    }
    void set_params(const vector<double> & params)
    {
        G_pi_pi_pi_pi = params[2];
        d_pi_pi = params[3];
    }
private:
    double G_pi_pi_pi_pi;
    double d_pi_pi;
    double m_pi = 134.976;
    double fm_times_MeV = 0.005067731;

    cdouble v_00(double k, double kp) const
    {
        return G_pi_pi_pi_pi/sqr(m_pi) * 1.0/sqr(1+sqr(d_pi_pi*k*fm_times_MeV)) *
            1.0/sqr(1+sqr(d_pi_pi*kp*fm_times_MeV));
    }
};

class ModelSpectrumDerived : public ModelSpectrumBase
{
public:
    ModelSpectrumDerived(uint num_levels) :
                             g({0,0,0,0}), v({0,0,0,0}),
                             m_bare({700.0}), m_tp({{134.976, 134.976}}),
                             fvss(100, 4.0 / hbar_c_over_fm, &(this->g), &(this->v), m_bare, m_tp)
    {
        this->num_levels = num_levels;
    }
    vector<double> operator() (const vector<double> & params)
    {
        g.set_params(params);
        v.set_params(params);
        const vector<double> spec = fvss.solve_spectrum();
        vector<double> resized_spec(num_levels);
        assert(spec.size() >= num_levels);
        for (uint i = 0; i < num_levels; i++)
            resized_spec[i] = spec[i];
        return resized_spec;
    }

private:
    uint num_levels;
    vector<double> m_bare;
    vector<pair<double, double> > m_tp;
    gDerived g;
    vDerived v;
    FVSpectrumSolver fvss;
};

int main()
{
    vector<double> initial_guess = {1,1,1,1};

    vector<vector<double> > spec_data =
        {{238.123,585.229,805.556,1019.33,1146.15,1279.1,1431.07,1554.46},
        {255.057193792,597.147840504,827.219541128,1026.09053447,1173.60910805,1303.74656926,1443.5481989,1560.36131749},
        {252.738957009,591.971920209,823.588316762,1024.78756813,1170.91422321,1282.71782832,1455.65880686,1562.89140883},
        {251.277972609,610.880999009,813.827670532,1039.37171301,1147.61229462,1293.90634237,1435.49832664,1555.0978255},
        {267.176933305,585.288050643,808.066368317,1048.52271245,1175.01277602,1293.55832325,1443.71685108,1571.75702411},
        {250.819011771,592.384884457,810.861085506,1048.63313811,1174.42653512,1301.34867463,1439.44789285,1565.49841845},
        {262.445800388,594.03815345,816.409573879,1022.59960291,1147.23214009,1296.05770383,1443.72095154,1555.63727455},
        {243.229664393,586.539980363,816.058137702,1036.53193278,1158.00920598,1305.80488564,1441.61752449,1575.12558267},
        {242.933968015,610.829028832,823.061061486,1024.4193765,1149.50657935,1289.52412762,1460.04344577,1578.81424938},
        {257.864802717,600.478546067,806.366667734,1028.01623189,1154.77765142,1290.72444052,1434.93270867,1564.84533377},
        {264.691026784,611.396397765,830.672735808,1046.77315017,1146.49881667,1292.11444102,1438.44610652,1577.35711676},
        {255.732683944,587.728652661,806.39419053,1030.16477223,1150.80902021,1299.646062,1461.06766509,1572.22300486},
        {251.873915322,589.542411035,819.169634306,1042.42281873,1170.99332855,1283.3950137,1455.20260882,1574.31269583},
        {256.637332439,589.485749813,822.417363489,1044.80559082,1151.45186744,1295.58678355,1447.16783778,1583.39064123},
        {250.814750388,612.364199883,814.138197672,1044.07895984,1154.05558635,1290.99719935,1458.97223842,1554.77017742},
        {261.927919186,607.694355063,830.048559841,1031.67288857,1172.93450522,1293.79847223,1449.59428598,1560.92456198},
        {238.409927482,611.419848049,824.753200934,1034.07374739,1150.40963734,1293.62354583,1450.85770247,1560.15716647},
        {254.95148913,594.509928068,805.795324151,1034.99582211,1167.55299112,1304.9648436,1443.81877902,1571.28133652},
        {258.052599408,601.288449715,830.011159426,1044.32645262,1160.29734941,1283.44682766,1452.38125254,1561.6387566},
        {246.3291879,598.630485425,820.925369631,1024.76035915,1149.36126135,1290.87938551,1440.82766595,1582.41813363},
        {245.522782095,587.691796898,815.05210258,1043.22187694,1151.95535027,1281.69473981,1447.24101818,1568.21828117}};

    ModelSpectrumDerived model_spectrum(spec_data[0].size());
    SpectrumFitter sf(spec_data, BOOTSTRAP, initial_guess, &model_spectrum);
    vector<double> result_params = sf.calc_result_params();
    for (auto a : result_params)
        cout << a << endl;
    return 0;
}