#include "fwdmodel_cvr.h"

#include <fabber_core/easylog.h>
#include <fabber_core/priors.h>

#include <miscmaths/miscprob.h>
#include <newimage/newimageall.h>

#include <newmatio.h>

#include <iostream>
#include <stdexcept>

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, CVRPETCO2Model> CVRPETCO2Model::registration("cvr_petco2");

std::string CVRPETCO2Model::GetDescription() const
{
    return "CVR-PETCO2 forward model";
}

static OptionSpec OPTIONS[] = {
    { "phys-data", OPT_MATRIX, "Physiological data file", OPT_REQ, "" },
    { "baseline", OPT_FLOAT, "Length of initial baseline block (s)", OPT_NONREQ, "60" },
    { "blocksize-on", OPT_FLOAT, "Length of ON block (s)", OPT_NONREQ, "120" },
    { "blocksize-off", OPT_FLOAT, "Length of OFF block (s)", OPT_NONREQ, "120" },
    { "samp-rate", OPT_FLOAT, "Powerlab sampling rate (Hz)", OPT_NONREQ, "100" },
    { "air-pressure", OPT_FLOAT, "Barometric pressure (mbar)", OPT_NONREQ, "1020" },
    { "threshold-trig", OPT_FLOAT, "Threshold to detect triggers", OPT_NONREQ, "3" },
    { "delay", OPT_FLOAT, "Mechanical delay", OPT_NONREQ, "15" },
    { "infer-sig0", OPT_BOOL, "Infer baseline signal", OPT_NONREQ, "" },
    { "infer-delay", OPT_BOOL, "Infer the delay parameter", OPT_NONREQ, "" },
    { "" },
};

void CVRPETCO2Model::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string CVRPETCO2Model::ModelVersion() const
{
    string version = "Fabber CVR models: ";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

void CVRPETCO2Model::Initialize(FabberRunData &rundata)
{
    // Mandatory parameters
    m_phys_data = read_ascii_matrix(rundata.GetString("phys-data"));
    m_baseline = rundata.GetDoubleDefault("baseline", 60);
    m_blocksize_on = rundata.GetDoubleDefault("blocksize-on", 120);
    m_blocksize_off = rundata.GetDoubleDefault("blocksize-off", 120);
    m_samp_rate = rundata.GetDoubleDefault("samp-rate", 100);
    m_air_pressure = rundata.GetDoubleDefault("air-pressure", 1020);
    m_threshold_trig = rundata.GetDoubleDefault("threshold-trig", 3);
    m_delay = rundata.GetDoubleDefault("delay", 15);

    m_infer_sig0 = rundata.ReadBool("infer-sig0");
    m_infer_delay = rundata.ReadBool("infer-delay");

    // TEMP until we have PCO2 preprocessing
    m_pco2 = m_phys_data.Column(1);
    m_diff_pco2 = 46.992615 - 41.151955;
    cerr << m_pco2 << endl;
    cerr << endl << m_diff_pco2 << endl;
}

void CVRPETCO2Model::GetParameterDefaults(std::vector<Parameter> &params) const
{
    unsigned int p = 0;
    params.push_back(Parameter(p++, "cvr", DistParams(1, 2000), DistParams(1, 10), PRIOR_NORMAL, TRANSFORM_ABS()));
    if (m_infer_sig0)
        params.push_back(Parameter(p++, "sig0", DistParams(0, 1e9), DistParams(1, 10), PRIOR_NORMAL, TRANSFORM_ABS()));
    if (m_infer_delay)
        params.push_back(Parameter(p++, "delay", DistParams(0, 100), DistParams(0, 10)));
}

void CVRPETCO2Model::InitVoxelPosterior(MVNDist &posterior) const
{
/*    int delay_idx = m_sig0_idx;
    if (m_infer_sig0) {
        double sig0_data = data(1);
        posterior.means(m_sig0_idx+1) = sig0_data;
        delay_idx++;
    }
    if (m_infer_delay && m_auto_init_delay) {
        posterior.means(delay_idx+1) = m_dt * fit_step(data);
    }*/
}

void CVRPETCO2Model::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    // Parameters that are inferred - extract and give sensible names
    int p = 1;
    double cvr = params(p++);
    double sig0 = 0;
    double delay = 0;

    // Optional parameters to infer
    if (m_infer_sig0) {
        sig0 = params(p++);   
    }

    if (m_infer_delay)
    {
        delay = params(p++);
    }

    // Convert CVR back to BOLD signal
    result.ReSize(data.Nrows());
    for (int i = 1; i <= data.Nrows(); i++)
    {
        double t = double(i) - 1 - delay;
        result(i) = sig0 * (1 + cvr * m_diff_pco2 * PCO2(t) / 100);
        //cerr << t << ", " << PCO2(t) << ", " << result(i) << endl;
    }
}

double CVRPETCO2Model::PCO2(double t) const
{
    if (t <= 0) 
    {
        return 0;
    }
    else 
    {
        // Linearly interpolate PETCO2 curve. Assume
        // constant before and after last point
        int t_idx_0 = int(floor(t));
        double t_idx_frac = t - double(t_idx_0);
        if (t_idx_0 < 0)
        {
            return m_pco2(1);
        }
        if (t_idx_0 < m_pco2.Nrows()-1)
        {
            return t_idx_frac*m_pco2(t_idx_0+2) + (1-t_idx_frac)*m_pco2(t_idx_0+1);
        }
        else 
        {
            return m_pco2(m_pco2.Nrows());
        }
    }   
}

FwdModel *CVRPETCO2Model::NewInstance()
{
    return new CVRPETCO2Model();
}
