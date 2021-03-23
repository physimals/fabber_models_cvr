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
    LOG << "CVRPETCO2Model::PCO2" << m_pco2 << endl;
}

/**
 * Preprocess CO2 data to extract end-tidal CO2 in mmHg at scan trigger times
 */
void CVRPETCO2Model::PreprocCO2()
{
#if 0
    // Physiological data stored in columns of data file
    ColumnVector timings = m_phys_data.Column(1);
    ColumnVector petco2 = m_phys_data.Column(2);
    ColumnVector peto2 = m_phys_data.Column(3);
    ColumnVector trig = m_phys_data.Column(4);

    // Determine number of volumes and TR from MR triggers
    // Note that we ignore first and last triggers?
    unsigned int first_trig_idx;
    ColumnVector trig_times;
    for (unsigned int idx=1; idx<=timings.nCols(); idx++)
    {
        if (trig(idx) > m_threshold_trig)
        {
            if (first_trig_idx < 0)
            {
                first_trig_idx = idx;
            }
            trig_times << timings(idx);
        }
    }

    ColumnVector trig_time_diff = trig_times.Rows(3, trig_times.Nrows()-1) - trig_times.Rows(2, trig_times.Nrows()-2);
    unsigned int vols = trig_time_diff.Nrows();
    double tr = trig_time_diff.Sum() / vols;
    cerr << "vols=" << vols << ", tr=" << tr << endl;

    // Temporal shift of end-tidal time courses by mechanical delay
    double trim_time_begin = trig_times(1) - self.delay // 1st trigger - delay (s)
    ColumnVector petco2_trim;
    for (unsigned int idx=1; idx<=timings.nCols(); idx++)
        if (timings(idx) >= trim_time_begin)
        {
            petco2_trim << petco2(idx);
        }
    }
    cerr << "len(petco2_trim)=" << petco2_trim.Nrows() << endl;
    
    // Determined respiratory frequency during baseline and use info to
    // determine size of end-tidal search window
    double samp_period = 1/m_samp_rate
    int baseline_vols = int(m_baseline * m_samp_rate)
    cerr << "baseline vols=" << baseline_vols << endl;
    ColumnVector petco2_trim_postbaseline = petco2_trim.Rows(1, baseline_vols);
    ColumnVector zeros = petco2_trim_postbaseline * 0;
    ColumnVector bfft_r, bfft_i;
    FFT(petco2_trim_postbaseline, zeros, bfft_r, bfft_i);
    ColumnVector p2 = (bfft_r.Square() + bfft_i.Square()).Sqrt()/baseline_vols;
    ColumnVector p1 = p2.Rows(1, int(baseline_vols/2)+1)
    ColumnVector f;
    for (unsigned int idx=1; idx<=p1.Ncols(); idx++)
    {
        if ((idx > 1) && (idx < pt.Ncols()))
        {
            p1(idx) = 2*p1(idx);
        }
        f << (idx-1)*(self.samp_rate/2)/(p1.Ncols()-1);
    }
    cerr << "f=" << f.t();
    cerr << "p1=" << p1.t();

    loc = np.argmax(p1[1:])
    pk = p1[loc+1]

    pkloc = loc+1
    harm = f[pkloc]
    resp_period = round(1/harm) # e.g. 8s

    # Search window = 1 second more than the respiratory period
    nsearch_vols = (resp_period+1)*self.samp_rate
    windows = int(np.floor(self.petco2_trim.shape[0]/nsearch_vols))

    # Find peak PETCO2 in each window - it's value and index position
    posmax = np.zeros(windows, dtype=np.int)
    winmax = np.zeros(windows)
    k=0
    for i in range(windows):
        for j in range(nsearch_vols):
            if j == 0 or self.petco2_trim[i*nsearch_vols+j] > winmax[i]:
                winmax[i] = self.petco2_trim[i*nsearch_vols+j]
                posmax[i] = i*nsearch_vols+j

    # Make new full sample ET time course where the PETCO2 changes linearly
    # between window maxima
    self.petco2_resamp = np.zeros((self.petco2_trim.shape[0], 1))
    for x in range(windows-1):
        dist_c = posmax[x+1] - posmax[x]
        step_c = winmax[x+1] - winmax[x]
        ramp_c = step_c / dist_c
        for g in range(dist_c+1):
            self.petco2_resamp[posmax[x]+g] = winmax[x] + (ramp_c * g)

    # Pad the start and end with repeats of first and last value to maintain
    # length and phase
    self.petco2_resamp[:posmax[0]] = self.petco2_resamp[posmax[0]]
    self.petco2_resamp[posmax[-1]:] = self.petco2_resamp[posmax[-1]]

    # Create a timecourse of the end tidal CO2 values at the TR's for use with CVR sigmoids
    # Make new time course at the TR resolution and normalise timecourse betwwen 0 and 1 to create EV
    block = round(tr*self.samp_rate)
    ev_co2 = np.zeros((vols,), dtype=np.float32)
    for i in range(vols):
        ev_co2[i] = self.petco2_resamp[block * i + block-1]

    # Calculate normo/hypercapnea in mmHg
    self.air_pressure_mmhg = self.air_pressure/1.33322387415 # pressure mbar
    self.co2_mmHg = (ev_co2 * self.air_pressure_mmhg) / 100 # div by 100 as values are in percent

    # Differences between timepoints for quick interpolation. Given a delay
    # time > 0 we can compute co2 = co2[int(delay)] + frac(delay) * diff[int(delay)]
    self.co2_diff = np.zeros(len(self.co2_mmHg), dtype=np.float32)
    self.co2_diff[:-1] = self.co2_mmHg[1:] - self.co2_mmHg[:-1]

    # Convert time periods to number of volumes
    baseline_vols = self.baseline/tr
    blocksize_on_vols = self.blocksize_on/tr
    blocksize_off_vols = self.blocksize_off/tr

    # Average all of first baseline block
    self.normocap = np.mean(self.co2_mmHg[:int(baseline_vols+self.delay)])

    s1 = (baseline_vols+self.delay+blocksize_on_vols/2)
    s2 = (baseline_vols+self.delay+blocksize_on_vols)
    s3 = (baseline_vols+self.delay+blocksize_on_vols+blocksize_off_vols+blocksize_on_vols/2)
    s4 = (baseline_vols+self.delay+blocksize_on_vols+blocksize_off_vols+blocksize_on_vols)
    s1, s2, s3, s4 = int(s1), int(s2), int(s3), int(s4)
    # Select 2nd half of each hypercapnic block to average
    hyperblock = np.concatenate([self.co2_mmHg[s1-1:s2], self.co2_mmHg[s3-1:s4]])
    self.hypercap = np.mean(hyperblock)
#endif
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
        result(i) = sig0 * (1 + cvr * PCO2(t) / 100);
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
