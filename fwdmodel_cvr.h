/**
 * fwdmodel_cvr.h
 *
 * CVR forward model
 *
 * Martin Craig
 *
 * Copyright (C) 2021 University of Nottingham  
 */

/*  CCOPYRIGHT */
#pragma once

#include "fabber_core/fwdmodel.h"

#include <newmat.h>

#include <string>
#include <vector>

/**
 * CVR-PETCO2 forward model
 */
class CVRPETCO2Model : public FwdModel
{
public:
    static FwdModel *NewInstance();

    virtual ~CVRPETCO2Model()
    {
    }

    std::string GetDescription() const;
    std::string ModelVersion() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;

    void Initialize(FabberRunData &rundata);
    void GetParameterDefaults(std::vector<Parameter> &params) const;
    
    void InitVoxelPosterior(MVNDist &posterior) const;
    void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;
    
protected:
    // Physiological data timecourse
    NEWMAT::Matrix m_phys_data;

    // CVR-PETCO2 configuration
    double m_baseline, m_blocksize_on, m_blocksize_off, m_air_pressure, m_threshold_trig, m_delay, m_samp_rate;
    
    // Inference flags
    bool m_infer_delay, m_infer_sig0;
    
    NEWMAT::ColumnVector m_pco2;

    // Hyper/normocapnia difference
    double m_diff_pco2;

    double PCO2(double t) const;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, CVRPETCO2Model> registration;
};
