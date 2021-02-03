/* cvr_models.cc 

Michael Chappell - IBME & FMRIB Analysis Group

Copyright (C) 2010-2011 University of Oxford */

/* CCOPYRIGHT  */

#include "fabber_core/fwdmodel.h"

#include <algorithm>

#include "fwdmodel_cvr.h"

#include "cvr_models.h"

extern "C" {

int CALL get_num_models()
{
    return 1;
}

const char *CALL get_model_name(int index)
{
    switch (index)
    {
    case 0:
        return "cvr_petco2";
        break;
    default:
        return NULL;
    }
}

NewInstanceFptr CALL get_new_instance_func(const char *name)
{
    if (string(name) == "cvr_petco2")
    {
        return CVRPETCO2Model::NewInstance;
    }
    else
    {
        return NULL;
    }
}
}
