#pragma once

#ifdef _WIN32
#ifdef fabber_cvr_EXPORTS
#define FABBER_CVR_API __declspec(dllexport)
#else
#define FABBER_CVR_API __declspec(dllimport)
#endif
#define CALL __stdcall
#else
#define FABBER_CVR_API
#define CALL
#endif

#include "fabber_core/fwdmodel.h"

extern "C" {
FABBER_CVR_API int CALL get_num_models();
FABBER_CVR_API const char *CALL get_model_name(int index);
FABBER_CVR_API NewInstanceFptr CALL get_new_instance_func(const char *name);
}
