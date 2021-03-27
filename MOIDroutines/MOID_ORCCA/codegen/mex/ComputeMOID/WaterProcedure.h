/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * WaterProcedure.h
 *
 * Code generation for function 'WaterProcedure'
 *
 */

#pragma once

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "ComputeMOID_types.h"

/* Function Declarations */
void WaterProcedure(const emlrtStack *sp, real_T incliB, real_T omegaB, real_T
                    argpB, real_T vtrueB[4], real_T vL[4], real_T vdis[4],
                    real_T *N);
void b_WaterProcedure(const emlrtStack *sp, real_T incliB, real_T omegaB, real_T
                      argpB, real_T *N, const real_T vtrueB[100], const real_T
                      vL[100], real_T vdis[4], real_T vtrueB_data[], int32_T
                      vtrueB_size[1], real_T vL_data[], int32_T vL_size[1]);

/* End of code generation (WaterProcedure.h) */
