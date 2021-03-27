/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ComputeMOID.h
 *
 * Code generation for function 'ComputeMOID'
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
void ComputeMOID(const emlrtStack *sp, const struct0_T *A, const struct0_T *B,
                 real_T *moid, real_T *vA_out, real_T *vB_out);

/* End of code generation (ComputeMOID.h) */
