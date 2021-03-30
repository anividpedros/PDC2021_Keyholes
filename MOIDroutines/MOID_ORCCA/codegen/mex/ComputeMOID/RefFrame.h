/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * RefFrame.h
 *
 * Code generation for function 'RefFrame'
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
void RefFrame(const emlrtStack *sp, real_T A_i, real_T A_Omega, real_T A_argp,
              real_T B_i, real_T B_Omega, real_T B_argp, real_T *incliB, real_T *
              omegaB, real_T *argpB);

/* End of code generation (RefFrame.h) */
