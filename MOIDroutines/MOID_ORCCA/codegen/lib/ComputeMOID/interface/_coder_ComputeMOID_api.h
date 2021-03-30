/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_ComputeMOID_api.h
 *
 * Code generation for function '_coder_ComputeMOID_api'
 *
 */

#ifndef _CODER_COMPUTEMOID_API_H
#define _CODER_COMPUTEMOID_API_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"

/* Type Definitions */
#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  real_T sma;
  real_T e;
  real_T i;
  real_T Omega;
  real_T argp;
} struct0_T;

#endif                                 /*typedef_struct0_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void ComputeMOID(struct0_T *A, struct0_T *B, real_T *moid, real_T *vA_out,
  real_T *vB_out);
extern void ComputeMOID_api(const mxArray * const prhs[2], int32_T nlhs, const
  mxArray *plhs[3]);
extern void ComputeMOID_atexit(void);
extern void ComputeMOID_initialize(void);
extern void ComputeMOID_terminate(void);
extern void ComputeMOID_xil_shutdown(void);
extern void ComputeMOID_xil_terminate(void);

#endif

/* End of code generation (_coder_ComputeMOID_api.h) */
