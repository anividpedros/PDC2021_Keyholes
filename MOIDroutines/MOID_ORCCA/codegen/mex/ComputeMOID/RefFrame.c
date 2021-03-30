/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * RefFrame.c
 *
 * Code generation for function 'RefFrame'
 *
 */

/* Include files */
#include "RefFrame.h"
#include "ComputeMOID.h"
#include "ComputeMOID_data.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo s_emlrtRSI = { 35,  /* lineNo */
  "RefFrame",                          /* fcnName */
  "C:\\Users\\Oscar\\Documents\\GitHub\\PDC2021_Keyholes\\MOIDroutines\\MOID_MATLAB\\RefFrame.m"/* pathName */
};

/* Function Definitions */
void RefFrame(const emlrtStack *sp, real_T A_i, real_T A_Omega, real_T A_argp,
              real_T B_i, real_T B_Omega, real_T B_argp, real_T *incliB, real_T *
              omegaB, real_T *argpB)
{
  real_T RA_tmp;
  real_T b_RA_tmp;
  real_T c_RA_tmp;
  real_T d_RA_tmp;
  real_T e_RA_tmp;
  real_T RA[9];
  real_T f_RA_tmp;
  real_T x_tmp;
  real_T b_x_tmp;
  real_T c_x_tmp;
  int32_T i;
  real_T zn[3];
  real_T b_RA[3];
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;

  /*  Reference frame rotation */
  /*  3-1-3 Tranformation (orbital elements) */
  /* from reference to A frame */
  RA_tmp = muDoubleScalarSin(A_Omega);
  b_RA_tmp = muDoubleScalarCos(A_argp);
  c_RA_tmp = muDoubleScalarCos(A_Omega);
  d_RA_tmp = muDoubleScalarCos(A_i);
  e_RA_tmp = muDoubleScalarSin(A_argp);
  RA[0] = b_RA_tmp * c_RA_tmp - e_RA_tmp * RA_tmp * d_RA_tmp;
  RA[3] = RA_tmp * b_RA_tmp + c_RA_tmp * d_RA_tmp * e_RA_tmp;
  f_RA_tmp = muDoubleScalarSin(A_i);
  RA[6] = f_RA_tmp * e_RA_tmp;
  RA[1] = -c_RA_tmp * e_RA_tmp - RA_tmp * d_RA_tmp * b_RA_tmp;
  RA[4] = -e_RA_tmp * RA_tmp + b_RA_tmp * d_RA_tmp * c_RA_tmp;
  RA[7] = f_RA_tmp * b_RA_tmp;
  RA[2] = f_RA_tmp * RA_tmp;
  RA[5] = -f_RA_tmp * c_RA_tmp;
  RA[8] = d_RA_tmp;

  /*  Preparing the orbit */
  d_RA_tmp = muDoubleScalarSin(B_Omega);
  e_RA_tmp = muDoubleScalarCos(B_argp);
  f_RA_tmp = muDoubleScalarCos(B_Omega);
  x_tmp = muDoubleScalarSin(B_argp);
  b_x_tmp = muDoubleScalarCos(B_i);
  c_x_tmp = muDoubleScalarSin(B_i);
  RA_tmp = c_x_tmp * d_RA_tmp;
  b_RA_tmp = -c_x_tmp * f_RA_tmp;
  for (i = 0; i < 3; i++) {
    zn[i] = (RA[i] * RA_tmp + RA[i + 3] * b_RA_tmp) + RA[i + 6] * b_x_tmp;
  }

  st.site = &s_emlrtRSI;
  RA_tmp = zn[0] * zn[0] + zn[1] * zn[1];
  if (RA_tmp < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  *incliB = muDoubleScalarAtan2(muDoubleScalarSqrt(RA_tmp), zn[2]);
  *omegaB = -muDoubleScalarAtan2(zn[0], -zn[1]);
  RA_tmp = f_RA_tmp * e_RA_tmp - x_tmp * b_x_tmp * d_RA_tmp;
  b_RA_tmp = d_RA_tmp * e_RA_tmp + muDoubleScalarSin(B_argp) * muDoubleScalarCos
    (B_i) * f_RA_tmp;
  c_RA_tmp = c_x_tmp * x_tmp;
  for (i = 0; i < 3; i++) {
    zn[i] = (RA[i] * RA_tmp + RA[i + 3] * b_RA_tmp) + RA[i + 6] * c_RA_tmp;
  }

  RA_tmp = -f_RA_tmp * x_tmp - e_RA_tmp * b_x_tmp * d_RA_tmp;
  b_RA_tmp = -d_RA_tmp * x_tmp + muDoubleScalarCos(B_argp) * muDoubleScalarCos
    (B_i) * f_RA_tmp;
  c_RA_tmp = c_x_tmp * e_RA_tmp;
  for (i = 0; i < 3; i++) {
    b_RA[i] = (RA[i] * RA_tmp + RA[i + 3] * b_RA_tmp) + RA[i + 6] * c_RA_tmp;
  }

  *argpB = -muDoubleScalarAtan2(zn[2], b_RA[2]);
}

/* End of code generation (RefFrame.c) */
