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
#include "ComputeMOID_rtwutil.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void RefFrame(double A_i, double A_Omega, double A_argp, double B_i, double
              B_Omega, double B_argp, double *incliB, double *omegaB, double
              *argpB)
{
  double RA_tmp;
  double b_RA_tmp;
  double c_RA_tmp;
  double d_RA_tmp;
  double e_RA_tmp;
  double RA[9];
  double f_RA_tmp;
  double x_tmp_idx_0;
  double x_tmp_idx_1;
  int i;
  double zn[3];
  double x_tmp_idx_2;
  double b_RA[3];

  /*  Reference frame rotation */
  /*  3-1-3 Tranformation (orbital elements) */
  /* from reference to A frame */
  RA_tmp = sin(A_Omega);
  b_RA_tmp = cos(A_argp);
  c_RA_tmp = cos(A_Omega);
  d_RA_tmp = cos(A_i);
  e_RA_tmp = sin(A_argp);
  RA[0] = b_RA_tmp * c_RA_tmp - e_RA_tmp * RA_tmp * d_RA_tmp;
  RA[3] = RA_tmp * b_RA_tmp + c_RA_tmp * d_RA_tmp * e_RA_tmp;
  f_RA_tmp = sin(A_i);
  RA[6] = f_RA_tmp * e_RA_tmp;
  RA[1] = -c_RA_tmp * e_RA_tmp - RA_tmp * d_RA_tmp * b_RA_tmp;
  RA[4] = -e_RA_tmp * RA_tmp + b_RA_tmp * d_RA_tmp * c_RA_tmp;
  RA[7] = f_RA_tmp * b_RA_tmp;
  RA[2] = f_RA_tmp * RA_tmp;
  RA[5] = -f_RA_tmp * c_RA_tmp;
  RA[8] = d_RA_tmp;

  /*  Preparing the orbit */
  RA_tmp = sin(B_Omega);
  b_RA_tmp = cos(B_argp);
  c_RA_tmp = cos(B_Omega);
  d_RA_tmp = sin(B_argp);
  e_RA_tmp = cos(B_i);
  f_RA_tmp = sin(B_i);
  x_tmp_idx_0 = f_RA_tmp * RA_tmp;
  x_tmp_idx_1 = -f_RA_tmp * c_RA_tmp;
  for (i = 0; i < 3; i++) {
    zn[i] = (RA[i] * x_tmp_idx_0 + RA[i + 3] * x_tmp_idx_1) + RA[i + 6] *
      e_RA_tmp;
  }

  *incliB = rt_atan2d_snf(sqrt(zn[0] * zn[0] + zn[1] * zn[1]), zn[2]);
  *omegaB = -rt_atan2d_snf(zn[0], -zn[1]);
  x_tmp_idx_0 = c_RA_tmp * b_RA_tmp - d_RA_tmp * e_RA_tmp * RA_tmp;
  x_tmp_idx_1 = RA_tmp * b_RA_tmp + sin(B_argp) * cos(B_i) * c_RA_tmp;
  x_tmp_idx_2 = f_RA_tmp * d_RA_tmp;
  for (i = 0; i < 3; i++) {
    zn[i] = (RA[i] * x_tmp_idx_0 + RA[i + 3] * x_tmp_idx_1) + RA[i + 6] *
      x_tmp_idx_2;
  }

  x_tmp_idx_0 = -c_RA_tmp * d_RA_tmp - b_RA_tmp * e_RA_tmp * RA_tmp;
  x_tmp_idx_1 = -RA_tmp * d_RA_tmp + cos(B_argp) * cos(B_i) * c_RA_tmp;
  x_tmp_idx_2 = f_RA_tmp * b_RA_tmp;
  for (i = 0; i < 3; i++) {
    b_RA[i] = (RA[i] * x_tmp_idx_0 + RA[i + 3] * x_tmp_idx_1) + RA[i + 6] *
      x_tmp_idx_2;
  }

  *argpB = -rt_atan2d_snf(zn[2], b_RA[2]);
}

/* End of code generation (RefFrame.c) */
