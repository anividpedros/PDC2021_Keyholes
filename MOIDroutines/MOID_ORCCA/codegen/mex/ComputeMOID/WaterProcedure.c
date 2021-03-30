/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * WaterProcedure.c
 *
 * Code generation for function 'WaterProcedure'
 *
 */

/* Include files */
#include "WaterProcedure.h"
#include "ComputeMOID.h"
#include "ComputeMOID_data.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void WaterProcedure(const emlrtStack *sp, real_T incliB, real_T omegaB, real_T
                    argpB, real_T vtrueB[4], real_T vL[4], real_T vdis[4],
                    real_T *N)
{
  real_T a_tmp;
  real_T b_a_tmp;
  real_T c_a_tmp;
  real_T d_a_tmp;
  real_T e_a_tmp;

  /*  "Water Procedure" */
  /*  To avoid missing minima . */
  /*  4 positions of the meridional plane, evenly distributed along the  */
  /*  inclined orbit. We move these points along as water droplets. */
  /* Special care when only 1 minima found */
  *N = 4.0;
  a_tmp = muDoubleScalarSin(omegaB);
  b_a_tmp = muDoubleScalarCos(omegaB);
  c_a_tmp = muDoubleScalarCos(incliB);
  vtrueB[0] = 2.3561944901923448;

  /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
  d_a_tmp = muDoubleScalarCos(argpB + 2.3561944901923448);
  e_a_tmp = muDoubleScalarSin(argpB + 2.3561944901923448);
  vL[0] = muDoubleScalarAtan2(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp * c_a_tmp,
    b_a_tmp * d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
  vdis[0] = 1.0E+6;

  /* set to something large */
  if (*emlrtBreakCheckR2012bFlagVar != 0) {
    emlrtBreakCheckR2012b(sp);
  }

  vtrueB[1] = 3.9269908169872414;

  /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
  d_a_tmp = muDoubleScalarCos(argpB + 3.9269908169872414);
  e_a_tmp = muDoubleScalarSin(argpB + 3.9269908169872414);
  vL[1] = muDoubleScalarAtan2(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp * c_a_tmp,
    b_a_tmp * d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
  vdis[1] = 1.0E+6;

  /* set to something large */
  if (*emlrtBreakCheckR2012bFlagVar != 0) {
    emlrtBreakCheckR2012b(sp);
  }

  vtrueB[2] = 5.497787143782138;

  /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
  d_a_tmp = muDoubleScalarCos(argpB + 5.497787143782138);
  e_a_tmp = muDoubleScalarSin(argpB + 5.497787143782138);
  vL[2] = muDoubleScalarAtan2(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp * c_a_tmp,
    b_a_tmp * d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
  vdis[2] = 1.0E+6;

  /* set to something large */
  if (*emlrtBreakCheckR2012bFlagVar != 0) {
    emlrtBreakCheckR2012b(sp);
  }

  vtrueB[3] = 7.0685834705770345;

  /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
  d_a_tmp = muDoubleScalarCos(argpB + 7.0685834705770345);
  e_a_tmp = muDoubleScalarSin(argpB + 7.0685834705770345);
  vL[3] = muDoubleScalarAtan2(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp * c_a_tmp,
    b_a_tmp * d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
  vdis[3] = 1.0E+6;

  /* set to something large */
  if (*emlrtBreakCheckR2012bFlagVar != 0) {
    emlrtBreakCheckR2012b(sp);
  }
}

void b_WaterProcedure(const emlrtStack *sp, real_T incliB, real_T omegaB, real_T
                      argpB, real_T *N, const real_T vtrueB[100], const real_T
                      vL[100], real_T vdis[4], real_T vtrueB_data[], int32_T
                      vtrueB_size[1], real_T vL_data[], int32_T vL_size[1])
{
  real_T a_tmp;
  real_T b_a_tmp;
  real_T c_a_tmp;
  real_T d_a_tmp;
  real_T e_a_tmp;
  vL_size[0] = 100;
  vtrueB_size[0] = 100;
  memcpy(&vL_data[0], &vL[0], 100U * sizeof(real_T));
  memcpy(&vtrueB_data[0], &vtrueB[0], 100U * sizeof(real_T));

  /*  "Water Procedure" */
  /*  To avoid missing minima . */
  /*  4 positions of the meridional plane, evenly distributed along the  */
  /*  inclined orbit. We move these points along as water droplets. */
  /* Special care when only 1 minima found */
  if (*N < 2.0) {
    vtrueB_size[0] = 4;
    vL_size[0] = 4;
    *N = 4.0;
    a_tmp = muDoubleScalarSin(omegaB);
    b_a_tmp = muDoubleScalarCos(omegaB);
    c_a_tmp = muDoubleScalarCos(incliB);
    vtrueB_data[0] = 2.3561944901923448;

    /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
    d_a_tmp = muDoubleScalarCos(argpB + 2.3561944901923448);
    e_a_tmp = muDoubleScalarSin(argpB + 2.3561944901923448);
    vL_data[0] = muDoubleScalarAtan2(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp *
      c_a_tmp, b_a_tmp * d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
    vdis[0] = 1.0E+6;

    /* set to something large */
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }

    vtrueB_data[1] = 3.9269908169872414;

    /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
    d_a_tmp = muDoubleScalarCos(argpB + 3.9269908169872414);
    e_a_tmp = muDoubleScalarSin(argpB + 3.9269908169872414);
    vL_data[1] = muDoubleScalarAtan2(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp *
      c_a_tmp, b_a_tmp * d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
    vdis[1] = 1.0E+6;

    /* set to something large */
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }

    vtrueB_data[2] = 5.497787143782138;

    /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
    d_a_tmp = muDoubleScalarCos(argpB + 5.497787143782138);
    e_a_tmp = muDoubleScalarSin(argpB + 5.497787143782138);
    vL_data[2] = muDoubleScalarAtan2(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp *
      c_a_tmp, b_a_tmp * d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
    vdis[2] = 1.0E+6;

    /* set to something large */
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }

    vtrueB_data[3] = 7.0685834705770345;

    /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
    d_a_tmp = muDoubleScalarCos(argpB + 7.0685834705770345);
    e_a_tmp = muDoubleScalarSin(argpB + 7.0685834705770345);
    vL_data[3] = muDoubleScalarAtan2(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp *
      c_a_tmp, b_a_tmp * d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
    vdis[3] = 1.0E+6;

    /* set to something large */
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }
}

/* End of code generation (WaterProcedure.c) */
