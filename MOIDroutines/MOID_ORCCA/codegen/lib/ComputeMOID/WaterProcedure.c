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
#include "ComputeMOID_rtwutil.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
void WaterProcedure(double incliB, double omegaB, double argpB, double vtrueB[4],
                    double vL[4], double vdis[4], double *N)
{
  double a_tmp;
  double b_a_tmp;
  double c_a_tmp;
  double d_a_tmp;
  double e_a_tmp;

  /*  "Water Procedure" */
  /*  To avoid missing minima . */
  /*  4 positions of the meridional plane, evenly distributed along the  */
  /*  inclined orbit. We move these points along as water droplets. */
  /* Special care when only 1 minima found */
  *N = 4.0;
  a_tmp = sin(omegaB);
  b_a_tmp = cos(omegaB);
  c_a_tmp = cos(incliB);
  vtrueB[0] = 2.3561944901923448;

  /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
  d_a_tmp = cos(argpB + 2.3561944901923448);
  e_a_tmp = sin(argpB + 2.3561944901923448);
  vL[0] = rt_atan2d_snf(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp * c_a_tmp, b_a_tmp *
                        d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
  vdis[0] = 1.0E+6;

  /* set to something large */
  vtrueB[1] = 3.9269908169872414;

  /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
  d_a_tmp = cos(argpB + 3.9269908169872414);
  e_a_tmp = sin(argpB + 3.9269908169872414);
  vL[1] = rt_atan2d_snf(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp * c_a_tmp, b_a_tmp *
                        d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
  vdis[1] = 1.0E+6;

  /* set to something large */
  vtrueB[2] = 5.497787143782138;

  /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
  d_a_tmp = cos(argpB + 5.497787143782138);
  e_a_tmp = sin(argpB + 5.497787143782138);
  vL[2] = rt_atan2d_snf(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp * c_a_tmp, b_a_tmp *
                        d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
  vdis[2] = 1.0E+6;

  /* set to something large */
  vtrueB[3] = 7.0685834705770345;

  /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
  d_a_tmp = cos(argpB + 7.0685834705770345);
  e_a_tmp = sin(argpB + 7.0685834705770345);
  vL[3] = rt_atan2d_snf(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp * c_a_tmp, b_a_tmp *
                        d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
  vdis[3] = 1.0E+6;

  /* set to something large */
}

void b_WaterProcedure(double incliB, double omegaB, double argpB, double *N,
                      const double vtrueB[100], const double vL[100], double
                      vdis[4], double vtrueB_data[], int vtrueB_size[1], double
                      vL_data[], int vL_size[1])
{
  double a_tmp;
  double b_a_tmp;
  double c_a_tmp;
  double d_a_tmp;
  double e_a_tmp;
  vL_size[0] = 100;
  vtrueB_size[0] = 100;
  memcpy(&vL_data[0], &vL[0], 100U * sizeof(double));
  memcpy(&vtrueB_data[0], &vtrueB[0], 100U * sizeof(double));

  /*  "Water Procedure" */
  /*  To avoid missing minima . */
  /*  4 positions of the meridional plane, evenly distributed along the  */
  /*  inclined orbit. We move these points along as water droplets. */
  /* Special care when only 1 minima found */
  if (*N < 2.0) {
    vtrueB_size[0] = 4;
    vL_size[0] = 4;
    *N = 4.0;
    a_tmp = sin(omegaB);
    b_a_tmp = cos(omegaB);
    c_a_tmp = cos(incliB);
    vtrueB_data[0] = 2.3561944901923448;

    /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
    d_a_tmp = cos(argpB + 2.3561944901923448);
    e_a_tmp = sin(argpB + 2.3561944901923448);
    vL_data[0] = rt_atan2d_snf(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp * c_a_tmp,
      b_a_tmp * d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
    vdis[0] = 1.0E+6;

    /* set to something large */
    vtrueB_data[1] = 3.9269908169872414;

    /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
    d_a_tmp = cos(argpB + 3.9269908169872414);
    e_a_tmp = sin(argpB + 3.9269908169872414);
    vL_data[1] = rt_atan2d_snf(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp * c_a_tmp,
      b_a_tmp * d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
    vdis[1] = 1.0E+6;

    /* set to something large */
    vtrueB_data[2] = 5.497787143782138;

    /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
    d_a_tmp = cos(argpB + 5.497787143782138);
    e_a_tmp = sin(argpB + 5.497787143782138);
    vL_data[2] = rt_atan2d_snf(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp * c_a_tmp,
      b_a_tmp * d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
    vdis[2] = 1.0E+6;

    /* set to something large */
    vtrueB_data[3] = 7.0685834705770345;

    /* evenly distributed 0.75 - 1.25 - 1.75 - 2.25 */
    d_a_tmp = cos(argpB + 7.0685834705770345);
    e_a_tmp = sin(argpB + 7.0685834705770345);
    vL_data[3] = rt_atan2d_snf(a_tmp * d_a_tmp + b_a_tmp * e_a_tmp * c_a_tmp,
      b_a_tmp * d_a_tmp - a_tmp * e_a_tmp * c_a_tmp);
    vdis[3] = 1.0E+6;

    /* set to something large */
  }
}

/* End of code generation (WaterProcedure.c) */
