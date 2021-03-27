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

#ifndef WATERPROCEDURE_H
#define WATERPROCEDURE_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "ComputeMOID_types.h"

/* Function Declarations */
extern void WaterProcedure(double incliB, double omegaB, double argpB, double
  vtrueB[4], double vL[4], double vdis[4], double *N);
extern void b_WaterProcedure(double incliB, double omegaB, double argpB, double *
  N, const double vtrueB[100], const double vL[100], double vdis[4], double
  vtrueB_data[], int vtrueB_size[1], double vL_data[], int vL_size[1]);

#endif

/* End of code generation (WaterProcedure.h) */
