/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ComputeMOID_initialize.c
 *
 * Code generation for function 'ComputeMOID_initialize'
 *
 */

/* Include files */
#include "ComputeMOID_initialize.h"
#include "ComputeMOID.h"
#include "ComputeMOID_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void ComputeMOID_initialize(void)
{
  rt_InitInfAndNaN();
  isInitialized_ComputeMOID = true;
}

/* End of code generation (ComputeMOID_initialize.c) */
