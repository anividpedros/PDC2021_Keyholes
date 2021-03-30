/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

/* Include files */
#include "main.h"
#include "ComputeMOID.h"
#include "ComputeMOID_terminate.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static double argInit_real_T(void);
static struct0_T argInit_struct0_T(void);
static void main_ComputeMOID(void);

/* Function Definitions */
static double argInit_real_T(void)
{
  return 0.0;
}

static struct0_T argInit_struct0_T(void)
{
  struct0_T result;
  double result_tmp;

  /* Set the value of each structure field.
     Change this value to the value that the application requires. */
  result_tmp = argInit_real_T();
  result.e = result_tmp;
  result.i = result_tmp;
  result.Omega = result_tmp;
  result.argp = result_tmp;
  result.sma = result_tmp;
  return result;
}

static void main_ComputeMOID(void)
{
  struct0_T A_tmp;
  double moid;
  double vA_out;
  double vB_out;

  /* Initialize function 'ComputeMOID' input arguments. */
  /* Initialize function input argument 'A'. */
  A_tmp = argInit_struct0_T();

  /* Initialize function input argument 'B'. */
  /* Call the entry-point 'ComputeMOID'. */
  ComputeMOID(&A_tmp, &A_tmp, &moid, &vA_out, &vB_out);
}

int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_ComputeMOID();

  /* Terminate the application.
     You do not need to do this more than one time. */
  ComputeMOID_terminate();
  return 0;
}

/* End of code generation (main.c) */
