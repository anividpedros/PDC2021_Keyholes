/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_ComputeMOID_mex.c
 *
 * Code generation for function '_coder_ComputeMOID_mex'
 *
 */

/* Include files */
#include "_coder_ComputeMOID_mex.h"
#include "ComputeMOID.h"
#include "ComputeMOID_data.h"
#include "ComputeMOID_initialize.h"
#include "ComputeMOID_terminate.h"
#include "_coder_ComputeMOID_api.h"

/* Function Declarations */
MEXFUNCTION_LINKAGE void ComputeMOID_mexFunction(int32_T nlhs, mxArray *plhs[3],
  int32_T nrhs, const mxArray *prhs[2]);

/* Function Definitions */
void ComputeMOID_mexFunction(int32_T nlhs, mxArray *plhs[3], int32_T nrhs, const
  mxArray *prhs[2])
{
  const mxArray *outputs[3];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 2, 4,
                        11, "ComputeMOID");
  }

  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 11,
                        "ComputeMOID");
  }

  /* Call the function. */
  ComputeMOID_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(&ComputeMOID_atexit);

  /* Module initialization. */
  ComputeMOID_initialize();

  /* Dispatch the entry-point. */
  ComputeMOID_mexFunction(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  ComputeMOID_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_ComputeMOID_mex.c) */
