/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ComputeMOID.c
 *
 * Code generation for function 'ComputeMOID'
 *
 */

/* Include files */
#include "ComputeMOID.h"
#include "ComputeMOID_data.h"
#include "RefFrame.h"
#include "WaterProcedure.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 4,     /* lineNo */
  "ComputeMOID",                       /* fcnName */
  "C:\\Users\\Oscar\\Documents\\GitHub\\PDC2021_Keyholes\\MOIDroutines\\MOID_MATLAB\\ComputeMOID.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 36,  /* lineNo */
  "ComputeMOID",                       /* fcnName */
  "C:\\Users\\Oscar\\Documents\\GitHub\\PDC2021_Keyholes\\MOIDroutines\\MOID_MATLAB\\ComputeMOID.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 74,  /* lineNo */
  "ComputeMOID",                       /* fcnName */
  "C:\\Users\\Oscar\\Documents\\GitHub\\PDC2021_Keyholes\\MOIDroutines\\MOID_MATLAB\\ComputeMOID.m"/* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 305, /* lineNo */
  "ComputeMOID",                       /* fcnName */
  "C:\\Users\\Oscar\\Documents\\GitHub\\PDC2021_Keyholes\\MOIDroutines\\MOID_MATLAB\\ComputeMOID.m"/* pathName */
};

static emlrtBCInfo emlrtBCI = { 1,     /* iFirst */
  4,                                   /* iLast */
  94,                                  /* lineNo */
  18,                                  /* colNo */
  "vdis",                              /* aName */
  "ComputeMOID",                       /* fName */
  "C:\\Users\\Oscar\\Documents\\GitHub\\PDC2021_Keyholes\\MOIDroutines\\MOID_MATLAB\\ComputeMOID.m",/* pName */
  3                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { 1,   /* iFirst */
  3,                                   /* iLast */
  249,                                 /* lineNo */
  26,                                  /* colNo */
  "xBt",                               /* aName */
  "ComputeMOID",                       /* fName */
  "C:\\Users\\Oscar\\Documents\\GitHub\\PDC2021_Keyholes\\MOIDroutines\\MOID_MATLAB\\ComputeMOID.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { 1,   /* iFirst */
  3,                                   /* iLast */
  249,                                 /* lineNo */
  34,                                  /* colNo */
  "xAt",                               /* aName */
  "ComputeMOID",                       /* fName */
  "C:\\Users\\Oscar\\Documents\\GitHub\\PDC2021_Keyholes\\MOIDroutines\\MOID_MATLAB\\ComputeMOID.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRSInfo u_emlrtRSI = { 132, /* lineNo */
  "ComputeMOID",                       /* fcnName */
  "C:\\Users\\Oscar\\Documents\\GitHub\\PDC2021_Keyholes\\MOIDroutines\\MOID_MATLAB\\ComputeMOID.m"/* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 107, /* lineNo */
  "ComputeMOID",                       /* fcnName */
  "C:\\Users\\Oscar\\Documents\\GitHub\\PDC2021_Keyholes\\MOIDroutines\\MOID_MATLAB\\ComputeMOID.m"/* pathName */
};

/* Function Definitions */
void ComputeMOID(const emlrtStack *sp, const struct0_T *A, const struct0_T *B,
                 real_T *moid, real_T *vA_out, real_T *vB_out)
{
  real_T incliB;
  real_T omegaB;
  real_T argpB;
  real_T vdis[4];
  real_T L_m;
  real_T trueB_m;
  real_T xB_tmp_tmp;
  real_T b_xB_tmp_tmp;
  real_T c_xB_tmp_tmp;
  real_T rB_tmp_tmp;
  real_T rB;
  real_T yB;
  real_T rAt_idx_1;
  real_T xB;
  real_T a_tmp;
  real_T k;
  real_T rhoB;
  real_T rA_tmp_tmp;
  real_T dist_oo;
  real_T dist_o;
  real_T trueB_o;
  real_T L_o;
  real_T trueB;
  real_T N;
  real_T vtrueB[100];
  real_T vL[100];
  real_T vtrueB_data[100];
  int32_T vtrueB_size[1];
  real_T vL_data[100];
  int32_T vL_size[1];
  real_T xBt[3];
  real_T yBt[3];
  real_T zBt[3];
  real_T xAt[3];
  real_T yAt[3];
  int32_T i;
  int32_T b_i;
  real_T b_vtrueB[4];
  real_T b_vL[4];
  boolean_T aleft;
  boolean_T aright;
  boolean_T bleft;
  boolean_T bright;
  int32_T lpoints;
  int32_T k1max;
  int32_T k1min;
  int32_T i1min;
  int32_T i1max;
  boolean_T calc1;
  boolean_T calc2;
  boolean_T calc3;
  boolean_T calc4;
  int32_T k1_t;
  int32_T i1_t;
  int32_T k1;
  int32_T i1;
  int32_T b_i1;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;

  /*     %% New orbital elements: A frame */
  /*  incliA = 0, omegaA = 0, argpA = 0 */
  st.site = &emlrtRSI;
  RefFrame(&st, A->i, A->Omega, A->argp, B->i, B->Omega, B->argp, &incliB,
           &omegaB, &argpB);

  /*     %% Scanning orbits */
  /*  Scanning one full revolution of meridional plane to find local minima */
  /*  First guess for MOID is set to be big ~1e6AU */
  /*  Angle sweeping steps */
  /* rad - based on Wisnioski paper ~ 0.12rad */
  /* rad - for initial tuning step */
  /* rad - for final step of first tuning */
  /* rad - threshold step for secondtuning */
  /*  Initial Guess */
  *moid = 1.0E+6;
  vdis[0] = 1.0E+6;
  vdis[1] = 1.0E+6;
  vdis[2] = 1.0E+6;
  vdis[3] = 1.0E+6;

  /*  Pre-define variables to know type in code generation */
  L_m = 0.0;
  trueB_m = 0.0;
  xB_tmp_tmp = muDoubleScalarSin(omegaB);
  b_xB_tmp_tmp = muDoubleScalarCos(omegaB);
  c_xB_tmp_tmp = muDoubleScalarCos(incliB);

  /*  First triplet */
  rB_tmp_tmp = B->sma * (1.0 - B->e * B->e);
  rB = rB_tmp_tmp / (B->e * 0.97133797485202966 + 1.0);

  /* compute the radius for B */
  yB = muDoubleScalarCos(argpB + -0.24);
  rAt_idx_1 = muDoubleScalarSin(argpB + -0.24);
  xB = rB * (b_xB_tmp_tmp * yB - xB_tmp_tmp * rAt_idx_1 * c_xB_tmp_tmp);
  yB = rB * (xB_tmp_tmp * yB + b_xB_tmp_tmp * rAt_idx_1 * c_xB_tmp_tmp);
  a_tmp = muDoubleScalarSin(incliB);
  k = rB * rAt_idx_1 * a_tmp;
  rhoB = xB * xB + yB * yB;
  st.site = &c_emlrtRSI;
  if (rhoB < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  rhoB = muDoubleScalarSqrt(rhoB);
  rA_tmp_tmp = A->sma * (1.0 - A->e * A->e);
  rAt_idx_1 = A->e * muDoubleScalarCos(muDoubleScalarAtan2(yB, xB));
  xB = rA_tmp_tmp / (rAt_idx_1 + 1.0);

  /* compute the radius for A */
  rB = rhoB - xB;
  if (muDoubleScalarAbs(rB) > muDoubleScalarAbs(rhoB + xB)) {
    rB = rhoB + rA_tmp_tmp / (1.0 - rAt_idx_1);
  }

  /*  square of the distance */
  /*  storing */
  dist_oo = k * k + rB * rB;
  if (*emlrtBreakCheckR2012bFlagVar != 0) {
    emlrtBreakCheckR2012b(sp);
  }

  /*  First triplet */
  rB = rB_tmp_tmp / (B->e * 0.99280863585386625 + 1.0);

  /* compute the radius for B */
  yB = muDoubleScalarCos(argpB + -0.12);
  rAt_idx_1 = muDoubleScalarSin(argpB + -0.12);
  xB = rB * (b_xB_tmp_tmp * yB - xB_tmp_tmp * rAt_idx_1 * c_xB_tmp_tmp);
  yB = rB * (xB_tmp_tmp * yB + b_xB_tmp_tmp * rAt_idx_1 * c_xB_tmp_tmp);
  k = rB * rAt_idx_1 * a_tmp;
  rhoB = xB * xB + yB * yB;
  st.site = &c_emlrtRSI;
  if (rhoB < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  rhoB = muDoubleScalarSqrt(rhoB);
  yB = muDoubleScalarAtan2(yB, xB);
  rAt_idx_1 = A->e * muDoubleScalarCos(yB);
  xB = rA_tmp_tmp / (rAt_idx_1 + 1.0);

  /* compute the radius for A */
  rB = rhoB - xB;
  if (muDoubleScalarAbs(rB) > muDoubleScalarAbs(rhoB + xB)) {
    yB -= 3.1415926535897931;
    rB = rhoB + rA_tmp_tmp / (1.0 - rAt_idx_1);
  }

  /*  square of the distance */
  /*  storing */
  dist_o = k * k + rB * rB;
  trueB_o = -0.12;
  L_o = yB;
  trueB = 0.0;
  if (*emlrtBreakCheckR2012bFlagVar != 0) {
    emlrtBreakCheckR2012b(sp);
  }

  /*  Full revolution */
  N = 0.0;

  /*  number of minima */
  memset(&vtrueB[0], 0, 100U * sizeof(real_T));
  memset(&vL[0], 0, 100U * sizeof(real_T));
  while (trueB < 6.4031853071795863) {
    rB = rB_tmp_tmp / (B->e * muDoubleScalarCos(trueB) + 1.0);

    /* compute the radius for B */
    rAt_idx_1 = argpB + trueB;
    yB = muDoubleScalarCos(rAt_idx_1);
    rAt_idx_1 = muDoubleScalarSin(rAt_idx_1);
    xB = rB * (b_xB_tmp_tmp * yB - xB_tmp_tmp * rAt_idx_1 * c_xB_tmp_tmp);
    yB = rB * (xB_tmp_tmp * yB + b_xB_tmp_tmp * rAt_idx_1 * c_xB_tmp_tmp);
    k = rB * rAt_idx_1 * a_tmp;
    rhoB = xB * xB + yB * yB;
    st.site = &h_emlrtRSI;
    if (rhoB < 0.0) {
      emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    rhoB = muDoubleScalarSqrt(rhoB);
    yB = muDoubleScalarAtan2(yB, xB);
    rAt_idx_1 = A->e * muDoubleScalarCos(yB);
    xB = rA_tmp_tmp / (rAt_idx_1 + 1.0);

    /* compute the radius for A */
    rB = rhoB - xB;
    if (muDoubleScalarAbs(rB) > muDoubleScalarAbs(rhoB + xB)) {
      yB -= 3.1415926535897931;
      rB = rhoB + rA_tmp_tmp / (1.0 - rAt_idx_1);
    }

    rAt_idx_1 = k * k + rB * rB;

    /*  square of the distance */
    if ((dist_o <= rAt_idx_1) && (dist_o <= dist_oo)) {
      N++;
      i = (int32_T)N;
      k1max = i - 1;
      vtrueB[k1max] = trueB_o;
      vL[k1max] = L_o;
      if ((i < 1) || (i > 4)) {
        emlrtDynamicBoundsCheckR2012b(i, 1, 4, &emlrtBCI, sp);
      }

      vdis[i - 1] = dist_o;
    }

    dist_oo = dist_o;
    trueB_o = trueB;
    L_o = yB;
    dist_o = rAt_idx_1;
    trueB += 0.12;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  /*     %% Water Procedure */
  st.site = &v_emlrtRSI;
  b_WaterProcedure(&st, incliB, omegaB, argpB, &N, vtrueB, vL, vdis, vtrueB_data,
                   vtrueB_size, vL_data, vL_size);

  /*     %% PARALLEL TUNING */
  /*  Move objects separately along their orbits */
  /*  Smallest possible distance no longer meridional distance */
  xBt[0] = rtNaN;
  yBt[0] = rtNaN;
  zBt[0] = rtNaN;
  xAt[0] = rtNaN;
  yAt[0] = rtNaN;
  xBt[2] = rtNaN;
  yBt[2] = rtNaN;
  zBt[2] = rtNaN;
  xAt[2] = rtNaN;
  yAt[2] = rtNaN;
  k = 1.0;
  while (k < N + 2.0) {
    if (k <= N) {
      i = (int32_T)k - 1;
      *moid = vdis[i];
      trueB_m = vtrueB_data[i];
      L_m = vL_data[i];
      rhoB = 0.07;
      dist_oo = 1.0E-5;

      /* maybe here problem */
    } else {
      if (N == 2.0) {
        /*  in case of two minima are very close to each other(<1E-4 a.u.) */
        /*  go to "water procedure" */
        if (muDoubleScalarAbs(vdis[0] - vdis[1]) < 0.0001) {
          st.site = &u_emlrtRSI;
          WaterProcedure(&st, incliB, omegaB, argpB, b_vtrueB, b_vL, vdis, &N);
          vtrueB_data[0] = b_vtrueB[0];
          vL_data[0] = b_vL[0];
          vtrueB_data[1] = b_vtrueB[1];
          vL_data[1] = b_vL[1];
          vtrueB_data[2] = b_vtrueB[2];
          vL_data[2] = b_vL[2];
          vtrueB_data[3] = b_vtrueB[3];
          vL_data[3] = b_vL[3];
          N = 4.0;
          k = 1.0;
        } else {
          if (vdis[0] < *moid) {
            *moid = vdis[0];
            trueB_m = vtrueB_data[0];
            L_m = vL_data[0];
          }
        }
      } else {
        /*  final tuning */
        b_i = (int32_T)(N - 1.0);
        for (i = 0; i < b_i; i++) {
          if (vdis[i] < *moid) {
            *moid = vdis[i];
            trueB_m = vtrueB_data[i];
            L_m = vL_data[i];
          }

          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b(sp);
          }
        }
      }

      rhoB = 0.14;

      /* inital state */
      dist_oo = 1.0E-14;

      /* terminal state */
    }

    yB = rB_tmp_tmp / (B->e * muDoubleScalarCos(trueB_m) + 1.0);

    /* compute the radius for B */
    rAt_idx_1 = argpB + trueB_m;
    xB = muDoubleScalarCos(rAt_idx_1);
    rB = muDoubleScalarSin(rAt_idx_1);
    xBt[1] = yB * (b_xB_tmp_tmp * xB - xB_tmp_tmp * rB * c_xB_tmp_tmp);
    yBt[1] = yB * (xB_tmp_tmp * xB + b_xB_tmp_tmp * rB * c_xB_tmp_tmp);
    zBt[1] = yB * rB * a_tmp;
    rB = muDoubleScalarCos(L_m);
    rAt_idx_1 = rA_tmp_tmp / (A->e * rB + 1.0);

    /* compute the radius for A */
    xAt[1] = rAt_idx_1 * rB;
    yAt[1] = rAt_idx_1 * muDoubleScalarSin(L_m);
    aleft = true;
    aright = true;
    bleft = true;
    bright = true;
    while (rhoB >= dist_oo) {
      lpoints = 0;
      k1min = 1;
      k1max = 3;
      i1min = 1;
      i1max = 3;
      calc1 = false;
      calc2 = false;
      calc3 = false;
      calc4 = false;
      if (bleft) {
        rAt_idx_1 = rB_tmp_tmp / (B->e * muDoubleScalarCos(trueB_m - rhoB) + 1.0);

        /* compute the radius for B */
        rB = (argpB + trueB_m) - rhoB;
        xB = muDoubleScalarCos(rB);
        rB = muDoubleScalarSin(rB);
        xBt[0] = rAt_idx_1 * (b_xB_tmp_tmp * xB - xB_tmp_tmp * rB * c_xB_tmp_tmp);
        yBt[0] = rAt_idx_1 * (xB_tmp_tmp * xB + b_xB_tmp_tmp * rB * c_xB_tmp_tmp);
        zBt[0] = rAt_idx_1 * rB * a_tmp;
        lpoints = 1;
      }

      if (bright) {
        rAt_idx_1 = rB_tmp_tmp / (B->e * muDoubleScalarCos(trueB_m + rhoB) + 1.0);

        /* compute the radius for B */
        rB = (argpB + trueB_m) + rhoB;
        xB = muDoubleScalarCos(rB);
        rB = muDoubleScalarSin(rB);
        xBt[2] = rAt_idx_1 * (b_xB_tmp_tmp * xB - xB_tmp_tmp * rB * c_xB_tmp_tmp);
        yBt[2] = rAt_idx_1 * (xB_tmp_tmp * xB + b_xB_tmp_tmp * rB * c_xB_tmp_tmp);
        zBt[2] = rAt_idx_1 * rB * a_tmp;
        lpoints++;
      }

      if (aleft) {
        yB = L_m - rhoB;
        rB = muDoubleScalarCos(yB);
        rAt_idx_1 = rA_tmp_tmp / (A->e * rB + 1.0);

        /* compute the radius for A */
        xAt[0] = rAt_idx_1 * rB;
        yAt[0] = rAt_idx_1 * muDoubleScalarSin(yB);
        lpoints++;
      }

      if (aright) {
        yB = L_m + rhoB;
        rB = muDoubleScalarCos(yB);
        rAt_idx_1 = rA_tmp_tmp / (A->e * rB + 1.0);

        /* compute the radius for A */
        xAt[2] = rAt_idx_1 * rB;
        yAt[2] = rAt_idx_1 * muDoubleScalarSin(yB);
        lpoints++;
      }

      k1_t = 2;
      i1_t = 2;
      if (lpoints == 1) {
        if (aleft) {
          i1max = 1;
        }

        if (aright) {
          i1min = 3;
        }

        if (bright) {
          k1min = 3;
        }

        if (bleft) {
          k1max = 1;
        }
      }

      if (lpoints == 2) {
        if (aleft && bright) {
          calc1 = true;
        }

        if (aleft && bleft) {
          calc2 = true;
        }

        if (aright && bright) {
          calc3 = true;
        }

        if (aright && bleft) {
          calc4 = true;
        }
      }

      b_i = k1max - k1min;
      for (k1 = 0; k1 <= b_i; k1++) {
        i = k1min + k1;
        k1max = i1max - i1min;
        for (i1 = 0; i1 <= k1max; i1++) {
          b_i1 = i1min + i1;
          aleft = true;
          if (lpoints == 2) {
            if ((b_i1 != 1) && (((i != 3) && calc1) || ((i != 1) && calc2))) {
              aleft = false;
            }

            if ((b_i1 != 3) && (((i != 3) && calc3) || ((i != 1) && calc4))) {
              aleft = false;
            }
          }

          if ((i == 2) && (b_i1 == 2)) {
            aleft = false;
          }

          if (aleft) {
            if (i > 3) {
              emlrtDynamicBoundsCheckR2012b(i, 1, 3, &b_emlrtBCI, sp);
            }

            if (b_i1 > 3) {
              emlrtDynamicBoundsCheckR2012b(b_i1, 1, 3, &c_emlrtBCI, sp);
            }

            rAt_idx_1 = xBt[i - 1] - xAt[b_i1 - 1];
            yB = yBt[i - 1] - yAt[b_i1 - 1];
            rB = zBt[i - 1];
            rAt_idx_1 = (rAt_idx_1 * rAt_idx_1 + yB * yB) + rB * rB;
            if (rAt_idx_1 < *moid) {
              *moid = rAt_idx_1;
              k1_t = i;
              i1_t = b_i1;
            }
          }

          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b(sp);
          }
        }

        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(sp);
        }
      }

      if ((k1_t != 2) || (i1_t != 2)) {
        aleft = false;
        aright = false;
        bleft = false;
        bright = false;
        if (i1_t != 2) {
          if (i1_t == 1) {
            aleft = true;
            L_m -= rhoB;
            xAt[2] = xAt[1];
            yAt[2] = yAt[1];
            xAt[1] = xAt[0];
            yAt[1] = yAt[0];
          } else {
            aright = true;
            L_m += rhoB;
            xAt[0] = xAt[1];
            yAt[0] = yAt[1];
            xAt[1] = xAt[2];
            yAt[1] = yAt[2];
          }
        }

        if (k1_t != 2) {
          if (k1_t == 1) {
            bleft = true;
            trueB_m -= rhoB;
            xBt[2] = xBt[1];
            yBt[2] = yBt[1];
            zBt[2] = zBt[1];
            xBt[1] = xBt[0];
            yBt[1] = yBt[0];
            zBt[1] = zBt[0];
          } else {
            bright = true;
            trueB_m += rhoB;
            xBt[0] = xBt[1];
            yBt[0] = yBt[1];
            zBt[0] = zBt[1];
            xBt[1] = xBt[2];
            yBt[1] = yBt[2];
            zBt[1] = zBt[2];
          }
        }
      } else {
        aleft = true;
        aright = true;
        bleft = true;
        bright = true;
        rhoB *= 0.15;
      }

      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }

    if (k <= N) {
      i = (int32_T)k - 1;
      vdis[i] = *moid;
      vtrueB_data[i] = trueB_m;
      vL_data[i] = L_m;
    }

    k++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  st.site = &r_emlrtRSI;
  if (*moid < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  *moid = muDoubleScalarSqrt(*moid);
  *vA_out = -L_m;
  *vB_out = -trueB_m;
}

/* End of code generation (ComputeMOID.c) */
