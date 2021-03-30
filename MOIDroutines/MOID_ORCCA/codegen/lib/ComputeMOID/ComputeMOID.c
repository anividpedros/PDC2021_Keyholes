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
#include "ComputeMOID_initialize.h"
#include "ComputeMOID_rtwutil.h"
#include "RefFrame.h"
#include "WaterProcedure.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
void ComputeMOID(const struct0_T *A, const struct0_T *B, double *moid, double
                 *vA_out, double *vB_out)
{
  double incliB;
  double omegaB;
  double argpB;
  double vdis[4];
  double L_m;
  double trueB_m;
  double xB_tmp_tmp;
  double b_xB_tmp_tmp;
  double c_xB_tmp_tmp;
  double rB_tmp_tmp;
  double rB;
  double yB;
  double rAt_idx_1;
  double xB;
  double a_tmp;
  double k;
  double rhoB;
  double rA_tmp_tmp;
  double dist_oo;
  double dist_o;
  double trueB_o;
  double L_o;
  double trueB;
  double N;
  double vtrueB[100];
  double vL[100];
  double vtrueB_data[100];
  int vtrueB_size[1];
  double vL_data[100];
  int vL_size[1];
  double xBt[3];
  double yBt[3];
  double zBt[3];
  double xAt[3];
  double yAt[3];
  int k1max;
  int i;
  double b_vtrueB[4];
  double b_vL[4];
  boolean_T aleft;
  boolean_T aright;
  boolean_T bleft;
  boolean_T bright;
  int lpoints;
  int k1min;
  int i1min;
  int i1max;
  boolean_T calc1;
  boolean_T calc2;
  boolean_T calc3;
  boolean_T calc4;
  int k1_t;
  int i1_t;
  int k1;
  int i1;
  int b_i1;
  int c_i1;
  if (!isInitialized_ComputeMOID) {
    ComputeMOID_initialize();
  }

  /*     %% New orbital elements: A frame */
  /*  incliA = 0, omegaA = 0, argpA = 0 */
  RefFrame(A->i, A->Omega, A->argp, B->i, B->Omega, B->argp, &incliB, &omegaB,
           &argpB);

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
  xB_tmp_tmp = sin(omegaB);
  b_xB_tmp_tmp = cos(omegaB);
  c_xB_tmp_tmp = cos(incliB);

  /*  First triplet */
  rB_tmp_tmp = B->sma * (1.0 - B->e * B->e);
  rB = rB_tmp_tmp / (B->e * 0.97133797485202966 + 1.0);

  /* compute the radius for B */
  yB = cos(argpB + -0.24);
  rAt_idx_1 = sin(argpB + -0.24);
  xB = rB * (b_xB_tmp_tmp * yB - xB_tmp_tmp * rAt_idx_1 * c_xB_tmp_tmp);
  yB = rB * (xB_tmp_tmp * yB + b_xB_tmp_tmp * rAt_idx_1 * c_xB_tmp_tmp);
  a_tmp = sin(incliB);
  k = rB * rAt_idx_1 * a_tmp;
  rhoB = sqrt(xB * xB + yB * yB);
  rA_tmp_tmp = A->sma * (1.0 - A->e * A->e);
  rAt_idx_1 = A->e * cos(rt_atan2d_snf(yB, xB));
  xB = rA_tmp_tmp / (rAt_idx_1 + 1.0);

  /* compute the radius for A */
  rB = rhoB - xB;
  if (fabs(rB) > fabs(rhoB + xB)) {
    rB = rhoB + rA_tmp_tmp / (1.0 - rAt_idx_1);
  }

  /*  square of the distance */
  /*  storing */
  dist_oo = k * k + rB * rB;

  /*  First triplet */
  rB = rB_tmp_tmp / (B->e * 0.99280863585386625 + 1.0);

  /* compute the radius for B */
  yB = cos(argpB + -0.12);
  rAt_idx_1 = sin(argpB + -0.12);
  xB = rB * (b_xB_tmp_tmp * yB - xB_tmp_tmp * rAt_idx_1 * c_xB_tmp_tmp);
  yB = rB * (xB_tmp_tmp * yB + b_xB_tmp_tmp * rAt_idx_1 * c_xB_tmp_tmp);
  k = rB * rAt_idx_1 * a_tmp;
  rhoB = sqrt(xB * xB + yB * yB);
  yB = rt_atan2d_snf(yB, xB);
  rAt_idx_1 = A->e * cos(yB);
  xB = rA_tmp_tmp / (rAt_idx_1 + 1.0);

  /* compute the radius for A */
  rB = rhoB - xB;
  if (fabs(rB) > fabs(rhoB + xB)) {
    yB -= 3.1415926535897931;
    rB = rhoB + rA_tmp_tmp / (1.0 - rAt_idx_1);
  }

  /*  square of the distance */
  /*  storing */
  dist_o = k * k + rB * rB;
  trueB_o = -0.12;
  L_o = yB;
  trueB = 0.0;

  /*  Full revolution */
  N = 0.0;

  /*  number of minima */
  memset(&vtrueB[0], 0, 100U * sizeof(double));
  memset(&vL[0], 0, 100U * sizeof(double));
  while (trueB < 6.4031853071795863) {
    rB = rB_tmp_tmp / (B->e * cos(trueB) + 1.0);

    /* compute the radius for B */
    rAt_idx_1 = argpB + trueB;
    yB = cos(rAt_idx_1);
    rAt_idx_1 = sin(rAt_idx_1);
    xB = rB * (b_xB_tmp_tmp * yB - xB_tmp_tmp * rAt_idx_1 * c_xB_tmp_tmp);
    yB = rB * (xB_tmp_tmp * yB + b_xB_tmp_tmp * rAt_idx_1 * c_xB_tmp_tmp);
    k = rB * rAt_idx_1 * a_tmp;
    rhoB = sqrt(xB * xB + yB * yB);
    yB = rt_atan2d_snf(yB, xB);
    rAt_idx_1 = A->e * cos(yB);
    xB = rA_tmp_tmp / (rAt_idx_1 + 1.0);

    /* compute the radius for A */
    rB = rhoB - xB;
    if (fabs(rB) > fabs(rhoB + xB)) {
      yB -= 3.1415926535897931;
      rB = rhoB + rA_tmp_tmp / (1.0 - rAt_idx_1);
    }

    rAt_idx_1 = k * k + rB * rB;

    /*  square of the distance */
    if ((dist_o <= rAt_idx_1) && (dist_o <= dist_oo)) {
      N++;
      k1max = (int)N - 1;
      vtrueB[k1max] = trueB_o;
      vL[k1max] = L_o;
      vdis[k1max] = dist_o;
    }

    dist_oo = dist_o;
    trueB_o = trueB;
    L_o = yB;
    dist_o = rAt_idx_1;
    trueB += 0.12;
  }

  /*     %% Water Procedure */
  b_WaterProcedure(incliB, omegaB, argpB, &N, vtrueB, vL, vdis, vtrueB_data,
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
      k1max = (int)k - 1;
      *moid = vdis[k1max];
      trueB_m = vtrueB_data[k1max];
      L_m = vL_data[k1max];
      rhoB = 0.07;
      dist_oo = 1.0E-5;

      /* maybe here problem */
    } else {
      if (N == 2.0) {
        /*  in case of two minima are very close to each other(<1E-4 a.u.) */
        /*  go to "water procedure" */
        if (fabs(vdis[0] - vdis[1]) < 0.0001) {
          WaterProcedure(incliB, omegaB, argpB, b_vtrueB, b_vL, vdis, &N);
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
        i = (int)(N - 1.0);
        for (k1max = 0; k1max < i; k1max++) {
          if (vdis[k1max] < *moid) {
            *moid = vdis[k1max];
            trueB_m = vtrueB_data[k1max];
            L_m = vL_data[k1max];
          }
        }
      }

      rhoB = 0.14;

      /* inital state */
      dist_oo = 1.0E-14;

      /* terminal state */
    }

    yB = rB_tmp_tmp / (B->e * cos(trueB_m) + 1.0);

    /* compute the radius for B */
    rAt_idx_1 = argpB + trueB_m;
    xB = cos(rAt_idx_1);
    rB = sin(rAt_idx_1);
    xBt[1] = yB * (b_xB_tmp_tmp * xB - xB_tmp_tmp * rB * c_xB_tmp_tmp);
    yBt[1] = yB * (xB_tmp_tmp * xB + b_xB_tmp_tmp * rB * c_xB_tmp_tmp);
    zBt[1] = yB * rB * a_tmp;
    rB = cos(L_m);
    rAt_idx_1 = rA_tmp_tmp / (A->e * rB + 1.0);

    /* compute the radius for A */
    xAt[1] = rAt_idx_1 * rB;
    yAt[1] = rAt_idx_1 * sin(L_m);
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
        rAt_idx_1 = rB_tmp_tmp / (B->e * cos(trueB_m - rhoB) + 1.0);

        /* compute the radius for B */
        rB = (argpB + trueB_m) - rhoB;
        xB = cos(rB);
        rB = sin(rB);
        xBt[0] = rAt_idx_1 * (b_xB_tmp_tmp * xB - xB_tmp_tmp * rB * c_xB_tmp_tmp);
        yBt[0] = rAt_idx_1 * (xB_tmp_tmp * xB + b_xB_tmp_tmp * rB * c_xB_tmp_tmp);
        zBt[0] = rAt_idx_1 * rB * a_tmp;
        lpoints = 1;
      }

      if (bright) {
        rAt_idx_1 = rB_tmp_tmp / (B->e * cos(trueB_m + rhoB) + 1.0);

        /* compute the radius for B */
        rB = (argpB + trueB_m) + rhoB;
        xB = cos(rB);
        rB = sin(rB);
        xBt[2] = rAt_idx_1 * (b_xB_tmp_tmp * xB - xB_tmp_tmp * rB * c_xB_tmp_tmp);
        yBt[2] = rAt_idx_1 * (xB_tmp_tmp * xB + b_xB_tmp_tmp * rB * c_xB_tmp_tmp);
        zBt[2] = rAt_idx_1 * rB * a_tmp;
        lpoints++;
      }

      if (aleft) {
        yB = L_m - rhoB;
        rB = cos(yB);
        rAt_idx_1 = rA_tmp_tmp / (A->e * rB + 1.0);

        /* compute the radius for A */
        xAt[0] = rAt_idx_1 * rB;
        yAt[0] = rAt_idx_1 * sin(yB);
        lpoints++;
      }

      if (aright) {
        yB = L_m + rhoB;
        rB = cos(yB);
        rAt_idx_1 = rA_tmp_tmp / (A->e * rB + 1.0);

        /* compute the radius for A */
        xAt[2] = rAt_idx_1 * rB;
        yAt[2] = rAt_idx_1 * sin(yB);
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

      i = k1max - k1min;
      for (k1max = 0; k1max <= i; k1max++) {
        k1 = k1min + k1max;
        i1 = i1max - i1min;
        for (b_i1 = 0; b_i1 <= i1; b_i1++) {
          c_i1 = i1min + b_i1;
          aleft = true;
          if (lpoints == 2) {
            if ((c_i1 != 1) && (((k1 != 3) && calc1) || ((k1 != 1) && calc2))) {
              aleft = false;
            }

            if ((c_i1 != 3) && (((k1 != 3) && calc3) || ((k1 != 1) && calc4))) {
              aleft = false;
            }
          }

          if ((k1 == 2) && (c_i1 == 2)) {
            aleft = false;
          }

          if (aleft) {
            rAt_idx_1 = xBt[k1 - 1] - xAt[c_i1 - 1];
            yB = yBt[k1 - 1] - yAt[c_i1 - 1];
            rB = zBt[k1 - 1];
            rAt_idx_1 = (rAt_idx_1 * rAt_idx_1 + yB * yB) + rB * rB;
            if (rAt_idx_1 < *moid) {
              *moid = rAt_idx_1;
              k1_t = k1;
              i1_t = c_i1;
            }
          }
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
    }

    if (k <= N) {
      k1max = (int)k - 1;
      vdis[k1max] = *moid;
      vtrueB_data[k1max] = trueB_m;
      vL_data[k1max] = L_m;
    }

    k++;
  }

  *moid = sqrt(*moid);
  *vA_out = -L_m;
  *vB_out = -trueB_m;
}

/* End of code generation (ComputeMOID.c) */
