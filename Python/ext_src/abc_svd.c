/*
 * svdcomp - SVD decomposition routine.
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 *
 * Input to dsvd is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   ROW = row dimension of a
 *   COL = column dimension of a
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
*/

//Regula experssion to replace array referencing:
//v\[([\w]+)\]\[([\w]+)\]
//v\[([\w]+)\]\[([\w]+)\]
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "abc_000_warning.h"

#include "abc_datatype.h"
#include "abc_ide_util.h"
 

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define MAX(x,y) ((x)>(y)?(x):(y))

static double PYTHAG(double a, double b) {
    double at = fabs(a), bt = fabs(b), ct, result;
    if (at > bt) { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}


int dsvd(double *a, int ROW, int COL, double* w, double *v)
{
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double* rv1;

    if (ROW < COL)
    {
        r_printf("#rows must be > #cols \n");
        return(0);
    }

    rv1 = (double*)malloc((unsigned int)COL * sizeof(double));

    /* Householder reduction to bidiagonal form */
    for (i = 0; i < COL; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < ROW)
        {
            for (k = i; k < ROW; k++)
                scale += fabs((double)a[i*ROW+k]);
            if (scale)
            {
                for (k = i; k < ROW; k++)
                {
                    a[i*ROW+k] = (double)((double)a[i*ROW+k] / scale);
                    s += ((double)a[i*ROW+k] * (double)a[i*ROW+k]);
                }
                f = (double)a[i*ROW+i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i*ROW+i] = (double)(f - g);
                if (i != COL - 1)
                {
                    for (j = l; j < COL; j++)
                    {
                        for (s = 0.0, k = i; k < ROW; k++)
                            s += ((double)a[i*ROW+k] * (double)a[j*ROW+k]);
                        f = s / h;
                        for (k = i; k < ROW; k++)
                            a[j*ROW+k] += (double)(f * (double)a[i*ROW+k]);
                    }
                }
                for (k = i; k < ROW; k++)
                    a[i*ROW+k] = (double)((double)a[i*ROW+k] * scale);
            }
        }
        w[i] = (double)(scale * g);

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < ROW && i != COL - 1)
        {
            for (k = l; k < COL; k++)
                scale += fabs((double)a[k*ROW+i]);
            if (scale)
            {
                for (k = l; k < COL; k++)
                {
                    a[k*ROW+i] = (double)((double)a[k*ROW+i] / scale);
                    s += ((double)a[k*ROW+i] * (double)a[k*ROW+i]);
                }
                f = (double)a[l*ROW+i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[l*ROW+i] = (double)(f - g);
                for (k = l; k < COL; k++)
                    rv1[k] = (double)a[k*ROW+i] / h;
                if (i != ROW - 1)
                {
                    for (j = l; j < ROW; j++)
                    {
                        for (s = 0.0, k = l; k < COL; k++)
                            s += ((double)a[k*ROW+j] * (double)a[k*ROW+i]);
                        for (k = l; k < COL; k++)
                            a[k*ROW+j] += (double)(s * rv1[k]);
                    }
                }
                for (k = l; k < COL; k++)
                    a[k*ROW+i] = (double)((double)a[k*ROW+i] * scale);
            }
        }
        anorm = MAX(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = COL - 1; i >= 0; i--)
    {
        if (i < COL - 1)
        {
            if (g)
            {
                for (j = l; j < COL; j++)
                    v[i*COL+j] = (double)(((double)a[j*ROW+i] / (double)a[l*ROW+i]) / g);
                /* double division to avoid underflow */
                for (j = l; j < COL; j++)
                {
                    for (s = 0.0, k = l; k < COL; k++)
                        s += ((double)a[k*ROW+i] * (double)v[j*COL+k]);
                    for (k = l; k < COL; k++)
                        v[j*COL+k] += (double)(s * (double)v[i*COL+k]);
                }
            }
            for (j = l; j < COL; j++)
                v[j*COL+i] = v[i*COL+j] = 0.0;
        }
        v[i*COL+i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = COL - 1; i >= 0; i--)
    {
        l = i + 1;
        g = (double)w[i];
        if (i < COL - 1)
            for (j = l; j < COL; j++)
                a[j*ROW+i] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != COL - 1)
            {
                for (j = l; j < COL; j++)
                {
                    for (s = 0.0, k = l; k < ROW; k++)
                        s += ((double)a[i*ROW+k] * (double)a[j*ROW+k]);
                    f = (s / (double)a[i*ROW+i]) * g;
                    for (k = i; k < ROW; k++)
                        a[j*ROW+k] += (double)(f * (double)a[i*ROW+k]);
                }
            }
            for (j = i; j < ROW; j++)
                a[i*ROW+j] = (double)((double)a[i*ROW+j] * g);
        }
        else
        {
            for (j = i; j < ROW; j++)
                a[i*ROW+j] = 0.0;
        }
        ++a[i*ROW+i];
    }

    /* diagonalize the bidiagonal form */
    for (k = COL - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs((double)w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = (double)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (double)h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (-f * h);
                        for (j = 0; j < ROW; j++)
                        {
                            y = (double)a[nm*ROW+j];
                            z = (double)a[i*ROW+j];
                            a[nm*ROW+j] = (double)(y * c + z * s);
                            a[i*ROW+j] = (double)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (double)(-z);
                    for (j = 0; j < COL; j++)
                        v[k*COL+j] = (-v[k*COL+j]);
                }
                break;
            }
            if (its >= 30) {
                free((void*)rv1);
                r_printf("No convergence after 30,000! iterations \n");
                return(0);
            }

            /* shift from bottom 2 x 2 minor */
            x = (double)w[l];
            nm = k - 1;
            y = (double)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < COL; jj++)
                {
                    x = (double)v[j*COL+jj];
                    z = (double)v[i*COL+jj];
                    v[j*COL+jj] = (double)(x * c + z * s);
                    v[i*COL+jj] = (double)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (double)z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < ROW; jj++)
                {
                    y = (double)a[j*ROW+jj];
                    z = (double)a[i*ROW+jj];
                    a[j*ROW+jj] = (double)(y * c + z * s);
                    a[i*ROW+jj] = (double)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (double)x;
        }
    }
    free((void*)rv1);
    return(1);
}

#include "abc_000_warning.h"

