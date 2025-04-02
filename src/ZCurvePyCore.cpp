#include"ZCurvePyCore.h"

/* 
 * Map for converting ASCII chars into one-hot vectors
 *
 * A = [1, 0, 0, 0] G = [0, 1, 0, 0]
 * C = [0, 0, 1, 0] T = [0, 0, 0, 1]
 * 
 * Degenerate symbols are calculated by probabilities
 */
static float ONE_HOT[][4] = 
{   
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0},  
    {1, 0, 0, 0}, {0, W, W, W}, {0, 0, 1, 0}, {W, W, 0, W},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 1, 0, 0}, {W, 0, W, W}, 
    {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, V, 0, V}, {0, 0, 0, 0},
    {V, 0, V, 0}, {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0}, {V, V, 0, 0}, {0, V, V, 0}, {0, 0, 0, 1},
    {0, 0, 0, 1}, {W, W, W, 0}, {V, 0, 0, V}, {0, 0, 0, 0}, 
    {0, 0, V, V}, {1, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {1, 0, 0, 0}, {0, W, W, W}, {0, 0, 1, 0}, {W, W, 0, W},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 1, 0, 0}, {W, 0, W, W}, 
    {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, V, 0, V}, {0, 0, 0, 0}, 
    {V, 0, V, 0}, {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {V, V, 0, 0}, {0, V, V, 0}, {0, 0, 0, 1},
    {0, 0, 0, 1}, {W, W, W, 0}, {V, 0, 0, V}, {0, 0, 0, 0},
    {0, 0, V, V}, {1, 0, 0, 0}
};

/* 
 * Map for converting ASCII chars into one-hot vectors
 * A or G = [1, 0] C or T = [0, 1]
 * 
 * Degenerate symbols are calculated by probabilities
 */
static float RY_HOT[][2] = 
{
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0},
    {1, 0}, {W, R}, {0, 1}, {R, W},
    {0, 0}, {0, 0}, {1, 0}, {W, R}, 
    {V, V}, {0, 0}, {V, V}, {0, 0},
    {V, V}, {V, V}, {0, 0}, {0, 0},
    {0, 0}, {1, 0}, {V, V}, {0, 1},
    {0, 1}, {R, W}, {V, V}, {0, 0}, 
    {0, 1}, {1, 0}, {0, 0}, {0, 0},
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {1, 0}, {W, R}, {0, 1}, {R, W},
    {0, 0}, {0, 0}, {1, 0}, {W, R}, 
    {V, V}, {0, 0}, {V, V}, {0, 0},
    {V, V}, {V, V}, {0, 0}, {0, 0},
    {0, 0}, {1, 0}, {V, V}, {0, 1},
    {0, 1}, {R, W}, {V, V}, {0, 0}, 
    {0, 1}, {1, 0}
};

/* 
 * Map for converting ASCII chars into one-hot vectors
 * A or C = [1, 0] G or T = [0, 1]
 * 
 * Degenerate symbols are calculated by probabilities
 */
static float MK_HOT[][2] = 
{
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0},
    {1, 0}, {W, R}, {1, 0}, {W, R},
    {0, 0}, {0, 0}, {0, 1}, {W, R}, 
    {V, V}, {0, 0}, {0, 1}, {0, 0},
    {1, 0}, {V, V}, {0, 0}, {0, 0},
    {0, 0}, {V, V}, {V, V}, {0, 1},
    {0, 1}, {R, W}, {V, V}, {0, 0}, 
    {V, V}, {1, 0}, {0, 0}, {0, 0},
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {1, 0}, {W, R}, {1, 0}, {W, R},
    {0, 0}, {0, 0}, {0, 1}, {W, R}, 
    {V, V}, {0, 0}, {0, 1}, {0, 0},
    {1, 0}, {V, V}, {0, 0}, {0, 0},
    {0, 0}, {V, V}, {V, V}, {0, 1},
    {0, 1}, {R, W}, {V, V}, {0, 0}, 
    {V, V}, {1, 0}
};

/* 
 * Map for converting ASCII chars into one-hot vectors
 * G or C = [1, 0] A or T = [0, 1]
 * 
 * Degenerate symbols are calculated by probabilities
 */
static float WS_HOT[][2] = 
{
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0},
    {0, 1}, {R, W}, {1, 0}, {W, R},
    {0, 0}, {0, 0}, {1, 0}, {W, R}, 
    {V, V}, {0, 0}, {V, V}, {0, 0},
    {V, V}, {V, V}, {0, 0}, {0, 0},
    {0, 0}, {V, V}, {1, 0}, {0, 1},
    {0, 1}, {R, W}, {0, 1}, {0, 0}, 
    {V, V}, {0, 1}, {0, 0}, {0, 0},
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 1}, {R, W}, {1, 0}, {W, R},
    {0, 0}, {0, 0}, {1, 0}, {W, R},
    {V, V}, {0, 0}, {V, V}, {0, 0}, 
    {V, V}, {V, V}, {0, 0}, {0, 0}, 
    {0, 0}, {V, V}, {1, 0}, {0, 1},
    {0, 1}, {R, W}, {0, 1}, {0, 0},
    {V, V}, {0, 1}
};

/* 
 * Map for converting ASCII chars into one-hot vectors
 * A = [1, 0] T = [0, 1]
 * 
 * Degenerate symbols are calculated by probabilities
 */
static float AT_HOT[][2] = 
{
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0},
    {1, 0}, {0, W}, {0, 0}, {W, W},
    {0, 0}, {0, 0}, {0, 0}, {W, W}, 
    {Q, Q}, {0, 0}, {0, V}, {0, 0},
    {V, 0}, {Q, Q}, {0, 0}, {0, 0},
    {0, 0}, {V, 0}, {0, 0}, {0, 1},
    {0, 1}, {W, 0}, {V, V}, {0, 0}, 
    {0, V}, {1, 0}, {0, 0}, {0, 0},
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {1, 0}, {0, W}, {0, 0}, {W, W},
    {0, 0}, {0, 0}, {0, 0}, {W, W}, 
    {Q, Q}, {0, 0}, {0, V}, {0, 0},
    {V, 0}, {Q, Q}, {0, 0}, {0, 0},
    {0, 0}, {V, 0}, {0, 0}, {0, 1},
    {0, 1}, {W, 0}, {V, V}, {0, 0}, 
    {0, V}, {1, 0}
};

/* 
 * Map for converting ASCII chars into one-hot vectors
 * G = [1, 0] C = [0, 1]
 * 
 * Degenerate symbols are calculated by probabilities
 */
static float GC_HOT[][2] = 
{
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0},
    {0, 0}, {W, W}, {0, 1}, {W, 0},
    {0, 0}, {0, 0}, {1, 0}, {0, W},
    {Q, Q}, {0, 0}, {V, 0}, {0, 0},
    {0, V}, {Q, Q}, {0, 0}, {0, 0},
    {0, 0}, {V, 0}, {V, V}, {0, 0},
    {0, 0}, {W, W}, {0, 0}, {0, 0}, 
    {0, V}, {0, 0}, {0, 0}, {0, 0},
    {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
    {0, 0}, {W, W}, {0, 1}, {W, 0},
    {0, 0}, {0, 0}, {1, 0}, {0, W},
    {Q, Q}, {0, 0}, {V, 0}, {0, 0},
    {0, V}, {Q, Q}, {0, 0}, {0, 0},
    {0, 0}, {V, 0}, {V, V}, {0, 0},
    {0, 0}, {W, W}, {0, 0}, {0, 0}, 
    {0, V}, {0, 0}
};

/* 
 * Map for converting ASCII chars into Z-curve coordinate
 *
 * A = [+1, +1, +1]  G = [+1, -1, -1]
 * C = [-1, +1, -1]  T = [-1, -1, +1]
 * 
 * Degenerate symbols are calculated by weighted vector sum
 */
static float Z_HOT[][3] = {
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0},
    {+1, +1, +1}, {-W, -W, -W}, {-1, +1, -1}, {+W, -W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+1, -1, -1}, {-W, +W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+0, -1, +0}, {+0, +0, +0},
    {+0, +1, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+1, +0, +0}, {+0, +0, -1}, {-1, -1, +1},
    {-1, -1, +1}, {+W, +W, -W}, {+0, +0, +1}, {+0, +0, +0},
    {+1, +0, +0}, {+1, +1, +1}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+1, +1, +1}, {-W, -W, -W}, {-1, +1, -1}, {+W, -W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+1, -1, -1}, {-W, +W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+0, -1, +0}, {+0, +0, +0},
    {+0, +1, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+1, +0, +0}, {+0, +0, -1}, {-1, -1, +1},
    {-1, -1, +1}, {+W, +W, -W}, {+0, +0, +1}, {+0, +0, +0},
    {+1, +0, +0}, {+1, +1, +1}
};
/* 
 * Mean-value smoothing progress for curves' y-values
 * Choose 'same' as the mode of convolution calculation
 * Using variable window so there is no padding method
 */
static void meanSmoothing(float *params, int len, int window) {
    // window = 0 means do nothing
    if (window) {
        float *buf = new float[len], ySum;
        int i, s, e, half = window / 2;

        // Init sum quene
        for (i = 0, ySum = 0; i < half - 1; i ++)
            ySum += params[i];

        /* 
         * In the following progress, every time a value added
         * from the front, a value substract from the back 
         */
        for (i = 0; i < len; i ++) {
            // s points the first element of the region
            s = (i - half) > 0 ? (i - half) : 0;
            // e points the last element of the region
            /* when e = len, it is to avoid adding the last element repeatedly */
            e = (i + half - 1) < len ? (i + half - 1) : len;

            if (s > 0) 
                ySum -= params[s - 1];

            if (e < len) 
                ySum += params[e];
            else if (e == len) 
                e -= 1; // rollback to calculate the length of region correctly

            buf[i] = ySum / (e - s + 1);
        }

        for (i = 0; i < len; i ++)
            params[i] = buf[i];

        delete[] buf;
    }
}

#ifdef __cplusplus
extern "C" {
#endif

void monoTrans(char *seq, int len, float *params, bool norml, bool local) {
    int i;

    for (i = 0; i < len; i ++) {
        params[X] += Z_HOT[seq[i]][X];
        params[Y] += Z_HOT[seq[i]][Y];
        params[Z] += Z_HOT[seq[i]][Z];
    }

    if (norml) {
        params[X] /= len;
        params[Y] /= len;
        params[Z] /= len;
    }
}

void diTrans(char *seq, int len, float *params, bool norml, bool local) {
    float counts[4][4] = {{0.0f}};
    int i, j, sublen = len - 1;
    float denom;

    for (i = 0; i < sublen; i ++)
        for (j = 0; j < 4; j ++) {
            counts[j][X] += ONE_HOT[seq[i]][j] * Z_HOT[seq[i + 1]][X];
            counts[j][Y] += ONE_HOT[seq[i]][j] * Z_HOT[seq[i + 1]][Y];
            counts[j][Z] += ONE_HOT[seq[i]][j] * Z_HOT[seq[i + 1]][Z];
            counts[j][L] += ONE_HOT[seq[i]][j];
        }
        
    for (i = 0; i < 4; i ++, params += 3) {
        params[X] = counts[i][X];
        params[Y] = counts[i][Y];
        params[Z] = counts[i][Z];

        if (norml) {
            denom = local ? counts[i][L] : sublen;

            if (denom > 0) {
                params[X] /= denom;
                params[Y] /= denom;
                params[Z] /= denom;
            }
        }
    }
}

void triTrans(char *seq, int len, float *params, bool norml, bool local) {
    float counts[4][4][4] = {{{0.0f}}};
    int i, j, l, sublen = len - 2;
    float denom;

    for (i = 0; i < sublen; i ++)
        for (j = 0; j < 4; j ++)
        for (l = 0; l < 4; l ++) {
            counts[j][l][X] += ONE_HOT[seq[i]][j] * ONE_HOT[seq[i + 1]][l] * Z_HOT[seq[i + 2]][X];
            counts[j][l][Y] += ONE_HOT[seq[i]][j] * ONE_HOT[seq[i + 1]][l] * Z_HOT[seq[i + 2]][Y];
            counts[j][l][Z] += ONE_HOT[seq[i]][j] * ONE_HOT[seq[i + 1]][l] * Z_HOT[seq[i + 2]][Z];
            counts[j][l][L] += ONE_HOT[seq[i]][j] * ONE_HOT[seq[i + 1]][l];
        }

    for (i = 0; i < 4; i ++)
    for (j = 0; j < 4; j ++, params += 3) {
        params[X] = counts[i][j][X];
        params[Y] = counts[i][j][Y];
        params[Z] = counts[i][j][Z];

        if (norml) {
            denom = local ? counts[i][j][L] : sublen;

            if (denom > 0) {
                params[X] /= denom;
                params[Y] /= denom;
                params[Z] /= denom;
            }
        }
    }
}

void monoPhaseTrans(char *seq, int len, float *params, int phase, bool norml, bool local) {
    float counts[6][3] = {{0.0f}};
    int phaseCount[6] = {0};
    int base = len / phase, residue = len % phase;
    int i, p;

    for (i = 0; i < phase; i ++)
        phaseCount[i] += base;

    for (i = 0; i < residue; i ++)
        phaseCount[i] ++;
        
    for (i = 0; i < len; i ++) {
        p = i % phase;

        counts[p][X] += Z_HOT[seq[i]][X];
        counts[p][Y] += Z_HOT[seq[i]][Y];
        counts[p][Z] += Z_HOT[seq[i]][Z];
    }
        
    for (p = 0; p < phase; p ++, params += 3) {
        params[X] = counts[p][X];
        params[Y] = counts[p][Y];
        params[Z] = counts[p][Z];

        if (norml) {
            params[X] /= phaseCount[p];
            params[Y] /= phaseCount[p];
            params[Z] /= phaseCount[p];
        }
    }
}

void diPhaseTrans(char *seq, int len, float *params, int phase, bool norml, bool local) {
    float counts[6][4][4] = {{{0.0f}}};
    int phaseCount[6] = {0};
    int sublen = len - 1;
    int base = sublen / phase, residue = sublen % phase;
    int i, p, b;
    float denom;

    for (i = 0; i < phase; i ++)
        phaseCount[i] += base;

    for (i = 0; i < residue; i ++)
        phaseCount[i] ++;

    for (i = 0; i < sublen; i++) {
        p = i % phase;
        for (b = 0; b < 4; b ++) {
            counts[p][b][X] += ONE_HOT[seq[i]][b] * Z_HOT[seq[i + 1]][X];
            counts[p][b][Y] += ONE_HOT[seq[i]][b] * Z_HOT[seq[i + 1]][X];
            counts[p][b][Z] += ONE_HOT[seq[i]][b] * Z_HOT[seq[i + 1]][Z];
            counts[p][b][L] += ONE_HOT[seq[i]][b];
        }
    }
        
    for (p = 0; p < phase; p ++)
        for (b = 0; b < 4; b ++, params += 3) {
            params[X] = counts[p][b][X];
            params[Y] = counts[p][b][Y];
            params[Z] = counts[p][b][Z];

            if (norml) {
                denom = local ? counts[p][b][L] : phaseCount[p];

                if (denom > 0) {
                    params[X] /= denom;
                    params[Y] /= denom;
                    params[Z] /= denom;
                }
            }
        }
}

void triPhaseTrans(char *seq, int len, float *params, int phase, bool norml, bool local) {
    float counts[6][4][4][4] = {{{{0.0f}}}};
    int phaseCount[6] = {0};
    int sublen = len - 2;
    int base = sublen / phase, residue = sublen % phase;
    int i, p, s, b;
    float denom;

    for (i = 0; i < phase; i ++)
        phaseCount[i] += base;

    for (i = 0; i < residue; i ++)
        phaseCount[i] ++;

    for (i = 0; i < sublen; i ++) {
        p = i % phase;
        for (s = 0; s < 4; s ++)
        for (b = 0; b < 4; b ++) {
            counts[p][s][b][X] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * Z_HOT[seq[i + 2]][X];
            counts[p][s][b][Y] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * Z_HOT[seq[i + 2]][Y];
            counts[p][s][b][Z] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * Z_HOT[seq[i + 2]][Z];
            counts[p][s][b][L] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b];
        }
    }

    for (p = 0; p < phase; p ++)
    for (s = 0; s < 4; s ++)
    for (b = 0; b < 4; b ++, params += 3) {
        params[X] = counts[p][s][b][X];
        params[Y] = counts[p][s][b][Y];
        params[Z] = counts[p][s][b][Z];

        if (norml) {
            denom = local ? counts[p][s][b][L] : phaseCount[p];

            if (denom > 0) {
                params[X] /= denom;
                params[Y] /= denom;
                params[Z] /= denom;
            }
        }
    }
}

void quartPhaseTrans(char *seq, int len, float *params, int phase, bool norml, bool local) {
    float counts[6][4][4][4][4] = {{{{{0.0f}}}}};
    int phaseCount[6] = {0};
    int sublen = len - 3;
    int base = sublen / phase, residue = sublen % phase;
    int i, p, s, b, t;
    float denom;

    for (i = 0; i < phase; i ++)
        phaseCount[i] += base;

    for (i = 0; i < residue; i ++)
        phaseCount[i] ++;

    for (i = 0; i < sublen; i ++) {
        p = i % phase;
        for (s = 0; s < 4; s ++)
        for (b = 0; b < 4; b ++)
        for (t = 0; t < 4; t ++) {
            counts[p][s][b][t][X] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * ONE_HOT[seq[i + 2]][t] * Z_HOT[seq[i + 3]][X];
            counts[p][s][b][t][Y] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * ONE_HOT[seq[i + 2]][t] * Z_HOT[seq[i + 3]][Y];
            counts[p][s][b][t][Z] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * ONE_HOT[seq[i + 2]][t] * Z_HOT[seq[i + 3]][Z];
            counts[p][s][b][t][L] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * ONE_HOT[seq[i + 2]][t];
        }
    }

    for (p = 0; p < phase; p ++)
    for (s = 0; s < 4; s ++)
    for (b = 0; b < 4; b ++)
    for (t = 0; t < 4; t ++, params += 3) {
        params[X] = counts[p][s][b][t][X];
        params[Y] = counts[p][s][b][t][Y];
        params[Z] = counts[p][s][b][t][Z];

        if (norml) {
            denom = local ? counts[p][s][b][t][L] : phaseCount[p];

            if (denom > 0) {
                params[X] /= phaseCount[p];
                params[Y] /= phaseCount[p];
                params[Z] /= phaseCount[p];
            }
        }
    }
}

void quintPhaseTrans(char *seq, int len, float *params, int phase, bool norml, bool local) {
    float counts[6][4][4][4][4][4] = {{{{{0.0f}}}}};
    int phaseCount[6] = {0};
    int sublen = len - 4;
    int base = sublen / phase, residue = sublen % phase;
    int i, p, s, b, t, l;
    float denom;

    for (i = 0; i < phase; i ++)
        phaseCount[i] += base;

    for (i = 0; i < residue; i ++)
        phaseCount[i] ++;

    for (i = 0; i < sublen; i ++) {
        p = i % phase;
        for (s = 0; s < 4; s ++)
        for (b = 0; b < 4; b ++)
        for (t = 0; t < 4; t ++)
        for (l = 0; l < 3; l ++) {
            counts[p][s][b][t][l][X] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * ONE_HOT[seq[i + 2]][t] *
                                        ONE_HOT[seq[i + 3]][l] * Z_HOT[seq[i + 4]][X];
            counts[p][s][b][t][l][Y] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * ONE_HOT[seq[i + 2]][t] *
                                        ONE_HOT[seq[i + 3]][l] * Z_HOT[seq[i + 4]][Y];
            counts[p][s][b][t][l][Z] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * ONE_HOT[seq[i + 2]][t] *
                                        ONE_HOT[seq[i + 3]][l] * Z_HOT[seq[i + 4]][Z];
            counts[p][s][b][t][l][L] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * 
                                        ONE_HOT[seq[i + 2]][t] * ONE_HOT[seq[i + 3]][l];
        }
    }

    for (p = 0; p < phase; p ++)
    for (s = 0; s < 4; s ++)
    for (b = 0; b < 4; b ++)
    for (t = 0; t < 4; t ++)
    for (l = 0; l < 4; l ++, params += 3) {
        params[X] = counts[p][s][b][t][l][X];
        params[Y] = counts[p][s][b][t][l][Y];
        params[Z] = counts[p][s][b][t][l][Z];

        if (norml) {
            denom = local ? counts[p][s][b][t][l][L] : phaseCount[p];

            if (denom > 0) {
                params[X] /= denom;
                params[Y] /= denom;
                params[Z] /= denom;
            }
        }
    }
}

void sexPhaseTrans(char *seq, int len, float *params, int phase, bool norml, bool local) {
    float counts[6][4][4][4][4][4][4] = {{{{{{0.0f}}}}}};
    int phaseCount[6] = {0};
    int sublen = len - 5;
    int base = sublen / phase, residue = sublen % phase;
    int i, p, s, b, t, l, m;
    float denom;

    for (i = 0; i < phase; i ++)
        phaseCount[i] += base;

    for (i = 0; i < residue; i ++)
        phaseCount[i] ++;

    for (i = 0; i < sublen; i ++) {
        p = i % phase;
        for (s = 0; s < 4; s ++)
        for (b = 0; b < 4; b ++)
        for (t = 0; t < 4; t ++)
        for (l = 0; l < 4; l ++)
        for (m = 0; m < 4; m ++) {
            counts[p][s][b][t][l][m][X] += ONE_HOT[seq[i + 0]][s] * ONE_HOT[seq[i + 1]][b] * ONE_HOT[seq[i + 2]][t] * 
                                           ONE_HOT[seq[i + 3]][l] * ONE_HOT[seq[i + 4]][m] * Z_HOT[seq[i + 5]][X];
            counts[p][s][b][t][l][m][Y] += ONE_HOT[seq[i + 0]][s] * ONE_HOT[seq[i + 1]][b] * ONE_HOT[seq[i + 2]][t] * 
                                           ONE_HOT[seq[i + 3]][l] * ONE_HOT[seq[i + 4]][m] * Z_HOT[seq[i + 5]][Y];
            counts[p][s][b][t][l][m][Z] += ONE_HOT[seq[i + 0]][s] * ONE_HOT[seq[i + 1]][b] * ONE_HOT[seq[i + 2]][t] * 
                                           ONE_HOT[seq[i + 3]][l] * ONE_HOT[seq[i + 4]][m] * Z_HOT[seq[i + 5]][Z];
            counts[p][s][b][t][l][m][L] += ONE_HOT[seq[i + 0]][s] * ONE_HOT[seq[i + 1]][b] * ONE_HOT[seq[i + 2]][t] * 
                                           ONE_HOT[seq[i + 3]][l] * ONE_HOT[seq[i + 4]][m];
        }
    }

    for (p = 0; p < phase; p ++)
    for (s = 0; s < 4; s ++)
    for (b = 0; b < 4; b ++)
    for (t = 0; t < 4; t ++)
    for (l = 0; l < 4; l ++)
    for (m = 0; m < 4; m ++, params += 3) {
        params[X] = counts[p][s][b][t][l][m][X];
        params[Y] = counts[p][s][b][t][l][m][Y];
        params[Z] = counts[p][s][b][t][l][m][Z];

        if (norml) {
            denom = local ? counts[p][s][b][t][l][m][L] : phaseCount[p];

            if (denom > 0) {
                params[X] /= denom;
                params[Y] /= denom;
                params[Z] /= denom;
            }
        }
    }
}

void (*PhaseTrans[6])(char *, int, float *, int, bool, bool) = {
    monoPhaseTrans,  diPhaseTrans,    triPhaseTrans,
    quartPhaseTrans, quintPhaseTrans, sexPhaseTrans
};

void kPhaseTrans(char *seq, int len, float *params, int k, int phase, bool norml, bool local) {
    if (k < 7) {
        (*PhaseTrans[k - 1])(seq, len, params, phase, norml, local);
    } else {
        /* Experimental function */
        /* Don't set k >= 7 in actual use */
        float phaseCount[6] = {0};
        int base = (len - k + 1) / phase, residue = (len - k + 1) % phase;
        int i, j, m, num, p;
        float denom;

        for (i = 0; i < phase; i ++)
            phaseCount[i] += base;

        for (i = 0; i < residue; i ++)
            phaseCount[i] ++;

        for (i = 0, num = 1; i < k; i ++) num *= 4;

        float *counts[6], weight;

        for (i = 0; i < phase; i ++)
            counts[i] = new float[num]();

        for (i = 0; i < len - k + 1; i ++) {
            p = i % phase;
            for (k = 0; k < num; k ++) {
                for (j = 0, m = k, weight = 1.0; j < k; j ++, m /= 4)
                    weight *= ONE_HOT[seq[i + k - j - 1]][m % 4];
            
                counts[p][k] += weight;
            }
        }

        for (i = 0; i < phase; i ++)
            for (j = 0; j < num; j += 4, params += 3) {
                params[X] = counts[i][j + 0] + counts[i][j + 1] - counts[i][j + 2] - counts[i][j + 3];
                params[Y] = counts[i][j + 0] + counts[i][j + 2] - counts[i][j + 1] - counts[i][j + 3];
                params[Z] = counts[i][j + 0] + counts[i][j + 3] - counts[i][j + 2] - counts[i][j + 1];

                if (norml) {
                    denom = local ? (counts[i][j + 0] + counts[i][j + 1] + counts[i][j + 2] + counts[i][j + 3]) : phaseCount[i];

                    if (denom > 0) {
                        params[X] /= denom;
                        params[Y] /= denom;
                        params[Z] /= denom;
                    }
                }
            }
            
        for (i = 0; i < phase; i ++)
            delete[] counts[i];
    }
}

float genomeOrderIndex(char *seq, int len) {
    float counts[4] = {0}, s;
    int i, j;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 4; j ++)
            counts[j] += ONE_HOT[seq[i]][j];
    
    for (i = 0, s = 0; i < 4; i ++)
        s += (counts[i] / len * counts[i] / len);

    return s;
}

float ryOrderIndex(char *seq, int len) {
    float counts[2] = {0}, s;
    int i, j;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 2; j ++)
            counts[j] += RY_HOT[seq[i]][j];
    
    for (i = 0, s = 0; i < 2; i ++)
        s += (counts[i] / len * counts[i] / len);

    return s;
}

float mkOrderIndex(char *seq, int len) {
    float counts[2] = {0}, s;
    int i, j;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 2; j ++)
            counts[j] += MK_HOT[seq[i]][j];
    
    for (i = 0, s = 0; i < 2; i ++)
        s += (counts[i] / len * counts[i] / len);

    return s;
}

float wsOrderIndex(char *seq, int len) {
    float counts[2] = {0}, s;
    int i, j;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 2; j ++)
            counts[j] += WS_HOT[seq[i]][j];
    
    for (i = 0, s = 0; i < 2; i ++)
        s += (counts[i] / len * counts[i] / len);

    return s;
}

float atOrderIndex(char *seq, int len) {
    float counts[2] = {0}, s;
    int i, j;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 2; j ++)
            counts[j] += AT_HOT[seq[i]][j];
    
    for (i = 0, s = 0; i < 2; i ++)
        s += (counts[i] / len * counts[i] / len);

    return s;
}

float gcOrderIndex(char *seq, int len) {
    float counts[2] = {0}, s;
    int i, j;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 2; j ++)
            counts[j] += GC_HOT[seq[i]][j];
    
    for (i = 0, s = 0; i < 2; i ++)
        s += (counts[i] / len * counts[i] / len);

    return s;
}

float CpGOrderIndex(char *seq, int len) {
    float counts[2] = {0}, s;
    int i, j, k;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 4; j ++)
            for (k = 0; k < 4; k ++)
                counts[j == 2 && k == 1] += ONE_HOT[seq[i]][j] * ONE_HOT[seq[i + 1]][k];
    
    for (i = 0, s = 0; i < 2; i ++)
        s += (counts[i] / (len - 1) * counts[i] / (len - 1));
    
    return s;
}

void zTrans(char *seq, int len, float **params, int window) {
    float counts[3] = {0.0f};
    int i, j;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 3; j ++) {
            counts[j] += Z_HOT[seq[i]][j];
            params[j][i] = counts[j];
        }

    if (window) for (i = 0; i < 3; i ++)
        meanSmoothing(params[i], len, window);
}

void ryTrans(char *seq, int len, float *params, int window) {
    float counts = 0.0f;
    int i;

    for (i = 0; i < len; i ++) {
        counts += Z_HOT[seq[i]][X];
        params[i] = counts;
    }

    meanSmoothing(params, len, window);
}

void mkTrans(char *seq, int len, float *params, int window) {
    float counts = 0.0f;
    int i;

    for (i = 0; i < len; i ++) {
        counts += Z_HOT[seq[i]][Y];
        params[i] = counts;
    }

    meanSmoothing(params, len, window);
}

void wsTrans(char *seq, int len, float *params, int window) {
    float counts = 0.0f;
    int i;

    for (i = 0; i < len; i ++) {
        counts += Z_HOT[seq[i]][Z];
        params[i] = counts;
    }

    meanSmoothing(params, len, window);
}

void atTrans(char *seq, int len, float *params, int window) {
    float counts[4] = {0.0f};
    int i, j;

    for (i = 0; i < len; i ++, params ++) {
        for (j = 0; j < 4; j += 3)
            counts[j] += ONE_HOT[seq[i]][j];

        *params = counts[0] - counts[3];
    }

    meanSmoothing(params - len, len, window);
}

void gcTrans(char *seq, int len, float *params, int window) {
    float counts[4] = {0.0f};
    int i, j;

    for (i = 0; i < len; i ++, params ++) {
        for (j = 1; j < 3; j ++)
            counts[j] += ONE_HOT[seq[i]][j];

        *params = counts[1] - counts[2];
    }

    meanSmoothing(params - len, len, window);
}

float xPrimeTrans(char *seq, int len, float *params, int window) {
    float counts = 0.0f;
    double xAvr, x2Sum, ySum, xySum, kp;
    int i;

    for (i = 0; i < len; i ++) {
        counts += Z_HOT[seq[i]][X];
        params[i] = counts;
    }

    xAvr = (len - 1) / 2.0f;
    x2Sum = (2 * len - 1) / 6.0f * len * (len - 1);
    
    for (i = 0, ySum = 0, xySum = 0; i < len; i ++)
        ySum += params[i], xySum += i * params[i];

    kp = (xySum - xAvr * ySum) / (x2Sum - xAvr * xAvr * len);
    
    for (i = 1; i < len; i++)
        params[i] -= ((float) kp) * i;
    
    meanSmoothing(params, len, window);
    
    return (float) kp;
}

float yPrimeTrans(char *seq, int len, float *params, int window) {
    float counts = 0.0f;
    double xAvr, x2Sum, ySum, xySum, kp;
    int i;

    for (i = 0; i < len; i ++) {
        counts += Z_HOT[seq[i]][Y];
        params[i] = counts;
    }

    xAvr = (len - 1) / 2.0f;
    x2Sum = (2 * len - 1) / 6.0f * len * (len - 1);
    
    for (i = 0, ySum = 0, xySum = 0; i < len; i ++)
        ySum += params[i], xySum += i * params[i];

    kp = (xySum - xAvr * ySum) / (x2Sum - xAvr * xAvr * len);
    
    for (i = 1; i < len; i++)
        params[i] -= ((float) kp) * i;
    
    meanSmoothing(params, len, window);
    
    return (float) kp;
}

float zPrimeTrans(char *seq, int len, float *params, int window) {
    float counts = 0.0f;
    double xAvr, x2Sum, ySum, xySum, kp;
    int i;

    for (i = 0; i < len; i ++) {
        counts += Z_HOT[seq[i]][Z];
        params[i] = counts;
    }

    xAvr = (len - 1) / 2.0f;
    x2Sum = (2 * len - 1) / 6.0f * len * (len - 1);

    for (i = 0, ySum = 0, xySum = 0; i < len; i ++)
        ySum += params[i], xySum += i * params[i];

    kp = (xySum - xAvr * ySum) / (x2Sum - xAvr * xAvr * len);
    
    for (i = 1; i < len; i++)
        params[i] -= ((float) kp) * i;
    
    meanSmoothing(params, len, window);
    
    return (float) kp;
}

float atPrimeTrans(char *seq, int len, float *params, int window) {
    float counts[4] = { 0.0f };
    double xAvr, x2Sum, ySum, xySum, kp;
    int i, j;

    for (i = 0; i < len; i ++) {
        for (j = 0; j < 4; j += 3)
            counts[j] += ONE_HOT[seq[i]][j];

        params[i] = counts[0] - counts[3];
    }

    xAvr = (len - 1) / 2.0f;
    x2Sum = (2 * len - 1) / 6.0f * len * (len - 1);

    for (i = 0, ySum = 0, xySum = 0; i < len; i ++)
        ySum += params[i], xySum += i * params[i];

    kp = (xySum - xAvr * ySum) / (x2Sum - xAvr * xAvr * len);
    
    for (i = 1; i < len; i++)
        params[i] -= ((float) kp) * i;
    
    meanSmoothing(params, len, window);
    
    return (float) kp;
}

float gcPrimeTrans(char *seq, int len, float *params, int window) {
    float counts[4] = { 0.0f };
    double xAvr, x2Sum, ySum, xySum, kp;
    int i, j;

    for (i = 0; i < len; i ++) {
        for (j = 1; j < 3; j ++)
            counts[j] += ONE_HOT[seq[i]][j];

        params[i] = counts[1] - counts[2];
    }

    xAvr = (len - 1) / 2.0f;
    x2Sum = (2 * len - 1) / 6.0f * len * (len - 1);

    for (i = 0, ySum = 0, xySum = 0; i < len; i ++)
        ySum += params[i], xySum += i * params[i];

    kp = (xySum - xAvr * ySum) / (x2Sum - xAvr * xAvr * len);
    
    for (i = 1; i < len; i++)
        params[i] -= ((float) kp) * i;
    
    meanSmoothing(params, len, window);
    
    return (float) kp;
}

float CpGPrimeTrans(char *seq, int len, float *params, int window) {
    float counts[2] = {0};
    double xAvr, x2Sum, ySum, xySum, kp;
    int i, j, k;

    for (i = 0; i < len - 1; i ++) {
        for (j = 0; j < 4; j ++)
            for (k = 0; k < 4; k ++)
                counts[j == 2 && k == 1] += ONE_HOT[seq[i]][j] * ONE_HOT[seq[i + 1]][k];
        
        params[i] = counts[1] - counts[0];
    }

    xAvr = len / 2.0f - 1;
    x2Sum = (2 * len - 3) / 6.0f * (len - 1) * (len - 2);

    for (i = 0, ySum = 0, xySum = 0; i < len - 1; i ++)
        ySum += params[i], xySum += i * params[i];

    kp = (xySum - xAvr * ySum) / (x2Sum - xAvr * xAvr * (len - 1));
    
    for (i = 1; i < len - 1; i ++)
        params[i] -= ((float) kp) * i;

    params[i] = params[i - 1];
    meanSmoothing(params, len - 1, window);

    return (float) kp;
}

int genomeDeltaSTrans(char *seq, int len, float *params, int window, float &max) {
    float p[4] = {0}, q[4] = {0}, dif[4] = {0}, s;
    int i, j, maxPoint = -1;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 4; j ++)
            q[j] += ONE_HOT[seq[i]][j];
    
    for (i = 0; i < len - 1; i ++) {
        for (j = 0, s = 0; j < 4; j ++) {
            p[j] += ONE_HOT[seq[i]][j];
            q[j] -= ONE_HOT[seq[i]][j];
            dif[j] = p[j] / (i + 1) - q[j] / (len - i - 1);
            s += dif[j] * dif[j];
        }

        params[i] = (s * (i + 1) / len * (len - i - 1));
    }

    params[len - 1] = 0;

    meanSmoothing(params, len, window);

    for (i = 0, max = -1; i < len; i ++)
        if (params[i] > max)
            max = params[i], maxPoint = i;
    
    return maxPoint;
}

int ryDeltaSTrans(char *seq, int len, float *params, int window, float &max) {
    float p[2] = {0}, q[2] = {0}, dif[2] = {0}, s;
    int i, j, maxPoint = -1;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 2; j ++)
            q[j] += RY_HOT[seq[i]][j];
    
    for (i = 0; i < len - 1; i ++) {
        for (j = 0, s = 0; j < 2; j ++) {
            p[j] += RY_HOT[seq[i]][j];
            q[j] -= RY_HOT[seq[i]][j];
            dif[j] = p[j] / (i + 1) - q[j] / (len - i - 1);
            s += dif[j] * dif[j];
        }

        params[i] = (s * (i + 1) / len * (len - i - 1));
    }

    params[len - 1] = 0;

    meanSmoothing(params, len, window);

    for (i = 0, max = -1; i < len; i ++)
        if (params[i] > max)
            max = params[i], maxPoint = i;
    
    return maxPoint;
}

int mkDeltaSTrans(char *seq, int len, float *params, int window, float &max) {
    float p[2] = {0}, q[2] = {0}, dif[2] = {0}, s;
    int i, j, maxPoint = -1;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 2; j ++)
            q[j] += MK_HOT[seq[i]][j];
    
    for (i = 0; i < len - 1; i ++) {
        for (j = 0, s = 0; j < 2; j ++) {
            p[j] += MK_HOT[seq[i]][j];
            q[j] -= MK_HOT[seq[i]][j];
            dif[j] = p[j] / (i + 1) - q[j] / (len - i - 1);
            s += dif[j] * dif[j];
        }

        params[i] = (s * (i + 1) / len * (len - i - 1));
    }

    params[len - 1] = 0;

    meanSmoothing(params, len, window);

    for (i = 0, max = -1; i < len; i ++)
        if (params[i] > max)
            max = params[i], maxPoint = i;
    
    return maxPoint;
}

int wsDeltaSTrans(char *seq, int len, float *params, int window, float &max) {
    float p[2] = {0}, q[2] = {0}, dif[2] = {0}, s;
    int i, j, maxPoint = -1;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 2; j ++)
            q[j] += WS_HOT[seq[i]][j];
    
    for (i = 0; i < len - 1; i ++) {
        for (j = 0, s = 0; j < 2; j ++) {
            p[j] += WS_HOT[seq[i]][j];
            q[j] -= WS_HOT[seq[i]][j];
            dif[j] = p[j] / (i + 1) - q[j] / (len - i - 1);
            s += dif[j] * dif[j];
        }

        params[i] = (s * (i + 1) / len * (len - i - 1));
    }

    params[len - 1] = 0;

    meanSmoothing(params, len, window);

    for (i = 0, max = -1; i < len; i ++)
        if (params[i] > max)
            max = params[i], maxPoint = i;
    
    return maxPoint;
}

int atDeltaSTrans(char *seq, int len, float *params, int window, float &max) {
    float p[2] = {0}, q[2] = {0}, dif[2] = {0}, s;
    int i, j, maxPoint = -1;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 2; j ++)
            q[j] += AT_HOT[seq[i]][j];
    
    for (i = 0; i < len - 1; i ++) {
        for (j = 0, s = 0; j < 2; j ++) {
            p[j] += AT_HOT[seq[i]][j];
            q[j] -= AT_HOT[seq[i]][j];
            dif[j] = p[j] / (i + 1) - q[j] / (len - i - 1);
            s += dif[j] * dif[j];
        }

        params[i] = (s * (i + 1) / len * (len - i - 1));
    }

    params[len - 1] = 0;

    meanSmoothing(params, len, window);

    for (i = 0, max = -1; i < len; i ++)
        if (params[i] > max)
            max = params[i], maxPoint = i;
    
    return maxPoint;
}

int gcDeltaSTrans(char *seq, int len, float *params, int window, float &max) {
    float p[2] = {0}, q[2] = {0}, dif[2] = {0}, s;
    int i, j, maxPoint = -1;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 2; j ++)
            q[j] += GC_HOT[seq[i]][j];
    
    for (i = 0; i < len - 1; i ++) {
        for (j = 0, s = 0; j < 2; j ++) {
            p[j] += GC_HOT[seq[i]][j];
            q[j] -= GC_HOT[seq[i]][j];
            dif[j] = p[j] / (i + 1) - q[j] / (len - i - 1);
            s += dif[j] * dif[j];
        }

        params[i] = (s * (i + 1) / len * (len - i - 1));
    }

    params[len - 1] = 0;

    meanSmoothing(params, len, window);

    for (i = 0, max = -1; i < len; i ++)
        if (params[i] > max)
            max = params[i], maxPoint = i;
    
    return maxPoint;
}

int CpGDeltaSTrans(char *seq, int len, float *params, int window, float &max) {
    float p[2] = {0}, q[2] = {0}, dif[2] = {0}, s;
    int i, j, k, maxPoint = -1;

    for (i = 0; i < len; i ++)
        for (j = 0; j < 4; j ++)
            for (k = 0; k < 4; k ++)
                q[j == 2 && k == 1] += ONE_HOT[seq[i]][j] * ONE_HOT[seq[i + 1]][k];
    
    for (i = 0; i < len - 1; i ++) {
        for (j = 0; j < 4; j ++)
            for (k = 0; k < 4; k ++) {
                p[j == 2 && k == 1] += ONE_HOT[seq[i]][j] * ONE_HOT[seq[i + 1]][k];
                q[j == 2 && k == 1] -= ONE_HOT[seq[i]][j] * ONE_HOT[seq[i + 1]][k];
            }
        
        for (j = 0, s = 0; j < 2; j ++) {
            dif[j] = p[j] / (i + 1) - q[j] / (len - i - 1);
            s += dif[j] * dif[j];
        }

        params[i] =  (s * (i + 1) / len * (len - i - 1));
    }

    params[len - 1] = 0;

    meanSmoothing(params, len, window);

    for (i = 0, max = -1; i < len; i ++)
        if (params[i] > max)
            max = params[i], maxPoint = i;
    
    return maxPoint;
}

void tetrahedron(char *seq, int len, float **params) {
    for (int i = 0; i < 3; i ++)
        for (int j = 0; j < len; j ++)
            params[i][j] = Z_HOT[seq[j]][i];
}


#ifdef __cplusplus
}
#endif
