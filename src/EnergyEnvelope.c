//#include "EnergyEnvelope.h"
//
///* if you rely on write() somewhere, keep extern; otherwise harmless */
//extern long write(int fd, const void *buf, unsigned long count);
//
///* --------------------
//   Math replacements (no math.h)
//   -------------------- */
//static double EnergyE_fabs(double x) {
//    return (x < 0.0) ? -x : x;
//}
//
//static double EnergyE_sqrt(double x) {
//    double r;
//    int i;
//    if (x <= 0.0) return 0.0;
//    r = x * 0.5;
//    if (r <= 0.0) r = 1.0;
//    for (i = 0; i < 20; ++i) {
//        r = 0.5 * (r + x / r);
//    }
//    return r;
//}
//
//static double EnergyE_cos(double x) {
//    double x2, x4, x6;
//    while (x > ENERGYE_PI) x -= 2.0 * ENERGYE_PI;
//    while (x < -ENERGYE_PI) x += 2.0 * ENERGYE_PI;
//    x2 = x * x;
//    x4 = x2 * x2;
//    x6 = x4 * x2;
//    return 1.0 - x2/2.0 + x4/24.0 - x6/720.0;
//}
//
//static double EnergyE_sin(double x) {
//    double x2, x3, x5, x7;
//    while (x > ENERGYE_PI) x -= 2.0 * ENERGYE_PI;
//    while (x < -ENERGYE_PI) x += 2.0 * ENERGYE_PI;
//    x2 = x * x;
//    x3 = x2 * x;
//    x5 = x3 * x2;
//    x7 = x5 * x2;
//    return x - x3/6.0 + x5/120.0 - x7/5040.0;
//}
//
///* Simple string copy without stdlib */
//void EnergyE_strcpy(char *dest, const char *src) {
//    int i;
//    if (!dest || !src) return;
//    i = 0;
//    while (src[i] != '\0' && i < 19) {
//        dest[i] = src[i];
//        i++;
//    }
//    dest[i] = '\0';
//}
//
///* --------------------
//   Helper buffers & state
//   -------------------- */
//
///* NOTE: changed to 6 columns so we can store Ax,Ay,Az,Gx,Gy,Gz per sample */
//static double E_window[WINDOW_SAMPLES][6];
//static int E_window_index = 0;
//static int E_cooldown = 0;
//static int E_windows_processed = 0;
//
///* FCFR stats (Welford) */
//static double fcfr_mean = 0.0;
//static double fcfr_M2 = 0.0;
//static int fcfr_count = 0;
//static double adaptive_threshold = FCFR_THRESHOLD;
//
///* last detection tracking */
//static int lastDetectedIdx = -100000;
//
///* --------------------
//   DFT (naive) - C90 style
//   -------------------- */
//static void simpleDFT(const double *in, double *out_real, double *out_imag, int N) {
//    int k, n;
//    for (k = 0; k < N; k++) {
//        double sum_r = 0.0;
//        double sum_i = 0.0;
//        for (n = 0; n < N; n++) {
//            double angle;
//            angle = -2.0 * ENERGYE_PI * (double)k * (double)n / (double)N;
//            sum_r += in[n] * EnergyE_cos(angle);
//            sum_i += in[n] * EnergyE_sin(angle);
//        }
//        out_real[k] = sum_r;
//        out_imag[k] = sum_i;
//    }
//}
//
///* Teager Energy Operator */
//static void computeTeagerEnergy(const double *x, double *E, int n) {
//    int i;
//    for (i = 1; i < n - 1; i++) {
//        E[i] = x[i] * x[i] - x[i - 1] * x[i + 1];
//    }
//    E[0] = 0.0;
//    E[n - 1] = 0.0;
//}
//
///* computeFCFR - returns 1 if fcfr_val > adaptive_threshold, also updates running stats */
//static int computeFCFR(const double *E2, int n, double fs, double *max_freq, double *fcfr_val) {
//    int i;
//    int upper;
//    double max_val;
//    double mean;
//    int max_idx;
//    int count;
//    double val;
//    double delta;
//    double fcfr_std;
//
//    max_val = 0.0;
//    mean = 0.0;
//    max_idx = 0;
//    count = 0;
//
//    upper = n / 2;
//    for (i = 1; i < upper; i++) {
//        val = EnergyE_fabs(E2[i]);
//        mean += val;
//        count++;
//        if (val > max_val) {
//            max_val = val;
//            max_idx = i;
//        }
//    }
//
//    if (count == 0) {
//        *fcfr_val = 0.0;
//        *max_freq = 0.0;
//        return 0;
//    }
//
//    mean /= (double)count;
//    if (mean < EPS) mean = EPS;
//
//    *fcfr_val = max_val / mean;
//    *max_freq = (double)max_idx * fs / (double)n;
//
//    /* update running stats (Welford) */
//    fcfr_count++;
//    delta = *fcfr_val - fcfr_mean;
//    fcfr_mean += delta / (double)fcfr_count;
//    fcfr_M2 += delta * (*fcfr_val - fcfr_mean);
//    if (fcfr_count > 1) {
//        fcfr_std = EnergyE_sqrt(fcfr_M2 / (double)(fcfr_count - 1));
//    } else {
//        fcfr_std = 0.0;
//    }
//    adaptive_threshold = fcfr_mean + FCFR_WEIGHT * fcfr_std;
//
//    return ((*fcfr_val) > adaptive_threshold) ? 1 : 0;
//}
//
///* --------------------
//   Public init
//   -------------------- */
//void EnergyEnvelope_init(void) {
//    int i, j;
//    E_window_index = 0;
//    E_cooldown = 0;
//    E_windows_processed = 0;
//    fcfr_mean = 0.0;
//    fcfr_M2 = 0.0;
//    fcfr_count = 0;
//    adaptive_threshold = FCFR_THRESHOLD;
//    lastDetectedIdx = -100000;
//
//    for (i = 0; i < WINDOW_SAMPLES; i++) {
//        for (j = 0; j < 6; j++) {
//            E_window[i][j] = 0.0;
//        }
//    }
//}
//
///* --------------------
//   Main processing (C90-style: all locals declared at top)
//   -------------------- */
//int EnergyEnvelope_processSample(double Ax, double Ay, double Az, double Gx, double Gy, double Gz) {
//    int i, j;
//    int result;
//    int minSeparation;
//    double amplitude;
//    double accel_magnitude;
//    double gyro_magnitude;
//    double x[WINDOW_SAMPLES];
//    double E1[WINDOW_SAMPLES];
//    double real[FFT_SIZE];
//    double imag[FFT_SIZE];
//    double mag[FFT_SIZE];
//    double E2[FFT_SIZE];
//    double fs;
//    double fcfr_val;
//    double max_freq;
//    int detected;
//    double maxE;
//    int maxIdx;
//
//    /* Always insert new sample */
//    if (E_window_index < WINDOW_SAMPLES) {
//        /* append */
//        E_window[E_window_index][0] = Ax;
//        E_window[E_window_index][1] = Ay;
//        E_window[E_window_index][2] = Az;
//        E_window[E_window_index][3] = Gx;
//        E_window[E_window_index][4] = Gy;
//        E_window[E_window_index][5] = Gz;
//        E_window_index++;
//        if (E_window_index < WINDOW_SAMPLES) {
//            return 0; /* still filling */
//        }
//        /* else continue to processing */
//    } else {
//        /* shift left by one and append */
//        for (i = 1; i < WINDOW_SAMPLES; i++) {
//            E_window[i - 1][0] = E_window[i][0];
//            E_window[i - 1][1] = E_window[i][1];
//            E_window[i - 1][2] = E_window[i][2];
//            E_window[i - 1][3] = E_window[i][3];
//            E_window[i - 1][4] = E_window[i][4];
//            E_window[i - 1][5] = E_window[i][5];
//        }
//        E_window[WINDOW_SAMPLES - 1][0] = Ax;
//        E_window[WINDOW_SAMPLES - 1][1] = Ay;
//        E_window[WINDOW_SAMPLES - 1][2] = Az;
//        E_window[WINDOW_SAMPLES - 1][3] = Gx;
//        E_window[WINDOW_SAMPLES - 1][4] = Gy;
//        E_window[WINDOW_SAMPLES - 1][5] = Gz;
//        /* increment windows processed (you had this earlier when shift happened) */
//        E_windows_processed++;
//    }
//
//    /* Update cooldown */
//    if (E_cooldown > 0) E_cooldown--;
//
//    /* Build input vector x from Ax,Ay,Az,Gx,Gy,Gz magnitudes */
//    for (i = 0; i < WINDOW_SAMPLES; i++) {
//        accel_magnitude = EnergyE_sqrt(E_window[i][0] * E_window[i][0] +
//                                       E_window[i][1] * E_window[i][1] +
//                                       E_window[i][2] * E_window[i][2]);
//        gyro_magnitude = EnergyE_sqrt(E_window[i][3] * E_window[i][3] +
//                                      E_window[i][4] * E_window[i][4] +
//                                      E_window[i][5] * E_window[i][5]);
//        x[i] = accel_magnitude + gyro_magnitude;
//    }
//
//    /* Teager Energy Operator */
//    computeTeagerEnergy(x, E1, WINDOW_SAMPLES);
//
//    /* DFT prep */
//    for (i = 0; i < FFT_SIZE; i++) {
//        real[i] = 0.0;
//        imag[i] = 0.0;
//        mag[i] = 0.0;
//    }
//    for (i = 0; i < WINDOW_SAMPLES && i < FFT_SIZE; i++) {
//        real[i] = E1[i];
//    }
//
//    /* naive DFT */
//    simpleDFT(real, real, imag, FFT_SIZE);
//
//    for (i = 0; i < FFT_SIZE; i++) {
//        double rv = real[i];
//        double iv = imag[i];
//        mag[i] = EnergyE_sqrt(rv * rv + iv * iv);
//    }
//
//    for (i = 0; i < FFT_SIZE; i++) E2[i] = 0.0;
//    for (i = 1; i < FFT_SIZE - 1; i++) {
//        double m_prev;
//        double m_next;
//        m_prev = (mag[i - 1] + mag[i]) / 2.0;
//        m_next = (mag[i + 1] + mag[i]) / 2.0;
//        E2[i] = mag[i] * mag[i] - m_prev * m_next;
//    }
//
//    /* Detection logic */
//    fs = 50.0;
//    detected = computeFCFR(E2, FFT_SIZE, fs, &max_freq, &fcfr_val);
//    result = 0;
//
//    if (detected && E_cooldown == 0) {
//        /* Peak detection on E1 to find event index */
//        maxE = -1e300;
//        maxIdx = 0;
//        for (i = 0; i < WINDOW_SAMPLES; i++) {
//            if (E1[i] > maxE) {
//                maxE = E1[i];
//                maxIdx = i;
//            }
//        }
//
//        /* amplitude (absolute acceleration amplitude; original used sqrt of Ax^2) */
//        amplitude = EnergyE_sqrt(E_window[maxIdx][0] * E_window[maxIdx][0]); /* this effectively abs(Ax) */
//        /* compute current global index approximated by windows processed count */
//        /* original: current_global_idx = E_windows_processed * WINDOW_SAMPLES + maxIdx; */
//        /* minSeparation in samples */
//        minSeparation = WINDOW_SAMPLES / 2;
//
//        if (amplitude >= 1.5 && ( (E_windows_processed * WINDOW_SAMPLES + maxIdx) - lastDetectedIdx > minSeparation)) {
//            result = 1;
//            lastDetectedIdx = E_windows_processed * WINDOW_SAMPLES + maxIdx;
//            E_cooldown = 100;
//        }
//    }
//
//    return result;
//}


#include "EnergyEnvelope.h"

/* if you rely on write() somewhere, keep extern; otherwise harmless */
extern long write(int fd, const void *buf, unsigned long count);

/* --------------------
   Math replacements (float instead of double)
   -------------------- */
static float EnergyE_fabs(float x) {
    return (x < 0.0f) ? -x : x;
}

static float EnergyE_sqrt(float x) {
    float r;
    int i;
    if (x <= 0.0f) return 0.0f;
    r = x * 0.5f;
    if (r <= 0.0f) r = 1.0f;
    for (i = 0; i < 20; ++i) {
        r = 0.5f * (r + x / r);
    }
    return r;
}

static float EnergyE_cos(float x) {
    float x2, x4, x6;
    while (x > ENERGYE_PI) x -= 2.0f * ENERGYE_PI;
    while (x < -ENERGYE_PI) x += 2.0f * ENERGYE_PI;
    x2 = x * x;
    x4 = x2 * x2;
    x6 = x4 * x2;
    return 1.0f - x2/2.0f + x4/24.0f - x6/720.0f;
}

static float EnergyE_sin(float x) {
    float x2, x3, x5, x7;
    while (x > ENERGYE_PI) x -= 2.0f * ENERGYE_PI;
    while (x < -ENERGYE_PI) x += 2.0f * ENERGYE_PI;
    x2 = x * x;
    x3 = x2 * x;
    x5 = x3 * x2;
    x7 = x5 * x2;
    return x - x3/6.0f + x5/120.0f - x7/5040.0f;
}

/* Simple string copy without stdlib */
void EnergyE_strcpy(char *dest, const char *src) {
    int i;
    if (!dest || !src) return;
    i = 0;
    while (src[i] != '\0' && i < 19) {
        dest[i] = src[i];
        i++;
    }
    dest[i] = '\0';
}

/* --------------------
   Helper buffers & state (Moved to static - Priority 2)
   -------------------- */

/* NOTE: changed to float (Priority 1) */
static float E_window[WINDOW_SAMPLES][6];  /* Reduced size (Priority 3) */
static int E_window_index = 0;
static int E_cooldown = 0;
static int E_windows_processed = 0;

/* FCFR stats (Welford) - changed to float */
static float fcfr_mean = 0.0f;
static float fcfr_M2 = 0.0f;
static int fcfr_count = 0;
static float adaptive_threshold = FCFR_THRESHOLD;

/* last detection tracking */
static int lastDetectedIdx = -100000;

/* --------------------
   Large arrays moved from stack to static (Priority 2)
   -------------------- */
static float x[WINDOW_SAMPLES];
static float E1[WINDOW_SAMPLES];
static float real[FFT_SIZE];
static float imag[FFT_SIZE];
static float mag[FFT_SIZE];
static float E2[FFT_SIZE];

/* --------------------
   DFT (naive) - C90 style (float version)
   -------------------- */
static void simpleDFT(const float *in, float *out_real, float *out_imag, int N) {
    int k, n;
    for (k = 0; k < N; k++) {
        float sum_r = 0.0f;
        float sum_i = 0.0f;
        for (n = 0; n < N; n++) {
            float angle;
            angle = -2.0f * ENERGYE_PI * (float)k * (float)n / (float)N;
            sum_r += in[n] * EnergyE_cos(angle);
            sum_i += in[n] * EnergyE_sin(angle);
        }
        out_real[k] = sum_r;
        out_imag[k] = sum_i;
    }
}

/* Teager Energy Operator (float version) */
static void computeTeagerEnergy(const float *x, float *E, int n) {
    int i;
    for (i = 1; i < n - 1; i++) {
        E[i] = x[i] * x[i] - x[i - 1] * x[i + 1];
    }
    E[0] = 0.0f;
    E[n - 1] = 0.0f;
}

/* computeFCFR - returns 1 if fcfr_val > adaptive_threshold */
static int computeFCFR(const float *E2, int n, float fs, float *max_freq, float *fcfr_val) {
    int i;
    int upper;
    float max_val;
    float mean;
    int max_idx;
    int count;
    float val;
    float delta;
    float fcfr_std;

    max_val = 0.0f;
    mean = 0.0f;
    max_idx = 0;
    count = 0;

    upper = n / 2;
    for (i = 1; i < upper; i++) {
        val = EnergyE_fabs(E2[i]);
        mean += val;
        count++;
        if (val > max_val) {
            max_val = val;
            max_idx = i;
        }
    }

    if (count == 0) {
        *fcfr_val = 0.0f;
        *max_freq = 0.0f;
        return 0;
    }

    mean /= (float)count;
    if (mean < EPS) mean = EPS;

    *fcfr_val = max_val / mean;
    *max_freq = (float)max_idx * fs / (float)n;

    /* update running stats (Welford) */
    fcfr_count++;
    delta = *fcfr_val - fcfr_mean;
    fcfr_mean += delta / (float)fcfr_count;
    fcfr_M2 += delta * (*fcfr_val - fcfr_mean);
    if (fcfr_count > 1) {
        fcfr_std = EnergyE_sqrt(fcfr_M2 / (float)(fcfr_count - 1));
    } else {
        fcfr_std = 0.0f;
    }
    adaptive_threshold = fcfr_mean + FCFR_WEIGHT * fcfr_std;

    return ((*fcfr_val) > adaptive_threshold) ? 1 : 0;
}

/* --------------------
   Public init
   -------------------- */
void EnergyEnvelope_init(void) {
    int i, j;
    E_window_index = 0;
    E_cooldown = 0;
    E_windows_processed = 0;
    fcfr_mean = 0.0f;
    fcfr_M2 = 0.0f;
    fcfr_count = 0;
    adaptive_threshold = FCFR_THRESHOLD;
    lastDetectedIdx = -100000;

    for (i = 0; i < WINDOW_SAMPLES; i++) {
        for (j = 0; j < 6; j++) {
            E_window[i][j] = 0.0f;
        }
    }
}

/* --------------------
   Main processing (updated with float and static arrays)
   -------------------- */
int EnergyEnvelope_processSample(float Ax, float Ay, float Az,
                                 float Gx, float Gy, float Gz) {
    int i, j;
    int result;
    int minSeparation;
    float amplitude;
    float accel_magnitude;
    float gyro_magnitude;
    float fs;
    float fcfr_val;
    float max_freq;
    int detected;
    float maxE;
    int maxIdx;

    /* Always insert new sample */
    if (E_window_index < WINDOW_SAMPLES) {
        /* append */
        E_window[E_window_index][0] = Ax;
        E_window[E_window_index][1] = Ay;
        E_window[E_window_index][2] = Az;
        E_window[E_window_index][3] = Gx;
        E_window[E_window_index][4] = Gy;
        E_window[E_window_index][5] = Gz;
        E_window_index++;
        if (E_window_index < WINDOW_SAMPLES) {
            return 0; /* still filling */
        }
        /* else continue to processing */
    } else {
        /* shift left by one and append */
        for (i = 1; i < WINDOW_SAMPLES; i++) {
            E_window[i - 1][0] = E_window[i][0];
            E_window[i - 1][1] = E_window[i][1];
            E_window[i - 1][2] = E_window[i][2];
            E_window[i - 1][3] = E_window[i][3];
            E_window[i - 1][4] = E_window[i][4];
            E_window[i - 1][5] = E_window[i][5];
        }
        E_window[WINDOW_SAMPLES - 1][0] = Ax;
        E_window[WINDOW_SAMPLES - 1][1] = Ay;
        E_window[WINDOW_SAMPLES - 1][2] = Az;
        E_window[WINDOW_SAMPLES - 1][3] = Gx;
        E_window[WINDOW_SAMPLES - 1][4] = Gy;
        E_window[WINDOW_SAMPLES - 1][5] = Gz;
        /* increment windows processed */
        E_windows_processed++;
    }

    /* Update cooldown */
    if (E_cooldown > 0) E_cooldown--;

    /* Build input vector x from Ax,Ay,Az,Gx,Gy,Gz magnitudes */
    for (i = 0; i < WINDOW_SAMPLES; i++) {
        accel_magnitude = EnergyE_sqrt(E_window[i][0] * E_window[i][0] +
                                       E_window[i][1] * E_window[i][1] +
                                       E_window[i][2] * E_window[i][2]);
        gyro_magnitude = EnergyE_sqrt(E_window[i][3] * E_window[i][3] +
                                      E_window[i][4] * E_window[i][4] +
                                      E_window[i][5] * E_window[i][5]);
        x[i] = accel_magnitude + gyro_magnitude;
    }

    /* Teager Energy Operator */
    computeTeagerEnergy(x, E1, WINDOW_SAMPLES);

    /* DFT prep */
    for (i = 0; i < FFT_SIZE; i++) {
        real[i] = 0.0f;
        imag[i] = 0.0f;
        mag[i] = 0.0f;
    }
    for (i = 0; i < WINDOW_SAMPLES && i < FFT_SIZE; i++) {
        real[i] = E1[i];
    }

    /* naive DFT */
    simpleDFT(real, real, imag, FFT_SIZE);

    for (i = 0; i < FFT_SIZE; i++) {
        float rv = real[i];
        float iv = imag[i];
        mag[i] = EnergyE_sqrt(rv * rv + iv * iv);
    }

    for (i = 0; i < FFT_SIZE; i++) E2[i] = 0.0f;
    for (i = 1; i < FFT_SIZE - 1; i++) {
        float m_prev;
        float m_next;
        m_prev = (mag[i - 1] + mag[i]) / 2.0f;
        m_next = (mag[i + 1] + mag[i]) / 2.0f;
        E2[i] = mag[i] * mag[i] - m_prev * m_next;
    }

    /* Detection logic */
    fs = 50.0f;
    detected = computeFCFR(E2, FFT_SIZE, fs, &max_freq, &fcfr_val);
    result = 0;

    if (detected && E_cooldown == 0) {
        /* Peak detection on E1 to find event index */
        maxE = -1e30f;
        maxIdx = 0;
        for (i = 0; i < WINDOW_SAMPLES; i++) {
            if (E1[i] > maxE) {
                maxE = E1[i];
                maxIdx = i;
            }
        }

        /* amplitude (absolute acceleration amplitude) */
        amplitude = EnergyE_sqrt(E_window[maxIdx][0] * E_window[maxIdx][0]);
        minSeparation = WINDOW_SAMPLES / 2;

        if (amplitude >= 1.5f &&
            ((E_windows_processed * WINDOW_SAMPLES + maxIdx) - lastDetectedIdx > minSeparation)) {
            result = 1;
            lastDetectedIdx = E_windows_processed * WINDOW_SAMPLES + maxIdx;
            E_cooldown = 100;
        }
    }

    return result;
}
