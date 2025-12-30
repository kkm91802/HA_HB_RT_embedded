//#include "DerivativeDetect.h"
//
///* -------------------- Custom boolean -------------------- */
//typedef enum {
//    DERIV_FALSE = 0,
//    DERIV_TRUE  = 1
//} DerivativeD_bool;
//
///* -------------------- Global state -------------------- */
//static DerivSample D_window[DERIV_WINDOW_SIZE];
//static int D_window_index = 0;
//static int D_cooldown = 0;
//static int D_windows_processed = 0;
//
///* -------------------- Math replacements -------------------- */
//static double DerivativeD_fabs(double x) {
//    return (x < 0.0) ? -x : x;
//}
//
///* C90-style sqrt replacement (Newton-Raphson) */
//static double DerivativeD_sqrt(double x) {
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
///* -------------------- Derivative window metrics -------------------- */
//
///* Compute derivative energy for current window */
//static double computeWindowDerivativeEnergy(const DerivSample *window, int window_size) {
//    double sum;
//    double prev;
//    int i;
//
//    if (!window || window_size <= 1) return 0.0;
//
//    sum = 0.0;
//    prev = window[0].Ax;
//
//    for (i = 1; i < window_size; ++i) {
//        double current = window[i].Ax;
//        double diff = current - prev;
//        sum += DerivativeD_fabs(diff);
//        prev = current;
//    }
//
//    return sum / (double)(window_size - 1);
//}
//
///* Compute maximum absolute Ax in window */
//static double computeWindowMaxAbs(const DerivSample *window, int window_size) {
//    double max_val;
//    int i;
//
//    if (!window || window_size <= 0) return 0.0;
//
//    max_val = DerivativeD_fabs(window[0].Ax);
//    for (i = 1; i < window_size; ++i) {
//        double val = DerivativeD_fabs(window[i].Ax);
//        if (val > max_val) {
//            max_val = val;
//        }
//    }
//    return max_val;
//}
//
///* -------------------- Energy statistics -------------------- */
//static double D_recent_energies[50];
//static int D_energy_index = 0;
//static int D_energy_count = 0;
//
///* Add one energy value into rolling history */
//static void addEnergyToHistory(double energy) {
//    D_recent_energies[D_energy_index] = energy;
//    D_energy_index = (D_energy_index + 1) % 50;
//    if (D_energy_count < 50) D_energy_count++;
//}
//
///* Compute mean and std of energy history */
//static void computeEnergyStats(double *mean, double *std) {
//    double mu;
//    double sigma;
//    int i;
//
//    if (D_energy_count == 0) {
//        *mean = 0.0;
//        *std = 0.0;
//        return;
//    }
//
//    mu = 0.0;
//    for (i = 0; i < D_energy_count; ++i) {
//        mu += D_recent_energies[i];
//    }
//    mu /= (double)D_energy_count;
//
//    sigma = 0.0;
//    for (i = 0; i < D_energy_count; ++i) {
//        double d = D_recent_energies[i] - mu;
//        sigma += d * d;
//    }
//    sigma = DerivativeD_sqrt(sigma / (double)D_energy_count);
//
//    *mean = mu;
//    *std = sigma;
//}
//
///* -------------------- Initialization -------------------- */
//
//void DerivativeD_init(void) {
//    int i;
//
//    D_window_index = 0;
//    D_cooldown = 0;
//    D_windows_processed = 0;
//    D_energy_index = 0;
//    D_energy_count = 0;
//
//    for (i = 0; i < DERIV_WINDOW_SIZE; ++i) {
//        D_window[i].Ax = 0.0;
//        D_window[i].Ay = 0.0;
//        D_window[i].Az = 0.0;
//        D_window[i].Gx = 0.0;
//        D_window[i].Gy = 0.0;
//        D_window[i].Gz = 0.0;
//    }
//}
//
///* -------------------- Real-time processing -------------------- */
//
///*
//    INPUT:  ax,ay,az,gx,gy,gz  (one sample per call)
//    OUTPUT:  1 = acceleration spike
//            -1 = braking spike
//             0 = none
//*/
//int DerivativeD_processSample(double ax, double ay, double az,
//                              double gx, double gy, double gz)
//{
//    /* All locals declared at top (C89/C90 requirement) */
//    double energy;
//    double maxabs;
//    int i;
//    double mu;
//    double sigma;
//    double pos_threshold;
//    double neg_threshold;
//    int detection;
//
//    /* Add sample into window */
//    if (D_window_index < DERIV_WINDOW_SIZE) {
//        D_window[D_window_index].Ax = ax;
//        D_window[D_window_index].Ay = ay;
//        D_window[D_window_index].Az = az;
//        D_window[D_window_index].Gx = gx;
//        D_window[D_window_index].Gy = gy;
//        D_window[D_window_index].Gz = gz;
//        D_window_index++;
//
//        if (D_window_index < DERIV_WINDOW_SIZE) {
//            return 0;   /* Not enough samples yet */
//        }
//    }
//
//    /* --- Window full → process --- */
//    energy = computeWindowDerivativeEnergy(D_window, DERIV_WINDOW_SIZE);
//    maxabs = computeWindowMaxAbs(D_window, DERIV_WINDOW_SIZE);
//
//    addEnergyToHistory(energy);
//    D_windows_processed++;
//
//    if (D_cooldown > 0)
//        D_cooldown--;
//
//    /* Not enough data for statistical thresholds yet */
//    if (D_windows_processed < DERIV_MIN_WINDOWS_FOR_DETECTION) {
//        /* slide window forward */
//        for (i = 1; i < DERIV_WINDOW_SIZE; ++i) {
//            D_window[i-1] = D_window[i];
//        }
//        D_window_index = DERIV_WINDOW_SIZE - 1;
//        return 0;
//    }
//
//    /* Compute adaptive statistical thresholds */
//    computeEnergyStats(&mu, &sigma);
//
//    pos_threshold = mu + 3.5 * sigma;
//    neg_threshold = -(mu + 3.0 * sigma);
//
//    detection = 0;
//
//    /* Check acceleration spike */
//    if (energy > pos_threshold && D_cooldown == 0) {
//        detection = 1;
//        D_cooldown = DERIV_COOLDOWN_WINDOWS;
//    }
//    /* Check braking spike */
//    else if (energy < neg_threshold && D_cooldown == 0) {
//        detection = -1;
//        D_cooldown = DERIV_COOLDOWN_WINDOWS;
//    }
//
//    /* Slide window */
//    for (i = 1; i < DERIV_WINDOW_SIZE; ++i) {
//        D_window[i-1] = D_window[i];
//    }
//    D_window_index = DERIV_WINDOW_SIZE - 1;
//
//    return detection;
//}


#include "DerivativeDetect.h"

/* -------------------- Custom boolean -------------------- */
typedef enum {
    DERIV_FALSE = 0,
    DERIV_TRUE  = 1
} DerivativeD_bool;

/* -------------------- Global state -------------------- */
static DerivSample D_window[DERIV_WINDOW_SIZE];
static int D_window_index = 0;
static int D_cooldown = 0;
static int D_windows_processed = 0;

/* -------------------- Math replacements -------------------- */
static float DerivativeD_fabs(float x) {
    return (x < 0.0f) ? -x : x;
}

static float DerivativeD_sqrt(float x) {
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

/* -------------------- Derivative window metrics -------------------- */
static float computeWindowDerivativeEnergy(const DerivSample *window, int window_size) {
    float sum;
    float prev;
    int i;

    if (!window || window_size <= 1) return 0.0f;

    sum = 0.0f;
    prev = window[0].Ax;

    for (i = 1; i < window_size; ++i) {
        float current = window[i].Ax;
        float diff = current - prev;
        sum += DerivativeD_fabs(diff);
        prev = current;
    }

    return sum / (float)(window_size - 1);
}

static float computeWindowMaxAbs(const DerivSample *window, int window_size) {
    float max_val;
    int i;

    if (!window || window_size <= 0) return 0.0f;

    max_val = DerivativeD_fabs(window[0].Ax);
    for (i = 1; i < window_size; ++i) {
        float val = DerivativeD_fabs(window[i].Ax);
        if (val > max_val) {
            max_val = val;
        }
    }
    return max_val;
}

/* -------------------- Energy statistics -------------------- */
static float D_recent_energies[50];
static int D_energy_index = 0;
static int D_energy_count = 0;

static void addEnergyToHistory(float energy) {
    D_recent_energies[D_energy_index] = energy;
    D_energy_index = (D_energy_index + 1) % 50;
    if (D_energy_count < 50) D_energy_count++;
}

static void computeEnergyStats(float *mean, float *std) {
    float mu;
    float sigma;
    int i;

    if (D_energy_count == 0) {
        *mean = 0.0f;
        *std = 0.0f;
        return;
    }

    mu = 0.0f;
    for (i = 0; i < D_energy_count; ++i) {
        mu += D_recent_energies[i];
    }
    mu /= (float)D_energy_count;

    sigma = 0.0f;
    for (i = 0; i < D_energy_count; ++i) {
        float d = D_recent_energies[i] - mu;
        sigma += d * d;
    }
    sigma = DerivativeD_sqrt(sigma / (float)D_energy_count);

    *mean = mu;
    *std = sigma;
}

/* -------------------- Initialization -------------------- */
void DerivativeD_init(void) {
    int i;

    D_window_index = 0;
    D_cooldown = 0;
    D_windows_processed = 0;
    D_energy_index = 0;
    D_energy_count = 0;

    for (i = 0; i < DERIV_WINDOW_SIZE; ++i) {
        D_window[i].Ax = 0.0f;
        D_window[i].Ay = 0.0f;
        D_window[i].Az = 0.0f;
        D_window[i].Gx = 0.0f;
        D_window[i].Gy = 0.0f;
        D_window[i].Gz = 0.0f;
    }
}

/* -------------------- Real-time processing -------------------- */
int DerivativeD_processSample(float ax, float ay, float az,
                              float gx, float gy, float gz)
{
    float energy;
    float maxabs;
    int i;
    float mu;
    float sigma;
    float pos_threshold;
    float neg_threshold;
    int detection;

    /* Add sample into window */
    if (D_window_index < DERIV_WINDOW_SIZE) {
        D_window[D_window_index].Ax = ax;
        D_window[D_window_index].Ay = ay;
        D_window[D_window_index].Az = az;
        D_window[D_window_index].Gx = gx;
        D_window[D_window_index].Gy = gy;
        D_window[D_window_index].Gz = gz;
        D_window_index++;

        if (D_window_index < DERIV_WINDOW_SIZE) {
            return 0;   /* Not enough samples yet */
        }
    }

    /* --- Window full- process --- */
    energy = computeWindowDerivativeEnergy(D_window, DERIV_WINDOW_SIZE);
    maxabs = computeWindowMaxAbs(D_window, DERIV_WINDOW_SIZE);

    addEnergyToHistory(energy);
    D_windows_processed++;

    if (D_cooldown > 0)
        D_cooldown--;

    if (D_windows_processed < DERIV_MIN_WINDOWS_FOR_DETECTION) {
        for (i = 1; i < DERIV_WINDOW_SIZE; ++i) {
            D_window[i-1] = D_window[i];
        }
        D_window_index = DERIV_WINDOW_SIZE - 1;
        return 0;
    }

    computeEnergyStats(&mu, &sigma);

    pos_threshold = mu + 3.5f * sigma;
    neg_threshold = -(mu + 3.0f * sigma);

    detection = 0;

    if (energy > pos_threshold && D_cooldown == 0) {
        detection = 1;
        D_cooldown = DERIV_COOLDOWN_WINDOWS;
    }
    else if (energy < neg_threshold && D_cooldown == 0) {
        detection = -1;
        D_cooldown = DERIV_COOLDOWN_WINDOWS;
    }

    for (i = 1; i < DERIV_WINDOW_SIZE; ++i) {
        D_window[i-1] = D_window[i];
    }
    D_window_index = DERIV_WINDOW_SIZE - 1;

    return detection;
}
