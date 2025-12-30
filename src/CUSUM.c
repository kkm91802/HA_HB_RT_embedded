//#include "CUSUM.h"
//
///* MATH.H ALTERNATIVE FUNCTIONS */
//static double CUSUM_fabs(double x) {
//    return (x < 0.0) ? -x : x;
//}
//
//static double CUSUM_fmax(double a, double b) {
//    return (a > b) ? a : b;
//}
//
///* CUSUM HELPER FUNCTIONS */
///* -------------------- Global state -------------------- */
//static CUSUMSample C_buffer[REF_SIZE + SMOOTH_WIN_CENTERED];
//static int C_buffer_index = 0;
//static int C_cooldown = 0;
//static int C_samples_processed = 0;
//
///* CUSUM algorithm state */
///* C90: no designated initializer; initialize nested array first, then idx and full */
//static Window ref_window = { {0.0}, 0, 0 };
//static double S_pos = 0.0;
//static double S_neg = 0.0;
//
///* Smoothing buffer */
//static double smoothed_values[SMOOTH_WIN_CENTERED];
//static int smooth_index = 0;
//
///* -------------------- Window management -------------------- */
//static void add_to_window(Window *w, double value) {
//    w->data[w->idx] = value;
//    w->idx = (w->idx + 1) % REF_SIZE;
//    if (w->idx == 0) {
//        w->full = 1;
//    }
//}
//
//static double window_mean(Window *w) {
//    int n;
//    int i;
//    double sum = 0.0;
//
//    n = w->full ? REF_SIZE : w->idx;
//    if (n == 0) return 0.0;
//
//    for (i = 0; i < n; i++) {
//        sum += w->data[i];
//    }
//    return sum / n;
//}
//
//static double window_var(Window *w, double mu) {
//    int n;
//    int i;
//    double s = 0.0;
//
//    n = w->full ? REF_SIZE : w->idx;
//    if (n < 2) return 0.0;
//
//    for (i = 0; i < n; i++) {
//        double diff = w->data[i] - mu;
//        s += diff * diff;
//    }
//    return s / (n - 1);
//}
//
//static double window_mean_exclude_recent(Window *w, int exclude) {
//    int n;
//    int limit;
//    int i;
//    double sum = 0.0;
//
//    n = w->full ? REF_SIZE : w->idx;
//    if (n == 0) return 0.0;
//
//    limit = n - exclude;
//    if (limit < 1) limit = n;
//
//    for (i = 0; i < limit; i++) {
//        sum += w->data[i];
//    }
//    return sum / limit;
//}
//
///* -------------------- Smoothing function -------------------- */
//static double compute_smoothed_value(void) {
//    /* Simple moving average of recent samples */
//    double sum = 0.0;
//    int count = 0;
//    int i;
//
//    for (i = 0; i < SMOOTH_WIN_CENTERED; i++) {
//        if (C_buffer_index > i) {
//            sum += C_buffer[C_buffer_index - 1 - i].Ax; /* Using Ax for smoothing */
//            count++;
//        }
//    }
//    return (count > 0) ? sum / count : 0.0;
//}
//
///* MAIN LOGIC */
//void CUSUM_init(void) {
//    int i, j;
//
//    C_buffer_index = 0;
//    C_cooldown = 0;
//    C_samples_processed = 0;
//    ref_window.idx = 0;
//    ref_window.full = 0;
//    S_pos = 0.0;
//    S_neg = 0.0;
//    smooth_index = 0;
//
//    /* Initialize buffers */
//    for (i = 0; i < REF_SIZE + SMOOTH_WIN_CENTERED; i++) {
//        C_buffer[i].Ax = 0.0;
//        C_buffer[i].Ay = 0.0;
//        C_buffer[i].Az = 0.0;
//        for (j = 0; j < 20; j++) {
//            /* C_buffer[i].time[j] = '\0'; optional - left commented */
//        }
//    }
//
//    for (i = 0; i < SMOOTH_WIN_CENTERED; i++) {
//        smoothed_values[i] = 0.0;
//    }
//}
//
//int CUSUM_processSample(double Ax, double Ay, double Az) {
//    /* Declare everything at the top for C90 */
//    double current_smoothed;
//    double mu0;
//    double sigma2;
//    double delta_in;
//    double mu1_in;
//    double s_in;
//    double delta_de;
//    double mu1_de;
//    double s_de;
//    double expected_per_sample;
//    double expected_S;
//    double alpha;
//    double best_val;
//    int detection_result;
//    int i, k;
//    int center_idx, lo, hi;
//
//    /* Buffer insertion */
//    if (C_buffer_index < REF_SIZE + SMOOTH_WIN_CENTERED) {
//        C_buffer[C_buffer_index].Ax = Ax;
//        C_buffer[C_buffer_index].Ay = Ay;
//        C_buffer[C_buffer_index].Az = Az;
//        C_buffer_index++;
//    } else {
//        for (i = 1; i < REF_SIZE + SMOOTH_WIN_CENTERED; i++) {
//            C_buffer[i - 1] = C_buffer[i];
//        }
//        C_buffer[REF_SIZE + SMOOTH_WIN_CENTERED - 1].Ax = Ax;
//        C_buffer[REF_SIZE + SMOOTH_WIN_CENTERED - 1].Ay = Ay;
//        C_buffer[REF_SIZE + SMOOTH_WIN_CENTERED - 1].Az = Az;
//    }
//
//    if (C_buffer_index < SMOOTH_WIN_CENTERED + REF_SIZE) {
//        return 0;
//    }
//
//    /* Processing */
//    current_smoothed = compute_smoothed_value();
//    add_to_window(&ref_window, current_smoothed);
//
//    if (C_cooldown > 0) C_cooldown--;
//
//    detection_result = 0;
//
//    if (ref_window.full && C_cooldown == 0) {
//        mu0 = window_mean_exclude_recent(&ref_window, EXCLUDE_RECENT);
//        sigma2 = window_var(&ref_window, mu0);
//        if (sigma2 < CUSUM_EPS) sigma2 = CUSUM_EPS;
//
//        delta_in = DELTA_MAG;
//        mu1_in = mu0 + delta_in;
//        s_in = (delta_in / sigma2) * (current_smoothed - (mu1_in + mu0) / 2.0);
//
//        delta_de = -DELTA_MAG;
//        mu1_de = mu0 + delta_de;
//        s_de = (delta_de / sigma2) * (current_smoothed - (mu1_de + mu0) / 2.0);
//
//        expected_per_sample = (delta_in * delta_in) / (2.0 * sigma2);
//        expected_S = expected_per_sample * (double)EVENT_LEN;
//        alpha = ALPHA_BASE + ALPHA_FACTOR * expected_S;
//
//        S_pos = CUSUM_fmax(0.0, S_pos + s_in);
//        S_neg = CUSUM_fmax(0.0, S_neg + s_de);
//
//        if (S_pos > alpha) {
//            center_idx = C_buffer_index - SMOOTH_WIN_CENTERED / 2;
//            lo = center_idx - PEAK_PRE;
//            hi = center_idx + PEAK_POST;
//            if (lo < 0) lo = 0;
//            if (hi >= C_buffer_index) hi = C_buffer_index - 1;
//
//            best_val = -1e9;
//            for (k = lo; k <= hi; k++) {
//                if (C_buffer[k].Ax > best_val) {
//                    best_val = C_buffer[k].Ax;
//                }
//            }
//
//            if (best_val > 0.5) {  /* Amplitude filter */
//                detection_result = 1;  /* Acceleration detected */
//                S_pos = 0.0;
//                S_neg = 0.0;
//                C_cooldown = RESET_DELAY;
//            }
//        }
//
//        if (S_neg > alpha) {
//            center_idx = C_buffer_index - SMOOTH_WIN_CENTERED / 2;
//            lo = center_idx - PEAK_PRE;
//            hi = center_idx + PEAK_POST;
//            if (lo < 0) lo = 0;
//            if (hi >= C_buffer_index) hi = C_buffer_index - 1;
//
//            best_val = 1e9;
//            for (k = lo; k <= hi; k++) {
//                if (C_buffer[k].Ax < best_val) {
//                    best_val = C_buffer[k].Ax;
//                }
//            }
//
//            if (best_val < -0.5) {  /* Amplitude filter */
//                detection_result = -1;  /* Braking detected */
//                S_pos = 0.0;
//                S_neg = 0.0;
//                C_cooldown = RESET_DELAY;
//            }
//        }
//    }
//
//    C_samples_processed++;
//    return detection_result;
//}

#include "CUSUM.h"

/* MATH.H ALTERNATIVE FUNCTIONS (float version) */
static float CUSUM_fabs(float x) {
    return (x < 0.0f) ? -x : x;
}

static float CUSUM_fmax(float a, float b) {
    return (a > b) ? a : b;
}

/* CUSUM HELPER FUNCTIONS */
/* -------------------- Global state -------------------- */
static CUSUMSample C_buffer[REF_SIZE + SMOOTH_WIN_CENTERED];
static int C_buffer_index = 0;
static int C_cooldown = 0;
static int C_samples_processed = 0;

/* CUSUM algorithm state */
static Window ref_window = { {0.0f}, 0, 0 };
static float S_pos = 0.0f;
static float S_neg = 0.0f;

/* Smoothing buffer */
static float smoothed_values[SMOOTH_WIN_CENTERED];
static int smooth_index = 0;

/* -------------------- Window management -------------------- */
static void add_to_window(Window *w, float value) {
    w->data[w->idx] = value;
    w->idx = (w->idx + 1) % REF_SIZE;
    if (w->idx == 0) {
        w->full = 1;
    }
}

static float window_mean(Window *w) {
    int n;
    int i;
    float sum = 0.0f;

    n = w->full ? REF_SIZE : w->idx;
    if (n == 0) return 0.0f;

    for (i = 0; i < n; i++) {
        sum += w->data[i];
    }
    return sum / (float)n;
}

static float window_var(Window *w, float mu) {
    int n;
    int i;
    float s = 0.0f;

    n = w->full ? REF_SIZE : w->idx;
    if (n < 2) return 0.0f;

    for (i = 0; i < n; i++) {
        float diff = w->data[i] - mu;
        s += diff * diff;
    }
    return s / (float)(n - 1);
}

static float window_mean_exclude_recent(Window *w, int exclude) {
    int n;
    int limit;
    int i;
    float sum = 0.0f;

    n = w->full ? REF_SIZE : w->idx;
    if (n == 0) return 0.0f;

    limit = n - exclude;
    if (limit < 1) limit = n;

    for (i = 0; i < limit; i++) {
        sum += w->data[i];
    }
    return sum / (float)limit;
}

/* -------------------- Smoothing function -------------------- */
static float compute_smoothed_value(void) {
    /* Simple moving average of recent samples */
    float sum = 0.0f;
    int count = 0;
    int i;

    for (i = 0; i < SMOOTH_WIN_CENTERED; i++) {
        if (C_buffer_index > i) {
            sum += C_buffer[C_buffer_index - 1 - i].Ax; /* Using Ax for smoothing */
            count++;
        }
    }
    return (count > 0) ? sum / (float)count : 0.0f;
}

/* MAIN LOGIC */
void CUSUM_init(void) {
    int i;

    C_buffer_index = 0;
    C_cooldown = 0;
    C_samples_processed = 0;
    ref_window.idx = 0;
    ref_window.full = 0;
    S_pos = 0.0f;
    S_neg = 0.0f;
    smooth_index = 0;

    /* Initialize buffers */
    for (i = 0; i < REF_SIZE + SMOOTH_WIN_CENTERED; i++) {
        C_buffer[i].Ax = 0.0f;
        C_buffer[i].Ay = 0.0f;
        C_buffer[i].Az = 0.0f;
    }

    for (i = 0; i < SMOOTH_WIN_CENTERED; i++) {
        smoothed_values[i] = 0.0f;
    }
}

int CUSUM_processSample(float Ax, float Ay, float Az) {
    /* Declare everything at the top for C90 */
    float current_smoothed;
    float mu0;
    float sigma2;
    float delta_in;
    float mu1_in;
    float s_in;
    float delta_de;
    float mu1_de;
    float s_de;
    float expected_per_sample;
    float expected_S;
    float alpha;
    float best_val;
    int detection_result;
    int i, k;
    int center_idx, lo, hi;

    /* Buffer insertion */
    if (C_buffer_index < REF_SIZE + SMOOTH_WIN_CENTERED) {
        C_buffer[C_buffer_index].Ax = Ax;
        C_buffer[C_buffer_index].Ay = Ay;
        C_buffer[C_buffer_index].Az = Az;
        C_buffer_index++;
    } else {
        for (i = 1; i < REF_SIZE + SMOOTH_WIN_CENTERED; i++) {
            C_buffer[i - 1] = C_buffer[i];
        }
        C_buffer[REF_SIZE + SMOOTH_WIN_CENTERED - 1].Ax = Ax;
        C_buffer[REF_SIZE + SMOOTH_WIN_CENTERED - 1].Ay = Ay;
        C_buffer[REF_SIZE + SMOOTH_WIN_CENTERED - 1].Az = Az;
    }

    if (C_buffer_index < SMOOTH_WIN_CENTERED + REF_SIZE) {
        return 0;
    }

    /* Processing */
    current_smoothed = compute_smoothed_value();
    add_to_window(&ref_window, current_smoothed);

    if (C_cooldown > 0) C_cooldown--;

    detection_result = 0;

    if (ref_window.full && C_cooldown == 0) {
        mu0 = window_mean_exclude_recent(&ref_window, EXCLUDE_RECENT);
        sigma2 = window_var(&ref_window, mu0);
        if (sigma2 < CUSUM_EPS) sigma2 = CUSUM_EPS;

        delta_in = DELTA_MAG;
        mu1_in = mu0 + delta_in;
        s_in = (delta_in / sigma2) * (current_smoothed - (mu1_in + mu0) / 2.0f);

        delta_de = -DELTA_MAG;
        mu1_de = mu0 + delta_de;
        s_de = (delta_de / sigma2) * (current_smoothed - (mu1_de + mu0) / 2.0f);

        expected_per_sample = (delta_in * delta_in) / (2.0f * sigma2);
        expected_S = expected_per_sample * (float)EVENT_LEN;
        alpha = ALPHA_BASE + ALPHA_FACTOR * expected_S;

        S_pos = CUSUM_fmax(0.0f, S_pos + s_in);
        S_neg = CUSUM_fmax(0.0f, S_neg + s_de);

        if (S_pos > alpha) {
            center_idx = C_buffer_index - SMOOTH_WIN_CENTERED / 2;
            lo = center_idx - PEAK_PRE;
            hi = center_idx + PEAK_POST;
            if (lo < 0) lo = 0;
            if (hi >= C_buffer_index) hi = C_buffer_index - 1;

            best_val = -1e9f;
            for (k = lo; k <= hi; k++) {
                if (C_buffer[k].Ax > best_val) {
                    best_val = C_buffer[k].Ax;
                }
            }

            if (best_val > 0.5f) {  /* Amplitude filter */
                detection_result = 1;  /* Acceleration detected */
                S_pos = 0.0f;
                S_neg = 0.0f;
                C_cooldown = RESET_DELAY;
            }
        }

        if (S_neg > alpha) {
            center_idx = C_buffer_index - SMOOTH_WIN_CENTERED / 2;
            lo = center_idx - PEAK_PRE;
            hi = center_idx + PEAK_POST;
            if (lo < 0) lo = 0;
            if (hi >= C_buffer_index) hi = C_buffer_index - 1;

            best_val = 1e9f;
            for (k = lo; k <= hi; k++) {
                if (C_buffer[k].Ax < best_val) {
                    best_val = C_buffer[k].Ax;
                }
            }

            if (best_val < -0.5f) {  /* Amplitude filter */
                detection_result = -1;  /* Braking detected */
                S_pos = 0.0f;
                S_neg = 0.0f;
                C_cooldown = RESET_DELAY;
            }
        }
    }

    C_samples_processed++;
    return detection_result;
}
