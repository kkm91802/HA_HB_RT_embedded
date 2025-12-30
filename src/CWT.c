//#include "CWT.h"
//
///* -------------------------
//   Small helpers (no math.h)
//   ------------------------- */
//
///* absolute value */
//static double CWT_fabs(double x) {
//    return (x < 0.0) ? -x : x;
//}
//
///* small-range exponential approx (kept from original) */
//static double CWT_exp(double x) {
//    double x2;
//    if (x > 5.0 || x < -5.0) return 0.0;
//    x2 = x * x;
//    return 1.0 - x2 + (x2 * x2) / 2.0 - (x2 * x2 * x2) / 6.0;
//}
//
///* small-angle cosine approximation */
//static double CWT_cos(double x) {
//    double x2, x4, x6;
//    while (x > CWT_PI) x -= 2.0 * CWT_PI;
//    while (x < -CWT_PI) x += 2.0 * CWT_PI;
//    x2 = x * x;
//    x4 = x2 * x2;
//    x6 = x4 * x2;
//    return 1.0 - x2 / 2.0 + x4 / 24.0 - x6 / 720.0;
//}
//
///* simple Newton-Raphson sqrt replacement for embedded use */
//static double CWT_sqrt(double x) {
//    double r;
//    int iter;
//    if (x <= 0.0) return 0.0;
//    /* initial guess */
//    r = x * 0.5;
//    if (r <= 0.0) r = 1.0;
//    /* iterate few times */
//    for (iter = 0; iter < 10; iter++) {
//        r = 0.5 * (r + x / r);
//    }
//    return r;
//}
//
///* -------------------------
//   Sorting / median helpers
//   ------------------------- */
//
///* simple in-place bubble-sort median (C90-safe) */
//static double CWT_median_double(double *arr, int n) {
//    int i, j;
//    double temp;
//    if (n <= 0) return 0.0;
//
//    for (i = 0; i < n - 1; i++) {
//        for (j = 0; j < n - i - 1; j++) {
//            if (arr[j] > arr[j + 1]) {
//                temp = arr[j];
//                arr[j] = arr[j + 1];
//                arr[j + 1] = temp;
//            }
//        }
//    }
//
//    if (n % 2) {
//        return arr[n / 2];
//    } else {
//        return 0.5 * (arr[n / 2 - 1] + arr[n / 2]);
//    }
//}
//
///* -------------------------
//   Static pool allocator (embedded)
//   ------------------------- */
//static void* CWT_malloc(int bytes) {
//    #define CWT_MEMORY_POOL_SIZE 10000           //16384--100000
//    static unsigned char memory_pool[CWT_MEMORY_POOL_SIZE];
//    static int memory_used = 0;
//    void *ptr;
//
//    if (bytes <= 0) return 0;
//    /* align to 8 bytes for safety */
//    bytes = (bytes + 7) & ~7;
//    if ((memory_used + bytes) > CWT_MEMORY_POOL_SIZE) return 0;
//    ptr = (void*)&memory_pool[memory_used];
//    memory_used += bytes;
//    return ptr;
//}
//
///* compute MAD. Uses static pool allocator for scratch; if pool unavailable, does destructive in-place fallback */
//static double CWT_compute_mad(double *x, int n) {
//    double *scratch;
//    double med;
//    double mad;
//    int i;
//
//    if (n <= 0) return 0.0;
//
//    /* Try allocate scratch from pool */
//    scratch = (double*)CWT_malloc(n * sizeof(double));
//    if (scratch) {
//        /* copy x to scratch for median calculation */
//        for (i = 0; i < n; i++) scratch[i] = x[i];
//        med = CWT_median_double(scratch, n);
//        for (i = 0; i < n; i++) scratch[i] = CWT_fabs(x[i] - med);
//        mad = CWT_median_double(scratch, n);
//        /* do not free pool memory */
//        return mad / MAD_SCALE_FACTOR;
//    }
//
//    /* Fallback: destructive on original array (embedded constraint) */
//    med = CWT_median_double(x, n);
//    for (i = 0; i < n; i++) {
//        x[i] = CWT_fabs(x[i] - med);
//    }
//    mad = CWT_median_double(x, n);
//    return mad / MAD_SCALE_FACTOR;
//}
//
///* -------------------------
//   Global state (original names kept)
//   ------------------------- */
//static CWTSample C_buffer[MAX_LINES];
//static int C_buffer_index = 0;
//static int C_cooldown = 0;
//static double C_fs = 50.0;
//
//static double scales[N_SCALES];
//static int wavelet_lengths[N_SCALES];
//static double *wavelet_kernels_storage[N_SCALES]; /* pointers into pool */
//
//static double S_pos = 0.0;
//static double S_neg = 0.0;
//
///* -------------------------
//   Morlet Wavelet maker (kept original style)
//   ------------------------- */
//static void CWT_make_morlet(double *w, int L, double scale_samples) {
//    int i;
//    int mid;
//    double sigma;
//    double k0;
//    double sumabs;
//    double t, val;
//
//    mid = L / 2;
//    sigma = scale_samples;
//    if (sigma < 1.0) sigma = 1.0;
//    k0 = 5.0;
//    sumabs = 0.0;
//
//    for (i = 0; i < L; i++) {
//        t = ((double)i - (double)mid) / sigma;
//        val = CWT_cos(2.0 * CWT_PI * k0 * t) * CWT_exp(-0.5 * t * t);
//        w[i] = val;
//        sumabs += CWT_fabs(val);
//    }
//    if (sumabs == 0.0) return;
//    for (i = 0; i < L; i++) {
//        w[i] /= sumabs;
//    }
//}
//
///* -------------------------
//   Convolution helper (your original implementation restored)
//   ------------------------- */
//static double CWT_convolve_at_point(const double *sig, int N, const double *ker, int L, int pos) {
//    int half = L / 2;
//    double sum = 0.0;
//    int k;
//    for (k = 0; k < L; k++) {
//        int j = pos - k + half;
//        if (j >= 0 && j < N) {
//            sum += sig[j] * ker[k];
//        }
//    }
//    return sum;
//}
//
///* -------------------------
//   Initialization
//   ------------------------- */
//void CWT_init(void) {
//    int s;
//    int i;
//    int L;
//    double min_s;
//    double max_s;
//    double *kern_ptr;
//
//    C_buffer_index = 0;
//    C_cooldown = 0;
//    C_fs = 50.0;
//    S_pos = 0.0;
//    S_neg = 0.0;
//
//    for (s = 0; s < N_SCALES; s++) {
//        wavelet_kernels_storage[s] = 0;
//    }
//
//    min_s = MIN_DUR_SEC * C_fs;
//    max_s = MAX_DUR_SEC * C_fs;
//
//    for (s = 0; s < N_SCALES; s++) {
//        scales[s] = min_s + (max_s - min_s) * ((double)s) / ((double)(N_SCALES - 1));
//        L = (int)(10.0 * scales[s]);
//        if (L < 21) L = 21;
//        if (L > MAX_WAVELET_LEN) L = MAX_WAVELET_LEN;
//        if ((L % 2) == 0) L++;
//
//        wavelet_lengths[s] = L;
//
//        kern_ptr = (double*)CWT_malloc(L * sizeof(double));
//        if (!kern_ptr) {
//            /* pool exhausted: leave kernel NULL (will be treated as zero kernel) */
//            wavelet_kernels_storage[s] = 0;
//        } else {
//            wavelet_kernels_storage[s] = kern_ptr;
//            CWT_make_morlet(wavelet_kernels_storage[s], L, scales[s]);
//        }
//    }
//
//    /* zero the circular buffer */
//    for (i = 0; i < MAX_LINES; i++) {
//        C_buffer[i].ax = 0.0;
//        C_buffer[i].ay = 0.0;
//        C_buffer[i].az = 0.0;
//        C_buffer[i].amp = 0.0;
//        C_buffer[i].t_seconds = 0.0;
//        /* time_str left untouched */
//    }
//}
//
///* -------------------------
//   Processing function (matches header)
//   ------------------------- */
//int CWT_processSample(double Ax,double Ay,double Az) {
//    /* All locals declared at top for C90 */
//    int detection_result;
////    double Ax, Ay, Az;
//    int min_samples_needed;
//    int start_idx;
//    int window_size;
//    int i;
//    int idx;
//    double amp;
//    /* local static scratch to avoid big stack usage */
//    static double amp_window_local[MAX_WAVELET_LEN * 2];
//    double *amp_window;
//    int scale_index;
//    double coeffs_arr[N_SCALES];
//    double sum_coeff;
//    double avg_coeff;
//    int recent_idx;
//    int recent_limit;
//    int found;
//
//    detection_result = 0;
////    if (!s) return 0;
////
////    Ax = s->ax;
////    Ay = s->ay;
////    Az = s->az;
//
//    /* append to buffer or shift if full (keeps original behavior) */
//    if (C_buffer_index < MAX_LINES) {
//        amp = CWT_sqrt(Ax * Ax + Ay * Ay + Az * Az);
//        C_buffer[C_buffer_index].ax = Ax;
//        C_buffer[C_buffer_index].ay = Ay;
//        C_buffer[C_buffer_index].az = Az;
//        C_buffer[C_buffer_index].amp = amp;
//        C_buffer_index++;
//    } else {
//        for (i = 1; i < MAX_LINES; i++) {
//            C_buffer[i - 1] = C_buffer[i];
//        }
//        amp = CWT_sqrt(Ax * Ax + Ay * Ay + Az * Az);
//        C_buffer[MAX_LINES - 1].ax = Ax;
//        C_buffer[MAX_LINES - 1].ay = Ay;
//        C_buffer[MAX_LINES - 1].az = Az;
//        C_buffer[MAX_LINES - 1].amp = amp;
//    }
//
//    /* cooldown update */
//    if (C_cooldown > 0) C_cooldown--;
//
//    /* minimum samples required */
//    min_samples_needed = MAX_WAVELET_LEN * 2;
//    if (C_buffer_index < min_samples_needed) return 0;
//
//    start_idx = C_buffer_index - min_samples_needed;
//    if (start_idx < 0) start_idx = 0;
//    window_size = C_buffer_index - start_idx;
//    if (window_size <= 0) return 0;
//
//    /* copy amplitudes into local scratch (bounded by MAX_WAVELET_LEN*2) */
//    for (i = 0; i < window_size && i < (MAX_WAVELET_LEN * 2); i++) {
//        amp_window_local[i] = C_buffer[start_idx + i].amp;
//    }
//    amp_window = amp_window_local;
//    idx = window_size - 1;
//
//    /* zero coeffs */
//    for (scale_index = 0; scale_index < N_SCALES; scale_index++) {
//        coeffs_arr[scale_index] = 0.0;
//    }
//
//    /* compute coefficients across scales using restored kernels */
//    for (scale_index = 0; scale_index < N_SCALES; scale_index++) {
//        if (wavelet_kernels_storage[scale_index]) {
//            coeffs_arr[scale_index] = CWT_convolve_at_point(amp_window, window_size, wavelet_kernels_storage[scale_index], wavelet_lengths[scale_index], idx);
//        } else {
//            coeffs_arr[scale_index] = 0.0;
//        }
//    }
//
//    /* average */
//    sum_coeff = 0.0;
//    for (scale_index = 0; scale_index < N_SCALES; scale_index++) {
//        sum_coeff += coeffs_arr[scale_index];
//    }
//    avg_coeff = sum_coeff / (double)N_SCALES;
//
//    /* positive spike detection (accel) */
//    if ((avg_coeff > 0.09) && (C_buffer[C_buffer_index - 1].amp > 0.24) && (C_cooldown == 0)) {
//        found = 0;
//        recent_limit = C_buffer_index - 10;
//        if (recent_limit < 0) recent_limit = 0;
//        for (recent_idx = recent_limit; recent_idx < C_buffer_index; recent_idx++) {
//            if (C_buffer[recent_idx].ax > 0.5) {
//                found = 1;
//                break;
//            }
//        }
//        if (found) {
//            C_cooldown = 100;
//            detection_result = 1;
//        }
//    }
//    /* negative spike detection (brake) */
//    else if ((avg_coeff < -0.09) && (C_buffer[C_buffer_index - 1].amp < -0.24) && (C_cooldown == 0)) {
//        found = 0;
//        recent_limit = C_buffer_index - 10;
//        if (recent_limit < 0) recent_limit = 0;
//        for (recent_idx = recent_limit; recent_idx < C_buffer_index; recent_idx++) {
//            if (C_buffer[recent_idx].ax < -0.5) {
//                found = 1;
//                break;
//            }
//        }
//        if (found) {
//            C_cooldown = 100;
//            detection_result = -1;
//        }
//    }
//
//    return detection_result;
//}

#include "CWT.h"

/* -------------------------
   Small helpers (float version)
   ------------------------- */

static float CWT_fabs(float x) {
    return (x < 0.0f) ? -x : x;
}

static float CWT_exp(float x) {
    float x2;
    if (x > 5.0f || x < -5.0f) return 0.0f;
    x2 = x * x;
    return 1.0f - x2 + (x2 * x2) / 2.0f - (x2 * x2 * x2) / 6.0f;
}

static float CWT_cos(float x) {
    float x2, x4, x6;
    while (x > CWT_PI) x -= 2.0f * CWT_PI;
    while (x < -CWT_PI) x += 2.0f * CWT_PI;
    x2 = x * x;
    x4 = x2 * x2;
    x6 = x4 * x2;
    return 1.0f - x2 / 2.0f + x4 / 24.0f - x6 / 720.0f;
}

static float CWT_sqrt(float x) {
    float r;
    int iter;
    if (x <= 0.0f) return 0.0f;
    r = x * 0.5f;
    if (r <= 0.0f) r = 1.0f;
    for (iter = 0; iter < 10; iter++) {
        r = 0.5f * (r + x / r);
    }
    return r;
}

/* -------------------------
   Sorting / median helpers (float version)
   ------------------------- */

static float CWT_median_float(float *arr, int n) {
    int i, j;
    float temp;
    if (n <= 0) return 0.0f;

    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
            }
        }
    }

    if (n % 2) {
        return arr[n / 2];
    } else {
        return 0.5f * (arr[n / 2 - 1] + arr[n / 2]);
    }
}

/* -------------------------
   Static pool allocator (reduced)
   ------------------------- */
static void* CWT_malloc(int bytes) {
    #define CWT_MEMORY_POOL_SIZE 4096  /* Reduced from 10000 */
    static unsigned char memory_pool[CWT_MEMORY_POOL_SIZE];
    static int memory_used = 0;
    void *ptr;

    if (bytes <= 0) return 0;
    bytes = (bytes + 7) & ~7;
    if ((memory_used + bytes) > CWT_MEMORY_POOL_SIZE) return 0;
    ptr = (void*)&memory_pool[memory_used];
    memory_used += bytes;
    return ptr;
}

static float CWT_compute_mad(float *x, int n) {
    float *scratch;
    float med;
    float mad;
    int i;

    if (n <= 0) return 0.0f;

    scratch = (float*)CWT_malloc(n * sizeof(float));
    if (scratch) {
        for (i = 0; i < n; i++) scratch[i] = x[i];
        med = CWT_median_float(scratch, n);
        for (i = 0; i < n; i++) scratch[i] = CWT_fabs(x[i] - med);
        mad = CWT_median_float(scratch, n);
        return mad / MAD_SCALE_FACTOR;
    }

    med = CWT_median_float(x, n);
    for (i = 0; i < n; i++) {
        x[i] = CWT_fabs(x[i] - med);
    }
    mad = CWT_median_float(x, n);
    return mad / MAD_SCALE_FACTOR;
}

/* -------------------------
   Global state (float version)
   ------------------------- */
static CWTSample C_buffer[MAX_LINES];  /* Now MAX_LINES = 50 */
static int C_buffer_index = 0;
static int C_cooldown = 0;
static float C_fs = 50.0f;  /* Changed to float */

static float scales[N_SCALES];  /* Changed to float */
static int wavelet_lengths[N_SCALES];
static float *wavelet_kernels_storage[N_SCALES];  /* Changed to float* */

static float S_pos = 0.0f;  /* Changed to float */
static float S_neg = 0.0f;  /* Changed to float */

/* -------------------------
   Morlet Wavelet maker (float version)
   ------------------------- */
static void CWT_make_morlet(float *w, int L, float scale_samples) {
    int i;
    int mid;
    float sigma;
    float k0;
    float sumabs;
    float t, val;

    mid = L / 2;
    sigma = scale_samples;
    if (sigma < 1.0f) sigma = 1.0f;
    k0 = 5.0f;
    sumabs = 0.0f;

    for (i = 0; i < L; i++) {
        t = ((float)i - (float)mid) / sigma;
        val = CWT_cos(2.0f * CWT_PI * k0 * t) * CWT_exp(-0.5f * t * t);
        w[i] = val;
        sumabs += CWT_fabs(val);
    }
    if (sumabs == 0.0f) return;
    for (i = 0; i < L; i++) {
        w[i] /= sumabs;
    }
}

/* -------------------------
   Convolution helper (float version)
   ------------------------- */
static float CWT_convolve_at_point(const float *sig, int N, const float *ker, int L, int pos) {
    int half = L / 2;
    float sum = 0.0f;
    int k;
    for (k = 0; k < L; k++) {
        int j = pos - k + half;
        if (j >= 0 && j < N) {
            sum += sig[j] * ker[k];
        }
    }
    return sum;
}

/* -------------------------
   Initialization
   ------------------------- */
void CWT_init(void) {
    int s;
    int i;
    int L;
    float min_s;
    float max_s;
    float *kern_ptr;

    C_buffer_index = 0;
    C_cooldown = 0;
    C_fs = 50.0f;
    S_pos = 0.0f;
    S_neg = 0.0f;

    for (s = 0; s < N_SCALES; s++) {
        wavelet_kernels_storage[s] = 0;
    }

    min_s = MIN_DUR_SEC * C_fs;
    max_s = MAX_DUR_SEC * C_fs;

    for (s = 0; s < N_SCALES; s++) {
        scales[s] = min_s + (max_s - min_s) * ((float)s) / ((float)(N_SCALES - 1));
        L = (int)(10.0f * scales[s]);
        if (L < 21) L = 21;
        if (L > MAX_WAVELET_LEN) L = MAX_WAVELET_LEN;
        if ((L % 2) == 0) L++;

        wavelet_lengths[s] = L;

        kern_ptr = (float*)CWT_malloc(L * sizeof(float));
        if (!kern_ptr) {
            wavelet_kernels_storage[s] = 0;
        } else {
            wavelet_kernels_storage[s] = kern_ptr;
            CWT_make_morlet(wavelet_kernels_storage[s], L, scales[s]);
        }
    }

    /* zero the circular buffer */
    for (i = 0; i < MAX_LINES; i++) {
        C_buffer[i].ax = 0.0f;
        C_buffer[i].ay = 0.0f;
        C_buffer[i].az = 0.0f;
        C_buffer[i].amp = 0.0f;
        C_buffer[i].t_seconds = 0.0f;
    }
}

/* -------------------------
   Processing function (float version)
   ------------------------- */
int CWT_processSample(float Ax, float Ay, float Az) {
    int detection_result;
    int min_samples_needed;
    int start_idx;
    int window_size;
    int i;
    int idx;
    float amp;
    static float amp_window_local[MAX_WAVELET_LEN * 2];
    float *amp_window;
    int scale_index;
    float coeffs_arr[N_SCALES];
    float sum_coeff;
    float avg_coeff;
    int recent_idx;
    int recent_limit;
    int found;

    detection_result = 0;

    /* append to buffer or shift if full */
    if (C_buffer_index < MAX_LINES) {
        amp = CWT_sqrt(Ax * Ax + Ay * Ay + Az * Az);
        C_buffer[C_buffer_index].ax = Ax;
        C_buffer[C_buffer_index].ay = Ay;
        C_buffer[C_buffer_index].az = Az;
        C_buffer[C_buffer_index].amp = amp;
        C_buffer_index++;
    } else {
        for (i = 1; i < MAX_LINES; i++) {
            C_buffer[i - 1] = C_buffer[i];
        }
        amp = CWT_sqrt(Ax * Ax + Ay * Ay + Az * Az);
        C_buffer[MAX_LINES - 1].ax = Ax;
        C_buffer[MAX_LINES - 1].ay = Ay;
        C_buffer[MAX_LINES - 1].az = Az;
        C_buffer[MAX_LINES - 1].amp = amp;
    }

    /* cooldown update */
    if (C_cooldown > 0) C_cooldown--;

    /* minimum samples required */
    min_samples_needed = MAX_WAVELET_LEN * 2;
    if (C_buffer_index < min_samples_needed) return 0;

    start_idx = C_buffer_index - min_samples_needed;
    if (start_idx < 0) start_idx = 0;
    window_size = C_buffer_index - start_idx;
    if (window_size <= 0) return 0;

    /* copy amplitudes into local scratch */
    for (i = 0; i < window_size && i < (MAX_WAVELET_LEN * 2); i++) {
        amp_window_local[i] = C_buffer[start_idx + i].amp;
    }
    amp_window = amp_window_local;
    idx = window_size - 1;

    /* zero coeffs */
    for (scale_index = 0; scale_index < N_SCALES; scale_index++) {
        coeffs_arr[scale_index] = 0.0f;
    }

    /* compute coefficients across scales */
    for (scale_index = 0; scale_index < N_SCALES; scale_index++) {
        if (wavelet_kernels_storage[scale_index]) {
            coeffs_arr[scale_index] = CWT_convolve_at_point(
                amp_window, window_size,
                wavelet_kernels_storage[scale_index],
                wavelet_lengths[scale_index],
                idx
            );
        } else {
            coeffs_arr[scale_index] = 0.0f;
        }
    }

    /* average */
    sum_coeff = 0.0f;
    for (scale_index = 0; scale_index < N_SCALES; scale_index++) {
        sum_coeff += coeffs_arr[scale_index];
    }
    avg_coeff = sum_coeff / (float)N_SCALES;

    /* positive spike detection (accel) */
    if ((avg_coeff > 0.09f) && (C_buffer[C_buffer_index - 1].amp > 0.24f) && (C_cooldown == 0)) {
        found = 0;
        recent_limit = C_buffer_index - 10;
        if (recent_limit < 0) recent_limit = 0;
        for (recent_idx = recent_limit; recent_idx < C_buffer_index; recent_idx++) {
            if (C_buffer[recent_idx].ax > 0.5f) {
                found = 1;
                break;
            }
        }
        if (found) {
            C_cooldown = 100;
            detection_result = 1;
        }
    }
    /* negative spike detection (brake) */
    else if ((avg_coeff < -0.09f) && (C_buffer[C_buffer_index - 1].amp < -0.24f) && (C_cooldown == 0)) {
        found = 0;
        recent_limit = C_buffer_index - 10;
        if (recent_limit < 0) recent_limit = 0;
        for (recent_idx = recent_limit; recent_idx < C_buffer_index; recent_idx++) {
            if (C_buffer[recent_idx].ax < -0.5f) {
                found = 1;
                break;
            }
        }
        if (found) {
            C_cooldown = 100;
            detection_result = -1;
        }
    }

    return detection_result;
}
