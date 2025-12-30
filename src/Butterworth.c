//#include "Butterworth.h"
//
//int i;
//double out;
//
///* COEFFICENTS FOR FS=50Hz, FC=10Hz, ORDER=4
//   Two biquad sections; each row = {b0, b1, b2, a1, a2} */
//static const double BUTTERWORTH_COEFFICIENTS[2][5] = {
//    {0.046582, 0.093164, 0.046582, -1.454243, 0.574062},
//    {0.046582, 0.093164, 0.046582, -1.012899, 0.309597}
//};
//
///* Initialize a single biquad state */
//static void Butterworth_InitBiquad(Butterworth_Biquad* bq, const double coeffs[5]) {
//    bq->b0 = coeffs[0];
//    bq->b1 = coeffs[1];
//    bq->b2 = coeffs[2];
//    bq->a1 = coeffs[3];
//    bq->a2 = coeffs[4];
//
//    /* Direct Form I state: previous inputs x[n-1], x[n-2] and previous outputs y[n-1], y[n-2] */
//    bq->x1 = 0.0;
//    bq->x2 = 0.0;
//    bq->y1 = 0.0;
//    bq->y2 = 0.0;
//}
//
///* Process one biquad using Direct Form I:
//   y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2] */
//static double Butterworth_ProcessBiquad(Butterworth_Biquad* bq, double x) {
//    double y = bq->b0 * x
//             + bq->b1 * bq->x1
//             + bq->b2 * bq->x2
//             - bq->a1 * bq->y1
//             - bq->a2 * bq->y2;
//
//    /* shift states */
//    bq->x2 = bq->x1;
//    bq->x1 = x;
//
//    bq->y2 = bq->y1;
//    bq->y1 = y;
//
//    return y;
//}
//
///* Initialize the filter structure */
//static void Butterworth_InitFilter(Butterworth_Filter* filter) {
//    if (!filter) return;
//    for (i = 0 ; i < 2 ; ++i) {
//        Butterworth_InitBiquad(&filter->sections[i], BUTTERWORTH_COEFFICIENTS[i]);
//    }
//
//    filter->initialized = 1;
//}
//
///* Process one filter (two biquad sections in series) */
//static double Butterworth_ProcessFilter(Butterworth_Filter* filter, double input) {
//
//    if (!filter || !filter->initialized) return input;
//    out = input;
//    out = Butterworth_ProcessBiquad(&filter->sections[0], out);
//    out = Butterworth_ProcessBiquad(&filter->sections[1], out);
//    return out;
//}
//
///* Public API: initialize whole IMU filter */
//void Butterworth_Init(Butterworth_IMU* butterworth) {
//    if (!butterworth) return;
//    Butterworth_InitFilter(&butterworth->ax_filter);
//    Butterworth_InitFilter(&butterworth->ay_filter);
//    Butterworth_InitFilter(&butterworth->az_filter);
//    Butterworth_InitFilter(&butterworth->gx_filter);
//    Butterworth_InitFilter(&butterworth->gy_filter);
//    Butterworth_InitFilter(&butterworth->gz_filter);
//}
//
///* Public API: filter scalar IMU sample (no arrays) */
//void Butterworth_FilterIMU(
//    Butterworth_IMU* butterworth,
//    double ax, double ay, double az,
//    double gx, double gy, double gz,
//    double* out_ax, double* out_ay, double* out_az,
//    double* out_gx, double* out_gy, double* out_gz)
//{
//    if (!out_ax || !out_ay || !out_az || !out_gx || !out_gy || !out_gz) return;
//
//    if (!butterworth) {
//        *out_ax = ax; *out_ay = ay; *out_az = az;
//        *out_gx = gx; *out_gy = gy; *out_gz = gz;
//        return;
//    }
//
//    *out_ax = Butterworth_ProcessFilter(&butterworth->ax_filter, ax);
//    *out_ay = Butterworth_ProcessFilter(&butterworth->ay_filter, ay);
//    *out_az = Butterworth_ProcessFilter(&butterworth->az_filter, az);
//
//    *out_gx = Butterworth_ProcessFilter(&butterworth->gx_filter, gx);
//    *out_gy = Butterworth_ProcessFilter(&butterworth->gy_filter, gy);
//    *out_gz = Butterworth_ProcessFilter(&butterworth->gz_filter, gz);
//}

#include "Butterworth.h"

 int i;
 float out;

/* COEFFICENTS FOR FS=50Hz, FC=10Hz, ORDER=4
   Two biquad sections; each row = {b0, b1, b2, a1, a2} */
static const float BUTTERWORTH_COEFFICIENTS[2][5] = {
    {0.046582, 0.093164, 0.046582, -1.454243, 0.574062},
    {0.046582, 0.093164, 0.046582, -1.012899, 0.309597}
};

/* Initialize a single biquad state */
static void Butterworth_InitBiquad(Butterworth_Biquad* bq, const float coeffs[5]) {
    bq->b0 = coeffs[0];
    bq->b1 = coeffs[1];
    bq->b2 = coeffs[2];
    bq->a1 = coeffs[3];
    bq->a2 = coeffs[4];

    /* Direct Form I state: previous inputs x[n-1], x[n-2] and previous outputs y[n-1], y[n-2] */
    bq->x1 = 0.0;
    bq->x2 = 0.0;
    bq->y1 = 0.0;
    bq->y2 = 0.0;
}

/* Process one biquad using Direct Form I:
   y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2] */
static float Butterworth_ProcessBiquad(Butterworth_Biquad* bq, float x) {
    float y = bq->b0 * x
             + bq->b1 * bq->x1
             + bq->b2 * bq->x2
             - bq->a1 * bq->y1
             - bq->a2 * bq->y2;

    /* shift states */
    bq->x2 = bq->x1;
    bq->x1 = x;

    bq->y2 = bq->y1;
    bq->y1 = y;

    return y;
}

/* Initialize the filter structure */
static void Butterworth_InitFilter(Butterworth_Filter* filter) {
    if (!filter) return;

    for (i = 0 ; i < 2 ; ++i) {
        Butterworth_InitBiquad(&filter->sections[i], BUTTERWORTH_COEFFICIENTS[i]);
    }

    filter->initialized = 1;
}

/* Process one filter (two biquad sections in series) */
static float Butterworth_ProcessFilter(Butterworth_Filter* filter, float input) {

    if (!filter || !filter->initialized) return input;

    out = input;
    out = Butterworth_ProcessBiquad(&filter->sections[0], out);
    out = Butterworth_ProcessBiquad(&filter->sections[1], out);
    return out;
}

/* Public API: initialize whole IMU filter */
void Butterworth_Init(Butterworth_IMU* butterworth) {
    if (!butterworth) return;
    Butterworth_InitFilter(&butterworth->ax_filter);
    Butterworth_InitFilter(&butterworth->ay_filter);
    Butterworth_InitFilter(&butterworth->az_filter);
    Butterworth_InitFilter(&butterworth->gx_filter);
    Butterworth_InitFilter(&butterworth->gy_filter);
    Butterworth_InitFilter(&butterworth->gz_filter);
}

/* Public API: filter scalar IMU sample (no arrays) */
void Butterworth_FilterIMU(
    Butterworth_IMU* butterworth,
    float ax, float ay, float az,
    float gx, float gy, float gz,
    float* out_ax, float* out_ay, float* out_az,
    float* out_gx, float* out_gy, float* out_gz)
{
    if (!out_ax || !out_ay || !out_az || !out_gx || !out_gy || !out_gz) return;

    if (!butterworth) {
        *out_ax = ax; *out_ay = ay; *out_az = az;
        *out_gx = gx; *out_gy = gy; *out_gz = gz;
        return;
    }

    *out_ax = Butterworth_ProcessFilter(&butterworth->ax_filter, ax);
    *out_ay = Butterworth_ProcessFilter(&butterworth->ay_filter, ay);
    *out_az = Butterworth_ProcessFilter(&butterworth->az_filter, az);

    *out_gx = Butterworth_ProcessFilter(&butterworth->gx_filter, gx);
    *out_gy = Butterworth_ProcessFilter(&butterworth->gy_filter, gy);
    *out_gz = Butterworth_ProcessFilter(&butterworth->gz_filter, gz);
}
