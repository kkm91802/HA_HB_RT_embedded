//#ifndef BUTTERWORTH_H
//#define BUTTERWORTH_H
//
//typedef struct {
//    double b0, b1, b2;
//    double a1, a2;
//    double x1, x2; /* previous inputs */
//    double y1, y2; /* previous outputs */
//} Butterworth_Biquad;
//
//typedef struct {
//    Butterworth_Biquad sections[2];
//    int initialized;
//} Butterworth_Filter;
//
//typedef struct {
//    Butterworth_Filter ax_filter;
//    Butterworth_Filter ay_filter;
//    Butterworth_Filter az_filter;
//    Butterworth_Filter gx_filter;
//    Butterworth_Filter gy_filter;
//    Butterworth_Filter gz_filter;
//} Butterworth_IMU;
//
///* Public API */
//void Butterworth_Init(Butterworth_IMU* butterworth);
//
//void Butterworth_FilterIMU(
//    Butterworth_IMU* butterworth,
//    double ax, double ay, double az,
//    double gx, double gy, double gz,
//    double* out_ax, double* out_ay, double* out_az,
//    double* out_gx, double* out_gy, double* out_gz);
//
//#endif

#ifndef BUTTERWORTH_H
#define BUTTERWORTH_H

typedef struct {
    float b0, b1, b2;
    float a1, a2;
    float x1, x2; /* previous inputs */
    float y1, y2; /* previous outputs */
} Butterworth_Biquad;

typedef struct {
    Butterworth_Biquad sections[2];
    int initialized;
} Butterworth_Filter;

typedef struct {
    Butterworth_Filter ax_filter;
    Butterworth_Filter ay_filter;
    Butterworth_Filter az_filter;
    Butterworth_Filter gx_filter;
    Butterworth_Filter gy_filter;
    Butterworth_Filter gz_filter;
} Butterworth_IMU;

/* Public API */
void Butterworth_Init(Butterworth_IMU* butterworth);

void Butterworth_FilterIMU(
    Butterworth_IMU* butterworth,
    float ax, float ay, float az,
    float gx, float gy, float gz,
    float* out_ax, float* out_ay, float* out_az,
    float* out_gx, float* out_gy, float* out_gz);

#endif
