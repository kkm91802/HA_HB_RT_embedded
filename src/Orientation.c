//#include "Orientation.h"
//
//// Absolute value helper
//static inline double Orientation_fabs(double x)
//{
//    return (x < 0.0) ? -x : x;
//}
//
//void transform_orientation_int(
//        int ax_raw, int ay_raw, int az_raw,
//        int gx_raw, int gy_raw, int gz_raw,
//        double* outAx, double* outAy, double* outAz,
//        double* outGx, double* outGy, double* outGz)
//{
//    // Convert to double
//    double ax = (double)ax_raw;
//    double ay = (double)ay_raw;
//    double az = (double)az_raw;
//
//    double gx = (double)gx_raw;
//    double gy = (double)gy_raw;
//    double gz = (double)gz_raw;
//
//    // Determine dominant accelerometer axis
//    double absAx = Orientation_fabs(ax);
//    double absAy = Orientation_fabs(ay);
//    double absAz = Orientation_fabs(az);
//    int orientation = 0;
//
//    if (absAx >= absAy && absAx >= absAz) {
//        orientation = (ax > 0) ? 2 : 4;
//    }
//    else if (absAy >= absAx && absAy >= absAz) {
//        orientation = (ay > 0) ? 3 : 1;
//    }
//    else {
//        orientation = (az < 0) ? 6 : 5;
//    }
//
//    // Apply rotation logic
//    switch (orientation)
//    {
//        case 1:
//            ax = -az; ay =  ax; az = -ay;
//            gx = -gz; gy =  gx; gz = -gy;
//            break;
//
//        case 2:
//            ax = -az; ay =  ay; az =  ax;
//            gx = -gz; gy =  gy; gz =  gx;
//            break;
//
//        case 3:
//            ax = -az; ay = -ax; az =  ay;
//            gx = -gz; gy = -gx; gz =  gy;
//            break;
//
//        case 4:
//            ax = -az; ay = -ay; az = -ax;
//            gx = -gz; gy = -gy; gz = -gx;
//            break;
//
//        case 5:
//            // no change
//            break;
//
//        case 6:
//            ax = -ax; ay =  ay; az = -az;
//            gx = -gx; gy =  gy; gz = -gz;
//            break;
//
//        default:
//            break;
//    }
//
//    // Output
//    *outAx = ax;
//    *outAy = ay;
//    *outAz = az;
//
//    *outGx = gx;
//    *outGy = gy;
//    *outGz = gz;
//}

#include "Orientation.h"

// Absolute value helper
static inline float Orientation_fabs(double x)
{
    return (x < 0.0) ? -x : x;
}

void transform_orientation_int(
        int ax_raw, int ay_raw, int az_raw,
        int gx_raw, int gy_raw, int gz_raw,
        float* outAx, float* outAy, float* outAz,
        float* outGx, float* outGy, float* outGz)
{
    // Convert to float
    float ax = (float)ax_raw;
    float ay = (float)ay_raw;
    float az = (float)az_raw;

    float gx = (float)gx_raw;
    float gy = (float)gy_raw;
    float gz = (float)gz_raw;

    // Determine dominant accelerometer axis
    float absAx = Orientation_fabs(ax);
    float absAy = Orientation_fabs(ay);
    float absAz = Orientation_fabs(az);
    int orientation = 0;

    if (absAx >= absAy && absAx >= absAz) {
        orientation = (ax > 0) ? 2 : 4;
    }
    else if (absAy >= absAx && absAy >= absAz) {
        orientation = (ay > 0) ? 3 : 1;
    }
    else {
        orientation = (az < 0) ? 6 : 5;
    }

    // Apply rotation logic
    switch (orientation)
    {
        case 1:
            ax = -az; ay =  ax; az = -ay;
            gx = -gz; gy =  gx; gz = -gy;
            break;

        case 2:
            ax = -az; ay =  ay; az =  ax;
            gx = -gz; gy =  gy; gz =  gx;
            break;

        case 3:
            ax = -az; ay = -ax; az =  ay;
            gx = -gz; gy = -gx; gz =  gy;
            break;

        case 4:
            ax = -az; ay = -ay; az = -ax;
            gx = -gz; gy = -gy; gz = -gx;
            break;

        case 5:
            // no change
            break;

        case 6:
            ax = -ax; ay =  ay; az = -az;
            gx = -gx; gy =  gy; gz = -gz;
            break;

        default:
            break;
    }

    // Output
    *outAx = ax;
    *outAy = ay;
    *outAz = az;

    *outGx = gx;
    *outGy = gy;
    *outGz = gz;
}
