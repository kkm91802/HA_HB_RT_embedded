#ifndef ORIENTATION_H
#define ORIENTATION_H

void transform_orientation_int(
        int ax_raw, int ay_raw, int az_raw,
        int gx_raw, int gy_raw, int gz_raw,
        float* outAx, float* outAy, float* outAz,
        float* outGx, float* outGy, float* outGz);

#endif

//#ifndef ORIENTATION_H
//#define ORIENTATION_H
//
//void transform_orientation_int(
//        int ax_raw, int ay_raw, int az_raw,
//        int gx_raw, int gy_raw, int gz_raw,
//        double* outAx, double* outAy, double* outAz,
//        double* outGx, double* outGy, double* outGz);
//
//#endif
