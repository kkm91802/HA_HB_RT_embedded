//#ifndef CUSUM_H
//#define CUSUM_H
//
//#define REF_SIZE 100
//#define SMOOTH_WIN_CENTERED 5
//#define EXCLUDE_RECENT 10
//#define DELTA_MAG 1.0
//#define CUSUM_EPS 1e-9
//#define EVENT_LEN 5
//#define ALPHA_BASE 5.0
//#define ALPHA_FACTOR 0.1
//#define PEAK_PRE 5
//#define PEAK_POST 5
//#define RESET_DELAY 50
//
//typedef struct {
//    char time[20];
//    double Ax, Ay, Az;
//    double Gx, Gy, Gz;
//} CUSUMSample;
//
//typedef struct {
//    double data[REF_SIZE];
//    int idx;
//    int full;
//} Window;
//
///* Real-time CUSUM functions */
//void CUSUM_init(void);
//int CUSUM_processSample(double Ax, double Ay, double Az);  /* Returns: 1=accel, -1=brake, 0=none */
//
//#endif

#ifndef CUSUM_H
#define CUSUM_H

/* Reduced sizes for memory optimization */
#define REF_SIZE 50                /* Reduced from 100 */
#define SMOOTH_WIN_CENTERED 3      /* Reduced from 5 */
#define EXCLUDE_RECENT 10
#define DELTA_MAG 1.0f             /* Changed to float */
#define CUSUM_EPS 1e-9f            /* Changed to float */
#define EVENT_LEN 5
#define ALPHA_BASE 5.0f            /* Changed to float */
#define ALPHA_FACTOR 0.1f          /* Changed to float */
#define PEAK_PRE 5
#define PEAK_POST 5
#define RESET_DELAY 50

/* Simplified structure - only store what we use */
typedef struct {
    float Ax, Ay, Az;  /* Only accelerometer data, no gyro, no time */
} CUSUMSample;

typedef struct {
    float data[REF_SIZE];  /* Changed from double to float */
    int idx;
    int full;
} Window;

/* Real-time CUSUM functions */
void CUSUM_init(void);
int CUSUM_processSample(float Ax, float Ay, float Az);  /* Changed to float */

#endif
