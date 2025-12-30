//#ifndef DERIVATIVE_DETECT_H
//#define DERIVATIVE_DETECT_H
//
///* -------- Configurable Parameters -------- */
//#ifndef DERIV_WINDOW_SIZE
//#define DERIV_WINDOW_SIZE 20      /* Sliding window size */
//#endif
//
//#ifndef DERIV_COOLDOWN_WINDOWS
//#define DERIV_COOLDOWN_WINDOWS 15 /* Minimum windows between events */
//#endif
//
//#ifndef DERIV_MIN_WINDOWS_FOR_DETECTION
//#define DERIV_MIN_WINDOWS_FOR_DETECTION 10  /* Wait this many windows before detecting */
//#endif
//
///* -------- Internal IMU sample struct -------- */
//typedef struct {
//    double Ax;
//    double Ay;
//    double Az;
//    double Gx;
//    double Gy;
//    double Gz;
//} DerivSample;
//
///* -------- External API -------- */
//
///* Call once at startup */
//void DerivativeD_init(void);
//
///*
//   Process ONE IMU sample in real-time.
//
//   Returns:
//       +1  = Harsh Acceleration
//       -1  = Harsh Braking
//        0  = No event
//*/
//int DerivativeD_processSample(double ax, double ay, double az,
//                              double gx, double gy, double gz);
//
//#endif


#ifndef DERIVATIVE_DETECT_H
#define DERIVATIVE_DETECT_H

/* -------- Configurable Parameters -------- */
#ifndef DERIV_WINDOW_SIZE
#define DERIV_WINDOW_SIZE 20      /* Sliding window size */
#endif

#ifndef DERIV_COOLDOWN_WINDOWS
#define DERIV_COOLDOWN_WINDOWS 15 /* Minimum windows between events */
#endif

#ifndef DERIV_MIN_WINDOWS_FOR_DETECTION
#define DERIV_MIN_WINDOWS_FOR_DETECTION 10  /* Wait this many windows before detecting */
#endif

/* -------- Change: Use float instead of double -------- */
typedef struct {
    float Ax;  // Changed from double
    float Ay;
    float Az;
    float Gx;
    float Gy;
    float Gz;
} DerivSample;

/* -------- External API -------- */
void DerivativeD_init(void);

/* Change: Use float parameters instead of double */
int DerivativeD_processSample(float ax, float ay, float az,
                              float gx, float gy, float gz);

#endif
