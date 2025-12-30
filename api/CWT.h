//#ifndef CWT_H
//#define CWT_H
//
//#define MAX_LINES 50          //100--2000000
//#define MAX_TIMESTR_LEN 32
//#define N_SCALES 5
//#define MAX_WAVELET_LEN 201
//#define MIN_DUR_SEC 0.1
//#define MAX_DUR_SEC 0.5
//#define THRESH_K 3.0
//#define MERGE_WINDOW_MS 500.0
//#define MAD_SCALE_FACTOR 0.6745
//#define CWT_PI 3.14159265358979323846
//
//typedef struct {
//    char time_str[MAX_TIMESTR_LEN];
//    double t_seconds;
//    double ax, ay, az;
//    double amp;
//} CWTSample;
//
///* Real-time CWT functions */
//void CWT_init(void);
//int CWT_processSample(double Ax,double Ay,double Az);  /* Returns: 1=accel, -1=brake, 0=none */
//double CWT_parse_time_to_seconds(const char *tstr);  /* ADD THIS LINE */
//
//#endif

#ifndef CWT_H
#define CWT_H

/* Reduced sizes */
#define MAX_LINES 50           /* Reduced from 100 */
#define N_SCALES 5
#define MAX_WAVELET_LEN 201
#define MIN_DUR_SEC 0.1f       /* Added 'f' */
#define MAX_DUR_SEC 0.5f       /* Added 'f' */
#define THRESH_K 3.0f          /* Added 'f' */
#define MERGE_WINDOW_MS 500.0f /* Added 'f' */
#define MAD_SCALE_FACTOR 0.6745f /* Added 'f' */
#define CWT_PI 3.14159265358979323846f /* Added 'f' */

/* Simplified structure - removed time_str */
typedef struct {
    float t_seconds;           /* Changed from double */
    float ax, ay, az;          /* Changed from double */
    float amp;                 /* Changed from double */
} CWTSample;

/* Real-time CWT functions */
void CWT_init(void);
int CWT_processSample(float Ax, float Ay, float Az);  /* Changed to float */

#endif
