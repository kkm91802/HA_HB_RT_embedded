//#ifndef INEKF_H
//#define INEKF_H
//
//// Constants
//#define NSTATE 15
//#define GGRAV 9.81
//#define DEG2RAD 0.017453292519943295
//#define RAD2DEG 57.29577951308232
//
//// Configuration structure
//typedef struct {
//    double gyro_noise_sigma;
//    double accel_noise_sigma;
//    double gyro_bias_rw;
//    double accel_bias_rw;
//    double zupt_velocity_noise;
//    double zupt_accel_threshold;
//    double zupt_gyro_threshold_deg;
//} InEKFConfig;
//
//// State structure
//typedef struct {
//    double R[9];        // Orientation matrix (3x3)
//    double v[3];        // Velocity vector
//    double p[3];        // Position vector
//    double bw[3];       // Gyroscope bias vector
//    double ba[3];       // Accelerometer bias vector
//    double P[NSTATE*NSTATE]; // Covariance matrix (15x15)
//} InEKFState;
//
//// Main filter structure
//typedef struct {
//    InEKFState state;
//    InEKFConfig config;
//    double Qc[NSTATE*NSTATE]; // Continuous process noise
//    double R_zupt[9];         // ZUPT measurement noise (3x3)
//} InEKF;
//
//
//void inekf_init(InEKF *filter, const InEKFConfig *config);
//void inekf_init_default(InEKF *filter);
//void inekf_predict(InEKF *filter, const double omega_meas[3], const double acc_meas[3], double dt);
//void inekf_zupt_update(InEKF *filter, const double velocity_meas[3]);
//int inekf_detect_zupt(const InEKF *filter, const double acc_meas[3], const double gyro_meas_deg[3]);
//void inekf_get_bias_corrected_imu(const InEKF *filter, const double acc_meas[3], const double gyro_meas_deg[3],
//                                  double acc_corrected[3], double gyro_corrected_deg[3]);
//
//void inekf_get_rotation(const InEKF *filter, double R[9]);
//void inekf_get_velocity(const InEKF *filter, double v[3]);
//void inekf_get_position(const InEKF *filter, double p[3]);
//void inekf_get_biases(const InEKF *filter, double bw[3], double ba[3]);
//
//void InEKF_ProcessIMU(
//    InEKF *filter,
//    double ax, double ay, double az,
//    double gx, double gy, double gz,
//    double dt,
//    double *out_ax, double *out_ay, double *out_az,
//    double *out_gx, double *out_gy, double *out_gz);
//
////void InEKF_ProcessIMUArray(InEKF *filter, const double input_array[7], double output_array[7], double dt);
////
////void InEKF_ProcessIMUArrayWithTimeString(InEKF *filter, const char* time_str,
////                                        double ax, double ay, double az,
////                                        double gx, double gy, double gz,
////                                        char* output_time_str,
////                                        double* output_ax, double* output_ay, double* output_az,
////                                        double* output_gx, double* output_gy, double* output_gz,
////                                        double dt);
//
//
//#endif


#ifndef INEKF_H
#define INEKF_H

// Constants
#define NSTATE 15
#define GGRAV 9.81f
#define DEG2RAD 0.017453292519943295f
#define RAD2DEG 57.29577951308232f

// Configuration structure
typedef struct {
    float gyro_noise_sigma;
    float accel_noise_sigma;
    float gyro_bias_rw;
    float accel_bias_rw;
    float zupt_velocity_noise;
    float zupt_accel_threshold;
    float zupt_gyro_threshold_deg;
} InEKFConfig;

// State structure
typedef struct {
    float R[9];        // Orientation matrix (3x3)
    float v[3];        // Velocity vector
    float p[3];        // Position vector
    float bw[3];       // Gyroscope bias vector
    float ba[3];       // Accelerometer bias vector
    float P[NSTATE*NSTATE]; // Covariance matrix (15x15)
} InEKFState;

// Main filter structure
typedef struct {
    InEKFState state;
    InEKFConfig config;
    float Qc[NSTATE*NSTATE]; // Continuous process noise
    float R_zupt[9];         // ZUPT measurement noise (3x3)
} InEKF;


void inekf_init(InEKF *filter, const InEKFConfig *config);
void inekf_init_default(InEKF *filter);
void inekf_predict(InEKF *filter, const float omega_meas[3], const float acc_meas[3], float dt);
void inekf_zupt_update(InEKF *filter, const float velocity_meas[3]);
int inekf_detect_zupt(const InEKF *filter, const float acc_meas[3], const float gyro_meas_deg[3]);
void inekf_get_bias_corrected_imu(const InEKF *filter, const float acc_meas[3], const float gyro_meas_deg[3],
                                  float acc_corrected[3], float gyro_corrected_deg[3]);

void inekf_get_rotation(const InEKF *filter, float R[9]);
void inekf_get_velocity(const InEKF *filter, float v[3]);
void inekf_get_position(const InEKF *filter, float p[3]);
void inekf_get_biases(const InEKF *filter, float bw[3], float ba[3]);

void InEKF_ProcessIMU(
    InEKF *filter,
    float ax, float ay, float az,
    float gx, float gy, float gz,
    float dt,
    float *out_ax, float *out_ay, float *out_az,
    float *out_gx, float *out_gy, float *out_gz);

#endif
