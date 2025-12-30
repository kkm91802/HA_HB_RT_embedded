//#ifndef ENERGY_ENVELOPE_H
//#define ENERGY_ENVELOPE_H
//
//#define WINDOW_SAMPLES 50 //100
//#define FFT_SIZE 64
//#define FCFR_THRESHOLD 0.0
//#define FCFR_WEIGHT 0.0
//#define EPS 1e-9
//#define ENERGYE_PI 3.14159265358979323846
//
///* Public functions */
//
///* Initializes the energy envelope processing system */
//void EnergyEnvelope_init(void);
//
///* Processes a new sample with accelerometer and gyroscope data */
//int EnergyEnvelope_processSample(double Ax, double Ay, double Az,double Gx,double Gy,double Gz);
//
///* A utility function to copy a string */
//void EnergyE_strcpy(char *dest, const char *src);
//
//#endif

#ifndef ENERGY_ENVELOPE_H
#define ENERGY_ENVELOPE_H

/* Reduced sizes for memory optimization */
#define WINDOW_SAMPLES 50      /* Reduced from 100 */
#define FFT_SIZE 32            /* Reduced from 64 */
#define FCFR_THRESHOLD 0.0f
#define FCFR_WEIGHT 0.0f
#define EPS 1e-9f
#define ENERGYE_PI 3.14159265358979323846f  /* Added 'f' suffix */

/* Public functions */
void EnergyEnvelope_init(void);

/* Changed all parameters from double to float */
int EnergyEnvelope_processSample(float Ax, float Ay, float Az,
                                 float Gx, float Gy, float Gz);

void EnergyE_strcpy(char *dest, const char *src);

#endif
