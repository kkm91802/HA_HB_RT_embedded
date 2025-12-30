#ifndef THRESHOLD_H
#define THRESHOLD_H

// Configuration constants
#define TIME_THRESHOLD 150
#define INIT_THRESHOLD_POS 0.0f
#define INIT_THRESHOLD_NEG 0.0f
#define WINDOW_SIZE 50  //100

// Adaptive Threshold structure
typedef struct {
    float posThreshold;
    float negThreshold;
    int timeThreshold;
    int lastSpikeIndex;
} AdaptiveThreshold;

typedef enum {
    THRESHOLD_FALSE = 0,
    THRESHOLD_TRUE = 1
} Threshold_Bool;

// Custom string and math functions
void Threshold_strcpy(char *dest, const char *src);
int Threshold_strcmp(const char *str1, const char *str2);
float Threshold_sqrtf(float x);
float Threshold_fabsf(float x);

// MAIN THRESHOLD API
void initThreshold(AdaptiveThreshold *th, float posInit, float negInit, int timeTh);

// NEW fixed function signatures
void updateAdaptiveThreshold(AdaptiveThreshold *th);

// Spike detection helpers (internal)
static Threshold_Bool detectPositiveSpike(int centerIdx, AdaptiveThreshold *th, int currentSampleIndex);
static Threshold_Bool detectNegativeSpike(int centerIdx, AdaptiveThreshold *th, int currentSampleIndex);

// Main real-time process function
Threshold_Bool processIMUSampleRealTime(
        AdaptiveThreshold *th,
        float ax, float ay, float az,
        float gx, float gy, float gz,
        int currentIndex,
        Threshold_Bool *positiveSpike,
        Threshold_Bool *negativeSpike);


#endif
