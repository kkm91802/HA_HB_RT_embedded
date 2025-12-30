#include "Threshold.h"

/* ---------------------------------- */
/* STRLIB ALTERNATIVES                */
/* ---------------------------------- */
void Threshold_strcpy(char *dest, const char *src)
{
    while (*src) {
        *dest++ = *src++;
    }
    *dest = '\0';
}

int Threshold_strcmp(const char *str1, const char *str2)
{
    while (*str1 && (*str1 == *str2)) {
        str1++;
        str2++;
    }
    return ((unsigned char)*str1 - (unsigned char)*str2);
}

/* ---------------------------------- */
/* MATH ALTERNATIVES                   */
/* ---------------------------------- */
float Threshold_sqrtf(float x)
{
    float g;
    int i;

    if (x <= 0.0f) return 0.0f;

    g = x * 0.5f;
    for (i = 0; i < 10; i++) {
        g = 0.5f * (g + x / g);
    }
    return g;
}

float Threshold_fabsf(float x)
{
    if (x < 0.0f) return -x;
    return x;
}

/* ---------------------------------- */
/* CIRCULAR BUFFER                     */
/* ---------------------------------- */
#define BUF_SIZE (WINDOW_SIZE * 2)

static float circBuf[BUF_SIZE];
static int writeIdx = 0;
static int totalSamples = 0;

/* C89 requires no inline */
static int circ(int idx)
{
    if (idx < 0) return idx + BUF_SIZE;
    if (idx >= BUF_SIZE) return idx - BUF_SIZE;
    return idx;
}

/* ---------------------------------- */
/* INIT THRESHOLD                     */
/* ---------------------------------- */
void initThreshold(AdaptiveThreshold *th, float posInit, float negInit, int timeTh)
{
    int i;

    th->posThreshold = posInit;
    th->negThreshold = negInit;
    th->timeThreshold = timeTh;
    th->lastSpikeIndex = -timeTh;

    writeIdx = 0;
    totalSamples = 0;

    for (i = 0; i < BUF_SIZE; i++) {
        circBuf[i] = 0.0f;
    }
}

/* ---------------------------------- */
/* UPDATE THRESHOLD                   */
/* ---------------------------------- */
void updateAdaptiveThreshold(AdaptiveThreshold *th)
{
    float minVal, maxVal, v;
    int i;

    if (totalSamples < WINDOW_SIZE)
        return;

    minVal = circBuf[circ(writeIdx)];
    maxVal = minVal;

    for (i = 1; i < WINDOW_SIZE; i++) {
        v = circBuf[circ(writeIdx - i)];
        if (v < minVal) minVal = v;
        if (v > maxVal) maxVal = v;
    }

    th->posThreshold = 0.95f * th->posThreshold + 0.05f * maxVal;
    th->negThreshold = 0.95f * th->negThreshold + 0.05f * minVal - 0.5f;
}

/* ---------------------------------- */
/* SPIKE DETECT HELPERS               */
/* (must be non-static in C-file only)*/
/* ---------------------------------- */
Threshold_Bool detectPositiveSpike(int centerIdx, AdaptiveThreshold *th, int currentSampleIndex)
{
    int p = circ(centerIdx - 1);
    int n = circ(centerIdx + 1);
    float c = circBuf[centerIdx];

    if (c > th->posThreshold &&
        circBuf[p] < c &&
        circBuf[n] < c)
    {
        if (currentSampleIndex - th->lastSpikeIndex > th->timeThreshold) {
            th->lastSpikeIndex = currentSampleIndex;
            return THRESHOLD_TRUE;
        }
    }
    return THRESHOLD_FALSE;
}

Threshold_Bool detectNegativeSpike(int centerIdx, AdaptiveThreshold *th, int currentSampleIndex)
{
    int p = circ(centerIdx - 1);
    int n = circ(centerIdx + 1);
    float c = circBuf[centerIdx];

    if (c < th->negThreshold &&
        circBuf[p] > c &&
        circBuf[n] > c)
    {
        if (currentSampleIndex - th->lastSpikeIndex > th->timeThreshold) {
            th->lastSpikeIndex = currentSampleIndex;
            return THRESHOLD_TRUE;
        }
    }
    return THRESHOLD_FALSE;
}

/* ---------------------------------- */
/* MAIN REAL-TIME FUNCTION            */
/* ---------------------------------- */
Threshold_Bool processIMUSampleRealTime(
        AdaptiveThreshold *th,
        float ax, float ay, float az,
        float gx, float gy, float gz,
        int currentIndex,
        Threshold_Bool *positiveSpike,
        Threshold_Bool *negativeSpike)
{
    float mag = ax;      /* You only use ax */

    circBuf[writeIdx] = mag;

    *positiveSpike = THRESHOLD_FALSE;
    *negativeSpike = THRESHOLD_FALSE;

    if (totalSamples >= 2) {
        int center = writeIdx;
        *positiveSpike = detectPositiveSpike(center, th, currentIndex);
        *negativeSpike = detectNegativeSpike(center, th, currentIndex);
    }

    writeIdx = circ(writeIdx + 1);
    totalSamples++;

    if ((totalSamples % WINDOW_SIZE) == 0)
        updateAdaptiveThreshold(th);

    return (Threshold_Bool)(*positiveSpike || *negativeSpike);
}
