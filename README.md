# HA_HB_RT_Embedded
**Embedded Vehicle Harsh Event Detection – Application Version**

## Overview
This repository contains the **embedded application implementation** of a vehicle driver behavior detection system.  
The system is designed to run **directly on an embedded ECU / device**, processing **real-time IMU sensor data** to detect:

- Harsh Acceleration
- Harsh Braking  
- (Extendable to Rash Turn Left / Right using gyroscope data)

This project is the **real-time, device-integrated counterpart** of the simulation-based pipeline and removes all offline, plotting, and library-dependent components.

---

## Key Characteristics
- Fully **real-time execution**
- No CSV parsing or file-based input
- No plotting or visualization
- No external library dependencies
- Designed for **embedded debugging using breakpoints**
- Lightweight, user-defined signal processing functions

---

## Input Data
The system operates on **live sensor data** acquired directly from the IMU mounted on the vehicle.

### Sensor Inputs
- Accelerometer: `ax, ay, az`
- Gyroscope: `gx, gy, gz`

There is:
- ❌ No timestamp handling
- ❌ No file-based input
- ❌ No offline data replay

Sensor values are consumed **as they are received in real time** from the device drivers.

---

## Processing Pipeline
The embedded application follows a deterministic, end-to-end signal processing pipeline:

1. **Sensor Acquisition**
   - Raw accelerometer and gyroscope readings are obtained directly from the IMU interface

2. **Sensor Orientation Handling**
   - Aligns IMU axes to the vehicle coordinate frame
   - Handles sensor mounting orientation differences

3. **Signal Calibration**
   - Bias and scaling correction implemented using user-defined logic
   - Replaces calibration routines previously dependent on external libraries

4. **State Estimation**
   - InEKF-inspired estimation logic
   - Implemented using lightweight mathematical routines suitable for embedded targets

5. **Filtering**
   - Butterworth-style filtering implemented manually
   - Suppresses high-frequency noise while preserving motion dynamics

6. **Event Detection Algorithms**
   The following detection methods are implemented:
   - Adaptive Thresholding
   - Derivative-based Detection
   - Energy Envelope Analysis
   - CUSUM
   - Continuous Wavelet Transform (CWT – embedded-friendly approximation)

7. **Decision Polling**
   - Outputs of multiple detection modules are combined
   - Final event decision is produced to improve robustness and reduce false triggers

---

## Debugging & Validation
- No logging or visualization is performed
- Validation is done using:
  - Breakpoints
  - Watch variables
  - On-device debugging tools

This reflects real-world embedded system constraints where console output and plotting are not available.

---

## Project Structure
HA_HB_RT_Embedded/
│
├── src/ # Core embedded application logic
└── api/ # Hardware / sensor interface abstractions


Internal utility `_metadata` directories are intentionally excluded from version control.

---

## Platform Information
- **Target:** Embedded ECU / device firmware
- **Language:** C
- **Compiler:** GCC (embedded toolchain compatible)
- **Dependencies:** None (only compiler built-ins)
- **Execution Model:** Real-time, continuous operation

---

## Background
This project was developed as part of an **internship at Bosch Global Software Technologies Private Limited, Bangalore**, addressing a real-world problem of **driver behavior detection** in vehicles.

- Original threshold-based detection was replaced with adaptive detection techniques
- This repository represents the **application-level embedded realization** of that work
- All implementation work in this repository is authored by me

---

## Future Enhancements
- Production-grade real-time buffer management
- Further optimization for low-power ECUs
- Gyroscope-based rash turn detection (Left / Right)
- Platform-specific hardware acceleration

---

## Ownership & Rights
© Kamal Kishore Majhi  
All rights reserved.

GitHub: **kkm91802**

This repository contains **no confidential logs, datasets, or proprietary information** and is shared for educational and demonstrative purposes only.
