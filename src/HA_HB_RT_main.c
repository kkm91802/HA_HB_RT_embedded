//#include "Orientation.h"
//#include "InEKF.h"
//#include "Butterworth.h"
//#include "Threshold.h"
//#include "DerivativeDetect.h"
//#include "EnergyEnvelope.h"
//#include "CUSUM.h"
//#include "CWT.h"
//#include "HA_HB_RT_main.h"
////#include "PoS_Appl.h"
//
//Butterworth_IMU butterworth_filter;
//InEKF inekf_filter;
//AdaptiveThreshold th;
//
//const double dt = 0.02; // 50 Hz
//int sampleIndex = 0;
//
//int harshaccl_counter = 0;
//int harshbrak_counter = 0;
//
////IMU Data
//void Pos_AccDataRead(int *Ax,int *Ay, int *Az);
//void Pos_GyroDataRead(int *Gx, int *Gy,int *Gz);
//
//   int GyroX_s32;
//   int GyroY_s32;
//   int GyroZ_s32;
//
//   int AccX_s16;
//   int AccY_s16;
//   int AccZ_s16;
//
////Orientation
//double ax, ay, az, gx, gy, gz;
//
////Butterworth LPF
//double fax, fay, faz, fgx, fgy, fgz;
//
////InEKF
//double ekf_ax, ekf_ay, ekf_az, ekf_gx, ekf_gy, ekf_gz;
//
////Threshold
//Threshold_Bool positiveSpike = THRESHOLD_FALSE;
//Threshold_Bool negativeSpike = THRESHOLD_FALSE;
//
//Threshold_Bool detected;
//
////Derivative detect
//int deriv_flag;
//
////Energy envelope
//int energy_detected;
//
////CUSUM
//int cusum_result;
//
////CWT
//int cwt_result;
//
//void HA_HB_RT_Init(void)
//{
//       //FILTER INITIALIZATION
//       Butterworth_Init(&butterworth_filter);
//       inekf_init_default(&inekf_filter);
//
//       //THRESHOLD INITIALIZATION
//       initThreshold(&th, INIT_THRESHOLD_POS, INIT_THRESHOLD_NEG, TIME_THRESHOLD);
//
//       //DERIVATIVE_DETECT INITIALIZATION
//       DerivativeD_init();
//
//       //ENERGY ENVELOPE INITIALIZATION
//       EnergyEnvelope_init();
//
//       //CUSUM INTIALIZATION
//       CUSUM_init();
//
//       //CWT INITIALIZATION
//       CWT_init();
//
//}
//
//void HA_HB_RT_MainFuntion(void)
//{
//    //RAW IMU DATA (static sint16)
//
//    Pos_AccDataRead(&AccX_s16, &AccY_s16, &AccZ_s16);
//    Pos_GyroDataRead(&GyroX_s32, &GyroY_s32, &GyroZ_s32);
//
//    //ORIENTATION
//    transform_orientation_int(
//        AccX_s16, AccY_s16, AccZ_s16,
//        GyroX_s32, GyroY_s32, GyroZ_s32,
//        &ax, &ay, &az,
//        &gx, &gy, &gz);
//
//
//    // BUTTERWORTH LPF
//    Butterworth_FilterIMU(&butterworth_filter, ax, ay, az, gx, gy, gz,
//                          &fax, &fay, &faz, &fgx, &fgy, &fgz);
//
//    //InEKF FILTER
//    InEKF_ProcessIMU(&inekf_filter, fax, fay, faz, fgx, fgy, fgz, dt,
//                     &ekf_ax, &ekf_ay, &ekf_az, &ekf_gx, &ekf_gy, &ekf_gz);
//
//    //DETECTORS
//
//    //Threshold Detection
//    detected = processIMUSampleRealTime(
//        &th,
//        ekf_ax, ekf_ay, ekf_az,     // filtered accelerometer
//        ekf_gx, ekf_gy, ekf_gz,     // filtered gyroscope
//        sampleIndex,
//        &positiveSpike,
//        &negativeSpike
//    );
//
//    //Derivative Detection
//    deriv_flag = DerivativeD_processSample(
//        ekf_ax, ekf_ay, ekf_az,
//        ekf_gx, ekf_gy, ekf_gz
//     );
//
//    // EnergyEnvelope Detection
//      energy_detected = EnergyEnvelope_processSample(
//              ekf_ax, ekf_ay, ekf_az,
//              ekf_gx, ekf_gy, ekf_gz
//      );
//
//    //CUSUM Detection
//      cusum_result= CUSUM_processSample(
//            ekf_ax, ekf_ay, ekf_az
//      );
//
//    //CWT Detection
//    cwt_result= CWT_processSample(
//            ekf_ax, ekf_ay, ekf_az
//    );
//
//    //DETECTION OUTPUTS
//
//               //Harsh Acceleration
//               if (deriv_flag == 1 && harshaccl_counter == 0)
//               {    harshaccl_counter += 1;
//                    }
//
//               if (positiveSpike == THRESHOLD_TRUE && harshaccl_counter == 1)
//               {   harshaccl_counter += 1;
//                   }
//
//               if (harshaccl_counter == 2 && (cusum_result == 1 || energy_detected==1 )){
//                    }
//
//               if (harshaccl_counter == 2 && cwt_result == 1){
//                   harshaccl_counter=0;
//
//               }
//
//               //Harsh Braking
//               if (deriv_flag == -1 && harshbrak_counter == 0)
//               {harshbrak_counter += 1;}
//
//               if (negativeSpike == THRESHOLD_FALSE && harshbrak_counter == 1)
//               {harshbrak_counter += 1;}
//
//               if (harshbrak_counter == 2 && (cusum_result == -1 )){
//
//                               }
//
//               if (harshbrak_counter == 2 && cwt_result == -1){
//                   harshbrak_counter=0;
//
//               }
//
//
//               sampleIndex++;
////               return 0;
//}

#include "Orientation.h"
#include "InEKF.h"
#include "Butterworth.h"
#include "Threshold.h"
#include "DerivativeDetect.h"
#include "EnergyEnvelope.h"
#include "CUSUM.h"
#include "CWT.h"
#include "HA_HB_RT_main.h"
//#include "PoS_Appl.h"

Butterworth_IMU butterworth_filter;
InEKF inekf_filter;
AdaptiveThreshold th;

const float dt = 0.02; // 50 Hz
int sampleIndex = 0;

int harshaccl_counter = 0;
int harshbrak_counter = 0;

//IMU Data
void Pos_AccDataRead(int *Ax,int *Ay, int *Az);
void Pos_GyroDataRead(int *Gx, int *Gy,int *Gz);

   int GyroX_s32;
   int GyroY_s32;
   int GyroZ_s32;

   int AccX_s16;
   int AccY_s16;
   int AccZ_s16;

//Orientation
float ax, ay, az, gx, gy, gz;

//Butterworth LPF
float fax, fay, faz, fgx, fgy, fgz;

//InEKF
float ekf_ax, ekf_ay, ekf_az, ekf_gx, ekf_gy, ekf_gz;

//Threshold
Threshold_Bool positiveSpike = THRESHOLD_FALSE;
Threshold_Bool negativeSpike = THRESHOLD_FALSE;

Threshold_Bool detected;

//Derivative detect
int deriv_flag;

//Energy envelope
int energy_detected;

//CUSUM
int cusum_result;

//CWT
int cwt_result;

void HA_HB_RT_Init(void)
{
       //FILTER INITIALIZATION
       Butterworth_Init(&butterworth_filter);
//       inekf_init_default(&inekf_filter);

       //THRESHOLD INITIALIZATION
       initThreshold(&th, INIT_THRESHOLD_POS, INIT_THRESHOLD_NEG, TIME_THRESHOLD);

       //DERIVATIVE_DETECT INITIALIZATION
       DerivativeD_init();

       //ENERGY ENVELOPE INITIALIZATION
       EnergyEnvelope_init();

       //CUSUM INTIALIZATION
       CUSUM_init();

       //CWT INITIALIZATION
       CWT_init();

}

void HA_HB_RT_MainFuntion(void)
{
    //RAW IMU DATA (static sint16)

    Pos_AccDataRead(&AccX_s16, &AccY_s16, &AccZ_s16);
    Pos_GyroDataRead(&GyroX_s32, &GyroY_s32, &GyroZ_s32);

    //ORIENTATION
    transform_orientation_int(
        AccX_s16, AccY_s16, AccZ_s16,
        GyroX_s32, GyroY_s32, GyroZ_s32,
        &ax, &ay, &az,
        &gx, &gy, &gz);


    // BUTTERWORTH LPF
    Butterworth_FilterIMU(&butterworth_filter, ax, ay, az, gx, gy, gz,
                          &fax, &fay, &faz, &fgx, &fgy, &fgz);

    //InEKF FILTER
    InEKF_ProcessIMU(&inekf_filter, fax, fay, faz, fgx, fgy, fgz, dt,
                     &ekf_ax, &ekf_ay, &ekf_az, &ekf_gx, &ekf_gy, &ekf_gz);

    //DETECTORS

    //Threshold Detection
    detected = processIMUSampleRealTime(
        &th,
        ekf_ax, ekf_ay, ekf_az, ekf_gx, ekf_gy, ekf_gz,     // filtered gyroscope
        sampleIndex,
        &positiveSpike,
        &negativeSpike
    );

    //Derivative Detection
    deriv_flag = DerivativeD_processSample(
            ekf_ax, ekf_ay, ekf_az, ekf_gx, ekf_gy, ekf_gz
     );

    // EnergyEnvelope Detection
      energy_detected = EnergyEnvelope_processSample(
              ekf_ax, ekf_ay, ekf_az, ekf_gx, ekf_gy, ekf_gz
      );

    //CUSUM Detection
      cusum_result= CUSUM_processSample(
              ekf_ax, ekf_ay, ekf_az
      );

    //CWT Detection
    cwt_result= CWT_processSample(
            ekf_ax, ekf_ay, ekf_az
    );

    //DETECTION OUTPUTS

               //Harsh Acceleration
               if (deriv_flag == 1 && harshaccl_counter == 0)
               {    harshaccl_counter += 1;
                    }

               if (positiveSpike == THRESHOLD_TRUE && harshaccl_counter == 1)
               {   harshaccl_counter += 1;
                   }

               if (harshaccl_counter == 2 && (cusum_result == 1 || energy_detected==1 )){
                    }

               if (harshaccl_counter == 2 && cwt_result == 1){
                   harshaccl_counter=0;

               }

               //Harsh Braking
               if (deriv_flag == -1 && harshbrak_counter == 0)
               {harshbrak_counter += 1;}

               if (negativeSpike == THRESHOLD_FALSE && harshbrak_counter == 1)
               {harshbrak_counter += 1;}

               if (harshbrak_counter == 2 && (cusum_result == -1 )){

                               }

               if (harshbrak_counter == 2 && cwt_result == -1){
                   harshbrak_counter=0;

               }


               sampleIndex++;
//               return 0;
}
