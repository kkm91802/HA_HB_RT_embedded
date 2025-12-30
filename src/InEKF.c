//#include "InEKF.h"
//
///* ============================
//   Small math helpers (kept)
//   ============================ */
//
//static double InEKF_sqrt(double x) {
//    double guess;
//    int i;
//    if (x < 0.0) return 0.0;
//    if (x == 0.0 || x == 1.0) return x;
//    guess = x / 2.0;
//    for (i = 0; i < 20; i++) guess = (guess + x / guess) / 2.0;
//    return guess;
//}
//
//static double InEKF_sin(double x) {
//    double x2, x3, x5, x7;
//    if (x > 3.141592653589793) x -= 6.283185307179586;
//    if (x < -3.141592653589793) x += 6.283185307179586;
//    x2 = x * x;
//    x3 = x * x2;
//    x5 = x3 * x2;
//    x7 = x5 * x2;
//    return x - x3/6.0 + x5/120.0 - x7/5040.0;
//}
//
//static double InEKF_cos(double x) {
//    double x2, x4, x6;
//    if (x > 3.141592653589793) x -= 6.283185307179586;
//    if (x < -3.141592653589793) x += 6.283185307179586;
//    x2 = x * x;
//    x4 = x2 * x2;
//    x6 = x4 * x2;
//    return 1.0 - x2/2.0 + x4/24.0 - x6/720.0;
//}
//
//static double InEKF_fabs(double x) {
//    if (x < 0.0) return -x;
//    return x;
//}
//
///* Matrix helpers */
//static void InEKF_mat_identity(double *A, int n) {
//    int i, j;
//    for (i = 0; i < n * n; ++i) A[i] = 0.0;
//    for (i = 0; i < n; ++i) A[i * n + i] = 1.0;
//}
//
//static void InEKF_mat_zero_elems(double *A, int elems) {
//    int i;
//    for (i = 0; i < elems; ++i) A[i] = 0.0;
//}
//
//static void InEKF_mat_copy(double *dst, const double *src, int n) {
//    int i;
//    for (i = 0; i < n * n; i++) dst[i] = src[i];
//}
//
//static void InEKF_mat3_mul(double *C, const double *A, const double *B) {
//    int i, j, k;
//    for (i = 0; i < 3; ++i) {
//        for (j = 0; j < 3; ++j) {
//            double s = 0.0;
//            for (k = 0; k < 3; ++k) s += A[i*3+k] * B[k*3+j];
//            C[i*3+j] = s;
//        }
//    }
//}
//
//static void InEKF_mat3_vec(double *y, const double *A, const double *x) {
//    int i;
//    for (i = 0; i < 3; ++i) {
//        y[i] = A[i*3+0]*x[0] + A[i*3+1]*x[1] + A[i*3+2]*x[2];
//    }
//}
//
//static void InEKF_mat_transpose(double *B, const double *A, int n) {
//    int i, j;
//    for (i = 0; i < n; ++i) {
//        for (j = 0; j < n; ++j) {
//            B[i*n + j] = A[j*n + i];
//        }
//    }
//}
//
//static double InEKF_norm3(const double *v) {
//    double t0 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
//    return InEKF_sqrt(t0);
//}
//
//static void InEKF_skew3(double *M, const double *w) {
//    M[0] = 0.0;      M[1] = -w[2]; M[2] = w[1];
//    M[3] = w[2];    M[4] = 0.0;    M[5] = -w[0];
//    M[6] = -w[1];   M[7] = w[0];   M[8] = 0.0;
//}
//
//static void InEKF_expSO3(double *R, const double *w, double dt) {
//    double phi[3];
//    double th;
//    double I3[9];
//    double K[9];
//    double u[3];
//    double K2[9];
//    double c, s_val;
//    int i, j, k;
//
//    phi[0] = w[0] * dt;
//    phi[1] = w[1] * dt;
//    phi[2] = w[2] * dt;
//    th = InEKF_norm3(phi);
//
//    /* identity */
//    I3[0]=1.0; I3[1]=0.0; I3[2]=0.0;
//    I3[3]=0.0; I3[4]=1.0; I3[5]=0.0;
//    I3[6]=0.0; I3[7]=0.0; I3[8]=1.0;
//
//    InEKF_skew3(K, phi);
//
//    if (th < 1e-12) {
//        for (i = 0; i < 9; i++) R[i] = I3[i] + K[i];
//        return;
//    }
//
//    u[0] = phi[0] / th; u[1] = phi[1] / th; u[2] = phi[2] / th;
//    InEKF_skew3(K, u);
//
//    /* compute K2 = K * K */
//    for (i = 0; i < 3; i++) {
//        for (j = 0; j < 3; j++) {
//            double s = 0.0;
//            for (k = 0; k < 3; k++) s += K[i*3+k] * K[k*3+j];
//            K2[i*3+j] = s;
//        }
//    }
//
//    c = InEKF_cos(th);
//    s_val = InEKF_sin(th);
//
//    for (i = 0; i < 9; i++) R[i] = I3[i] + s_val * K[i] + (1.0 - c) * K2[i];
//}
//
///* ============================
//   Initialization / config
//   ============================ */
//
//void inekf_init(InEKF *filter, const InEKFConfig *config) {
//    int i;
//    if (!filter || !config) return;
//
//    filter->config = *config;
//
//    InEKF_mat_identity(filter->state.R, 3);
//    for (i = 0; i < 3; i++) {
//        filter->state.v[i] = 0.0;
//        filter->state.p[i] = 0.0;
//        filter->state.bw[i] = 0.0;
//        filter->state.ba[i] = 0.0;
//    }
//
//    InEKF_mat_zero_elems(filter->state.P, NSTATE * NSTATE);
//    for (i = 0; i < 3; i++) filter->state.P[i * NSTATE + i] = 1e-4;
//    for (i = 3; i < 6; i++) filter->state.P[i * NSTATE + i] = 1e-3;
//    for (i = 6; i < 9; i++) filter->state.P[i * NSTATE + i] = 1e-3;
//    for (i = 9; i < 12; i++) filter->state.P[i * NSTATE + i] = 1e-6;
//    for (i = 12; i < 15; i++) filter->state.P[i * NSTATE + i] = 1e-6;
//
//    InEKF_mat_zero_elems(filter->Qc, NSTATE * NSTATE);
//    for (i = 0; i < 3; i++) filter->Qc[i * NSTATE + i] = config->gyro_noise_sigma * config->gyro_noise_sigma;
//    for (i = 3; i < 6; i++) filter->Qc[i * NSTATE + i] = config->accel_noise_sigma * config->accel_noise_sigma;
//    for (i = 9; i < 12; i++) filter->Qc[i * NSTATE + i] = config->gyro_bias_rw * config->gyro_bias_rw;
//    for (i = 12; i < 15; i++) filter->Qc[i * NSTATE + i] = config->accel_bias_rw * config->accel_bias_rw;
//
//    InEKF_mat_zero_elems(filter->R_zupt, 9);
//    {
//        double noise = config->zupt_velocity_noise;
//        filter->R_zupt[0] = noise * noise;
//        filter->R_zupt[4] = noise * noise;
//        filter->R_zupt[8] = noise * noise;
//    }
//}
//
//void inekf_init_default(InEKF *filter) {
//    InEKFConfig cfg;
//    cfg.gyro_noise_sigma = 0.003;
//    cfg.accel_noise_sigma = 0.05;
//    cfg.gyro_bias_rw = 1e-5;
//    cfg.accel_bias_rw = 1e-3;
//    cfg.zupt_velocity_noise = 0.01;
//    cfg.zupt_accel_threshold = 0.1;
//    cfg.zupt_gyro_threshold_deg = 0.5;
//    inekf_init(filter, &cfg);
//}
//
///* ============================
//   F/G, propagate, update
//   ============================ */
//
//static void InEKF_build_F_G(double *F, double *G, const double *R, const double *v, const double *p) {
//    int i, j, k;
//    double gskew[9];
//    double vsk[9];
//    double psk[9];
//
//    InEKF_mat_zero_elems(F, NSTATE * NSTATE);
//    InEKF_mat_zero_elems(G, NSTATE * NSTATE);
//
//    for (i = 0; i < 3; i++) {
//        for (j = 0; j < 3; j++) {
//            F[i * NSTATE + (9 + j)] = -R[i*3 + j];
//        }
//    }
//
//    /* gskew */
//    gskew[0] = 0.0; gskew[1] = -GGRAV; gskew[2] = 0.0;
//    gskew[3] = GGRAV; gskew[4] = 0.0; gskew[5] = 0.0;
//    gskew[6] = 0.0; gskew[7] = 0.0; gskew[8] = 0.0;
//
//    for (i = 0; i < 3; i++) {
//        for (j = 0; j < 3; j++) {
//            F[(3 + i) * NSTATE + j] = gskew[i*3 + j];
//        }
//    }
//
//    InEKF_skew3(vsk, v);
//    for (i = 0; i < 3; i++) {
//        for (j = 0; j < 3; j++) {
//            double s = 0.0;
//            for (k = 0; k < 3; k++) s += vsk[i*3 + k] * R[k*3 + j];
//            F[(3 + i) * NSTATE + (9 + j)] = -s;
//        }
//    }
//
//    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) F[(3 + i) * NSTATE + (12 + j)] = -R[i*3 + j];
//
//    for (i = 0; i < 3; i++) F[(6 + i) * NSTATE + (3 + i)] = 1.0;
//
//    InEKF_skew3(psk, p);
//    for (i = 0; i < 3; i++) {
//        for (j = 0; j < 3; j++) {
//            double s = 0.0;
//            for (k = 0; k < 3; k++) s += psk[i*3 + k] * R[k*3 + j];
//            F[(6 + i) * NSTATE + (9 + j)] = -s;
//        }
//    }
//
//    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) G[i * NSTATE + (0 + j)] = R[i*3 + j];
//
//    for (i = 0; i < 3; i++) {
//        for (j = 0; j < 3; j++) {
//            double s = 0.0;
//            for (k = 0; k < 3; k++) s += vsk[i*3 + k] * R[k*3 + j];
//            G[(3 + i) * NSTATE + (0 + j)] = s;
//            G[(3 + i) * NSTATE + (3 + j)] = R[i*3 + j];
//        }
//    }
//
//    for (i = 0; i < 3; i++) {
//        for (j = 0; j < 3; j++) {
//            double s = 0.0;
//            for (k = 0; k < 3; k++) s += psk[i*3 + k] * R[k*3 + j];
//            G[(6 + i) * NSTATE + (0 + j)] = s;
//            G[(6 + i) * NSTATE + (6 + j)] = R[i*3 + j];
//        }
//    }
//
//    for (i = 0; i < 3; i++) G[(9 + i) * NSTATE + (9 + i)] = -1.0;
//    for (i = 0; i < 3; i++) G[(12 + i) * NSTATE + (12 + i)] = -1.0;
//}
//
//static void InEKF_compute_Phi_firstorder(double *Phi, const double *F, double dt) {
//    int i, j;
//    InEKF_mat_identity(Phi, NSTATE);
//    for (i = 0; i < NSTATE; i++) for (j = 0; j < NSTATE; j++) Phi[i*NSTATE + j] += F[i*NSTATE + j] * dt;
//}
//
//static void InEKF_propagate_covariance(double *Sigma, const double *F, const double *G, const double *Qc, double dt) {
//    int i, j, k;
//    static double Phi[NSTATE * NSTATE];
//    static double tmp1[NSTATE * NSTATE];
//    static double tmp2[NSTATE * NSTATE];
//    static double tmp3[NSTATE * NSTATE];
//    static double GQ[NSTATE * NSTATE];
//    static double GQGT[NSTATE * NSTATE];
//
//    InEKF_compute_Phi_firstorder(Phi, F, dt);
//
//    for (i = 0; i < NSTATE; i++) {
//        for (j = 0; j < NSTATE; j++) {
//            double s = 0.0;
//            for (k = 0; k < NSTATE; k++) s += Phi[i*NSTATE + k] * Sigma[k*NSTATE + j];
//            tmp1[i*NSTATE + j] = s;
//        }
//    }
//
//    for (i = 0; i < NSTATE; i++) {
//        for (j = 0; j < NSTATE; j++) {
//            double s = 0.0;
//            for (k = 0; k < NSTATE; k++) s += tmp1[i*NSTATE + k] * Phi[j*NSTATE + k];
//            tmp2[i*NSTATE + j] = s;
//        }
//    }
//
//    for (i = 0; i < NSTATE; i++) {
//        for (j = 0; j < NSTATE; j++) {
//            double s = 0.0;
//            for (k = 0; k < NSTATE; k++) s += G[i*NSTATE + k] * Qc[k*NSTATE + j];
//            GQ[i*NSTATE + j] = s;
//        }
//    }
//
//    for (i = 0; i < NSTATE; i++) {
//        for (j = 0; j < NSTATE; j++) {
//            double s = 0.0;
//            for (k = 0; k < NSTATE; k++) s += GQ[i*NSTATE + k] * G[j*NSTATE + k];
//            GQGT[i*NSTATE + j] = s;
//        }
//    }
//
//    for (i = 0; i < NSTATE; i++) {
//        for (j = 0; j < NSTATE; j++) {
//            double s = 0.0;
//            for (k = 0; k < NSTATE; k++) s += Phi[i*NSTATE + k] * GQGT[k*NSTATE + j];
//            tmp1[i*NSTATE + j] = s;
//        }
//    }
//
//    for (i = 0; i < NSTATE; i++) {
//        for (j = 0; j < NSTATE; j++) {
//            double s = 0.0;
//            for (k = 0; k < NSTATE; k++) s += tmp1[i*NSTATE + k] * Phi[j*NSTATE + k];
//            tmp3[i*NSTATE + j] = s * dt;
//        }
//    }
//
//    for (i = 0; i < NSTATE; i++) {
//        for (j = 0; j < NSTATE; j++) Sigma[i*NSTATE + j] = tmp2[i*NSTATE + j] + tmp3[i*NSTATE + j];
//    }
//}
//
//static void InEKF_cov_update_joseph(double *P, const double *K, const double *H, const double *Rmeas, int m) {
//    int i, j, r, k, c;
//    double I_N[NSTATE * NSTATE];
//    double KH[NSTATE * NSTATE];
//    double IminusKH[NSTATE * NSTATE];
//    double tmp[NSTATE * NSTATE];
//    double Pnew[NSTATE * NSTATE];
//    double KR[NSTATE * 3];
//
//    InEKF_mat_identity(I_N, NSTATE);
//
//    for (i = 0; i < NSTATE; i++) {
//        for (j = 0; j < NSTATE; j++) {
//            double s = 0.0;
//            for (r = 0; r < m; ++r) {
//                s += K[i*NSTATE + r] * H[r * NSTATE + j];
//            }
//            KH[i*NSTATE + j] = s;
//        }
//    }
//
//    for (i = 0; i < NSTATE; i++) for (j = 0; j < NSTATE; j++) IminusKH[i*NSTATE + j] = I_N[i*NSTATE + j] - KH[i*NSTATE + j];
//
//    for (i = 0; i < NSTATE; i++) {
//        for (j = 0; j < NSTATE; j++) {
//            double s = 0.0;
//            for (k = 0; k < NSTATE; k++) s += IminusKH[i*NSTATE + k] * P[k*NSTATE + j];
//            tmp[i*NSTATE + j] = s;
//        }
//    }
//
//    for (i = 0; i < NSTATE; i++) {
//        for (j = 0; j < NSTATE; j++) {
//            double s = 0.0;
//            for (k = 0; k < NSTATE; k++) s += tmp[i*NSTATE + k] * IminusKH[j*NSTATE + k];
//            Pnew[i*NSTATE + j] = s;
//        }
//    }
//
//    InEKF_mat_zero_elems(KR, NSTATE * 3);
//    for (i = 0; i < NSTATE; i++) {
//        for (r = 0; r < m; r++) {
//            double s = 0.0;
//            for (c = 0; c < m; c++) s += K[i*NSTATE + c] * Rmeas[c*3 + r];
//            KR[i*3 + r] = s;
//        }
//    }
//
//    for (i = 0; i < NSTATE; i++) {
//        for (j = 0; j < NSTATE; j++) {
//            double s = 0.0;
//            for (r = 0; r < m; r++) s += KR[i*3 + r] * K[j*NSTATE + r];
//            Pnew[i*NSTATE + j] += s;
//        }
//    }
//
//    for (i = 0; i < NSTATE * NSTATE; i++) P[i] = Pnew[i];
//}
//
///* ============================
//   Predict / update
//   ============================ */
//
//void inekf_predict(InEKF *filter, const double omega_meas[3], const double acc_meas[3], double dt) {
//    int i, j, k;
//    double omega_corr[3];
//    double acc_corr[3];
//    double Rinc[9];
//    double Rnew[9];
//    double tmpvec[3];
//    double a_world[3];
//    double F[NSTATE * NSTATE];
//    double G[NSTATE * NSTATE];
//
//    if (!filter || dt <= 0.0) return;
//
//    for (i = 0; i < 3; i++) {
//        omega_corr[i] = omega_meas[i] - filter->state.bw[i];
//        acc_corr[i] = acc_meas[i] - filter->state.ba[i];
//    }
//
//    InEKF_expSO3(Rinc, omega_corr, dt);
//
//    InEKF_mat3_mul(Rnew, filter->state.R, Rinc);
//    InEKF_mat_copy(filter->state.R, Rnew, 3);
//
//    InEKF_mat3_vec(tmpvec, filter->state.R, acc_corr);
//    a_world[0] = tmpvec[0];
//    a_world[1] = tmpvec[1];
//    a_world[2] = tmpvec[2] - GGRAV;
//
//    for (i = 0; i < 3; i++) filter->state.v[i] += a_world[i] * dt;
//    for (i = 0; i < 3; i++) filter->state.p[i] += filter->state.v[i] * dt;
//
//    InEKF_build_F_G(F, G, filter->state.R, filter->state.v, filter->state.p);
//    InEKF_propagate_covariance(filter->state.P, F, G, filter->Qc, dt);
//}
//
//static void InEKF_measurement_update_velocity(InEKF *filter, const double z[3]) {
//    int i, j, k, r;
//    double Rt[9];
//    double v_in_imu[3];
//    double y[3];
//    double H[3 * NSTATE];
//    double HP[3 * NSTATE];
//    double S[9];
//    double detS;
//    double Sinv[9];
//    double PHT[NSTATE * 3];
//    double K[NSTATE * 3];
//    double delta[NSTATE];
//    double dphi[3];
//    double Rcor[9];
//    double Rupdated[9];
//
//    if (!filter) return;
//
//    InEKF_mat_transpose(Rt, filter->state.R, 3);
//    InEKF_mat3_vec(v_in_imu, Rt, filter->state.v);
//
//    for (i = 0; i < 3; i++) y[i] = z[i] - v_in_imu[i];
//
//    InEKF_mat_zero_elems(H, 3 * NSTATE);
//    for (r = 0; r < 3; r++) H[r * NSTATE + (3 + r)] = 1.0;
//
//    for (i = 0; i < 3; i++) for (j = 0; j < NSTATE; j++) {
//        double s = 0.0;
//        for (k = 0; k < NSTATE; k++) s += H[i * NSTATE + k] * filter->state.P[k * NSTATE + j];
//        HP[i * NSTATE + j] = s;
//    }
//
//    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
//        double s = 0.0;
//        for (k = 0; k < NSTATE; k++) s += HP[i * NSTATE + k] * H[j * NSTATE + k];
//        S[i*3+j] = s + filter->R_zupt[i*3+j];
//    }
//
//    detS = S[0]*(S[4]*S[8]-S[5]*S[7]) - S[1]*(S[3]*S[8]-S[5]*S[6]) + S[2]*(S[3]*S[7]-S[4]*S[6]);
//
//    if (InEKF_fabs(detS) < 1e-12) {
//        InEKF_mat_zero_elems(Sinv, 9);
//        for (i = 0; i < 3; i++) if (InEKF_fabs(S[i*3 + i]) > 1e-12) Sinv[i*3 + i] = 1.0 / S[i*3 + i];
//    } else {
//        Sinv[0] =  (S[4]*S[8]-S[5]*S[7]) / detS;
//        Sinv[1] = -(S[1]*S[8]-S[2]*S[7]) / detS;
//        Sinv[2] =  (S[1]*S[5]-S[2]*S[4]) / detS;
//        Sinv[3] = -(S[3]*S[8]-S[5]*S[6]) / detS;
//        Sinv[4] =  (S[0]*S[8]-S[2]*S[6]) / detS;
//        Sinv[5] = -(S[0]*S[5]-S[2]*S[3]) / detS;
//        Sinv[6] =  (S[3]*S[7]-S[4]*S[6]) / detS;
//        Sinv[7] = -(S[0]*S[7]-S[1]*S[6]) / detS;
//        Sinv[8] =  (S[0]*S[4]-S[1]*S[3]) / detS;
//    }
//
//    for (i = 0; i < NSTATE; i++) {
//        for (j = 0; j < 3; j++) {
//            double s = 0.0;
//            for (k = 0; k < NSTATE; k++) s += filter->state.P[i * NSTATE + k] * H[j * NSTATE + k];
//            PHT[i*3 + j] = s;
//        }
//    }
//
//    for (i = 0; i < NSTATE; i++) {
//        for (j = 0; j < 3; j++) {
//            double s = 0.0;
//            for (k = 0; k < 3; k++) s += PHT[i*3 + k] * Sinv[k*3 + j];
//            K[i*3 + j] = s;
//        }
//    }
//
//    for (i = 0; i < NSTATE; i++) {
//        double s = 0.0;
//        for (j = 0; j < 3; j++) s += K[i*3 + j] * y[j];
//        delta[i] = s;
//    }
//
//    dphi[0] = delta[0]; dphi[1] = delta[1]; dphi[2] = delta[2];
//    InEKF_expSO3(Rcor, dphi, 1.0);
//    InEKF_mat3_mul(Rupdated, Rcor, filter->state.R);
//    InEKF_mat_copy(filter->state.R, Rupdated, 3);
//
//    for (i = 0; i < 3; i++) {
//        filter->state.v[i] += delta[3 + i];
//        filter->state.p[i] += delta[6 + i];
//        filter->state.bw[i] += delta[9 + i];
//        filter->state.ba[i] += delta[12 + i];
//    }
//
//    InEKF_cov_update_joseph(filter->state.P, K, H, filter->R_zupt, 3);
//}
//
//void inekf_zupt_update(InEKF *filter, const double velocity_meas[3]) {
//    if (!filter) return;
//    InEKF_measurement_update_velocity(filter, velocity_meas);
//}
//
//int inekf_detect_zupt(const InEKF *filter, const double acc_meas[3], const double gyro_meas_deg[3]) {
//    double acc_norm;
//    double gyro_norm_deg;
//    if (!filter) return 0;
//    acc_norm = InEKF_norm3(acc_meas);
//    gyro_norm_deg = InEKF_norm3(gyro_meas_deg);
//    if (InEKF_fabs(acc_norm - GGRAV) < filter->config.zupt_accel_threshold &&
//        gyro_norm_deg < filter->config.zupt_gyro_threshold_deg) return 1;
//    return 0;
//}
//
//void inekf_get_bias_corrected_imu(const InEKF *filter, const double acc_meas[3], const double gyro_meas_deg[3],
//                                  double acc_corrected[3], double gyro_corrected_deg[3]) {
//    int i;
//    double gyro_rad;
//    double corrected_rad;
//    if (!filter) return;
//    if (acc_corrected) {
//        for (i = 0; i < 3; i++) acc_corrected[i] = acc_meas[i] - filter->state.ba[i];
//    }
//    if (gyro_corrected_deg) {
//        for (i = 0; i < 3; i++) {
//            gyro_rad = gyro_meas_deg[i] * DEG2RAD;
//            corrected_rad = gyro_rad - filter->state.bw[i];
//            gyro_corrected_deg[i] = corrected_rad * RAD2DEG;
//        }
//    }
//}
//
///* ============================
//   Real-time scalar API
//   ============================ */
//
//void InEKF_ProcessIMU(
//    InEKF *filter,
//    double ax, double ay, double az,
//    double gx, double gy, double gz,
//    double dt,
//    double *out_ax, double *out_ay, double *out_az,
//    double *out_gx, double *out_gy, double *out_gz)
//{
//    double gyro_rad[3];
//    double acc_meas[3];
//    double gyro_deg[3];
//    double acc_corr[3];
//    double gyro_corr[3];
//
//    if (!filter || dt <= 0.0) return;
//
//    gyro_rad[0] = gx * DEG2RAD; gyro_rad[1] = gy * DEG2RAD; gyro_rad[2] = gz * DEG2RAD;
//    acc_meas[0] = ax; acc_meas[1] = ay; acc_meas[2] = az;
//    gyro_deg[0] = gx; gyro_deg[1] = gy; gyro_deg[2] = gz;
//
//    /* Predict */
//    inekf_predict(filter, gyro_rad, acc_meas, dt);
//
//    /* ZUPT detection and update */
//    if (inekf_detect_zupt(filter, acc_meas, gyro_deg)) {
//        double zero_v[3];
//        zero_v[0] = 0.0; zero_v[1] = 0.0; zero_v[2] = 0.0;
//        inekf_zupt_update(filter, zero_v);
//    }
//
//    /* Get bias-corrected data */
//    inekf_get_bias_corrected_imu(filter, acc_meas, gyro_deg, acc_corr, gyro_corr);
//
//    if (out_ax) *out_ax = acc_corr[0];
//    if (out_ay) *out_ay = acc_corr[1];
//    if (out_az) *out_az = acc_corr[2];
//    if (out_gx) *out_gx = gyro_corr[0];
//    if (out_gy) *out_gy = gyro_corr[1];
//    if (out_gz) *out_gz = gyro_corr[2];
//}


#include "InEKF.h"

/* ============================
   Small math helpers (kept)
   ============================ */

static float InEKF_sqrt(float x) {
    float guess;
    int i;
    if (x < 0.0f) return 0.0f;
    if (x == 0.0f || x == 1.0f) return x;
    guess = x / 2.0f;
    for (i = 0; i < 20; i++) guess = (guess + x / guess) / 2.0f;
    return guess;
}

static float InEKF_sin(float x) {
    float x2, x3, x5, x7;
    if (x > 3.141592653589793f) x -= 6.283185307179586f;
    if (x < -3.141592653589793f) x += 6.283185307179586f;
    x2 = x * x;
    x3 = x * x2;
    x5 = x3 * x2;
    x7 = x5 * x2;
    return x - x3/6.0f + x5/120.0f - x7/5040.0f;
}

static float InEKF_cos(float x) {
    float x2, x4, x6;
    if (x > 3.141592653589793f) x -= 6.283185307179586f;
    if (x < -3.141592653589793f) x += 6.283185307179586f;
    x2 = x * x;
    x4 = x2 * x2;
    x6 = x4 * x2;
    return 1.0f - x2/2.0f + x4/24.0f - x6/720.0f;
}

static float InEKF_fabs(float x) {
    if (x < 0.0f) return -x;
    return x;
}

/* Matrix helpers */
static void InEKF_mat_identity(float *A, int n) {
    int i, j;
    for (i = 0; i < n * n; ++i) A[i] = 0.0f;
    for (i = 0; i < n; ++i) A[i * n + i] = 1.0f;
}

static void InEKF_mat_zero_elems(float *A, int elems) {
    int i;
    for (i = 0; i < elems; ++i) A[i] = 0.0f;
}

static void InEKF_mat_copy(float *dst, const float *src, int n) {
    int i;
    for (i = 0; i < n * n; i++) dst[i] = src[i];
}

static void InEKF_mat3_mul(float *C, const float *A, const float *B) {
    int i, j, k;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            float s = 0.0f;
            for (k = 0; k < 3; ++k) s += A[i*3+k] * B[k*3+j];
            C[i*3+j] = s;
        }
    }
}

static void InEKF_mat3_vec(float *y, const float *A, const float *x) {
    int i;
    for (i = 0; i < 3; ++i) {
        y[i] = A[i*3+0]*x[0] + A[i*3+1]*x[1] + A[i*3+2]*x[2];
    }
}

static void InEKF_mat_transpose(float *B, const float *A, int n) {
    int i, j;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            B[i*n + j] = A[j*n + i];
        }
    }
}

static float InEKF_norm3(const float *v) {
    float t0 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    return InEKF_sqrt(t0);
}

static void InEKF_skew3(float *M, const float *w) {
    M[0] = 0.0f;      M[1] = -w[2]; M[2] = w[1];
    M[3] = w[2];    M[4] = 0.0f;    M[5] = -w[0];
    M[6] = -w[1];   M[7] = w[0];   M[8] = 0.0f;
}

static void InEKF_expSO3(float *R, const float *w, float dt) {
    float phi[3];
    float th;
    float I3[9];
    float K[9];
    float u[3];
    float K2[9];
    float c, s_val;
    int i, j, k;

    phi[0] = w[0] * dt;
    phi[1] = w[1] * dt;
    phi[2] = w[2] * dt;
    th = InEKF_norm3(phi);

    /* identity */
    I3[0]=1.0f; I3[1]=0.0f; I3[2]=0.0f;
    I3[3]=0.0f; I3[4]=1.0f; I3[5]=0.0f;
    I3[6]=0.0f; I3[7]=0.0f; I3[8]=1.0f;

    InEKF_skew3(K, phi);

    if (th < 1e-12f) {
        for (i = 0; i < 9; i++) R[i] = I3[i] + K[i];
        return;
    }

    u[0] = phi[0] / th; u[1] = phi[1] / th; u[2] = phi[2] / th;
    InEKF_skew3(K, u);

    /* compute K2 = K * K */
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            float s = 0.0f;
            for (k = 0; k < 3; k++) s += K[i*3+k] * K[k*3+j];
            K2[i*3+j] = s;
        }
    }

    c = InEKF_cos(th);
    s_val = InEKF_sin(th);

    for (i = 0; i < 9; i++) R[i] = I3[i] + s_val * K[i] + (1.0f - c) * K2[i];
}

/* ============================
   Initialization / config
   ============================ */

void inekf_init(InEKF *filter, const InEKFConfig *config) {
    int i;
    if (!filter || !config) return;

    filter->config = *config;

    InEKF_mat_identity(filter->state.R, 3);
    for (i = 0; i < 3; i++) {
        filter->state.v[i] = 0.0f;
        filter->state.p[i] = 0.0f;
        filter->state.bw[i] = 0.0f;
        filter->state.ba[i] = 0.0f;
    }

    InEKF_mat_zero_elems(filter->state.P, NSTATE * NSTATE);
    for (i = 0; i < 3; i++) filter->state.P[i * NSTATE + i] = 1e-4f;
    for (i = 3; i < 6; i++) filter->state.P[i * NSTATE + i] = 1e-3f;
    for (i = 6; i < 9; i++) filter->state.P[i * NSTATE + i] = 1e-3f;
    for (i = 9; i < 12; i++) filter->state.P[i * NSTATE + i] = 1e-6f;
    for (i = 12; i < 15; i++) filter->state.P[i * NSTATE + i] = 1e-6f;

    InEKF_mat_zero_elems(filter->Qc, NSTATE * NSTATE);
    for (i = 0; i < 3; i++) filter->Qc[i * NSTATE + i] = config->gyro_noise_sigma * config->gyro_noise_sigma;
    for (i = 3; i < 6; i++) filter->Qc[i * NSTATE + i] = config->accel_noise_sigma * config->accel_noise_sigma;
    for (i = 9; i < 12; i++) filter->Qc[i * NSTATE + i] = config->gyro_bias_rw * config->gyro_bias_rw;
    for (i = 12; i < 15; i++) filter->Qc[i * NSTATE + i] = config->accel_bias_rw * config->accel_bias_rw;

    InEKF_mat_zero_elems(filter->R_zupt, 9);
    {
        float noise = config->zupt_velocity_noise;
        filter->R_zupt[0] = noise * noise;
        filter->R_zupt[4] = noise * noise;
        filter->R_zupt[8] = noise * noise;
    }
}

void inekf_init_default(InEKF *filter) {
    InEKFConfig cfg;
    cfg.gyro_noise_sigma = 0.003f;
    cfg.accel_noise_sigma = 0.05f;
    cfg.gyro_bias_rw = 1e-5f;
    cfg.accel_bias_rw = 1e-3f;
    cfg.zupt_velocity_noise = 0.01f;
    cfg.zupt_accel_threshold = 0.1f;
    cfg.zupt_gyro_threshold_deg = 0.5f;
    inekf_init(filter, &cfg);
}

/* ============================
   F/G, propagate, update
   ============================ */

static void InEKF_build_F_G(float *F, float *G, const float *R, const float *v, const float *p) {
    int i, j, k;
    float gskew[9];
    float vsk[9];
    float psk[9];

    InEKF_mat_zero_elems(F, NSTATE * NSTATE);
    InEKF_mat_zero_elems(G, NSTATE * NSTATE);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            F[i * NSTATE + (9 + j)] = -R[i*3 + j];
        }
    }

    /* gskew */
    gskew[0] = 0.0f; gskew[1] = -GGRAV; gskew[2] = 0.0f;
    gskew[3] = GGRAV; gskew[4] = 0.0f; gskew[5] = 0.0f;
    gskew[6] = 0.0f; gskew[7] = 0.0f; gskew[8] = 0.0f;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            F[(3 + i) * NSTATE + j] = gskew[i*3 + j];
        }
    }

    InEKF_skew3(vsk, v);
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            float s = 0.0f;
            for (k = 0; k < 3; k++) s += vsk[i*3 + k] * R[k*3 + j];
            F[(3 + i) * NSTATE + (9 + j)] = -s;
        }
    }

    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) F[(3 + i) * NSTATE + (12 + j)] = -R[i*3 + j];

    for (i = 0; i < 3; i++) F[(6 + i) * NSTATE + (3 + i)] = 1.0f;

    InEKF_skew3(psk, p);
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            float s = 0.0f;
            for (k = 0; k < 3; k++) s += psk[i*3 + k] * R[k*3 + j];
            F[(6 + i) * NSTATE + (9 + j)] = -s;
        }
    }

    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) G[i * NSTATE + (0 + j)] = R[i*3 + j];

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            float s = 0.0f;
            for (k = 0; k < 3; k++) s += vsk[i*3 + k] * R[k*3 + j];
            G[(3 + i) * NSTATE + (0 + j)] = s;
            G[(3 + i) * NSTATE + (3 + j)] = R[i*3 + j];
        }
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            float s = 0.0f;
            for (k = 0; k < 3; k++) s += psk[i*3 + k] * R[k*3 + j];
            G[(6 + i) * NSTATE + (0 + j)] = s;
            G[(6 + i) * NSTATE + (6 + j)] = R[i*3 + j];
        }
    }

    for (i = 0; i < 3; i++) G[(9 + i) * NSTATE + (9 + i)] = -1.0f;
    for (i = 0; i < 3; i++) G[(12 + i) * NSTATE + (12 + i)] = -1.0f;
}

static void InEKF_compute_Phi_firstorder(float *Phi, const float *F, float dt) {
    int i, j;
    InEKF_mat_identity(Phi, NSTATE);
    for (i = 0; i < NSTATE; i++) for (j = 0; j < NSTATE; j++) Phi[i*NSTATE + j] += F[i*NSTATE + j] * dt;
}

static void InEKF_propagate_covariance(float *Sigma, const float *F, const float *G, const float *Qc, float dt) {
    int i, j, k;
    static float Phi[NSTATE * NSTATE];
    static float tmp1[NSTATE * NSTATE];
    static float tmp2[NSTATE * NSTATE];
    static float tmp3[NSTATE * NSTATE];
    static float GQ[NSTATE * NSTATE];
    static float GQGT[NSTATE * NSTATE];

    InEKF_compute_Phi_firstorder(Phi, F, dt);

    for (i = 0; i < NSTATE; i++) {
        for (j = 0; j < NSTATE; j++) {
            float s = 0.0f;
            for (k = 0; k < NSTATE; k++) s += Phi[i*NSTATE + k] * Sigma[k*NSTATE + j];
            tmp1[i*NSTATE + j] = s;
        }
    }

    for (i = 0; i < NSTATE; i++) {
        for (j = 0; j < NSTATE; j++) {
            float s = 0.0f;
            for (k = 0; k < NSTATE; k++) s += tmp1[i*NSTATE + k] * Phi[j*NSTATE + k];
            tmp2[i*NSTATE + j] = s;
        }
    }

    for (i = 0; i < NSTATE; i++) {
        for (j = 0; j < NSTATE; j++) {
            float s = 0.0f;
            for (k = 0; k < NSTATE; k++) s += G[i*NSTATE + k] * Qc[k*NSTATE + j];
            GQ[i*NSTATE + j] = s;
        }
    }

    for (i = 0; i < NSTATE; i++) {
        for (j = 0; j < NSTATE; j++) {
            float s = 0.0f;
            for (k = 0; k < NSTATE; k++) s += GQ[i*NSTATE + k] * G[j*NSTATE + k];
            GQGT[i*NSTATE + j] = s;
        }
    }

    for (i = 0; i < NSTATE; i++) {
        for (j = 0; j < NSTATE; j++) {
            float s = 0.0f;
            for (k = 0; k < NSTATE; k++) s += Phi[i*NSTATE + k] * GQGT[k*NSTATE + j];
            tmp1[i*NSTATE + j] = s;
        }
    }

    for (i = 0; i < NSTATE; i++) {
        for (j = 0; j < NSTATE; j++) {
            float s = 0.0f;
            for (k = 0; k < NSTATE; k++) s += tmp1[i*NSTATE + k] * Phi[j*NSTATE + k];
            tmp3[i*NSTATE + j] = s * dt;
        }
    }

    for (i = 0; i < NSTATE; i++) {
        for (j = 0; j < NSTATE; j++) Sigma[i*NSTATE + j] = tmp2[i*NSTATE + j] + tmp3[i*NSTATE + j];
    }
}

static void InEKF_cov_update_joseph(float *P, const float *K, const float *H, const float *Rmeas, int m) {
    int i, j, r, k, c;
    float I_N[NSTATE * NSTATE];
    float KH[NSTATE * NSTATE];
    float IminusKH[NSTATE * NSTATE];
    float tmp[NSTATE * NSTATE];
    float Pnew[NSTATE * NSTATE];
    float KR[NSTATE * 3];

    InEKF_mat_identity(I_N, NSTATE);

    for (i = 0; i < NSTATE; i++) {
        for (j = 0; j < NSTATE; j++) {
            float s = 0.0f;
            for (r = 0; r < m; ++r) {
                s += K[i*NSTATE + r] * H[r * NSTATE + j];
            }
            KH[i*NSTATE + j] = s;
        }
    }

    for (i = 0; i < NSTATE; i++) for (j = 0; j < NSTATE; j++) IminusKH[i*NSTATE + j] = I_N[i*NSTATE + j] - KH[i*NSTATE + j];

    for (i = 0; i < NSTATE; i++) {
        for (j = 0; j < NSTATE; j++) {
            float s = 0.0f;
            for (k = 0; k < NSTATE; k++) s += IminusKH[i*NSTATE + k] * P[k*NSTATE + j];
            tmp[i*NSTATE + j] = s;
        }
    }

    for (i = 0; i < NSTATE; i++) {
        for (j = 0; j < NSTATE; j++) {
            float s = 0.0f;
            for (k = 0; k < NSTATE; k++) s += tmp[i*NSTATE + k] * IminusKH[j*NSTATE + k];
            Pnew[i*NSTATE + j] = s;
        }
    }

    InEKF_mat_zero_elems(KR, NSTATE * 3);
    for (i = 0; i < NSTATE; i++) {
        for (r = 0; r < m; r++) {
            float s = 0.0f;
            for (c = 0; c < m; c++) s += K[i*NSTATE + c] * Rmeas[c*3 + r];
            KR[i*3 + r] = s;
        }
    }

    for (i = 0; i < NSTATE; i++) {
        for (j = 0; j < NSTATE; j++) {
            float s = 0.0f;
            for (r = 0; r < m; r++) s += KR[i*3 + r] * K[j*NSTATE + r];
            Pnew[i*NSTATE + j] += s;
        }
    }

    for (i = 0; i < NSTATE * NSTATE; i++) P[i] = Pnew[i];
}

/* ============================
   Predict / update
   ============================ */

void inekf_predict(InEKF *filter, const float omega_meas[3], const float acc_meas[3], float dt) {
    int i, j, k;
    float omega_corr[3];
    float acc_corr[3];
    float Rinc[9];
    float Rnew[9];
    float tmpvec[3];
    float a_world[3];
    float F[NSTATE * NSTATE];
    float G[NSTATE * NSTATE];

    if (!filter || dt <= 0.0f) return;

    for (i = 0; i < 3; i++) {
        omega_corr[i] = omega_meas[i] - filter->state.bw[i];
        acc_corr[i] = acc_meas[i] - filter->state.ba[i];
    }

    InEKF_expSO3(Rinc, omega_corr, dt);

    InEKF_mat3_mul(Rnew, filter->state.R, Rinc);
    InEKF_mat_copy(filter->state.R, Rnew, 3);

    InEKF_mat3_vec(tmpvec, filter->state.R, acc_corr);
    a_world[0] = tmpvec[0];
    a_world[1] = tmpvec[1];
    a_world[2] = tmpvec[2] - GGRAV;

    for (i = 0; i < 3; i++) filter->state.v[i] += a_world[i] * dt;
    for (i = 0; i < 3; i++) filter->state.p[i] += filter->state.v[i] * dt;

    InEKF_build_F_G(F, G, filter->state.R, filter->state.v, filter->state.p);
    InEKF_propagate_covariance(filter->state.P, F, G, filter->Qc, dt);
}

static void InEKF_measurement_update_velocity(InEKF *filter, const float z[3]) {
    int i, j, k, r;
    float Rt[9];
    float v_in_imu[3];
    float y[3];
    float H[3 * NSTATE];
    float HP[3 * NSTATE];
    float S[9];
    float detS;
    float Sinv[9];
    float PHT[NSTATE * 3];
    float K[NSTATE * 3];
    float delta[NSTATE];
    float dphi[3];
    float Rcor[9];
    float Rupdated[9];

    if (!filter) return;

    InEKF_mat_transpose(Rt, filter->state.R, 3);
    InEKF_mat3_vec(v_in_imu, Rt, filter->state.v);

    for (i = 0; i < 3; i++) y[i] = z[i] - v_in_imu[i];

    InEKF_mat_zero_elems(H, 3 * NSTATE);
    for (r = 0; r < 3; r++) H[r * NSTATE + (3 + r)] = 1.0f;

    for (i = 0; i < 3; i++) for (j = 0; j < NSTATE; j++) {
        float s = 0.0f;
        for (k = 0; k < NSTATE; k++) s += H[i * NSTATE + k] * filter->state.P[k * NSTATE + j];
        HP[i * NSTATE + j] = s;
    }

    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
        float s = 0.0f;
        for (k = 0; k < NSTATE; k++) s += HP[i * NSTATE + k] * H[j * NSTATE + k];
        S[i*3+j] = s + filter->R_zupt[i*3+j];
    }

    detS = S[0]*(S[4]*S[8]-S[5]*S[7]) - S[1]*(S[3]*S[8]-S[5]*S[6]) + S[2]*(S[3]*S[7]-S[4]*S[6]);

    if (InEKF_fabs(detS) < 1e-12f) {
        InEKF_mat_zero_elems(Sinv, 9);
        for (i = 0; i < 3; i++) if (InEKF_fabs(S[i*3 + i]) > 1e-12f) Sinv[i*3 + i] = 1.0f / S[i*3 + i];
    } else {
        Sinv[0] =  (S[4]*S[8]-S[5]*S[7]) / detS;
        Sinv[1] = -(S[1]*S[8]-S[2]*S[7]) / detS;
        Sinv[2] =  (S[1]*S[5]-S[2]*S[4]) / detS;
        Sinv[3] = -(S[3]*S[8]-S[5]*S[6]) / detS;
        Sinv[4] =  (S[0]*S[8]-S[2]*S[6]) / detS;
        Sinv[5] = -(S[0]*S[5]-S[2]*S[3]) / detS;
        Sinv[6] =  (S[3]*S[7]-S[4]*S[6]) / detS;
        Sinv[7] = -(S[0]*S[7]-S[1]*S[6]) / detS;
        Sinv[8] =  (S[0]*S[4]-S[1]*S[3]) / detS;
    }

    for (i = 0; i < NSTATE; i++) {
        for (j = 0; j < 3; j++) {
            float s = 0.0f;
            for (k = 0; k < NSTATE; k++) s += filter->state.P[i * NSTATE + k] * H[j * NSTATE + k];
            PHT[i*3 + j] = s;
        }
    }

    for (i = 0; i < NSTATE; i++) {
        for (j = 0; j < 3; j++) {
            float s = 0.0f;
            for (k = 0; k < 3; k++) s += PHT[i*3 + k] * Sinv[k*3 + j];
            K[i*3 + j] = s;
        }
    }

    for (i = 0; i < NSTATE; i++) {
        float s = 0.0f;
        for (j = 0; j < 3; j++) s += K[i*3 + j] * y[j];
        delta[i] = s;
    }

    dphi[0] = delta[0]; dphi[1] = delta[1]; dphi[2] = delta[2];
    InEKF_expSO3(Rcor, dphi, 1.0f);
    InEKF_mat3_mul(Rupdated, Rcor, filter->state.R);
    InEKF_mat_copy(filter->state.R, Rupdated, 3);

    for (i = 0; i < 3; i++) {
        filter->state.v[i] += delta[3 + i];
        filter->state.p[i] += delta[6 + i];
        filter->state.bw[i] += delta[9 + i];
        filter->state.ba[i] += delta[12 + i];
    }

    InEKF_cov_update_joseph(filter->state.P, K, H, filter->R_zupt, 3);
}

void inekf_zupt_update(InEKF *filter, const float velocity_meas[3]) {
    if (!filter) return;
    InEKF_measurement_update_velocity(filter, velocity_meas);
}

int inekf_detect_zupt(const InEKF *filter, const float acc_meas[3], const float gyro_meas_deg[3]) {
    float acc_norm;
    float gyro_norm_deg;
    if (!filter) return 0;
    acc_norm = InEKF_norm3(acc_meas);
    gyro_norm_deg = InEKF_norm3(gyro_meas_deg);
    if (InEKF_fabs(acc_norm - GGRAV) < filter->config.zupt_accel_threshold &&
        gyro_norm_deg < filter->config.zupt_gyro_threshold_deg) return 1;
    return 0;
}

void inekf_get_bias_corrected_imu(const InEKF *filter, const float acc_meas[3], const float gyro_meas_deg[3],
                                  float acc_corrected[3], float gyro_corrected_deg[3]) {
    int i;
    float gyro_rad;
    float corrected_rad;
    if (!filter) return;
    if (acc_corrected) {
        for (i = 0; i < 3; i++) acc_corrected[i] = acc_meas[i] - filter->state.ba[i];
    }
    if (gyro_corrected_deg) {
        for (i = 0; i < 3; i++) {
            gyro_rad = gyro_meas_deg[i] * DEG2RAD;
            corrected_rad = gyro_rad - filter->state.bw[i];
            gyro_corrected_deg[i] = corrected_rad * RAD2DEG;
        }
    }
}

/* ============================
   Real-time scalar API
   ============================ */

void InEKF_ProcessIMU(
    InEKF *filter,
    float ax, float ay, float az,
    float gx, float gy, float gz,
    float dt,
    float *out_ax, float *out_ay, float *out_az,
    float *out_gx, float *out_gy, float *out_gz)
{
    float gyro_rad[3];
    float acc_meas[3];
    float gyro_deg[3];
    float acc_corr[3];
    float gyro_corr[3];

    if (!filter || dt <= 0.0f) return;

    gyro_rad[0] = gx * DEG2RAD; gyro_rad[1] = gy * DEG2RAD; gyro_rad[2] = gz * DEG2RAD;
    acc_meas[0] = ax; acc_meas[1] = ay; acc_meas[2] = az;
    gyro_deg[0] = gx; gyro_deg[1] = gy; gyro_deg[2] = gz;

    /* Predict */
    inekf_predict(filter, gyro_rad, acc_meas, dt);

    /* ZUPT detection and update */
    if (inekf_detect_zupt(filter, acc_meas, gyro_deg)) {
        float zero_v[3];
        zero_v[0] = 0.0f; zero_v[1] = 0.0f; zero_v[2] = 0.0f;
        inekf_zupt_update(filter, zero_v);
    }

    /* Get bias-corrected data */
    inekf_get_bias_corrected_imu(filter, acc_meas, gyro_deg, acc_corr, gyro_corr);

    if (out_ax) *out_ax = acc_corr[0];
    if (out_ay) *out_ay = acc_corr[1];
    if (out_az) *out_az = acc_corr[2];
    if (out_gx) *out_gx = gyro_corr[0];
    if (out_gy) *out_gy = gyro_corr[1];
    if (out_gz) *out_gz = gyro_corr[2];
}
