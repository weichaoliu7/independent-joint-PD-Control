#include <stdio.h>
#include <stdlib.h>
#include "step.h"

int main()
{

    double Ts = 0.01; // sampling period
    double t0 = 0.0;
    double t1 = 10.0;
    double theta = 0.0; // threshold

    int ARRAY_SIZE = (t1 - t0) / Ts;

    Data *qd1 = malloc(sizeof(Data));
    Data *qd2 = malloc(sizeof(Data));

    step(qd1, theta, Ts, t0, t1);  // target position of link 1
    step(qd2, theta, Ts, t0, t1);  // target position of link 2

    double t;
    double ctrl_u1[ARRAY_SIZE], ctrl_u2[ARRAY_SIZE], ctrl_u3[ARRAY_SIZE];
    double ctrl_u4[ARRAY_SIZE], ctrl_u5[ARRAY_SIZE], ctrl_u6[ARRAY_SIZE];

    ctrl_u1[0] = qd1->y[0];
    ctrl_u2[0] = qd2->y[0];

    double R1[ARRAY_SIZE], R2[ARRAY_SIZE];

    R1[0] = ctrl_u1[0];
    R2[0] = ctrl_u2[0];

    double dr1 = 0;
    double dr2 = 0;

    ctrl_u3[0] = 0;
    ctrl_u4[0] = 0;
    ctrl_u5[0] = 0;
    ctrl_u6[0] = 0;

    double x1[ARRAY_SIZE], x2[ARRAY_SIZE], x3[ARRAY_SIZE], x4[ARRAY_SIZE];

    x1[0] = ctrl_u3[0];
    x2[0] = ctrl_u4[0];
    x3[0] = ctrl_u5[0];
    x4[0] = ctrl_u6[0];

    double e1[ARRAY_SIZE], e2[ARRAY_SIZE], e[2][ARRAY_SIZE];

    e1[0] = R1[0] - x1[0];  // position error of link 1
    e2[0] = R1[0] - x3[0];  // position error of link 2

    e[0][0] = e1[0];
    e[1][0] = e2[0];

    double de1[ARRAY_SIZE], de2[ARRAY_SIZE], de[2][ARRAY_SIZE];
    de1[0] = dr1 - x2[0];  // position error derivative of link 1
    de2[0] = dr2 - x4[0];  // position error derivative of link 2

    de[0][0] = de1[0];
    de[1][0] = de2[0];

    double Kp[2][2] = {{30, 0}, {0, 30}};  // proportional gain 
    double Kd[2][2] = {{30, 0}, {0, 30}};  // differential gain

    double tol[2][ARRAY_SIZE];  // control torque

    for (int j = 0; j < 2; j++){
        tol[j][0] = 0.0;

        for (int k = 0; k < 2; k++){
            tol[j][0] += Kp[j][k] * e[k][0] + Kd[j][k] * de[k][0];  // PD control law
        }
    }

    double ctrl_y1[ARRAY_SIZE], ctrl_y2[ARRAY_SIZE];

    ctrl_y1[0] = tol[0][0];
    ctrl_y2[0] = tol[1][0];

    double p[5] = {2.9, 0.76, 0.87, 3.04, 0.87};
    double D0[2][2][ARRAY_SIZE], C0[2][2][ARRAY_SIZE];

    D0[0][0][0] = p[0] + p[1] + 2 * p[2] * cos(x3[0]);  // inertia matrix
    D0[0][1][0] = p[1] + p[2] * cos(x3[0]);
    D0[1][0][0] = p[1] + p[2] * cos(x3[0]);
    D0[1][1][0] = p[1];

    C0[0][0][0] = -p[2] * x4[0] * sin(x3[0]) - p[2] * (x2[0] + x4[0]) * sin(x3[0]);  // coriolis force
    C0[0][1][0] = p[2] * x2[0] * sin(x3[0]);
    C0[1][0][0] = p[2] * x2[0] * sin(x3[0]);
    C0[1][1][0] = 0;

    double plant_u1[ARRAY_SIZE], plant_u2[ARRAY_SIZE];

    plant_u1[0] = tol[0][0];
    plant_u2[0] = tol[1][0];

    double dq[2][ARRAY_SIZE];

    dq[0][0] = x2[0];  // position derivative of link 1
    dq[1][0] = x4[0];  // position derivative of link 2

    double a = D0[0][0][0], b = D0[0][1][0], c = D0[1][0][0], d = D0[1][1][0];

    double D0_inv[2][2][ARRAY_SIZE], tol_C0dq[2][2][ARRAY_SIZE], S[2][ARRAY_SIZE];
    double tol_C0dq_1[ARRAY_SIZE], tol_C0dq_2[ARRAY_SIZE];

    D0_inv[0][0][0] = d / (a*d - b*c);
    D0_inv[0][1][0] = -b / (a*d - b *c);
    D0_inv[1][0][0] = -c / (a*d - b *c);
    D0_inv[1][1][0] = a / (a*d - b *c);

    tol_C0dq_1[0] = tol[0][0] - ( C0[0][0][0] * dq[0][0] + C0[0][1][0] * dq[1][0]);
    tol_C0dq_2[0] = tol[1][0] - ( C0[1][0][0] * dq[0][0] + C0[1][1][0] * dq[1][0] );
    S[0][0] = D0_inv[0][0][0] * tol_C0dq_1[0] + D0_inv[0][1][0] * tol_C0dq_2[0];
    S[1][0] = D0_inv[1][0][0] * tol_C0dq_1[0] + D0_inv[1][1][0] * tol_C0dq_2[0];

    double ddq[2][ARRAY_SIZE], q[2][ARRAY_SIZE];

    ddq[0][0]=S[0][0];  // second order position derivative of link 1
    ddq[1][0]=S[1][0];  // second order position derivative of link 2

    dq[0][1]=dq[0][0]+ddq[0][0]*Ts;
    dq[1][1]=dq[1][0]+ddq[1][0]*Ts;

    q[0][1]=q[0][0]+dq[0][1]*Ts;
    q[1][1]=q[1][0]+dq[1][1]*Ts;

    double plant_y1[ARRAY_SIZE], plant_y2[ARRAY_SIZE], plant_y3[ARRAY_SIZE], plant_y4[ARRAY_SIZE];

    plant_y1[0]=q[0][1];
    plant_y2[0]=dq[0][1];
    plant_y3[0]=q[1][1];
    plant_y4[0]=dq[1][1];


    for (int i = 1; i < ARRAY_SIZE; i++){

        ctrl_u1[i] = qd1->y[i];
        ctrl_u2[i] = qd2->y[i];

        R1[i] = ctrl_u1[i];
        R2[i] = ctrl_u2[i];

        ctrl_u3[i] = plant_y1[i-1];
        ctrl_u4[i] = plant_y2[i-1];
        ctrl_u5[i] = plant_y3[i-1];
        ctrl_u6[i] = plant_y4[i-1];

        x1[i] = ctrl_u3[i];
        x2[i] = ctrl_u4[i];
        x3[i] = ctrl_u5[i];
        x4[i] = ctrl_u6[i];

        e1[i] = R1[i] - x1[i];
        e2[i] = R1[i] - x3[i];

        e[0][i] = e1[i];
        e[1][i] = e2[i];

        de1[i] = dr1 - x2[i];
        de2[i] = dr2 - x4[i];

        de[0][i] = de1[i];
        de[1][i] = de2[i];

        for (int j = 0; j < 2; j++){
            tol[j][i] = 0.0;

            for (int k = 0; k < 2; k++){
                tol[j][i] += Kp[j][k] * e[k][i] + Kd[j][k] * de[k][i];
            }
        }

        ctrl_y1[i] = tol[0][i];
        ctrl_y2[i] = tol[1][i];

        D0[0][0][i] = p[0] + p[1] + 2 * p[2] * cos(x3[i]);
        D0[0][1][i] = p[1] + p[2] * cos(x3[i]);
        D0[1][0][i] = p[1] + p[2] * cos(x3[i]);
        D0[1][1][i] = p[1];

        C0[0][0][i] = -p[2] * x4[0] * sin(x3[i]) - p[2] * (x2[i] + x4[i]) * sin(x3[i]);
        C0[0][1][i] = p[2] * x2[0] * sin(x3[i]);
        C0[1][0][i] = p[2] * x2[0] * sin(x3[i]);
        C0[1][1][i] = 0;

        plant_u1[i] = tol[0][i];
        plant_u2[i] = tol[1][i];

        dq[0][i] = x2[i];
        dq[1][i] = x4[i];

        double a = D0[0][0][i], b = D0[0][1][i], c = D0[1][0][i], d = D0[1][1][i];

        D0_inv[0][0][i] = d / (a*d - b*c);
        D0_inv[0][1][i] = -b / (a*d - b *c);
        D0_inv[1][0][i] = -c / (a*d - b *c);
        D0_inv[1][1][i] = a / (a*d - b *c);

        tol_C0dq_1[i] = tol[0][i] - ( C0[0][0][i] * dq[0][i] + C0[0][1][i] * dq[1][i]);
        tol_C0dq_2[i] = tol[1][i] - ( C0[1][0][i] * dq[0][i] + C0[1][1][i] * dq[1][i] );
        S[0][i] = D0_inv[0][0][i] * tol_C0dq_1[i] + D0_inv[0][1][i] * tol_C0dq_2[i];
        S[1][i] = D0_inv[1][0][i] * tol_C0dq_1[i] + D0_inv[1][1][i] * tol_C0dq_2[i];

        ddq[0][i]=S[0][i];
        ddq[1][i]=S[1][i];

        dq[0][i+1]=dq[0][i]+ddq[0][i]*Ts;
        dq[1][i+1]=dq[1][i]+ddq[1][i]*Ts;

        q[0][i+1]=q[0][i]+dq[0][i+1]*Ts;
        q[1][i+1]=q[1][i]+dq[1][i+1]*Ts;

        plant_y1[i]=q[0][i+1];
        plant_y2[i]=dq[0][i+1];
        plant_y3[i]=q[1][i+1];
        plant_y4[i]=dq[1][i+1];

    }

    FILE *fp1 = fopen("q1.txt", "w");
    if (fp1 == NULL) {
        printf("Failed to open file q1.txt\n");
        return -1;
    }
    for (int i = 0; i < ARRAY_SIZE; i++) {
        fprintf(fp1, "%.6f\n", plant_y1[i]);
    }
    fclose(fp1);
    printf("q1.txt successfully saved\n");

    FILE *fp2 = fopen("q2.txt", "w");
    if (fp2 == NULL) {
        printf("Failed to open file q2.txt\n");
        return -1;
    }
    for (int i = 0; i < ARRAY_SIZE; i++) {
        fprintf(fp1, "%.6f\n", plant_y3[i]);
    }
    fclose(fp2);
    printf("q2.txt successfully saved\n");

    FILE *fp3 = fopen("tol1.txt", "w");
    if (fp3 == NULL) {
        printf("Failed to open file tol1.txt\n");
        return -1;
    }
    for (int i = 0; i < ARRAY_SIZE; i++) {
        fprintf(fp3, "%.6f\n", tol[0][i]);
    }
    fclose(fp3);
    printf("tol1.txt successfully saved\n");

    FILE *fp4 = fopen("tol2.txt", "w");
    if (fp4 == NULL) {
        printf("Failed to open file tol2.txt\n");
        return -1;
    }
    for (int i = 0; i < ARRAY_SIZE; i++) {
        fprintf(fp4, "%.6f\n", tol[1][i]);
    }
    fclose(fp4);
    printf("to12.txt successfully saved\n");

    return 0;

    // set title "plot of q1"
    // plot "q1.txt" with lines

}
