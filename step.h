#ifndef STEP_H
#define STEP_H


typedef struct Data {
    double* y;
    int size;
} Data;

void step(Data *output, double theta, double Ts, double t0, double t1) {
    int ARRAY_SIZE = (t1 - t0) / Ts;
    output->size = ARRAY_SIZE;
    output->y = malloc(sizeof(double) * ARRAY_SIZE);
    double t = t0;
    int i = 0;
    while (t0 <= t1 && i < ARRAY_SIZE) {
        double value = (t >= theta) ? 1.0 : 0.0;
        output->y[i] = value;
        t += Ts;
        i++;
    }
}

#endif
