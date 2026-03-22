#include "../qc2.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Parameters
#define N_QUBITS 5
#define FEATURE_DIM 20
#define K_NEIGHBORS 5

typedef struct {
    int id;
    float vector[FEATURE_DIM];
} TrainSample;

typedef struct {
    int true_label;
    float vector[FEATURE_DIM];
} TestSample;

typedef struct {
    int id;
    float score;
} Candidate;

// Globals
TrainSample *train_data = NULL;
int num_train_samples = 0;

// Helper: Encode vector into Quantum Register using Amplitude Encoding
void encode_data(QuantumRegister *reg, const float *vector) {
    if (!reg) return;

    double norm = 0.0;
    for (int i = 0; i < FEATURE_DIM; i++) {
        norm += vector[i] * vector[i];
    }
    norm = sqrt(norm);
    if (norm < 1e-9) norm = 1.0;

    // Direct Amplitude Loading (Simulation only)
    memset(reg->amplitudes, 0, reg->dim * sizeof(cfloat));
    for (size_t i = 0; i < FEATURE_DIM; i++) {
        if (i < reg->dim) {
            reg->amplitudes[i] = (qfloat)(vector[i] / norm);
        }
    }
}

// Calculate Fidelity |<psi|phi>|^2
float calculate_fidelity_registers_fast(const QuantumRegister *reg1, const QuantumRegister *reg2) {
    if (reg1->dim != reg2->dim) return 0.0f;
    cfloat dot = 0;
    // Sequential loop here because this function is called inside a parallel loop
    for (size_t i = 0; i < reg1->dim; i++) {
        // conj(a) * b
        qfloat re1 = Q_CREAL(reg1->amplitudes[i]);
        qfloat im1 = Q_CIMAG(reg1->amplitudes[i]);
        qfloat re2 = Q_CREAL(reg2->amplitudes[i]);
        qfloat im2 = Q_CIMAG(reg2->amplitudes[i]);
        qfloat re = re1*re2 + im1*im2;
        qfloat im = re1*im2 - im1*re2;
        dot += re + Q_I * im;
    }
    qfloat mag = Q_CABS(dot);
    return (float)(mag * mag);
}

int load_training_set(const char *filename) {
    FILE *f = fopen(filename, "rb");
    if (!f) return 0;
    if (fread(&num_train_samples, sizeof(int), 1, f) != 1) { fclose(f); return 0; }

    printf("Loading %d training samples...\n", num_train_samples);
    train_data = (TrainSample *)malloc(sizeof(TrainSample) * (size_t)num_train_samples);

    for (int i = 0; i < num_train_samples; i++) {
        if (fread(&train_data[i].id, sizeof(int), 1, f) != 1) break;
        if (fread(train_data[i].vector, sizeof(float), FEATURE_DIM, f) != FEATURE_DIM) break;
    }
    fclose(f);
    return 1;
}

TestSample* load_test_samples(const char *filename, int *count) {
    FILE *f = fopen(filename, "rb");
    if (!f) return NULL;
    int n;
    if (fread(&n, sizeof(int), 1, f) != 1) { fclose(f); return NULL; }
    TestSample *samples = (TestSample *)malloc(sizeof(TestSample) * (size_t)n);
    for (int i = 0; i < n; i++) {
        if (fread(&samples[i].true_label, sizeof(int), 1, f) != 1) break;
        if (fread(samples[i].vector, sizeof(float), FEATURE_DIM, f) != FEATURE_DIM) break;
    }
    *count = n;
    fclose(f);
    printf("Loaded %d test samples.\n", n);
    return samples;
}

int compare_candidates(const void *a, const void *b) {
    const Candidate *cA = (const Candidate *)a;
    const Candidate *cB = (const Candidate *)b;
    if (cB->score > cA->score) return 1;
    if (cB->score < cA->score) return -1;
    return 0;
}

// Full Quantum k-NN Classification
int classify_knn(const QuantumRegister *reg_input) {
    Candidate *results = (Candidate*)malloc(sizeof(Candidate) * (size_t)num_train_samples);

    // Parallel loop over all training samples
    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel
    #endif
    {
        // Each thread gets its own temp register to avoid reallocation
        QuantumRegister *reg_temp = create_register(N_QUBITS);

        #if defined(QC2_USE_OPENMP)
        #pragma omp for
        #endif
        for (int i = 0; i < num_train_samples; i++) {
            // Encode training sample
            encode_data(reg_temp, train_data[i].vector);

            // Calculate fidelity
            float fid = calculate_fidelity_registers_fast(reg_input, reg_temp);

            results[i].id = train_data[i].id;
            results[i].score = fid;
        }

        destroy_register(reg_temp);
    }

    // Sort by Fidelity (Descending)
    qsort(results, (size_t)num_train_samples, sizeof(Candidate), compare_candidates);

    // Majority Vote
    // Map ID -> Vote Count
    // IDs are ASCII (0-255). Simple array is enough.
    int votes[256] = {0};
    float scores[256] = {0.0f}; // Tie-breaker

    int k = K_NEIGHBORS;
    if (k > num_train_samples) k = num_train_samples;

    for (int i = 0; i < k; i++) {
        int id = results[i].id;
        votes[id]++;
        scores[id] += results[i].score;
    }

    // Find Winner
    int winner = -1;
    int max_votes = -1;
    float max_score = -1.0f;

    for (int i = 0; i < 256; i++) {
        if (votes[i] > max_votes) {
            max_votes = votes[i];
            winner = i;
            max_score = scores[i];
        } else if (votes[i] == max_votes && max_votes > 0) {
            // Tie-break by total fidelity score
            if (scores[i] > max_score) {
                winner = i;
                max_score = scores[i];
            }
        }
    }

    free(results);
    return winner;
}

int main(void) {
    printf("=== Real Dataset QML (Using QC Library: Quantum k-NN / Amplitude Encoding) ===\n");

    if (!qc2_init()) {
        fprintf(stderr, "Failed to init QC lib.\n");
        return 1;
    }

    if (!load_training_set("train_samples.bin")) {
        fprintf(stderr, "Failed to load training set. Run preprocess_qml.py first.\n");
        return 1;
    }

    int num_samples = 0;
    TestSample *samples = load_test_samples("test_samples.bin", &num_samples);
    if (!samples) return 1;

    QuantumRegister *reg_input = create_register(N_QUBITS);
    if (!reg_input) return 1;

    int correct = 0;
    int limit_test = 20; // Test more samples now that k-NN is robust
    if (limit_test > num_samples) limit_test = num_samples;

    printf("Running classification on first %d samples...\n", limit_test);

    for (int s = 0; s < limit_test; s++) {
        printf("\n--- Sample %d (True: '%c') ---\n", s+1, (char)samples[s].true_label);

        // Encode Input
        encode_data(reg_input, samples[s].vector);

        // Classify
        int predicted_id = classify_knn(reg_input);

        printf("Predicted: '%c'\n", (char)predicted_id);
        if (predicted_id == samples[s].true_label) {
            printf("[PASS]\n");
            correct++;
        } else {
            printf("[FAIL]\n");
        }
    }

    printf("\nAccuracy: %d / %d (%.2f%%)\n", correct, limit_test, ((float)correct)/((float)limit_test) * 100.0f);

    destroy_register(reg_input);
    free(train_data);
    free(samples);
    qc2_destroy();

    return 0;
}
