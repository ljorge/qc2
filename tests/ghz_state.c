#include "../qc2.h"

// Default to 4 qubits if not specified
#define DEFAULT_QUBITS 4

int main(int argc, const char *argv[]) {
    if (!qc2_init()) {
        fprintf(stderr, "Failed to initialize QC library.\n");
        return 1;
    }

    int n_qubits = DEFAULT_QUBITS;
    if (argc > 1) {
        n_qubits = atoi(argv[1]);
        if (n_qubits < 2) n_qubits = 2;
        if (n_qubits > 30) {
            printf("Warning: %d qubits requested. This requires significant RAM. Capping at 30? No, attempting as requested.\n", n_qubits);
            // We'll let it try, but typical hardware limits exist.
        }
    }

    printf("GHZ State Simulation for %d qubits\n", n_qubits);
    printf("Expected Entanglement: |00...0> + |11...1>\n");

    clock_t start_time = clock();
    QuantumRegister *reg = create_register(n_qubits);
    if (!reg) {
        fprintf(stderr, "Failed to allocate register for %d qubits.\n", n_qubits);
        qc2_destroy();
        return 1;
    }

    // Algorithm:
    // 1. H on Q0 -> |+00...0>
    // 2. CNOT(0, 1) -> |00...0> + |110...0>
    // 3. CNOT(1, 2) -> |000...0> + |111...0>
    // ...
    // n. CNOT(n-2, n-1)

    printf("1. Applying Hadamard to Q0...\n");
    q_hadamard(reg, 0);

    printf("2. Entangling chain (CNOTs 0->1, 1->2...)\n");
    for (int i = 0; i < n_qubits - 1; i++) {
        q_cnot(reg, i, i + 1);
    }

    clock_t end_comp = clock();
    double comp_time = (double)(end_comp - start_time) / CLOCKS_PER_SEC;
    printf("State preparation took: %.4f seconds\n", comp_time);

    // Verification
    // We will measure all qubits.
    // They should ALL be 0 or ALL be 1.

    printf("3. Measuring all qubits...\n");
    int first_measure = measure(reg, 0);
    int correlated = 1;

    // We just print a summary line unless it's small
    char *result_str = malloc((size_t)n_qubits + 1);
    if(result_str) {
        result_str[0] = (char)('0' + first_measure);
        for (int i = 1; i < n_qubits; i++) {
            int m = measure(reg, i);
            result_str[i] = (char)('0' + m);
            if (m != first_measure) correlated = 0;
        }
        result_str[n_qubits] = '\0';

        if (n_qubits <= 64) {
            printf("Measurement Result: |%s>\n", result_str);
        } else {
            printf("Measurement Result: (Too long to print)\n");
        }
        free(result_str);
    } else {
         // Fallback if malloc fails
        for (int i = 1; i < n_qubits; i++) {
            int m = measure(reg, i);
            if (m != first_measure) correlated = 0;
        }
    }

    if (correlated) {
        printf("\nSUCCESS: All qubits collapsed to %d (Perfect Correlation).\n", first_measure);
        // Additional probabilistic check?
        // For GHZ, Prob(All-0) = 0.5, Prob(All-1) = 0.5.
        // But since we measured, we collapsed it.
    } else {
        printf("\nFAILURE: Qubits are NOT correlated!\n");
    }

    destroy_register(reg);
    qc2_destroy();
    return correlated ? 0 : 1;
}
