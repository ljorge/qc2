#include "../qc2.h"

// Default to 10 qubits (9 data + 1 ancilla usually, or N data for BV where ancilla is optional depending on implementation?
// Standard BV uses N query qubits + 1 ancilla. Total N+1.
// Let's implement the standard N+1 version.
#define DEFAULT_QUBITS 10

int main(int argc, const char *argv[]) {
    if (!qc2_init()) {
        fprintf(stderr, "Failed to initialize QC library.\n");
        return 1;
    }

    int n_total = DEFAULT_QUBITS;
    if (argc > 1) {
        n_total = atoi(argv[1]);
        if (n_total < 2) n_total = 2;
    }

    int n_data = n_total - 1; // Last qubit is ancilla

    // Generate a secret string 's'
    // For simplicity, let's make it alternating 1s and 0s: 101010... or random?
    // Let's make it random based on time.
    srand((unsigned int)time(NULL));

    // We'll store s as a integer mask if n_data <= 64, or char array.
    // Since we support up to 30, a 32-bit int is fine, but wait 30 is total.
    // If user asks for 30, n_data=29.
    // Let's allocate an array of ints (0 or 1) for the secret.
    int *secret_s = malloc((size_t)n_data * sizeof(int));
    printf("Bernstein-Vazirani Algorithm for %d qubits (%d data + 1 ancilla)\n", n_total, n_data);
    printf("Hidden Secret String s: ");
    for (int i = 0; i < n_data; i++) {
        secret_s[i] = rand() % 2;
        printf("%d", secret_s[i]);
    }
    printf("\n");

    clock_t start_time = clock();
    QuantumRegister *reg = create_register(n_total);
    if (!reg) {
        fprintf(stderr, "Failed to allocate register.\n");
        free(secret_s);
        qc2_destroy();
        return 1;
    }

    // 1. Initialize data qubits to |0>, ancilla to |->
    // Ancilla is at index n_data (last one).
    // |-> = H|1> = H(X|0>)
    printf("1. Initialization (Ancilla to |->)...\n");
    q_pauli_x(reg, n_data); // Ancilla -> |1>
    q_hadamard(reg, n_data); // Ancilla -> |->

    // Hadamard on all data qubits
    printf("2. H on all data qubits...\n");
    for (int i = 0; i < n_data; i++) {
        q_hadamard(reg, i);
    }

    // 2. Oracle
    // For each bit i in s, if s[i] == 1, apply CNOT(control=i, target=ancilla)
    printf("3. Oracle (Querying secret)...\n");
    for (int i = 0; i < n_data; i++) {
        if (secret_s[i] == 1) {
            q_cnot(reg, i, n_data);
        }
    }

    // 3. Interference (H on all data qubits)
    printf("4. Final H on data qubits...\n");
    for (int i = 0; i < n_data; i++) {
        q_hadamard(reg, i);
    }

    clock_t end_comp = clock();
    double comp_time = (double)(end_comp - start_time) / CLOCKS_PER_SEC;
    printf("Computation took: %.4f seconds\n", comp_time);

    // 4. Measure data qubits
    printf("5. Measurement...\n");
    int success = 1;
    printf("Measured String:      ");
    for (int i = 0; i < n_data; i++) {
        int m = measure(reg, i); // Note: measure(0) corresponds to LSB or first qubit?
        // In our print above we printed index 0 first. We should match order.
        printf("%d", m);
        if (m != secret_s[i]) {
            success = 0;
        }
    }
    printf("\n");

    if (success) {
        printf("\nSUCCESS: The hidden string was found correctly!\n");
    } else {
        printf("\nFAILURE: Measured string does not match secret!\n");
    }

    destroy_register(reg);
    free(secret_s);
    qc2_destroy();
    return success ? 0 : 1;
}
