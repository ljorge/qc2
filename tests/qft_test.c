#include "../qc2.h"

// Default to 8 qubits. QFT is N^2 deep.
#define DEFAULT_QUBITS 8

// Helper: Apply QFT on 'n' qubits starting at 'start_idx'
// (Though typically we apply to whole register)
void apply_qft(QuantumRegister *reg, int n) {
    for (int i = 0; i < n; i++) {
        // H on qubit i
        q_hadamard(reg, i);

        // Controlled Phase Rotations from j > i
        // R_k where k = j - i + 1
        // R_k adds phase e^(i * 2pi / 2^k) if both 1.
        // In our library, `q_cp(control, target, theta)` adds phase theta/2 to 11?
        // Wait, standard CP(theta) adds phase e^(i theta) to |11>.
        // We need phase 2*pi / 2^k = pi / 2^(k-1).

        for (int j = i + 1; j < n; j++) {
            int k = j - i + 1;
            qfloat angle = Q_PI / (qfloat)(1 << (k - 1));
            // e.g. k=2 (j=i+1): angle = PI/2.

            // Note: Standard QFT diagrams often have Control on j, Target on i, or vice versa.
            // CP is symmetric.
            q_cp(reg, j, i, angle);
        }
    }

    // Swap qubits to reverse order (standard QFT output is bit-reversed)
    // We'll skip the SWAP for now if we just want to verify QFT^-1 * QFT = I
    // IF we apply QFT^-1 also with skipped swaps?
    // Actually, let's do the swaps to be "Textbook QFT".
    for (int i = 0; i < n / 2; i++) {
        q_swap(reg, i, n - 1 - i);
    }
}

void apply_inverse_qft(QuantumRegister *reg, int n) {
    // Reverse Swaps
    for (int i = 0; i < n / 2; i++) {
        q_swap(reg, i, n - 1 - i);
    }

    // Reverse Gates in Reverse Order
    for (int i = n - 1; i >= 0; i--) {
        for (int j = n - 1; j > i; j--) {
            int k = j - i + 1;
            qfloat angle = -Q_PI / (qfloat)(1 << (k - 1)); // Negative angle
            q_cp(reg, j, i, angle);
        }
        q_hadamard(reg, i);
    }
}

int main(int argc, const char *argv[]) {
    if (!qc2_init()) {
        fprintf(stderr, "Failed to initialize QC library.\n");
        return 1;
    }

    int n_qubits = DEFAULT_QUBITS;
    if (argc > 1) {
        n_qubits = atoi(argv[1]);
        if (n_qubits < 2) n_qubits = 2;
    }

    printf("QFT Identity Test (QFT^-1 * QFT) for %d qubits\n", n_qubits);

    QuantumRegister *reg = create_register(n_qubits);

    // 1. Prepare a random state or a specific superposition
    // Let's create a superposition to ensure we aren't just testing |00..0>
    printf("1. initializing: H on all qubits...\n");
    for (int i = 0; i < n_qubits; i++) {
        q_hadamard(reg, i);
    }
    // And maybe some random phases? T gate on even qubits
    for (int i = 0; i < n_qubits; i+=2) {
        q_t_gate(reg, i);
    }

    // Snapshot: We can't easily snapshot full state without copying.
    // But since this is a simulation, we trust valid math.
    // We could measure probability of |0> vs |1>?
    // Instead, let's trust that if we get back to initial, it works.
    // Wait, testing "Identity" on a known state is better.
    // Let's settle for: Start |0>, Apply X to Q0. State |100...0> (index 1).
    // Apply QFT -> Fourier Basis.
    // Apply InvQFT -> Should assume |100...0>.
    // Measure -> Should be 100...0 with prob 1.

    qc2_reset(reg);
    printf("1. Set Input State: |101...0> (Alternating)\n");
    // Let's set an alternating pattern: 1010...
    for(int i=0; i<n_qubits; i+=2) {
        q_pauli_x(reg, i);
    }
    // Store expected bits
    int *expected = malloc((size_t)n_qubits * sizeof(int));
    for(int i=0; i<n_qubits; i++) expected[i] = (i % 2 == 0) ? 1 : 0;

    clock_t start_time = clock();

    printf("2. Applying QFT...\n");
    apply_qft(reg, n_qubits);

    printf("3. Applying Inverse QFT...\n");
    apply_inverse_qft(reg, n_qubits);

    clock_t end_comp = clock();
    double comp_time = (double)(end_comp - start_time) / CLOCKS_PER_SEC;
    printf("Transformations took: %.4f seconds\n", comp_time);

    // 4. Measure
    printf("4. Measurement...\n");
    int success = 1;
    printf("Result:   ");
    for(int i=0; i<n_qubits; i++) {
        int m = measure(reg, i);
        printf("%d", m);
        if (m != expected[i]) success = 0;
    }
    printf("\nExpected: ");
    for(int i=0; i<n_qubits; i++) printf("%d", expected[i]);
    printf("\n");

    if (success) {
        printf("\nSUCCESS: QFT^-1 * QFT = Identity.\n");
    } else {
        printf("\nFAILURE: Did not return to initial state.\n");
    }

    free(expected);
    destroy_register(reg);
    qc2_destroy();
    return success ? 0 : 1;
}
