#include "../qc2.h"

// Default to 10 qubits (N data + 1 ancilla)
// If N=10, we have 9 data qubits.
#define DEFAULT_QUBITS 10

// Oracle Type
typedef enum {
    ORACLE_CONSTANT_ZERO = 0,
    ORACLE_CONSTANT_ONE = 1,
    ORACLE_BALANCED = 2
} OracleType;

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

    // Choose oracle type randomly
    srand((unsigned int)time(NULL));
    int rt = rand() % 2; // 0 or 1
    OracleType type = (rt == 0) ? ORACLE_CONSTANT_ZERO : ORACLE_BALANCED;
    // Let's sometimes make it Constant One too?
    // Deutsch-Jozsa distinguishes Constant vs Balanced.
    // Which constant it is (0 or 1) shouldn't matter for the final measurement being |0>.

    // For rigorous testing, let's force alternating runs if we could, but here random is fine.
    // Let's actually make BALANCED more likely since it's more interesting? No, 50/50.

    printf("Deutsch-Jozsa Algorithm for %d qubits (%d data + 1 ancilla)\n", n_total, n_data);
    printf("Oracle Type: %s\n", (type == ORACLE_BALANCED) ? "BALANCED" : "CONSTANT");

    clock_t start_time = clock();
    QuantumRegister *reg = create_register(n_total);
    if (!reg) {
        fprintf(stderr, "Failed to allocate register.\n");
        qc2_destroy();
        return 1;
    }

    // 1. Initialization
    // Data qubits |0>, Ancilla |1> -> H -> |->
    printf("1. Initialization (Ancilla to |->)...\n");
    q_pauli_x(reg, n_data);
    q_hadamard(reg, n_data);

    // H on data qubits
    printf("2. H on all data qubits...\n");
    for (int i = 0; i < n_data; i++) {
        q_hadamard(reg, i);
    }

    // 2. Oracle
    printf("3. Oracle Application...\n");
    if (type == ORACLE_CONSTANT_ZERO) {
        // f(x) = 0. Do nothing?
        // Identity.
    } else if (type == ORACLE_CONSTANT_ONE) {
        // f(x) = 1.
        // Target qubit (ancilla) flips: y XOR 1.
        q_pauli_x(reg, n_data);
    } else if (type == ORACLE_BALANCED) {
        // Balanced: f(x) is 0 for half inputs, 1 for half.
        // Simplest balanced function: f(x) = x_0 (LSB).
        // CNOT(0, ancilla).
        // Or f(x) = x_0 XOR x_1 ...
        // Let's do CNOT(0, ancilla) -> Balanced.
        q_cnot(reg, 0, n_data);
    }

    // 3. Final H on data qubits
    printf("4. Final H on data qubits...\n");
    for (int i = 0; i < n_data; i++) {
        q_hadamard(reg, i);
    }

    clock_t end_comp = clock();
    double comp_time = (double)(end_comp - start_time) / CLOCKS_PER_SEC;
    printf("Computation took: %.4f seconds\n", comp_time);

    // 4. Measure data qubits
    printf("5. Measurement...\n");
    int all_zeros = 1;
    printf("Measured String: ");
    for (int i = 0; i < n_data; i++) {
        int m = measure(reg, i);
        printf("%d", m);
        if (m != 0) all_zeros = 0;
    }
    printf("\n");

    // Interpretation
    // If Constant -> All zeros
    // If Balanced -> At least one non-zero
    int is_constant_prediction = all_zeros;

    printf("Prediction: %s\n", is_constant_prediction ? "CONSTANT" : "BALANCED");

    int success = 0;
    if (type == ORACLE_BALANCED && !is_constant_prediction) success = 1;
    if ((type == ORACLE_CONSTANT_ZERO || type == ORACLE_CONSTANT_ONE) && is_constant_prediction) success = 1;

    if (success) {
        printf("\nSUCCESS: Correctly identified the function type.\n");
    } else {
        printf("\nFAILURE: Wrong identification!\n");
    }

    destroy_register(reg);
    qc2_destroy();
    return success ? 0 : 1;
}
