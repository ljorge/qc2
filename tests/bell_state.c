#include "../qc2.h"

int main(void) {
    if (!qc2_init()) {
        fprintf(stderr, "Failed to initialize QC library.\n");
        return 1;
    }
    printf("Bell State Simulation (Correlated Qubits)\n");

    // Create a register with 2 qubits
    // Initial state: |00>
    QuantumRegister *reg = create_register(2);

    // Apply Hadamard to Qubit 0
    // State: (|00> + |10>) / sqrt(2) (Bit ordering dependent, let's assume standard tensor product)
    printf("Applying Hadamard to Q0...\n");
    q_hadamard(reg, 0);

    // Apply CNOT with Control Q0, Target Q1
    // If Q0 is 0, Q1 stays 0 -> |00>
    // If Q0 is 1, Q1 flips to 1 -> |11>
    // State: (|00> + |11>) / sqrt(2)
    printf("Applying CNOT (Control Q0, Target Q1)...\n");
    q_cnot(reg, 0, 1);

    print_state(reg);

    printf("\nProbabilities before measurement:\n");
    printf("P(Q0=1): %.2" Q_FMT "\n", get_prob_one(reg, 0));
    printf("P(Q1=1): %.2" Q_FMT "\n", get_prob_one(reg, 1));

    // Measure Q0
    printf("\nMeasuring Q0...\n");
    int m0 = measure(reg, 0);
    printf("Result Q0: %d\n", m0);

    print_state(reg);

    // Measure Q1
    printf("Measuring Q1...\n");
    int m1 = measure(reg, 1);
    printf("Result Q1: %d\n", m1);

    if (m0 == m1) {
        printf("\nSUCCESS: Measurements are correlated!\n");
    } else {
        printf("\nFAILURE: Measurements are NOT correlated (this should not happen for Bell state)!\n");
    }

    destroy_register(reg);
    qc2_destroy();
    return 0;
}
