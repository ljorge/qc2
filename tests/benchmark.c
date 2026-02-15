#include "../qc2.h"

// Simple timer
double get_time_sec(void) {
    return (double)clock() / CLOCKS_PER_SEC;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage: %s <num_qubits>\n", argv[0]);
        return 1;
    }

    int n_qubits = atoi(argv[1]);
    if (n_qubits < 1 || n_qubits > 30) {
        printf("Invalid number of qubits: %d (Must be 1-30)\n", n_qubits);
        return 1;
    }

    printf("========================================\n");
    printf(" QC Library Benchmark\n");
    printf("========================================\n");
    printf("Qubits: %d\n", n_qubits);
    size_t dim = 1ULL << n_qubits;
    printf("State Vector Size: %zu complex amplitudes\n", dim);
    printf("Memory: %.2f MB\n", (double)(dim * sizeof(cfloat)) / (1024.0 * 1024.0));

    // Config Info
    #if defined(QC2_USE_OPENCL)
        printf("Backend: OpenCL (GPU)\n");
    #elif defined(QC2_USE_OPENMP)
        printf("Backend: OpenMP (CPU Parallel)\n");
    #else
        printf("Backend: Single Thread (CPU)\n");
    #endif

    printf("Initializing Register... ");
    fflush(stdout);

    if (!qc2_init()) {
        fprintf(stderr, "Failed to init lib.\n");
        return 1;
    }

    double t0 = get_time_sec();
    QuantumRegister *reg = create_register(n_qubits);
    double t_init = get_time_sec() - t0;
    if (!reg) {
        fprintf(stderr, "Failed to allocate memory.\n");
        return 1;
    }
    printf("Done (%.4fs)\n", t_init);

    // --- BENCHMARK ---
    // Apply a sequence of gates.
    // For fair benchmarking, we want O(N) gates or O(1) multi-qubit gates.
    // Let's do a layer of Hadamards (superposition) followed by some entangling gates.

    printf("Running Gates...\n");
    t0 = get_time_sec();

    // 1. Layer of Hadamards (N gates)
    for (int i = 0; i < n_qubits; i++) {
        q_hadamard(reg, i);
    }

    // 2. Chain of CNOTs (N-1 gates)
    for (int i = 0; i < n_qubits - 1; i++) {
        q_cnot(reg, i, i + 1);
    }

    // 3. Some rotations (N gates)
    for (int i = 0; i < n_qubits; i++) {
        q_rx(reg, i, 0.5f);
    }

    // Force synchronization if on GPU to include transfer time/execution completion
    #if defined(QC2_USE_OPENCL)
    // We can force a read to ensure sync
    get_prob_one(reg, 0);
    #endif

    double t_run = get_time_sec() - t0;
    int total_gates = n_qubits + (n_qubits - 1) + n_qubits;

    printf("----------------------------------------\n");
    printf("Total Gates: %d\n", total_gates);
    printf("Execution Time: %.4f seconds\n", t_run);
    printf("Throughput: %.2f gates/sec\n", (double)total_gates / t_run);
    printf("========================================\n");

    qc2_destroy();
    return 0;
}
