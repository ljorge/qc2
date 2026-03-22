#include "../qc2.h"

int tests_passed = 0;
int tests_failed = 0;

void assert_prob(const QuantumRegister *reg, int qubit, int target_val, qfloat expected_prob, const char *test_name) {
    qfloat prob = (target_val == 1) ? get_prob_one(reg, qubit) : get_prob_zero(reg, qubit);
    if (Q_ABS(prob - expected_prob) < Q_TEST_TOLERANCE) {
        printf("[PASS] %s: Prob(Q%d=%d) = %.2" Q_FMT "\n", test_name, qubit, target_val, prob);
        tests_passed++;
    } else {
        printf("[FAIL] %s: Prob(Q%d=%d) = %.2" Q_FMT " (Expected: %.2" Q_FMT ")\n", test_name, qubit, target_val, prob, expected_prob);
        tests_failed++;
    }
}

void test_single_qubit_gates(void) {
    printf("\n--- Single Qubit Gates ---\n");

    // Pauli X
    QuantumRegister *r = create_register(1);
    q_pauli_x(r, 0);
    // Pauli Y (|0> -> i|1>)
    r = create_register(1);
    q_pauli_y(r, 0);
    assert_prob(r, 0, 1, Q_1_0, "Pauli-Y (|0> -> |1>)");
    destroy_register(r);

    // Pauli Z
    // Z|0> = |0>, Z|1> = -|1>. Magnitudes shouldn't change.
    r = create_register(1);
    q_pauli_x(r, 0); // |1>
    q_pauli_z(r, 0); // -|1>
    assert_prob(r, 0, 1, Q_1_0, "Pauli-Z (|1> -> -|1>) [Check Prob]");
    destroy_register(r);

    // Hadamard
    // H|0> = |+> -> Prob(0)=Q_0_5, Prob(1)=Q_0_5
    r = create_register(1);
    q_hadamard(r, 0);
    assert_prob(r, 0, 1, Q_0_5, "Hadamard (|0> -> |+>)");
    destroy_register(r);

    // Phase (S) Gate (|0> -> |0>, |1> -> i|1>)
    // Test S on |+>: S(|+>) = (|0> + i|1>)/sqrt(2) -> Prob(1) = Q_0_5
    r = create_register(1);
    q_hadamard(r, 0);
    q_phase(r, 0);
    assert_prob(r, 0, 1, Q_0_5, "Phase (S) on |+>");
    destroy_register(r);

    // T Gate (T^4 = Z, T^2 = S)
    // T on |1> adds pi/4 phase. Probabilities shouldn't change.
    r = create_register(1);
    q_pauli_x(r, 0); // |1>
    q_t_gate(r, 0);
    assert_prob(r, 0, 1, Q_1_0, "T Gate on |1> (Check Prob)");
    destroy_register(r);
}

void test_rotations(void) {
    printf("\n--- Rotation Gates ---\n");

    // RX(PI) == X
    QuantumRegister *r = create_register(1);
    q_rx(r, 0, Q_PI);
    assert_prob(r, 0, 1, Q_1_0, "RX(PI) similar to X");
    destroy_register(r);

    r = create_register(1);
    q_ry(r, 0, Q_PI_HALF);
    assert_prob(r, 0, 1, Q_0_5, "RY(PI/2) superposition");
    destroy_register(r);

    // RZ(PI) on |+> -> |- i.e. |0> - |1> (ignoring global phase)
    r = create_register(1);
    q_hadamard(r, 0);
    q_rz(r, 0, Q_PI);
    // Apply H again: H( |0> - |1> ) = |1>
    q_hadamard(r, 0);
    assert_prob(r, 0, 1, Q_1_0, "RZ(PI) equivalent to Z (verified via H-Z-H = X)");
    destroy_register(r);

    int steps = 360 / QC2_U3_TEST_STEP_DEG;
    int step_deg = QC2_U3_TEST_STEP_DEG;

    printf("\n[INFO] U3 Sweep: Steps of %d degrees (%d steps per loop -> %d total tests)\n",
            step_deg, steps, steps*steps*steps);

    int step_count = 0;
    // Optimize: Alloc once
    r = create_register(1);

    for (int t = 0; t < steps; t++) {
        for (int p = 0; p < steps; p++) {
            for (int l = 0; l < steps; l++) {
                // Calculate angle based on index and step size
                qfloat theta = (qfloat)(t * step_deg) * (Q_PI / Q_180_0);
                qfloat phi = (qfloat)(p * step_deg) * (Q_PI / Q_180_0);
                qfloat lambda = (qfloat)(l * step_deg) * (Q_PI / Q_180_0);

                // Reset state to |0> using API to ensure sync
                qc2_reset(r);


                q_u3(r, 0, theta, phi, lambda);

                qfloat expected_prob_1 = Q_SIN(theta / Q_2_0) * Q_SIN(theta / Q_2_0);
                qfloat prob = get_prob_one(r, 0);

                if (Q_ABS(prob - expected_prob_1) > Q_TEST_TOLERANCE) {
                    printf("[FAIL] U3(%.2" Q_FMT ", %.2" Q_FMT ", %.2" Q_FMT "): P(1)=%.2" Q_FMT " (Exp: %.2" Q_FMT ")\n",
                            theta, phi, lambda, prob, expected_prob_1);
                    tests_failed++;
                } else {
                    step_count++;
                }
            }
        }
    }
    // Free once
    destroy_register(r);
    tests_passed += step_count;
    printf("[INFO] U3 Parametric Sweep: %d combinations verified.\n", step_count);
}

void test_multi_qubit_gates(void) {
    printf("\n--- Multi Qubit Gates ---\n");

    // CNOT
    QuantumRegister *r = create_register(2);
    q_pauli_x(r, 0); // Control = 1
    q_cnot(r, 0, 1); // Target flips to 1
    assert_prob(r, 1, 1, Q_1_0, "CNOT (10 -> 11)");
    destroy_register(r);

    // CZ Gate
    // Control |1>, Target |+> -> |->
    // Initial: |1> |+> = |10> + |11>
    // After CZ: |10> - |11>
    // Apply H to target: |1> |1>
    r = create_register(2);
    q_pauli_x(r, 0); // Control = 1
    q_hadamard(r, 1); // Target = +
    q_cz(r, 0, 1);
    q_hadamard(r, 1); // Should become 1
    assert_prob(r, 1, 1, Q_1_0, "CZ (1+ -> 1-) checked via H");
    destroy_register(r);

    // SWAP
    r = create_register(2);
    // |01> (Q0=1, Q1=0 in our indexing? Wait, internal is bitmasks. Let's assume index 0 is first bit)
    q_pauli_x(r, 0);
    // Actually in my code: bit 0 is 1ULL << 0. So index 0 is LSB usually.
    // Let's verify SWAP by putting Q0=1, Q1=0.
    q_swap(r, 0, 1);
    // Now Q0 should be 0, Q1 should be 1.
    assert_prob(r, 0, 0, Q_1_0, "SWAP Q0 (1->0)");
    assert_prob(r, 1, 1, Q_1_0, "SWAP Q1 (0->1)");
    destroy_register(r);

    // Controlled-Y
    // Control |1>, Target |0> -> Target becomes Y|0> = i|1>
    r = create_register(2);
    q_pauli_x(r, 0); // Control = 1
    q_cy(r, 0, 1);
    assert_prob(r, 1, 1, Q_1_0, "CY (10 -> 11) [i|1>]");
    destroy_register(r);

    // CCNOT (Toffoli)
    r = create_register(3);
    q_pauli_x(r, 0);
    q_pauli_x(r, 1);
    q_ccnot(r, 0, 1, 2); // 110 -> 111
    assert_prob(r, 2, 1, Q_1_0, "Toffoli (110 -> 111)");
    destroy_register(r);
}

// Helper to print metrics and configuration
void print_test_metrics(clock_t start_time, clock_t end_time) {
    double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("\n--- Configuration Flags ---\n");
    #if defined(QC2_COMPILER_NAME)

        printf("Compiler: %s\n", QC2_COMPILER_NAME);
    #else
        printf("Compiler: Unknown\n");
    #endif

    // Precision
    #if defined(QC2_USE_FLOAT)
        printf("QC2_USE: FLOAT\n");
    #elif defined(QC2_USE_DOUBLE)
        printf("QC2_USE: DOUBLE\n");
    #elif defined(QC2_USE_LONG_DOUBLE)
        printf("QC2_USE: LONG_DOUBLE\n");
    #endif
    printf("sizeof(qfloat): %zu\n", sizeof(qfloat));

    // OpenMP
    #if defined(QC2_USE_OPENMP)
        printf("QC2_USE_OPENMP: Yes\n");
    #else
        printf("QC2_USE_OPENMP: No\n");
    #endif

    // OpenCL
    #if defined(QC2_USE_OPENCL)
        printf("QC2_USE_OPENCL: Yes\n");
        #if defined(CL_TARGET_OPENCL_VERSION)
            printf("CL_TARGET_OPENCL_VERSION: %d\n", CL_TARGET_OPENCL_VERSION);
        #endif
    #else
        printf("QC2_USE_OPENCL: No\n");
    #endif


    // Aligned Memory
    #if defined(QC2_ALIGNED_MEMORY)
        printf("QC2_ALIGNED_MEMORY: Yes\n");
    #else
        printf("QC2_ALIGNED_MEMORY: No\n");
    #endif

    // Trig Tables
    #if defined(QC2_USE_TRIG_TABLES)

        printf("QC2_USE_TRIG_TABLES: Yes\n");
    #else
        printf("QC2_USE_TRIG_TABLES: No\n");
    #endif

    // U3 Step
    printf("QC2_U3_TEST_STEP_DEG: %d\n", QC2_U3_TEST_STEP_DEG);

    printf("\n--- Test Suite Metrics ---\n");
    printf("Execution Time: %.4f seconds\n", time_spent);

    printf("--------------------------\n");
}

// ============================================================================
// ============================================================================
// New Tests - Extended Gates
// ============================================================================
// ============================================================================

void test_new_single_qubit_gates(void) {
    printf("\n--- New Single Qubit Gates ---\n");

    // Identity gate
    QuantumRegister *r = create_register(1);
    q_identity(r, 0);
    assert_prob(r, 0, 0, Q_1_0, "Identity on |0>");
    q_pauli_x(r, 0);
    q_identity(r, 0);
    assert_prob(r, 0, 1, Q_1_0, "Identity on |1>");
    destroy_register(r);

    // sqrt(X) - Square root of NOT
    r = create_register(1);
    q_sqrt_x(r, 0);
    q_sqrt_x(r, 0);
    assert_prob(r, 0, 1, Q_1_0, "sqrt(X)^2 = X");
    destroy_register(r);

    // sqrt(Z) - Square root of Z
    r = create_register(1);
    q_sqrt_z(r, 0);
    q_sqrt_z(r, 0);
    assert_prob(r, 0, 0, Q_1_0, "sqrt(Z)^2 = Z");
    destroy_register(r);

    // U1 gate
    r = create_register(1);
    q_u1(r, 0, Q_PI);
    assert_prob(r, 0, 0, Q_1_0, "U1(PI) = Z");
    destroy_register(r);

    // P and P dagger gates
    r = create_register(1);
    q_hadamard(r, 0);
    q_phase(r, 0);
    assert_prob(r, 0, 1, Q_0_5, "P on |+>");
    // q_p_dagger
    r = create_register(1);
    q_hadamard(r, 0);
    q_p_dagger(r, 0);
    assert_prob(r, 0, 1, Q_0_5, "P_dagger on |+>");
    destroy_register(r);

    // q_u2 - testing U2(0, 0)
    r = create_register(1);
    q_u2(r, 0, Q_0_0, Q_0_0);
    assert_prob(r, 0, 1, Q_0_5, "U2(0,0) acts as Hadamard-like superposition");
    destroy_register(r);

    // Walsh-Hadamard - verify uniform superposition |++>
    r = create_register(2);
    q_walsh(r);
    assert_prob(r, 0, 1, Q_0_5, "Walsh(2) - Q0 prob");
    assert_prob(r, 1, 1, Q_0_5, "Walsh(2) - Q1 prob");
    destroy_register(r);
}

void test_controlled_gates_new(void) {
    printf("\n--- New Controlled Gates ---\n");

    // CH - Controlled Hadamard
    QuantumRegister *r = create_register(2);
    q_pauli_x(r, 0);
    q_ch(r, 0, 1);
    assert_prob(r, 1, 1, Q_0_5, "CH on |10>");
    destroy_register(r);

    // CS - Controlled Phase
    r = create_register(2);
    q_pauli_x(r, 0);
    q_pauli_x(r, 1);
    q_cs(r, 0, 1);
    q_cs(r, 0, 1);
    q_cs(r, 0, 1);
    q_cs(r, 0, 1);
    assert_prob(r, 1, 1, Q_1_0, "CS^4 = I");
    destroy_register(r);

    // CT - Controlled T
    r = create_register(2);
    q_pauli_x(r, 0);
    q_pauli_x(r, 1);
    q_ct(r, 0, 1);
    q_ct(r, 0, 1);
    q_ct(r, 0, 1);
    q_ct(r, 0, 1);
    assert_prob(r, 1, 1, Q_1_0, "CT^4 = I");
    destroy_register(r);

    // CRX, CRY, CRZ
    r = create_register(2);
    q_pauli_x(r, 0);
    q_crx(r, 0, 1, Q_PI);
    assert_prob(r, 1, 1, Q_1_0, "CRX(PI)");
    destroy_register(r);

    r = create_register(2);
    q_pauli_x(r, 0);
    q_cry(r, 0, 1, Q_PI);
    assert_prob(r, 1, 1, Q_1_0, "CRY(PI)");
    destroy_register(r);

    // CRZ
    r = create_register(2);
    q_pauli_x(r, 0);
    q_hadamard(r, 1);
    q_crz(r, 0, 1, Q_PI);
    q_hadamard(r, 1);
    assert_prob(r, 1, 1, Q_1_0, "CRZ(PI)");
    destroy_register(r);

    // CU1
    r = create_register(2);
    q_pauli_x(r, 0); // Control = 1
    q_hadamard(r, 1); // Target = |+>
    q_cu1(r, 0, 1, Q_PI); // Phase flip -> Target = |->
    q_hadamard(r, 1); // Target = |1>
    assert_prob(r, 1, 1, Q_1_0, "CU1(PI) acts as CZ");
    destroy_register(r);

    // CZ_EQUIV
    r = create_register(2);
    q_pauli_x(r, 0); // Control = 1
    q_hadamard(r, 1);
    q_cz_equiv(r, 0, 1);
    q_hadamard(r, 1);
    assert_prob(r, 1, 1, Q_1_0, "CZ_EQUIV acts as CZ");
    destroy_register(r);
}

void test_two_qubit_rotation_gates(void) {
    printf("\n--- Two-Qubit Rotation Gates ---\n");

    // iSWAP
    QuantumRegister *r = create_register(2);
    q_pauli_x(r, 0); // |10>
    q_iswap(r, 0, 1); // -> i|01>
    assert_prob(r, 0, 0, Q_1_0, "iSWAP Q0 (1->0)");
    assert_prob(r, 1, 1, Q_1_0, "iSWAP Q1 (0->1)");
    destroy_register(r);

    // sqrt(CNOT)
    r = create_register(2);
    q_pauli_x(r, 0);
    q_sqrt_cnot(r, 0, 1);
    q_sqrt_cnot(r, 0, 1);
    assert_prob(r, 1, 1, Q_1_0, "sqrt(CNOT)^2 = CNOT");
    destroy_register(r);

    // RXX - test that RXX(PI) flips |00> to |11>
    r = create_register(2);
    q_rxx(r, 0, 1, Q_PI);
    assert_prob(r, 0, 1, Q_1_0, "RXX(PI) on |00> flips to |11>");
    destroy_register(r);

    // RYY - test that RYY(PI) flips |00> to |11>
    r = create_register(2);
    q_ryy(r, 0, 1, Q_PI);
    assert_prob(r, 0, 1, Q_1_0, "RYY(PI) on |00> flips to |11>");
    destroy_register(r);

    // RZZ - test that RZZ(PI) doesn't change |00>
    r = create_register(2);
    q_rzz(r, 0, 1, Q_PI);
    assert_prob(r, 0, 0, Q_1_0, "RZZ(PI) on |00>");
    destroy_register(r);

    // ECR - test
    r = create_register(2);
    q_ecr(r, 0, 1);
    assert_prob(r, 0, 1, Q_0_5, "ECR on |00>");
    destroy_register(r);
}

void test_three_qubit_gates_new(void) {
    printf("\n--- New Three-Qubit Gates ---\n");

    // Fredkin (CSWAP) - test with control=0 (no change)
    QuantumRegister *r = create_register(3);
    q_pauli_x(r, 2);
    q_fredkin(r, 0, 1, 2);
    assert_prob(r, 2, 1, Q_1_0, "Fredkin control=0");
    destroy_register(r);

    // CCZ - verify phase flip on |111>
    r = create_register(3);
    q_pauli_x(r, 0);
    q_pauli_x(r, 1);
    q_pauli_x(r, 2);
    q_ccz(r, 0, 1, 2);
    q_ccz(r, 0, 1, 2);
    assert_prob(r, 2, 1, Q_1_0, "CCZ^2 = I");
    destroy_register(r);

    // CCP - verify phase flip on |111> with PI
    r = create_register(3);
    q_pauli_x(r, 0);
    q_pauli_x(r, 1);
    q_hadamard(r, 2);
    q_ccp(r, 0, 1, 2, Q_PI);
    q_hadamard(r, 2);
    assert_prob(r, 2, 1, Q_1_0, "CCP(PI) acts as CCZ");
    destroy_register(r);

    // Toffoli alias
    r = create_register(3);
    q_pauli_x(r, 0);
    q_pauli_x(r, 1);
    q_ccnot(r, 0, 1, 2);
    assert_prob(r, 2, 1, Q_1_0, "Toffoli alias");
    destroy_register(r);

    // CCH - verify H applied when control is 11
    r = create_register(3);
    q_pauli_x(r, 0);
    q_pauli_x(r, 1);
    q_cch(r, 0, 1, 2);
    q_hadamard(r, 2);
    assert_prob(r, 2, 0, Q_1_0, "CCH on |110> maps 0 to + and then H maps to 0");
    destroy_register(r);
}

void test_utilities_new(void) {
    printf("\n--- Utility Functions ---\n");

    // q_normalize
    QuantumRegister *r = create_register(1);
    r->amplitudes[0] = Q_2_0;
    q_normalize(r);
    assert_prob(r, 0, 0, Q_1_0, "Normalize");
    destroy_register(r);

    // q_hadamard_row
    r = create_register(2);
    q_hadamard_row(r, 0, 2);
    assert_prob(r, 0, 1, Q_0_5, "Hadamard Row Q0");
    assert_prob(r, 1, 1, Q_0_5, "Hadamard Row Q1");
    destroy_register(r);

    // q_copy
    r = create_register(1);
    QuantumRegister *r_copy = q_copy(r);
    assert_prob(r_copy, 0, 0, Q_1_0, "Copy register");
    destroy_register(r);
    destroy_register(r_copy);

    // q_fidelity
    r = create_register(1);
    q_hadamard(r, 0);
    QuantumRegister *r2 = create_register(1);
    q_hadamard(r2, 0);
    qfloat fid = q_fidelity(r, r2);
    if (Q_ABS(fid - Q_1_0) < Q_TEST_TOLERANCE) {
        printf("[PASS] Fidelity of identical states: %.2" Q_FMT "\n", fid);
        tests_passed++;
    } else {
        printf("[FAIL] Fidelity: %.2" Q_FMT " (Expected: 1.0)\n", fid);
        tests_failed++;
    }
    destroy_register(r);
    destroy_register(r2);

    // q_measure_partial
    r = create_register(2);
    q_hadamard(r, 0);
    int qubits[1] = {0};
    qfloat *probs = q_measure_partial(r, qubits, 1);
    if (probs && Q_ABS(probs[0] - Q_0_5) < Q_TEST_TOLERANCE) {
        printf("[PASS] Measure partial P(0): %.2" Q_FMT "\n", probs[0]);
        tests_passed++;
    } else {
        printf("[FAIL] Measure partial\n");
        tests_failed++;
    }
    free(probs);
    destroy_register(r);
}

int main(void) {

    clock_t start_time = clock();

    if (!qc2_init()) {

        fprintf(stderr, "Failed to initialize QC library.\n");
        return 1;
    }
    printf("Running QC Test Suite...");
    printf("\n--------------------------\n");

    test_single_qubit_gates();
    test_rotations();
    test_multi_qubit_gates();

    // New tests
    test_new_single_qubit_gates();
    test_controlled_gates_new();
    test_two_qubit_rotation_gates();
    test_three_qubit_gates_new();
    test_utilities_new();

    printf("\n--------------------------\n");
    printf("Tests Passed: %d\n", tests_passed);
    printf("Tests Failed: %d\n", tests_failed);

    qc2_destroy();

    clock_t end_time = clock();
    print_test_metrics(start_time, end_time);

    return (tests_failed == 0) ? 0 : 1;
}
