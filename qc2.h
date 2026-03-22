// QC v2.0
// This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public License along with this library; if not, see <https://www.gnu.org/licenses/>.
#if !defined(_QC2_H_)
#define _QC2_H_

#if defined(__clang__)
#pragma clang diagnostic ignored "-Wgnu-imaginary-constant"
#endif

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#if defined(_WIN32)
#include <malloc.h>
#endif
#if defined(QC2_USE_OPENMP)
#include <omp.h>
#endif
#if defined(QC2_USE_OPENCL)
#include <CL/cl.h>
#endif

// ***************************************************************************
// Floating point precision settings
#if !defined(QC2_USE_FLOAT) && !defined(QC2_USE_DOUBLE) && !defined(QC2_USE_LONG_DOUBLE)
    #error "Precision not defined. Please define QC2_USE_FLOAT, QC2_USE_DOUBLE, or QC2_USE_LONG_DOUBLE."
#endif

// Use OpenMP for parallel processing
//#define QC2_USE_OPENMP

// Use OpenCL for GPU acceleration
//#define QC2_USE_OPENCL

#if defined(QC2_USE_OPENCL) && defined(QC2_USE_LONG_DOUBLE)
    #error "OpenCL does not support 'long double' precision. Please use FLOAT or DOUBLE."
#endif

// Use Aligned Memory (default for AVX/SIMD)
#if !defined(QC2_NO_ALIGNED_MEMORY)
    #define QC2_ALIGNED_MEMORY
#endif

// Use Trig Tables
//#define QC2_USE_TRIG_TABLES

// U3 Nested Sweep
#if !defined(QC2_U3_TEST_STEP_DEG)
#define QC2_U3_TEST_STEP_DEG 15
#endif
// ***************************************************************************

#if defined(QC2_USE_FLOAT)
    typedef float qfloat;
    typedef float _Complex cfloat;

    #define Q_0_0 0.0f
    #define Q_0_25 0.25f
    #define Q_0_5 0.5f
    #define Q_0_75 0.75f
    #define Q_1_0 1.0f
    #define Q_2_0 2.0f
    #define Q_4_0 4.0f
    #define Q_180_0 180.0f
    #define Q_TEST_TOLERANCE 0.01f
    #define Q_PRINT_THRESHOLD 0.0001f

    #define Q_SIN sinf
    #define Q_COS cosf
    #define Q_SQRT sqrtf
    #define Q_FMOD fmodf
    #define Q_CABS cabsf
    #define Q_ABS fabsf
    #define Q_CREAL crealf
    #define Q_CIMAG cimagf
    #define Q_CEXP cexpf
    #define Q_PI 3.14159265358979323846f
    #define Q_E 2.71828182845904523536f
    #define Q_FMT "f"
    #define Q_I _Complex_I

#elif defined(QC2_USE_DOUBLE)
    typedef double qfloat;
    typedef double _Complex cfloat;

    #define Q_0_0 0.0
    #define Q_0_25 0.25
    #define Q_0_5 0.5
    #define Q_0_75 0.75
    #define Q_1_0 1.0
    #define Q_2_0 2.0
    #define Q_4_0 4.0
    #define Q_180_0 180.0
    #define Q_TEST_TOLERANCE 0.01
    #define Q_PRINT_THRESHOLD 0.0001

    #define Q_SIN sin
    #define Q_COS cos
    #define Q_SQRT sqrt
    #define Q_FMOD fmod
    #define Q_CABS cabs
    #define Q_ABS fabs
    #define Q_CREAL creal
    #define Q_CIMAG cimag
    #define Q_CEXP cexp
    #define Q_PI 3.14159265358979323846
    #define Q_E 2.71828182845904523536
    #define Q_FMT "lf"
    #define Q_I _Complex_I

#elif defined(QC2_USE_LONG_DOUBLE)
    typedef long double qfloat;
    typedef long double _Complex cfloat;

    #define Q_0_0 0.0L
    #define Q_0_25 0.25L
    #define Q_0_5 0.5L
    #define Q_0_75 0.75L
    #define Q_1_0 1.0L
    #define Q_2_0 2.0L
    #define Q_4_0 4.0L
    #define Q_180_0 180.0L
    #define Q_TEST_TOLERANCE 0.01L
    #define Q_PRINT_THRESHOLD 0.0001L

    #define Q_SIN sinl
    #define Q_COS cosl
    #define Q_SQRT sqrtl
    #define Q_FMOD fmodl
    #define Q_CABS cabsl
    #define Q_ABS fabsl
    #define Q_CREAL creall
    #define Q_CIMAG cimagl
    #define Q_CEXP cexpl
    #define Q_PI 3.14159265358979323846264338327950288419716939937510L
    #define Q_E 2.71828182845904523536028747135266249775724709369995L
    #define Q_FMT "Lf"
    #define Q_I _Complex_I
#else
    #error "No precision defined. Please define QC2_USE_FLOAT, QC2_USE_DOUBLE, or QC2_USE_LONG_DOUBLE."
#endif

// Constants
#define Q_PI_HALF (Q_PI / Q_2_0)
#define Q_PI_QUARTER (Q_PI / Q_4_0)
#define Q_PI_DOUBLE (Q_PI * Q_2_0)

// If defined QC2_USE_TRIG_TABLES, use tables instead of Q_SIN and Q_COS
#if defined(QC2_USE_TRIG_TABLES)
#define Q_TRIG_TABLE_SIZE 360000
#endif

/**
 * @brief Represents the state of a quantum system with N qubits.
 *
 * The state is stored as a vector of 2^N complex amplitudes.
 * When OpenCL is enabled, it also manages the device memory and synchronization state.
 */
typedef struct {
    int n_qubits;       /**< Number of qubits in the register */
    size_t dim;         /**< Dimension of the state vector (2^n_qubits) */
    cfloat *amplitudes; /**< Array of complex amplitudes on the Host */
#if defined(QC2_USE_OPENCL)
    cl_mem device_amplitudes; /**< Array of amplitudes on the Device (OpenCL) */
    int device_dirty;   /**< Flag indicating if device memory has newer data than host */
#endif
} QuantumRegister;

/**
 * @brief Seeds the random number generator.
 *
 * @param seed The seed value.
 */
void qc2_seed(unsigned int seed);

/**
 * @brief Initializes the QC library.
 *
 * Ideally called once at the start of the program.
 * Initializes internal tables (trig tables) and OpenCL context if enabled.
 *
 * @return 1 on success, 0 on failure.
 */
int qc2_init(void);

/**
 * @brief Cleans up resources used by the QC library.
 *
 * Frees internal tables and OpenCL context.
 */
void qc2_destroy(void);

/**
 * @brief Creates a new quantum register.
 *
 * Allocates memory for the state vector and initializes it to |0...0>.
 *
 * @param n_qubits Number of qubits in the register.
 * @return Pointer to the new QuantumRegister, or NULL on failure.
 */
QuantumRegister* create_register(int n_qubits);

/**
 * @brief Destroys a quantum register.
 *
 * Frees the associated memory on host and device.
 *
 * @param reg Pointer to the QuantumRegister to destroy.
 */
void destroy_register(QuantumRegister *reg);

/**
 * @brief Resets the quantum register to the state |0...0>.
 *
 * @param reg Pointer to the QuantumRegister.
 */
void qc2_reset(QuantumRegister *reg);

/**
 * @brief Prints the current state vector to stdout.
 *
 * Only prints non-zero amplitudes (above Q_PRINT_THRESHOLD).
 *
 * @param reg Pointer to the QuantumRegister.
 */
void print_state(QuantumRegister *reg);

// Gates

/**
 * @brief Applies the Pauli-X gate (NOT) to a target qubit.
 *
 * Matrix:
 * │ 0  1 │
 * │ 1  0 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @return 1 on success.
 */
int q_pauli_x(QuantumRegister *reg, int target_qubit);

/**
 * @brief Applies the Pauli-Y gate to a target qubit.
 *
 * Matrix:
 * │ 0  -i │
 * │ i   0 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @return 1 on success.
 */
int q_pauli_y(QuantumRegister *reg, int target_qubit);

/**
 * @brief Applies the Pauli-Z gate to a target qubit.
 *
 * Matrix:
 * │ 1   0 │
 * │ 0  -1 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @return 1 on success.
 */
int q_pauli_z(QuantumRegister *reg, int target_qubit);

/**
 * @brief Applies the Hadamard gate to a target qubit.
 *
 * Creates superposition.
 * Matrix:
 * │ 1/√2   1/√2 │
 * │ 1/√2  -1/√2 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @return 1 on success.
 */
int q_hadamard(QuantumRegister *reg, int target_qubit);

/**
 * @brief Applies the Phase gate (S gate) to a target qubit.
 *
 * Rotates phase by PI/2 around Z axis.
 * Matrix:
 * │ 1  0 │
 * │ 0  i │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @return 1 on success.
 */
int q_phase(QuantumRegister *reg, int target_qubit); // S gate

/**
 * @brief Applies the T gate to a target qubit.
 *
 * Rotates phase by PI/4 around Z axis.
 * Matrix:
 * │ 1      0    │
 * │ 0  e^(iπ/4) │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @return 1 on success.
 */
int q_t_gate(QuantumRegister *reg, int target_qubit); // T gate

// Rotation Gates

/**
 * @brief Applies a rotation around the X axis.
 *
 * Matrix:
 * │   cos(θ/2)   -i·sin(θ/2) │
 * │ -i·sin(θ/2)    cos(θ/2)  │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @param theta Rotation angle in radians.
 * @return 1 on success.
 */
int q_rx(QuantumRegister *reg, int target_qubit, qfloat theta);

/**
 * @brief Applies a rotation around the Y axis.
 *
 * Matrix:
 * │ cos(θ/2)  -sin(θ/2) │
 * │ sin(θ/2)   cos(θ/2) │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @param theta Rotation angle in radians.
 * @return 1 on success.
 */
int q_ry(QuantumRegister *reg, int target_qubit, qfloat theta);

/**
 * @brief Applies a rotation around the Z axis.
 *
 * Matrix:
 * │ e^(-iθ/2)      0    │
 * │     0      e^(iθ/2) │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @param theta Rotation angle in radians.
 * @return 1 on success.
 */
int q_rz(QuantumRegister *reg, int target_qubit, qfloat theta);

// Generic Gate

/**
 * @brief Applies a universal single-qubit rotation gate U3.
 *
 * Matrix:
 * │    cos(θ/2)        -e^(iλ)·sin(θ/2)  │
 * │ e^(iφ)·sin(θ/2)  e^(i(φ+λ))·cos(θ/2) │
 *
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @param theta Rotation angle.
 * @param phi Phase parameter.
 * @param lambda Phase parameter.
 * @return 1 on success.
 */
int q_u3(QuantumRegister *reg, int target_qubit, qfloat theta, qfloat phi, qfloat lambda);

// Two Qubit Gates

/**
 * @brief Applies the Controlled-NOT (CNOT) gate.
 *
 * Flips the target qubit if the control qubit is |1>.
 * Matrix:
 * │ 1  0  0  0 │
 * │ 0  1  0  0 │
 * │ 0  0  0  1 │
 * │ 0  0  1  0 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target Index of the target qubit.
 * @return 1 on success.
 */
int q_cnot(QuantumRegister *reg, int control, int target);

/**
 * @brief Applies the Controlled-Y gate.
 *
 * Applies Y to target if control is |1>.
 * Matrix:
 * │ 1  0  0   0 │
 * │ 0  1  0   0 │
 * │ 0  0  0  -i │
 * │ 0  0  i   0 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target Index of the target qubit.
 * @return 1 on success.
 */
int q_cy(QuantumRegister *reg, int control, int target);

/**
 * @brief Applies the Controlled-Z gate.
 *
 * Applies Z to target if control is |1>. Equivalent to phase flip if both are |1>.
 * Matrix:
 * │ 1  0  0  0 │
 * │ 0  1  0  0 │
 * │ 0  0  1  0 │
 * │ 0  0  0 -1 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target Index of the target qubit.
 * @return 1 on success.
 */
int q_cz(QuantumRegister *reg, int control, int target);

/**
 * @brief Applies the Controlled-Phase gate.
 *
 * Applies Phase(theta) to target if control is |1>.
 * Matrix:
 * │ 1  0  0    0    │
 * │ 0  1  0    0    │
 * │ 0  0  1    0    │
 * │ 0  0  0  e^(iθ) │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target Index of the target qubit.
 * @param theta Phase angle.
 * @return 1 on success.
 */
int q_cp(QuantumRegister *reg, int control, int target, qfloat theta);

/**
 * @brief Swaps the states of two qubits.
 *
 * Matrix:
 * │ 1  0  0  0 │
 * │ 0  0  1  0 │
 * │ 0  1  0  0 │
 * │ 0  0  0  1 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param qubit1 Index of the first qubit.
 * @param qubit2 Index of the second qubit.
 * @return 1 on success.
 */
int q_swap(QuantumRegister *reg, int qubit1, int qubit2);

// Three Qubit Gates

/**
 * @brief Applies the Toffoli gate (CCNOT).
 *
 * Flips the target qubit if both control qubits are |1>.
 * Matrix (8x8):
 * │ 1  0  0  0  0  0  0  0 │
 * │ 0  1  0  0  0  0  0  0 │
 * │ 0  0  1  0  0  0  0  0 │
 * │ 0  0  0  1  0  0  0  0 │
 * │ 0  0  0  0  1  0  0  0 │
 * │ 0  0  0  0  0  1  0  0 │
 * │ 0  0  0  0  0  0  0  1 │
 * │ 0  0  0  0  0  0  1  0 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control1 Index of the first control qubit.
 * @param control2 Index of the second control qubit.
 * @param target Index of the target qubit.
 * @return 1 on success.
 */
int q_ccnot(QuantumRegister *reg, int control1, int control2, int target); // Toffoli

// Measurement

/**
 * @brief Measures a specific qubit, collapsing the state vector.
 *
 * The state becomes consistent with the measurement result.
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the qubit to measure.
 * @return The result of the measurement (0 or 1).
 */
int measure(QuantumRegister *reg, int target_qubit);

// Utilities

/**
 * @brief Calculates the probability of measuring 0 for a qubit.
 *
 * Does NOT collapse the state vector.
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the qubit.
 * @return Probability [0.0, 1.0].
 */
qfloat get_prob_zero(const QuantumRegister *reg, int target_qubit);

/**
 * @brief Calculates the probability of measuring 1 for a qubit.
 *
 * Does NOT collapse the state vector.
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the qubit.
 * @return Probability [0.0, 1.0].
 */
qfloat get_prob_one(const QuantumRegister *reg, int target_qubit);

/**
 * @brief Applies the Identity gate (I) to a target qubit.
 *
 * Matrix:
 * │ 1  0 │
 * │ 0  1 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @return 1 on success.
 */
int q_identity(QuantumRegister *reg, int target_qubit);

/**
 * @brief Applies the Square-Root of X gate (√X) to a target qubit.
 *
 * Also known as sqrt(NOT). Creates a superposition between |0> and |1>
 * with a relative phase.
 * Matrix:
 * │ (1+i)/2  (1-i)/2 │
 * │ (1-i)/2  (1+i)/2 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @return 1 on success.
 */
int q_sqrt_x(QuantumRegister *reg, int target_qubit);

/**
 * @brief Applies the Square-Root of Z gate (√Z) to a target qubit.
 *
 * Matrix:
 * │ 1        0    │
 * │ 0  e^(iπ/4) │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @return 1 on success.
 */
int q_sqrt_z(QuantumRegister *reg, int target_qubit);

/**
 * @brief Applies the U1 gate (phase rotation) to a target qubit.
 *
 * Also known as Rz(λ) up to a global phase.
 * Matrix:
 * │ 1     0   │
 * │ 0  e^(iλ) │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @param lambda Phase angle in radians.
 * @return 1 on success.
 */
int q_u1(QuantumRegister *reg, int target_qubit, qfloat lambda);

/**
 * @brief Applies the U2 gate to a target qubit.
 *
 * Matrix:
 * │   1/√2         -e^(iλ)/√2  │
 * │ e^(iφ)/√2    e^(i(φ+λ))/√2 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @param phi Phase parameter in radians.
 * @param lambda Phase parameter in radians.
 * @return 1 on success.
 */
int q_u2(QuantumRegister *reg, int target_qubit, qfloat phi, qfloat lambda);

/**
 * @brief Applies the Phase gate (P) to a target qubit.
 *
 * Equivalent to T^n where n is an integer. Rotates phase by π/2 around Z axis.
 * Matrix:
 * │ 1  0 │
 * │ 0  i │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @return 1 on success.
 */
int q_p(QuantumRegister *reg, int target_qubit);

/**
 * @brief Applies the Phase dagger gate (P) to a target qubit.
 *
 * Inverse of Phase gate. Rotates phase by -π/2 around Z axis.
 * Matrix:
 * │ 1   0 │
 * │ 0  -i │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param target_qubit Index of the target qubit.
 * @return 1 on success.
 */
int q_p_dagger(QuantumRegister *reg, int target_qubit);

/**
 * @brief Applies the Hadamard gate to a row of qubits.
 *
 * Applies H to all qubits from first_qubit to first_qubit+num_qubits-1.
 * Useful for creating uniform superposition states.
 *
 * @param reg Pointer to the QuantumRegister.
 * @param first_qubit First qubit index.
 * @param num_qubits Number of qubits to apply Hadamard.
 * @return 1 on success.
 */
int q_hadamard_row(QuantumRegister *reg, int first_qubit, int num_qubits);

/**
 * @brief Applies the Walsh-Hadamard transform.
 *
 * Applies Hadamard to ALL qubits in the register.
 * Creates the uniform superposition state |+...+>.
 *
 * @param reg Pointer to the QuantumRegister.
 * @return 1 on success.
 */
int q_walsh(QuantumRegister *reg);

/**
 * @brief Applies the Controlled-Hadamard (CH) gate.
 *
 * Applies Hadamard to target if control is |1>.
 * Matrix:
 * │ 1    0    0      0  │
 * │ 0    1    0      0  │
 * │ 0    0  1/√2   1/√2 │
 * │ 0    0  1/√2  -1/√2 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target Index of the target qubit.
 * @return 1 on success.
 */
int q_ch(QuantumRegister *reg, int control, int target);

/**
 * @brief Applies the Controlled-Phase (CS) gate.
 *
 * Applies S gate to target if control is |1>.
 * Matrix:
 * │ 1  0  0  0 │
 * │ 0  1  0  0 │
 * │ 0  0  1  0 │
 * │ 0  0  0  i │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target Index of the target qubit.
 * @return 1 on success.
 */
int q_cs(QuantumRegister *reg, int control, int target);

/**
 * @brief Applies the Controlled-T (CT) gate.
 *
 * Applies T gate to target if control is |1>.
 * Matrix:
 * │ 1  0  0      0    │
 * │ 0  1  0      0    │
 * │ 0  0  1      0    │
 * │ 0  0  0  e^(iπ/4) │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target Index of the target qubit.
 * @return 1 on success.
 */
int q_ct(QuantumRegister *reg, int control, int target);

/**
 * @brief Applies the Controlled-RX (CRX) gate.
 *
 * Applies RX(theta) to target if control is |1>.
 * Matrix:
 * │ 1  0       0           0      │
 * │ 0  1       0           0      │
 * │ 0  0    cos(θ/2)  -i·sin(θ/2) │
 * │ 0  0  -i·sin(θ/2)   cos(θ/2)  │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target Index of the target qubit.
 * @param theta Rotation angle in radians.
 * @return 1 on success.
 */
int q_crx(QuantumRegister *reg, int control, int target, qfloat theta);

/**
 * @brief Applies the Controlled-RY (CRY) gate.
 *
 * Applies RY(theta) to target if control is |1>.
 * Matrix:
 * │ 1  0     0          0     │
 * │ 0  1     0          0     │
 * │ 0  0  cos(θ/2)  -sin(θ/2) │
 * │ 0  0  sin(θ/2)   cos(θ/2) │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target Index of the target qubit.
 * @param theta Rotation angle in radians.
 * @return 1 on success.
 */
int q_cry(QuantumRegister *reg, int control, int target, qfloat theta);

/**
 * @brief Applies the Controlled-RZ (CRZ) gate.
 *
 * Applies RZ(theta) to target if control is |1>.
 * Matrix:
 * │ 1  0      0          0    │
 * │ 0  1      0          0    │
 * │ 0  0  e^(-iθ/2)      0    │
 * │ 0  0      0      e^(iθ/2) │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target Index of the target qubit.
 * @param theta Rotation angle in radians.
 * @return 1 on success.
 */
int q_crz(QuantumRegister *reg, int control, int target, qfloat theta);

/**
 * @brief Applies the Controlled-U1 (CU1) gate.
 *
 * Applies U1(lambda) to target if control is |1>.
 * Matrix:
 * │ 1  0  0    0    │
 * │ 0  1  0    0    │
 * │ 0  0  1    0    │
 * │ 0  0  0  e^(iλ) │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target Index of the target qubit.
 * @param lambda Phase angle in radians.
 * @return 1 on success.
 */
int q_cu1(QuantumRegister *reg, int control, int target, qfloat lambda);

/**
 * @brief Applies the iSWAP gate.
 *
 * Swaps qubits with a phase factor of i.
 * Matrix:
 * │ 1  0  0  0 │
 * │ 0  0  i  0 │
 * │ 0  i  0  0 │
 * │ 0  0  0  1 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param qubit1 Index of the first qubit.
 * @param qubit2 Index of the second qubit.
 * @return 1 on success.
 */
int q_iswap(QuantumRegister *reg, int qubit1, int qubit2);

/**
 * @brief Applies the Square-Root of CNOT (CNOT) gate.
 *
 * Also known as sqrt(CX). Useful for creating entangled states
 * from product states.
 * Matrix:
 * │ 1  0     0        0    │
 * │ 0  1     0        0    │
 * │ 0  0  (1+i)/2  (1-i)/2 │
 * │ 0  0  (1-i)/2  (1+i)/2 │
 *
 * (Matriz completa: √CNOT = (I ⊗ H) · CZ · (I ⊗ H))
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target Index of the target qubit.
 * @return 1 on success.
 */
int q_sqrt_cnot(QuantumRegister *reg, int control, int target);

/**
 * @brief Applies the RXX gate (rotation around XX axis).
 *
 * Implements exp(-i*theta/2 * X⊗X).
 * Matrix:
 * │   cos(θ/2)        0            0       -i·sin(θ/2) │
 * │      0         cos(θ/2)   -i·sin(θ/2)       0      │
 * │      0        -i·sin(θ/2)   cos(θ/2)        0      │
 * │ -i·sin(θ/2)       0            0         cos(θ/2)  │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param qubit1 Index of the first qubit.
 * @param qubit2 Index of the second qubit.
 * @param theta Rotation angle in radians.
 * @return 1 on success.
 */
int q_rxx(QuantumRegister *reg, int qubit1, int qubit2, qfloat theta);

/**
 * @brief Applies the RYY gate (rotation around YY axis).
 *
 * Implements exp(-i*theta/2 * Y⊗Y).
 * Matrix:
 * │  cos(θ/2)       0           0        i·sin(θ/2) │
 * │    0         cos(θ/2)   -i·sin(θ/2)       0     │
 * │    0        -i·sin(θ/2)   cos(θ/2)        0     │
 * │ i·sin(θ/2)      0           0         cos(θ/2)  │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param qubit1 Index of the first qubit.
 * @param qubit2 Index of the second qubit.
 * @param theta Rotation angle in radians.
 * @return 1 on success.
 */
int q_ryy(QuantumRegister *reg, int qubit1, int qubit2, qfloat theta);

/**
 * @brief Applies the RZZ gate (rotation around ZZ axis).
 *
 * Implements exp(-i*theta/2 * Z⊗Z).
 * Matrix:
 * │ e^(-iθ/2)    0        0         0     │
 * │    0      e^(iθ/2)    0         0     │
 * │    0         0     e^(iθ/2)     0     │
 * │    0         0        0     e^(-iθ/2) │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param qubit1 Index of the first qubit.
 * @param qubit2 Index of the second qubit.
 * @param theta Rotation angle in radians.
 * @return 1 on success.
 */
int q_rzz(QuantumRegister *reg, int qubit1, int qubit2, qfloat theta);

/**
 * @brief Applies the ECR (echoed RZX) gate.
 *
 * Equivalent to RZX(π/2) up to single-qubit rotations.
 * Matrix:
 * │      0           0      1/√2   i/√2 │
 * │      0           0      i/√2   1/√2 │
 * │    1/√2      -i/√2         0      0 │
 * │   -i/√2       1/√2         0      0 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param qubit1 Index of the first qubit.
 * @param qubit2 Index of the second qubit.
 * @return 1 on success.
 */
int q_ecr(QuantumRegister *reg, int qubit1, int qubit2);

/**
 * @brief Applies CZ equivalent gate.
 *
 * Applies CZ gate (equivalent to q_cz).
 * Matrix:
 * │ 1  0  0  0 │
 * │ 0  1  0  0 │
 * │ 0  0  1  0 │
 * │ 0  0  0 -1 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target Index of the target qubit.
 * @return 1 on success.
 */
int q_cz_equiv(QuantumRegister *reg, int control, int target);

/**
 * @brief Applies the Fredkin gate (CSWAP).
 *
 * Swaps target1 and target2 if control is |1>.
 * Matrix (8x8):
 * │ 1  0  0  0  0  0  0  0 │
 * │ 0  1  0  0  0  0  0  0 │
 * │ 0  0  1  0  0  0  0  0 │
 * │ 0  0  0  1  0  0  0  0 │
 * │ 0  0  0  0  1  0  0  0 │
 * │ 0  0  0  0  0  0  1  0 │
 * │ 0  0  0  0  0  1  0  0 │
 * │ 0  0  0  0  0  0  0  1 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control Index of the control qubit.
 * @param target1 Index of the first target qubit.
 * @param target2 Index of the second target qubit.
 * @return 1 on success.
 */
int q_fredkin(QuantumRegister *reg, int control, int target1, int target2);

/**
 * @brief Applies the Controlled-Controlled-Hadamard (CCH) gate.
 *
 * Applies Hadamard to target if both controls are |1>.
 * Matrix (8x8):
 * │ 1  0  0  0  0  0    0    0   │
 * │ 0  1  0  0  0  0    0    0   │
 * │ 0  0  1  0  0  0    0    0   │
 * │ 0  0  0  1  0  0    0    0   │
 * │ 0  0  0  0  1  0    0    0   │
 * │ 0  0  0  0  0  1    0    0   │
 * │ 0  0  0  0  0  0  1/√2  1/√2 │
 * │ 0  0  0  0  0  0  1/√2 -1/√2 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control1 Index of the first control qubit.
 * @param control2 Index of the second control qubit.
 * @param target Index of the target qubit.
 * @return 1 on success.
 */
int q_cch(QuantumRegister *reg, int control1, int control2, int target);

/**
 * @brief Applies the Controlled-Controlled-Z (CCZ) gate.
 *
 * Applies Z to target if both controls are |1>.
 * Applies a phase of -1 to the |111> state.
 * Matrix (8x8):
 * │ 1  0  0  0  0  0  0  0 │
 * │ 0  1  0  0  0  0  0  0 │
 * │ 0  0  1  0  0  0  0  0 │
 * │ 0  0  0  1  0  0  0  0 │
 * │ 0  0  0  0  1  0  0  0 │
 * │ 0  0  0  0  0  1  0  0 │
 * │ 0  0  0  0  0  0  1  0 │
 * │ 0  0  0  0  0  0  0 -1 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control1 Index of the first control qubit.
 * @param control2 Index of the second control qubit.
 * @param target Index of the target qubit.
 * @return 1 on success.
 */
int q_ccz(QuantumRegister *reg, int control1, int control2, int target);

/**
 * @brief Applies the Controlled-Controlled-Phase (CCP) gate.
 *
 * Applies phase exp(i*theta) to target if both controls are |1>.
 * Matrix (8x8):
 * │ 1  0  0  0  0  0  0   0    │
 * │ 0  1  0  0  0  0  0   0    │
 * │ 0  0  1  0  0  0  0   0    │
 * │ 0  0  0  1  0  0  0   0    │
 * │ 0  0  0  0  1  0  0   0    │
 * │ 0  0  0  0  0  1  0   0    │
 * │ 0  0  0  0  0  0  1   0    │
 * │ 0  0  0  0  0  0  0 e^(iθ) │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control1 Index of the first control qubit.
 * @param control2 Index of the second control qubit.
 * @param target Index of the target qubit.
 * @param theta Phase angle in radians.
 * @return 1 on success.
 */
int q_ccp(QuantumRegister *reg, int control1, int control2, int target, qfloat theta);

/**
 * @brief Applies the Toffoli gate - alias for compatibility.
 *
 * Flips the target qubit if both control qubits are |1>.
 * Matrix (8x8):
 * │ 1  0  0  0  0  0  0  0 │
 * │ 0  1  0  0  0  0  0  0 │
 * │ 0  0  1  0  0  0  0  0 │
 * │ 0  0  0  1  0  0  0  0 │
 * │ 0  0  0  0  1  0  0  0 │
 * │ 0  0  0  0  0  1  0  0 │
 * │ 0  0  0  0  0  0  0  1 │
 * │ 0  0  0  0  0  0  1  0 │
 *
 * @param reg Pointer to the QuantumRegister.
 * @param control1 Index of the first control qubit.
 * @param control2 Index of the second control qubit.
 * @param target Index of the target qubit.
 * @return 1 on success.
 */
int q_toffoli(QuantumRegister *reg, int control1, int control2, int target);

// Advanced Measurements

/**
 * @brief Measures multiple qubits.
 *
 * Collapses the state vector for all specified qubits.
 *
 * @param reg Pointer to the QuantumRegister.
 * @param qubits Array of qubit indices to measure.
 * @param num_qubits Number of qubits to measure.
 * @return Pointer to array of measurement results (0 or 1), or NULL on failure.
 *         Caller is responsible for freeing the returned array.
 */
int* q_measure_multiple(QuantumRegister *reg, const int *qubits, int num_qubits);

/**
 * @brief Gets probabilities for multiple qubits without collapsing.
 *
 * Returns the probability of measuring 0 for each specified qubit.
 *
 * @param reg Pointer to the QuantumRegister.
 * @param qubits Array of qubit indices.
 * @param num_qubits Number of qubits.
 * @return Pointer to array of probabilities, or NULL on failure.
 *         Caller is responsible for freeing the returned array.
 */
qfloat* q_measure_partial(const QuantumRegister *reg, const int *qubits, int num_qubits);

/**
 * @brief Measures a qubit in an arbitrary basis.
 *
 * Measures in basis X, Y, or Z.
 * Basis: 0=X, 1=Y, 2=Z
 *
 * @param reg Pointer to the QuantumRegister.
 * @param qubit Index of the qubit to measure.
 * @param basis Basis to measure in (0=X, 1=Y, 2=Z).
 * @return Measurement result (0 or 1).
 */
int q_measure_basis(QuantumRegister *reg, int qubit, int basis);

// Utilities

/**
 * @brief Normalizes the quantum register state vector.
 *
 * Ensures the state vector has unit norm.
 *
 * @param reg Pointer to the QuantumRegister.
 * @return 1 on success, 0 on failure.
 */
int q_normalize(QuantumRegister *reg);

/**
 * @brief Calculates the fidelity between two quantum states.
 *
 * Fidelity F = |<ψ|φ>|^2 for pure states.
 *
 * @param reg1 Pointer to the first QuantumRegister.
 * @param reg2 Pointer to the second QuantumRegister.
 * @return Fidelity value [0.0, 1.0].
 */
qfloat q_fidelity(const QuantumRegister *reg1, const QuantumRegister *reg2);

/**
 * @brief Computes the partial trace of a quantum register.
 *
 * Traces out specified qubits to get reduced density matrix probabilities.
 *
 * @param reg Pointer to the QuantumRegister.
 * @param qubits_to_trace Array of qubit indices to trace out.
 * @param num_qubits Number of qubits to trace out.
 * @param remaining_qubits Output array for remaining qubit indices.
 * @param num_remaining Number of remaining qubits.
 * @return Pointer to array of probabilities for each basis state of remaining qubits.
 *         Caller is responsible for freeing the returned array.
 */
qfloat* q_partial_trace(const QuantumRegister *reg, const int *qubits_to_trace, 
                        int num_qubits, const int *remaining_qubits, int num_remaining);

/**
 * @brief Creates a copy of a quantum register.
 *
 * @param reg Pointer to the QuantumRegister to copy.
 * @return Pointer to new QuantumRegister, or NULL on failure.
 */
QuantumRegister* q_copy(const QuantumRegister *reg);

// Aliases

#define q_i                 q_identity
#define q_swap_not          q_sqrt_x
#define q_cswap             q_fredkin
#define q_sqrt_cx           q_sqrt_cnot
#define q_toffoli_gate      q_toffoli
#define q_p                 q_phase
#define q_cz_equiv          q_cz
#define q_toffoli           q_ccnot

#endif
