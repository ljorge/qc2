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
    #define Q_0_5 0.5f
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
    #define Q_FMT "f"
    #define Q_I _Complex_I

#elif defined(QC2_USE_DOUBLE)
    typedef double qfloat;
    typedef double _Complex cfloat;

    #define Q_0_0 0.0
    #define Q_0_5 0.5
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
    #define Q_FMT "lf"
    #define Q_I _Complex_I

#elif defined(QC2_USE_LONG_DOUBLE)
    typedef long double qfloat;
    typedef long double _Complex cfloat;

    #define Q_0_0 0.0L
    #define Q_0_5 0.5L
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

// Single Qubit Gates

/**
 * @brief Applies the Pauli-X gate (NOT) to a target qubit.
 *
 * Matrix:
 * | 0 1 |
 * | 1 0 |
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
 * | 0 -i |
 * | i  0 |
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
 * | 1  0 |
 * | 0 -1 |
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
 *        | 1  1 |
 * 1/√2 * | 1 -1 |
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
 * | 1 0 |
 * | 0 i |
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
 * | 1 0           |
 * | 0 exp(i*pi/4) |
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
 * | cos(theta/2)     -i*sin(theta/2) |
 * | -i*sin(theta/2)  cos(theta/2)    |
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
 * | cos(theta/2)   -sin(theta/2) |
 * | sin(theta/2)    cos(theta/2) |
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
 * | exp(-i*theta/2)      0             |
 * | 0                    exp(i*theta/2)|
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
 * | cos(theta/2)            -exp(i*lambda)*sin(theta/2)       |
 * | exp(i*phi)*sin(theta/2)  exp(i*(phi+lambda))*cos(theta/2) |
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
 * | 1 0 0 0 |
 * | 0 1 0 0 |
 * | 0 0 0 1 |
 * | 0 0 1 0 |
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
 * | 1 0 0  0 |
 * | 0 1 0  0 |
 * | 0 0 0 -i |
 * | 0 0 i  0 |
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
 * | 1 0 0  0 |
 * | 0 1 0  0 |
 * | 0 0 1  0 |
 * | 0 0 0 -1 |
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
 * | 1 0 0 0            |
 * | 0 1 0 0            |
 * | 0 0 1 0            |
 * | 0 0 0 exp(i*theta) |
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
 * | 1 0 0 0 |
 * | 0 0 1 0 |
 * | 0 1 0 0 |
 * | 0 0 0 1 |
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
 * | 1 0 0 0 0 0 0 0 |
 * | 0 1 0 0 0 0 0 0 |
 * | 0 0 1 0 0 0 0 0 |
 * | 0 0 0 1 0 0 0 0 |
 * | 0 0 0 0 1 0 0 0 |
 * | 0 0 0 0 0 1 0 0 |
 * | 0 0 0 0 0 0 0 1 |
 * | 0 0 0 0 0 0 1 0 |
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

#endif
