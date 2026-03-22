// QC v2.0
// This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public License along with this library; if not, see <https://www.gnu.org/licenses/>.
#include "qc2_internal.h"

#ifdef QC2_DEBUG
static int is_unitary(cfloat u00, cfloat u01, cfloat u10, cfloat u11) {
    // Check U * U^dagger = I
    qfloat mag0 = Q_CABS(u00)*Q_CABS(u00) + Q_CABS(u01)*Q_CABS(u01);
    qfloat mag1 = Q_CABS(u10)*Q_CABS(u10) + Q_CABS(u11)*Q_CABS(u11);

    if (Q_ABS(mag0 - Q_1_0) > Q_TEST_TOLERANCE) return 0;
    if (Q_ABS(mag1 - Q_1_0) > Q_TEST_TOLERANCE) return 0;

    // Check orthogonality: row0 . row1* = 0
    // row0 = (u00, u01), row1 = (u10, u11)
    // dot = u00 * conj(u10) + u01 * conj(u11)
    qfloat dot_real = Q_CREAL(u00) * Q_CREAL(u10) + Q_CIMAG(u00) * Q_CIMAG(u10) 
                    + Q_CREAL(u01) * Q_CREAL(u11) + Q_CIMAG(u01) * Q_CIMAG(u11);
    qfloat dot_imag = -Q_CREAL(u00) * Q_CIMAG(u10) + Q_CIMAG(u00) * Q_CREAL(u10) 
                    - Q_CREAL(u01) * Q_CIMAG(u11) + Q_CIMAG(u01) * Q_CREAL(u11);

    if (Q_ABS(dot_real) > Q_TEST_TOLERANCE) return 0;
    if (Q_ABS(dot_imag) > Q_TEST_TOLERANCE) return 0;

    return 1;
}
#endif

/**
 * @brief Applies a general single-qubit gate defined by a 2x2 matrix.
 *
 * Optimized with single loop and bit-insertion for O(N) parallelism.
 *
 * @param reg Quantum Register
 * @param target_qubit Target qubit index
 * @param u00 Matrix element (0,0)
 * @param u01 Matrix element (0,1)
 * @param u10 Matrix element (1,0)
 * @param u11 Matrix element (1,1)
 * @return 1 on success
 */
static int apply_1q_gate(QuantumRegister *reg, int target_qubit, cfloat u00, cfloat u01, cfloat u10, cfloat u11) {
    #ifdef QC2_DEBUG
    assert(is_unitary(u00, u01, u10, u11) && "Gate matrix is not unitary");
    #endif

    if (!reg) return 0;
    if (target_qubit >= reg->n_qubits) return 0;

    #if defined(QC2_USE_OPENCL)
    if (reg->device_amplitudes) {
        return apply_gate_opencl(reg, target_qubit, u00, u01, u10, u11);
    }
    #endif

    size_t dim = reg->dim;
    size_t half_dim = dim >> 1;
    size_t target_bit = 1ULL << target_qubit;
    size_t mask_low = target_bit - 1;
    size_t mask_high = ~mask_low;

    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < half_dim; i++) {
        // Insert zero at target_qubit position.
        size_t idx0 = (i & mask_low) | ((i & mask_high) << 1);
        size_t idx1 = idx0 | target_bit;

        cfloat a = reg->amplitudes[idx0];
        cfloat b = reg->amplitudes[idx1];

        reg->amplitudes[idx0] = u00 * a + u01 * b;
        reg->amplitudes[idx1] = u10 * a + u11 * b;
    }
    return 1;
}

static int apply_2q_unitary(QuantumRegister *reg, int qubit1, int qubit2,
                            cfloat u00, cfloat u01, cfloat u02, cfloat u03,
                            cfloat u10, cfloat u11, cfloat u12, cfloat u13,
                            cfloat u20, cfloat u21, cfloat u22, cfloat u23,
                            cfloat u30, cfloat u31, cfloat u32, cfloat u33) {
    #ifdef QC2_DEBUG
    assert(0 && "apply_2q_unitary not fully implemented for debug");
    #endif

    if (!reg) return 0;
    if (qubit1 >= reg->n_qubits || qubit2 >= reg->n_qubits) return 0;
    if (qubit1 == qubit2) return 0;

    #if defined(QC2_USE_OPENCL)
    if (reg->device_amplitudes) {
        sync_to_host(reg);
    }
    #endif

    size_t dim = reg->dim;
    size_t quarter_dim = dim >> 2;

    int min_q = (qubit1 < qubit2) ? qubit1 : qubit2;
    int max_q = (qubit1 > qubit2) ? qubit1 : qubit2;

    size_t min_bit = 1ULL << min_q;
    size_t max_bit = 1ULL << max_q;
    size_t mask_low = min_bit - 1;

    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < quarter_dim; i++) {
        size_t t1 = (i & mask_low) | ((i & ~mask_low) << 1);
        size_t idx_base = (t1 & (max_bit - 1)) | ((t1 & ~(max_bit - 1)) << 1);

        size_t idx0 = idx_base;
        size_t idx1 = idx_base | min_bit;
        size_t idx2 = idx_base | max_bit;
        size_t idx3 = idx_base | min_bit | max_bit;

        cfloat a0 = reg->amplitudes[idx0];
        cfloat a1 = reg->amplitudes[idx1];
        cfloat a2 = reg->amplitudes[idx2];
        cfloat a3 = reg->amplitudes[idx3];

        reg->amplitudes[idx0] = u00 * a0 + u01 * a1 + u02 * a2 + u03 * a3;
        reg->amplitudes[idx1] = u10 * a0 + u11 * a1 + u12 * a2 + u13 * a3;
        reg->amplitudes[idx2] = u20 * a0 + u21 * a1 + u22 * a2 + u23 * a3;
        reg->amplitudes[idx3] = u30 * a0 + u31 * a1 + u32 * a2 + u33 * a3;
    }

    #if defined(QC2_USE_OPENCL)
    if (reg->device_amplitudes) {
        sync_to_device(reg);
    }
    #endif

    return 1;
}

/**
 * @brief Applies a controlled 2x2 matrix.
 *
 * Iterates only over states where control bit is 1.
 *
 * @param reg Quantum Register
 * @param control Control qubit index
 * @param target Target qubit index
 * @param u00 Matrix element (0,0)
 * @param u01 Matrix element (0,1)
 * @param u10 Matrix element (1,0)
 * @param u11 Matrix element (1,1)
 * @return 1 on success
 */
static int apply_controlled_gate(QuantumRegister *reg, int control, int target, cfloat u00, cfloat u01, cfloat u10, cfloat u11) {
    #ifdef QC2_DEBUG
    assert(is_unitary(u00, u01, u10, u11) && "Controlled Gate matrix is not unitary");
    #endif
    if (!reg) return 0;
    if (control >= reg->n_qubits || target >= reg->n_qubits) return 0;
    if (control == target) return 0;

    #if defined(QC2_USE_OPENCL)
    if (reg->device_amplitudes) {
        return apply_controlled_gate_opencl(reg, control, target, u00, u01, u10, u11);
    }
    #endif

    size_t dim = reg->dim;
    size_t quarter_dim = dim >> 2; // Total iterations = N/4 (states where C=1, T=0)

    // Calculate masks for double bit insertion
    int min_q = (control < target) ? control : target;
    int max_q = (control > target) ? control : target;

    size_t min_bit = 1ULL << min_q;
    size_t max_bit = 1ULL << max_q;

    // Masks to insert gaps at min_q and max_q
    size_t mask_low = min_bit - 1;
    // size_t mask_mid = (max_bit - 1) & ~mask_low;
    // size_t mask_high = ~(max_bit - 1);

    // Because we iterate only where Control=1, we need to know if Control is the Min or Max bit to determine insertion value.
    size_t control_mask = (1ULL << control);
    size_t target_mask  = (1ULL << target);

    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < quarter_dim; i++) {
        // We take 'i' (which has N-2 bits) and insert:
        // 1 at 'control' position
        // 0 at 'target' position

        // 1. Split i at min_q
        size_t t1 = (i & mask_low) | ((i & ~mask_low) << 1);
        // 2. Split result at max_q
        size_t t2 = (t1 & (max_bit - 1)) | ((t1 & ~(max_bit - 1)) << 1);

        // t2 has gaps at min_q and max_q (both 0)
        // Now fill gaps:
        // If min_q is control, add min_bit. Else if max_q is control, add max_bit.

        size_t idx0 = t2 | control_mask; // Sets Control=1, Target=0 (since target bit in t2 is 0)
        size_t idx1 = idx0 | target_mask; // Sets Control=1, Target=1

        cfloat a = reg->amplitudes[idx0];
        cfloat b = reg->amplitudes[idx1];

        reg->amplitudes[idx0] = u00 * a + u01 * b;
        reg->amplitudes[idx1] = u10 * a + u11 * b;
    }
    return 1;
}

// -- Gates --

int q_pauli_x(QuantumRegister *reg, int target_qubit) {
    return apply_1q_gate(reg, target_qubit, Q_0_0, Q_1_0, Q_1_0, Q_0_0);
}

int q_pauli_y(QuantumRegister *reg, int target_qubit) {
    return apply_1q_gate(reg, target_qubit, Q_0_0, -Q_I, Q_I, Q_0_0);
}

int q_pauli_z(QuantumRegister *reg, int target_qubit) {
    return apply_1q_gate(reg, target_qubit, Q_1_0, Q_0_0, Q_0_0, -1);
}

int q_hadamard(QuantumRegister *reg, int target_qubit) {
    qfloat s = Q_1_0 / Q_SQRT(Q_2_0);
    return apply_1q_gate(reg, target_qubit, s, s, s, -s);
}

int q_phase(QuantumRegister *reg, int target_qubit) {
    return apply_1q_gate(reg, target_qubit, Q_1_0, Q_0_0, Q_0_0, Q_I);
}

int q_t_gate(QuantumRegister *reg, int target_qubit) {
    return apply_1q_gate(reg, target_qubit, Q_1_0, Q_0_0, Q_0_0, Q_CEXP(Q_I * Q_PI_QUARTER));
}

int q_rx(QuantumRegister *reg, int target_qubit, qfloat theta) {
    cfloat half_theta = theta / Q_2_0;
    cfloat cos_ht = fast_cos(Q_CREAL(half_theta));
    cfloat sin_ht = fast_sin(Q_CREAL(half_theta));
    return apply_1q_gate(reg, target_qubit, cos_ht, -Q_I * sin_ht, -Q_I * sin_ht, cos_ht);
}

int q_ry(QuantumRegister *reg, int target_qubit, qfloat theta) {
    cfloat half_theta = theta / Q_2_0;
    cfloat cos_ht = fast_cos(Q_CREAL(half_theta));
    cfloat sin_ht = fast_sin(Q_CREAL(half_theta));
    return apply_1q_gate(reg, target_qubit, cos_ht, -sin_ht, sin_ht, cos_ht);
}

int q_rz(QuantumRegister *reg, int target_qubit, qfloat theta) {
    cfloat half_theta = theta / Q_2_0;
    cfloat cos_ht = fast_cos(Q_CREAL(half_theta));
    cfloat sin_ht = fast_sin(Q_CREAL(half_theta));
    return apply_1q_gate(reg, target_qubit, cos_ht - Q_I * sin_ht, Q_0_0, Q_0_0, cos_ht + Q_I * sin_ht);
}

int q_u3(QuantumRegister *reg, int target_qubit, qfloat theta, qfloat phi, qfloat lambda) {
    cfloat half_theta = theta / Q_2_0;
    cfloat cos_ht = fast_cos(Q_CREAL(half_theta));
    cfloat sin_ht = fast_sin(Q_CREAL(half_theta));

    cfloat exp_phi = fast_cos(phi) + Q_I * fast_sin(phi);
    cfloat exp_lambda = fast_cos(lambda) + Q_I * fast_sin(lambda);
    cfloat exp_phi_lambda = fast_cos(phi + lambda) + Q_I * fast_sin(phi + lambda);

    cfloat u00 = cos_ht;
    cfloat u01 = -exp_lambda * sin_ht;
    cfloat u10 = exp_phi * sin_ht;
    cfloat u11 = exp_phi_lambda * cos_ht;

    return apply_1q_gate(reg, target_qubit, u00, u01, u10, u11);
}

int q_cnot(QuantumRegister *reg, int control, int target) {
    return apply_controlled_gate(reg, control, target, Q_0_0, Q_1_0, Q_1_0, Q_0_0);
}

int q_cy(QuantumRegister *reg, int control, int target) {
    return apply_controlled_gate(reg, control, target, Q_0_0, -Q_I, Q_I, Q_0_0);
}

int q_cz(QuantumRegister *reg, int control, int target) {
    return apply_controlled_gate(reg, control, target, Q_1_0, Q_0_0, Q_0_0, -1);
}

int q_cp(QuantumRegister *reg, int control, int target, qfloat theta) {
    cfloat phase = fast_cos(theta) + Q_I * fast_sin(theta);
    return apply_controlled_gate(reg, control, target, Q_1_0, Q_0_0, Q_0_0, phase);
}

// OPTIMIZED SWAP
int q_swap(QuantumRegister *reg, int qubit1, int qubit2) {
    if (!reg) return 0;

    // If using OpenCL, fallback to CNOT decomposition for now unless we implement SWAP kernel.
    // For this refactor, let's keep it simple: check OpenCL status.
    // To properly optimize on OpenCL, we'd need a kernel.
    // BUT the task is to optimize SWAP on CPU mostly (assuming report context).
    // If OpenCL is active, we can just call CNOTs or pull data back (slow).
    // Since we didn't add a SWAP kernel to qc2_opencl.c, let's fallback to CNOT-based SWAP if on device.
    #if defined(QC2_USE_OPENCL)
    if (reg->device_amplitudes) {
        if (!q_cnot(reg, qubit1, qubit2)) return 0;
        if (!q_cnot(reg, qubit2, qubit1)) return 0;
        return q_cnot(reg, qubit1, qubit2);
    }
    #endif

    if (qubit1 == qubit2) return 1;
    if (qubit1 >= reg->n_qubits || qubit2 >= reg->n_qubits) return 0;

    size_t dim = reg->dim;
    size_t quarter_dim = dim >> 2;

    int min_q = (qubit1 < qubit2) ? qubit1 : qubit2;
    int max_q = (qubit1 > qubit2) ? qubit1 : qubit2;
    // We iterate states where q1 != q2.
    // Specifically: q1=0, q2=1 AND q1=1, q2=0.
    // Actually we can iterate over indices with 0 at q1 and 0 at q2, then construct the two swap variants.

    size_t min_bit = 1ULL << min_q;
    size_t max_bit = 1ULL << max_q;

    size_t mask_low = min_bit - 1;
    // size_t mask_mid = (max_bit - 1) & ~mask_low;
    // size_t mask_high = ~(max_bit - 1);

    // Actually, simpler logic for masks:
    // We just insert two gaps at min_q and max_q.
    // Loop i goes 0..dim/4.

    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < quarter_dim; i++) {
        // Insert gap at min_q
        size_t t1 = (i & mask_low) | ((i & ~mask_low) << 1);
        // Insert gap at max_q (which is shifted by 1 relative to i if max_q > min_q)
        // Note: min_q is bit index. max_q is bit index.
        // t1 has a 0 at min_q.
        // We need to insert a 0 at max_q.
        // Since t1 has expanded, the bit index for max_q in t1 is unchanged IF max_q > min_q?
        // No, indices > min_q correspond to bits >= min_q+1.
        // So max_q corresponds to position max_q.
        // So we mask t1 at max_bit.

        size_t idx_base = (t1 & (max_bit - 1)) | ((t1 & ~(max_bit - 1)) << 1);

        // idx_base has 0 at min_q and 0 at max_q.

        // State 01: min=0, max=1 (if q1=min, q2=max => q1=0, q2=1)
        // OR q1=max, q2=min => q1=1, q2=0

        // We want to swap state (...0...1...) with (...1...0...)
        // idxA = idx_base | max_bit; // 0 at min, 1 at max
        // idxB = idx_base | min_bit; // 1 at min, 0 at max

        size_t idxA = idx_base | max_bit;
        size_t idxB = idx_base | min_bit;

        cfloat tmp = reg->amplitudes[idxA];
        reg->amplitudes[idxA] = reg->amplitudes[idxB];
        reg->amplitudes[idxB] = tmp;
    }
    return 1;
}

int q_ccnot(QuantumRegister *reg, int control1, int control2, int target) {
    if (!reg) return 0;
    if (control1 >= reg->n_qubits || control2 >= reg->n_qubits || target >= reg->n_qubits) return 0;
    if (control1 == control2 || control1 == target || control2 == target) return 0;

    #if defined(QC2_USE_OPENCL)
    if (reg->device_amplitudes) {
        return apply_toffoli_gate_opencl(reg, control1, control2, target, Q_0_0, Q_1_0, Q_1_0, Q_0_0);
    }
    #endif

    size_t dim = reg->dim;
    size_t eighth_dim = dim >> 3;

    // Sort bits for insertion: q0 < q1 < q2
    int q_sorted[3] = {control1, control2, target};
    if (q_sorted[0] > q_sorted[1]) { int temp = q_sorted[0]; q_sorted[0] = q_sorted[1]; q_sorted[1] = temp; }
    if (q_sorted[1] > q_sorted[2]) { int temp = q_sorted[1]; q_sorted[1] = q_sorted[2]; q_sorted[2] = temp; }
    if (q_sorted[0] > q_sorted[1]) { int temp = q_sorted[0]; q_sorted[0] = q_sorted[1]; q_sorted[1] = temp; }

    int q0 = q_sorted[0];
    int q1 = q_sorted[1];
    int q2 = q_sorted[2];

    size_t mask0 = (1ULL << q0) - 1;
    size_t mask1 = (1ULL << q1) - 1;
    size_t mask2 = (1ULL << q2) - 1;

    size_t bit_c1 = 1ULL << control1;
    size_t bit_c2 = 1ULL << control2;
    size_t bit_t  = 1ULL << target;
    size_t fixed_bits = bit_c1 | bit_c2; // We want c1=1, c2=1, t=0 (from insertion)

    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < eighth_dim; i++) {
        // Insert gap at q0
        size_t t0 = (i & mask0) | ((i & ~mask0) << 1);
        // Insert gap at q1
        size_t t1 = (t0 & mask1) | ((t0 & ~mask1) << 1);
        // Insert gap at q2
        size_t t2 = (t1 & mask2) | ((t1 & ~mask2) << 1);

        size_t idx0 = t2 | fixed_bits; // c1=1, c2=1, t=0
        size_t idx1 = idx0 | bit_t;    // c1=1, c2=1, t=1

        cfloat tmp = reg->amplitudes[idx0];
        reg->amplitudes[idx0] = reg->amplitudes[idx1];
        reg->amplitudes[idx1] = tmp;
    }

    return 1;
}

int q_identity(QuantumRegister *reg, int target_qubit) {
    return apply_1q_gate(reg, target_qubit, Q_1_0, Q_0_0, Q_0_0, Q_1_0);
}

int q_sqrt_x(QuantumRegister *reg, int target_qubit) {
    qfloat half = Q_1_0 / Q_2_0;
    cfloat sqrt_half = half + Q_I * half;
    cfloat sqrt_half_conj = half - Q_I * half;
    return apply_1q_gate(reg, target_qubit, sqrt_half, sqrt_half_conj, sqrt_half_conj, sqrt_half);
}

int q_sqrt_z(QuantumRegister *reg, int target_qubit) {
    return apply_1q_gate(reg, target_qubit, Q_1_0, Q_0_0, Q_0_0, Q_CEXP(Q_I * Q_PI_QUARTER));
}

int q_u1(QuantumRegister *reg, int target_qubit, qfloat lambda) {
    cfloat phase = fast_cos(lambda) + Q_I * fast_sin(lambda);
    return apply_1q_gate(reg, target_qubit, Q_1_0, Q_0_0, Q_0_0, phase);
}

int q_u2(QuantumRegister *reg, int target_qubit, qfloat phi, qfloat lambda) {
    qfloat inv_sqrt_2 = Q_1_0 / Q_SQRT(Q_2_0);
    cfloat exp_phi = fast_cos(phi) + Q_I * fast_sin(phi);
    cfloat exp_lambda = fast_cos(lambda) + Q_I * fast_sin(lambda);
    cfloat exp_phi_lambda = fast_cos(phi + lambda) + Q_I * fast_sin(phi + lambda);
    
    cfloat u00 = inv_sqrt_2;
    cfloat u01 = -inv_sqrt_2 * exp_lambda;
    cfloat u10 = inv_sqrt_2 * exp_phi;
    cfloat u11 = inv_sqrt_2 * exp_phi_lambda;
    
    return apply_1q_gate(reg, target_qubit, u00, u01, u10, u11);
}

int q_p_dagger(QuantumRegister *reg, int target_qubit) {
    return apply_1q_gate(reg, target_qubit, Q_1_0, Q_0_0, Q_0_0, -Q_I);
}

int q_hadamard_row(QuantumRegister *reg, int first_qubit, int num_qubits) {
    if (!reg) return 0;
    if (first_qubit < 0 || num_qubits < 1) return 0;
    if (first_qubit + num_qubits > reg->n_qubits) return 0;
    
    for (int i = 0; i < num_qubits; i++) {
        if (!q_hadamard(reg, first_qubit + i)) return 0;
    }
    return 1;
}

int q_walsh(QuantumRegister *reg) {
    if (!reg) return 0;
    return q_hadamard_row(reg, 0, reg->n_qubits);
}

int q_ch(QuantumRegister *reg, int control, int target) {
    qfloat inv_sqrt_2 = Q_1_0 / Q_SQRT(Q_2_0);
    return apply_controlled_gate(reg, control, target, inv_sqrt_2, inv_sqrt_2, inv_sqrt_2, -inv_sqrt_2);
}

int q_cs(QuantumRegister *reg, int control, int target) {
    return q_phase(reg, target) ? apply_controlled_gate(reg, control, target, Q_1_0, Q_0_0, Q_0_0, Q_I) : 0;
}

int q_ct(QuantumRegister *reg, int control, int target) {
    return apply_controlled_gate(reg, control, target, Q_1_0, Q_0_0, Q_0_0, Q_CEXP(Q_I * Q_PI_QUARTER));
}

int q_crx(QuantumRegister *reg, int control, int target, qfloat theta) {
    cfloat half_theta = theta / Q_2_0;
    cfloat cos_ht = fast_cos(Q_CREAL(half_theta));
    cfloat sin_ht = fast_sin(Q_CREAL(half_theta));
    return apply_controlled_gate(reg, control, target, cos_ht, -Q_I * sin_ht, -Q_I * sin_ht, cos_ht);
}

int q_cry(QuantumRegister *reg, int control, int target, qfloat theta) {
    cfloat half_theta = theta / Q_2_0;
    cfloat cos_ht = fast_cos(Q_CREAL(half_theta));
    cfloat sin_ht = fast_sin(Q_CREAL(half_theta));
    return apply_controlled_gate(reg, control, target, cos_ht, -sin_ht, sin_ht, cos_ht);
}

int q_crz(QuantumRegister *reg, int control, int target, qfloat theta) {
    cfloat half_theta = theta / Q_2_0;
    cfloat cos_ht = fast_cos(Q_CREAL(half_theta));
    cfloat sin_ht = fast_sin(Q_CREAL(half_theta));
    return apply_controlled_gate(reg, control, target, cos_ht - Q_I * sin_ht, Q_0_0, Q_0_0, cos_ht + Q_I * sin_ht);
}

int q_cu1(QuantumRegister *reg, int control, int target, qfloat lambda) {
    cfloat phase = fast_cos(lambda) + Q_I * fast_sin(lambda);
    return apply_controlled_gate(reg, control, target, Q_1_0, Q_0_0, Q_0_0, phase);
}

int q_iswap(QuantumRegister *reg, int qubit1, int qubit2) {
    if (!reg) return 0;
    if (qubit1 == qubit2) return 1;
    if (qubit1 >= reg->n_qubits || qubit2 >= reg->n_qubits) return 0;

    #if defined(QC2_USE_OPENCL)
    if (reg->device_amplitudes) {
        if (!q_swap(reg, qubit1, qubit2)) return 0;
        if (!q_phase(reg, qubit1)) return 0;
        if (!q_phase(reg, qubit2)) return 0;
        return q_cz(reg, qubit1, qubit2);
    }
    #endif

    size_t dim = reg->dim;
    size_t quarter_dim = dim >> 2;

    int min_q = (qubit1 < qubit2) ? qubit1 : qubit2;
    int max_q = (qubit1 > qubit2) ? qubit1 : qubit2;

    size_t min_bit = 1ULL << min_q;
    size_t max_bit = 1ULL << max_q;
    size_t mask_low = min_bit - 1;

    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < quarter_dim; i++) {
        size_t t1 = (i & mask_low) | ((i & ~mask_low) << 1);
        size_t idx_base = (t1 & (max_bit - 1)) | ((t1 & ~(max_bit - 1)) << 1);

        //size_t idx00 = idx_base;
        size_t idx01 = idx_base | min_bit;
        size_t idx10 = idx_base | max_bit;
        //size_t idx11 = idx_base | min_bit | max_bit;

        cfloat amp01 = reg->amplitudes[idx01];
        cfloat amp10 = reg->amplitudes[idx10];

        reg->amplitudes[idx01] = Q_I * amp10;
        reg->amplitudes[idx10] = Q_I * amp01;
    }
    return 1;
}

int q_sqrt_cnot(QuantumRegister *reg, int control, int target) {
    qfloat half = Q_1_0 / Q_2_0;
    cfloat sqrt_half = half + Q_I * half;
    cfloat sqrt_half_conj = half - Q_I * half;
    return apply_controlled_gate(reg, control, target, sqrt_half, sqrt_half_conj, sqrt_half_conj, sqrt_half);
}

int q_rxx(QuantumRegister *reg, int qubit1, int qubit2, qfloat theta) {
    if (!reg) return 0;
    if (qubit1 == qubit2) return 1;
    if (qubit1 >= reg->n_qubits || qubit2 >= reg->n_qubits) return 0;

    cfloat half_theta = theta / Q_2_0;
    cfloat cos_ht = fast_cos(Q_CREAL(half_theta));
    cfloat sin_ht = fast_sin(Q_CREAL(half_theta));

    return apply_2q_unitary(reg, qubit1, qubit2,
                            cos_ht, Q_0_0, Q_0_0, -Q_I * sin_ht,
                            Q_0_0, cos_ht, Q_I * sin_ht, Q_0_0,
                            Q_0_0, Q_I * sin_ht, cos_ht, Q_0_0,
                            -Q_I * sin_ht, Q_0_0, Q_0_0, cos_ht);
}

int q_ryy(QuantumRegister *reg, int qubit1, int qubit2, qfloat theta) {
    if (!reg) return 0;
    if (qubit1 == qubit2) return 1;
    if (qubit1 >= reg->n_qubits || qubit2 >= reg->n_qubits) return 0;

    cfloat half_theta = theta / Q_2_0;
    cfloat cos_ht = fast_cos(Q_CREAL(half_theta));
    cfloat sin_ht = fast_sin(Q_CREAL(half_theta));

    return apply_2q_unitary(reg, qubit1, qubit2,
                            cos_ht, Q_0_0, Q_0_0, Q_I * sin_ht,
                            Q_0_0, cos_ht, -Q_I * sin_ht, Q_0_0,
                            Q_0_0, -Q_I * sin_ht, cos_ht, Q_0_0,
                            Q_I * sin_ht, Q_0_0, Q_0_0, cos_ht);
}

int q_rzz(QuantumRegister *reg, int qubit1, int qubit2, qfloat theta) {
    if (!reg) return 0;
    if (qubit1 == qubit2) return 1;
    if (qubit1 >= reg->n_qubits || qubit2 >= reg->n_qubits) return 0;

    cfloat half_theta = theta / Q_2_0;
    cfloat exp_minus = fast_cos(Q_CREAL(half_theta)) - Q_I * fast_sin(Q_CREAL(half_theta));
    cfloat exp_plus = fast_cos(Q_CREAL(half_theta)) + Q_I * fast_sin(Q_CREAL(half_theta));

    return apply_2q_unitary(reg, qubit1, qubit2,
                            exp_minus, Q_0_0, Q_0_0, Q_0_0,
                            Q_0_0, exp_plus, Q_0_0, Q_0_0,
                            Q_0_0, Q_0_0, exp_plus, Q_0_0,
                            Q_0_0, Q_0_0, Q_0_0, exp_minus);
}

int q_ecr(QuantumRegister *reg, int qubit1, int qubit2) {
    if (!reg) return 0;
    q_hadamard(reg, qubit2);
    q_crz(reg, qubit1, qubit2, Q_PI_HALF);
    q_hadamard(reg, qubit1);
    return 1;
}

int q_fredkin(QuantumRegister *reg, int control, int target1, int target2) {
    if (!reg) return 0;
    if (control >= reg->n_qubits || target1 >= reg->n_qubits || target2 >= reg->n_qubits) return 0;
    if (control == target1 || control == target2 || target1 == target2) return 0;

    #if defined(QC2_USE_OPENCL)
    if (reg->device_amplitudes) {
        return apply_fredkin_gate_opencl(reg, control, target1, target2);
    }
    #endif

    size_t dim = reg->dim;
    size_t sixth_dim = dim >> 3;

    int q_sorted[3] = {control, target1, target2};
    if (q_sorted[0] > q_sorted[1]) { int temp = q_sorted[0]; q_sorted[0] = q_sorted[1]; q_sorted[1] = temp; }
    if (q_sorted[1] > q_sorted[2]) { int temp = q_sorted[1]; q_sorted[1] = q_sorted[2]; q_sorted[2] = temp; }
    if (q_sorted[0] > q_sorted[1]) { int temp = q_sorted[0]; q_sorted[0] = q_sorted[1]; q_sorted[1] = temp; }

    int q0 = q_sorted[0];
    int q1 = q_sorted[1];
    int q2 = q_sorted[2];

    size_t mask0 = (1ULL << q0) - 1;
    size_t mask1 = (1ULL << q1) - 1;
    size_t mask2 = (1ULL << q2) - 1;

    size_t bit_c = 1ULL << control;
    size_t bit_t1 = 1ULL << target1;
    size_t bit_t2 = 1ULL << target2;
    size_t fixed_bits = bit_c;

    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < sixth_dim; i++) {
        size_t t0 = (i & mask0) | ((i & ~mask0) << 1);
        size_t t1 = (t0 & mask1) | ((t0 & ~mask1) << 1);
        size_t t2 = (t1 & mask2) | ((t1 & ~mask2) << 1);

        size_t idx_base = t2 | fixed_bits;
        size_t idx_t1_0_t2_1 = idx_base | bit_t2;
        size_t idx_t1_1_t2_0 = idx_base | bit_t1;
        size_t idx_t1_0_t2_0 = idx_base;
        size_t idx_t1_1_t2_1 = idx_base | bit_t1 | bit_t2;

        cfloat tmp0 = reg->amplitudes[idx_t1_0_t2_1];
        cfloat tmp1 = reg->amplitudes[idx_t1_1_t2_0];

        reg->amplitudes[idx_t1_0_t2_1] = reg->amplitudes[idx_t1_1_t2_1];
        reg->amplitudes[idx_t1_1_t2_0] = reg->amplitudes[idx_t1_0_t2_0];
        reg->amplitudes[idx_t1_1_t2_1] = tmp0;
        reg->amplitudes[idx_t1_0_t2_0] = tmp1;
    }

    return 1;
}

int q_cch(QuantumRegister *reg, int control1, int control2, int target) {
    if (!reg) return 0;
    if (control1 >= reg->n_qubits || control2 >= reg->n_qubits || target >= reg->n_qubits) return 0;
    if (control1 == control2 || control1 == target || control2 == target) return 0;

    #if defined(QC2_USE_OPENCL)
    if (reg->device_amplitudes) {
        sync_to_host(reg);
    }
    #endif

    size_t dim = reg->dim;
    size_t eighth_dim = dim >> 3;

    int q_sorted[3] = {control1, control2, target};
    if (q_sorted[0] > q_sorted[1]) { int temp = q_sorted[0]; q_sorted[0] = q_sorted[1]; q_sorted[1] = temp; }
    if (q_sorted[1] > q_sorted[2]) { int temp = q_sorted[1]; q_sorted[1] = q_sorted[2]; q_sorted[2] = temp; }
    if (q_sorted[0] > q_sorted[1]) { int temp = q_sorted[0]; q_sorted[0] = q_sorted[1]; q_sorted[1] = temp; }

    int q0 = q_sorted[0];
    int q1 = q_sorted[1];
    int q2 = q_sorted[2];

    size_t mask0 = (1ULL << q0) - 1;
    size_t mask1 = (1ULL << q1) - 1;
    size_t mask2 = (1ULL << q2) - 1;

    size_t bit_c1 = 1ULL << control1;
    size_t bit_c2 = 1ULL << control2;
    size_t bit_t = 1ULL << target;
    size_t fixed_bits = bit_c1 | bit_c2;
    
    qfloat inv_sqrt_2 = Q_1_0 / Q_SQRT(Q_2_0);

    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < eighth_dim; i++) {
        size_t t0 = (i & mask0) | ((i & ~mask0) << 1);
        size_t t1 = (t0 & mask1) | ((t0 & ~mask1) << 1);
        size_t t2 = (t1 & mask2) | ((t1 & ~mask2) << 1);

        size_t idx0 = t2 | fixed_bits;
        size_t idx1 = idx0 | bit_t;
        
        cfloat a0 = reg->amplitudes[idx0];
        cfloat a1 = reg->amplitudes[idx1];
        
        reg->amplitudes[idx0] = (a0 + a1) * inv_sqrt_2;
        reg->amplitudes[idx1] = (a0 - a1) * inv_sqrt_2;
    }

    #if defined(QC2_USE_OPENCL)
    if (reg->device_amplitudes) {
        sync_to_device(reg);
    }
    #endif

    return 1;
}

int q_ccz(QuantumRegister *reg, int control1, int control2, int target) {
    if (!reg) return 0;
    if (control1 >= reg->n_qubits || control2 >= reg->n_qubits || target >= reg->n_qubits) return 0;
    if (control1 == control2 || control1 == target || control2 == target) return 0;

    #if defined(QC2_USE_OPENCL)
    if (reg->device_amplitudes) {
        return apply_ccz_gate_opencl(reg, control1, control2, target);
    }
    #endif

    size_t dim = reg->dim;
    size_t eighth_dim = dim >> 3;

    int q_sorted[3] = {control1, control2, target};
    if (q_sorted[0] > q_sorted[1]) { int temp = q_sorted[0]; q_sorted[0] = q_sorted[1]; q_sorted[1] = temp; }
    if (q_sorted[1] > q_sorted[2]) { int temp = q_sorted[1]; q_sorted[1] = q_sorted[2]; q_sorted[2] = temp; }
    if (q_sorted[0] > q_sorted[1]) { int temp = q_sorted[0]; q_sorted[0] = q_sorted[1]; q_sorted[1] = temp; }

    int q0 = q_sorted[0];
    int q1 = q_sorted[1];
    int q2 = q_sorted[2];

    size_t mask0 = (1ULL << q0) - 1;
    size_t mask1 = (1ULL << q1) - 1;
    size_t mask2 = (1ULL << q2) - 1;

    size_t bit_c1 = 1ULL << control1;
    size_t bit_c2 = 1ULL << control2;
    size_t bit_t = 1ULL << target;
    size_t fixed_bits = bit_c1 | bit_c2 | bit_t;

    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < eighth_dim; i++) {
        size_t t0 = (i & mask0) | ((i & ~mask0) << 1);
        size_t t1 = (t0 & mask1) | ((t0 & ~mask1) << 1);
        size_t t2 = (t1 & mask2) | ((t1 & ~mask2) << 1);

        size_t idx = t2 | fixed_bits;
        reg->amplitudes[idx] = -reg->amplitudes[idx];
    }

    return 1;
}

int q_ccp(QuantumRegister *reg, int control1, int control2, int target, qfloat theta) {
    if (!reg) return 0;
    if (control1 >= reg->n_qubits || control2 >= reg->n_qubits || target >= reg->n_qubits) return 0;
    if (control1 == control2 || control1 == target || control2 == target) return 0;

    size_t dim = reg->dim;
    size_t eighth_dim = dim >> 3;

    #if defined(QC2_USE_OPENCL)
    if (reg->device_amplitudes) sync_to_host(reg);
    #endif

    int q_sorted[3] = {control1, control2, target};
    if (q_sorted[0] > q_sorted[1]) { int temp = q_sorted[0]; q_sorted[0] = q_sorted[1]; q_sorted[1] = temp; }
    if (q_sorted[1] > q_sorted[2]) { int temp = q_sorted[1]; q_sorted[1] = q_sorted[2]; q_sorted[2] = temp; }
    if (q_sorted[0] > q_sorted[1]) { int temp = q_sorted[0]; q_sorted[0] = q_sorted[1]; q_sorted[1] = temp; }

    int q0 = q_sorted[0];
    int q1 = q_sorted[1];
    int q2 = q_sorted[2];

    size_t mask0 = (1ULL << q0) - 1;
    size_t mask1 = (1ULL << q1) - 1;
    size_t mask2 = (1ULL << q2) - 1;

    size_t bit_c1 = 1ULL << control1;
    size_t bit_c2 = 1ULL << control2;
    size_t bit_t = 1ULL << target;
    size_t fixed_bits = bit_c1 | bit_c2 | bit_t;

    cfloat phase = fast_cos(theta) + Q_I * fast_sin(theta);

    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < eighth_dim; i++) {
        size_t t0 = (i & mask0) | ((i & ~mask0) << 1);
        size_t t1 = (t0 & mask1) | ((t0 & ~mask1) << 1);
        size_t t2 = (t1 & mask2) | ((t1 & ~mask2) << 1);

        size_t idx = t2 | fixed_bits;
        reg->amplitudes[idx] = phase * reg->amplitudes[idx];
    }

    #if defined(QC2_USE_OPENCL)
    if (reg->device_amplitudes) sync_to_device(reg);
    #endif

    return 1;
}

