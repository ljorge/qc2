// QC v2.0
// This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public License along with this library; if not, see <https://www.gnu.org/licenses/>.
#include "qc2_internal.h"

// Type used for random number generation
typedef unsigned long long qlong;

// -- Randomness --
static _Thread_local qlong qc2_rng_state = 88172645463325252LL; // Default non-zero seed (Thread Local)

// Xorshift64 implementation
static qlong qc2_xorshift64(void) {
    qlong x = qc2_rng_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    return qc2_rng_state = x;
}

// Generate qfloat in [0, 1)
static qfloat qc2_rand_float(void) {
    qlong r = qc2_xorshift64();

    #if defined(QC2_USE_FLOAT)
    // 23 bits of mantissa, 32 bits is plenty.
    // Divide by 2^32
    return (qfloat)(r & 0xFFFFFFFF) / (qfloat)4294967296.0;

    #elif defined(QC2_USE_DOUBLE)
    // 52 bits of mantissa.
    // Standard trick: r >> 11 (keep 53 bits) multiplied by 2^-53
    // 9007199254740992.0 = 2^53
    return (qfloat)(r >> 11) * (1.0 / 9007199254740992.0);

    #elif defined(QC2_USE_LONG_DOUBLE)
    // Use full 64 bits.
    // 1.8446744073709551616e19L = 2^64
    return (qfloat)r / 18446744073709551616.0L;
    #endif
}

void qc2_seed(unsigned int seed) {
    if (seed == 0) seed = 123456789; // Xorshift cannot handle 0
    qc2_rng_state = (qlong)seed;
    // Warm up
    qc2_xorshift64();
}

// -- Memory Abstraction --

/**
 * @brief internal wrapper for aligned memory allocation.
 *
 * Uses posix_memalign on Posix or _aligned_malloc on Windows.
 *
 * @param size bytes to allocate.
 * @return void* pointer to memory or NULL.
 */
void* qc2_malloc(size_t size) {
    #if defined(QC2_ALIGNED_MEMORY)
        #if defined(_WIN32)
        return _aligned_malloc(size, 64);
        #else
        void *ptr = NULL;
        if (posix_memalign(&ptr, 64, size) != 0) return NULL;
        return ptr;
        #endif
    #else
    return malloc(size);
    #endif
}

/**
 * @brief Internal wrapper for freeing aligned memory.
 *
 * @param ptr Pointer to free.
 */
void qc2_free(void *ptr) {
    #if defined(QC2_ALIGNED_MEMORY) && defined(_WIN32)
    _aligned_free(ptr);
    #else
    free(ptr);
    #endif
}

// -- Trig Tables --
#if defined(QC2_USE_TRIG_TABLES)
static qfloat *sin_table = NULL;

/**
 * @brief Initializes the pre-computed sine table.
 *
 * Allocated on heap to avoid large BSS capability.
 *
 * @return 1 on success.
 */
int init_trig_tables(void) {
    if (sin_table) return 1;
    sin_table = qc2_malloc(Q_TRIG_TABLE_SIZE * sizeof(qfloat));
    if (!sin_table) return 0;

    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel for
    #endif
    for (int i = 0; i < Q_TRIG_TABLE_SIZE; i++) {
        qfloat angle = (qfloat)i * Q_PI_DOUBLE / (qfloat)Q_TRIG_TABLE_SIZE;
        sin_table[i] = Q_SIN(angle);
    }
    return 1;
}

/**
 * @brief Frees the trig table.
 */
void cleanup_trig_tables(void) {
    if (sin_table) {
        qc2_free(sin_table);
        sin_table = NULL;
    }
}

/**
 * @brief Fast sine lookup.
 *
 * @param theta Angle in radians.
 * @return Sine value from table.
 */
qfloat fast_sin(qfloat theta) {
    theta = Q_FMOD(theta, Q_PI_DOUBLE);
    if (theta < 0) theta += Q_PI_DOUBLE;
    int index = (int)(theta * (qfloat)Q_TRIG_TABLE_SIZE / Q_PI_DOUBLE);
    if (index >= Q_TRIG_TABLE_SIZE) index = 0;
    return sin_table[index];
}

/**
 * @brief Fast cosine lookup using sine table.
 *
 * cos(x) = sin(x + PI/2)
 *
 * @param theta Angle in radians.
 * @return Cosine value.
 */
qfloat fast_cos(qfloat theta) {
    theta = Q_FMOD(theta, Q_PI_DOUBLE);
    if (theta < 0) theta += Q_PI_DOUBLE;
    int index = (int)(theta * (qfloat)Q_TRIG_TABLE_SIZE / Q_PI_DOUBLE);
    // cos(x) = sin(x + pi/2)
    index += Q_TRIG_TABLE_SIZE / 4;
    // Wrap around
    if (index >= Q_TRIG_TABLE_SIZE) index -= Q_TRIG_TABLE_SIZE;
    return sin_table[index];
}
#endif

// -- Init/Destroy --
int qc2_init(void) {
    qc2_seed((unsigned)time(NULL));

    return  1
    #if defined(QC2_USE_TRIG_TABLES)
            && init_trig_tables()
    #endif
    #if defined(QC2_USE_OPENCL)
            && init_opencl()
    #endif
    ;
}

void qc2_destroy(void) {
    #if defined(QC2_USE_TRIG_TABLES)
    cleanup_trig_tables();
    #endif

    #if defined(QC2_USE_OPENCL)
    cleanup_opencl();
    #endif
}

// -- Register Management --
QuantumRegister* create_register(int n_qubits) {
    if (n_qubits < 1) return NULL;

    QuantumRegister *reg = qc2_malloc(sizeof(QuantumRegister));
    if (!reg) return NULL;

    reg->n_qubits = n_qubits;
    reg->dim = 1ULL << n_qubits;

    reg->amplitudes = qc2_malloc(reg->dim * sizeof(cfloat));
    if (!reg->amplitudes) {
        qc2_free(reg);
        return NULL;
    }

    // Initialize state to |0...0>
    memset(reg->amplitudes, 0, reg->dim * sizeof(cfloat));
    reg->amplitudes[0] = Q_1_0;

    #if defined(QC2_USE_OPENCL)
    // Initialization of device buffer is complex (fill vs copy).
    // Let's call a helper: qc2_opencl_init_register(reg);
    // For this refactor, I will declare `int qc2_opencl_init_register(QuantumRegister *reg)` in header.
    if (!qc2_opencl_init_register(reg)) {
        // If opencl init fails, we might just fail entirely or fallback?
        // Original code returned NULL if buffer creation failed.
        qc2_free(reg->amplitudes);
        qc2_free(reg);
        return NULL;
    }
    #endif

    return reg;
}

void destroy_register(QuantumRegister *reg) {
    if (!reg) return;

    if (reg->amplitudes) qc2_free(reg->amplitudes);

    #if defined(QC2_USE_OPENCL)
    // Helper: qc2_opencl_free_register(reg);
    qc2_opencl_free_register(reg);
    #endif

    qc2_free(reg);
}

void qc2_reset(QuantumRegister *reg) {
    if (!reg) return;

    memset(reg->amplitudes, 0, reg->dim * sizeof(cfloat));
    reg->amplitudes[0] = Q_1_0;

    #if defined(QC2_USE_OPENCL)
    qc2_opencl_reset_register(reg);
    #endif
}

void print_state(QuantumRegister *reg) {
    #if defined(QC2_USE_OPENCL)
    sync_to_host(reg);
    #endif
    printf("Quantum Register State (N=%d):\n", reg->n_qubits);

    for (size_t i = 0; i < reg->dim; i++) {
        cfloat amp = reg->amplitudes[i];
        if (Q_CABS(amp) > Q_PRINT_THRESHOLD) {
            printf("|%zu>: %.3" Q_FMT " + %.3" Q_FMT "i (prob: %.3" Q_FMT ")\n",
                    i, Q_CREAL(amp), Q_CIMAG(amp), Q_CABS(amp)*Q_CABS(amp));
        }
    }
}

// -- Measurement --
qfloat get_prob_one(const QuantumRegister *reg, int target_qubit) {
    #if defined(QC2_USE_OPENCL)
    sync_to_host((QuantumRegister*)reg);
    #endif
    if (target_qubit >= reg->n_qubits) return Q_0_0;

    qfloat prob = Q_0_0;
    size_t bit = 1ULL << target_qubit;

    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel for reduction(+:prob)
    #endif
    for (size_t i = 0; i < reg->dim; i++) {
        if ((i & bit) != 0) {
            prob += Q_CABS(reg->amplitudes[i]) * Q_CABS(reg->amplitudes[i]);
        }
    }
    return prob;
}

qfloat get_prob_zero(const QuantumRegister *reg, int target_qubit) {
    return Q_1_0 - get_prob_one(reg, target_qubit);
}

int measure(QuantumRegister *reg, int target_qubit) {
    #if defined(QC2_USE_OPENCL)
    sync_to_host(reg);
    #endif
    qfloat p1 = get_prob_one(reg, target_qubit);

    // Use internal RNG instead of rand()
    qfloat r = qc2_rand_float();
    int result = (r < p1) ? 1 : 0;

    size_t bit = 1ULL << target_qubit;
    qfloat norm_factor = Q_0_0;
    if (result == 1) norm_factor = Q_1_0 / Q_SQRT(p1);
    else             norm_factor = Q_1_0 / Q_SQRT(Q_1_0 - p1);

    #if defined(QC2_USE_OPENMP)
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < reg->dim; i++) {
        int is_one = (i & bit) != 0;
        if (is_one == result) {
            reg->amplitudes[i] *= norm_factor;
        } else {
            reg->amplitudes[i] = Q_0_0;
        }
    }

    // Sync back to device
    #if defined(QC2_USE_OPENCL)
    sync_to_device(reg);
    #endif

    return result;
}
