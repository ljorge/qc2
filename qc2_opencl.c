// QC v2.0
// This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public License along with this library; if not, see <https://www.gnu.org/licenses/>.
#include "qc2_internal.h"

#if defined(QC2_USE_OPENCL)
// -- OpenCL Globals --
static cl_context cl_ctx = NULL;
static cl_command_queue cl_queue = NULL;
static cl_program cl_prog = NULL;
static cl_kernel k_apply_gate = NULL;
static cl_kernel k_apply_controlled_gate = NULL;
static cl_kernel k_apply_toffoli_gate = NULL;

// -- Kernel Source --
static const char *CL_KERNEL_SOURCE =
"// Complex multiplication helpers \n"
"inline qcomplex cmul(qcomplex a, qcomplex b) { \n"
"    return (qcomplex)(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x); \n"
"} \n"
"inline qcomplex cadd(qcomplex a, qcomplex b) { \n"
"    return (qcomplex)(a.x+b.x, a.y+b.y); \n"
"} \n"
"\n"
"__kernel void apply_gate(__global qcomplex *state, \n"
"                         unsigned int target_qubit, \n"
"                         qreal u00_r, qreal u00_i, \n"
"                         qreal u01_r, qreal u01_i, \n"
"                         qreal u10_r, qreal u10_i, \n"
"                         qreal u11_r, qreal u11_i) \n"
"{ \n"
"    size_t i = get_global_id(0); \n"
"    size_t target_bit = 1UL << target_qubit; \n"
"    size_t mask_low = target_bit - 1; \n"
"    size_t mask_high = ~mask_low; \n"
"    // Insert 0 at target \n"
"    size_t idx0 = (i & mask_low) | ((i & mask_high) << 1); \n"
"    size_t idx1 = idx0 | target_bit; \n"
"    \n"
"    qcomplex a = state[idx0]; \n"
"    qcomplex b = state[idx1]; \n"
"    \n"
"    qcomplex u00 = (qcomplex)(u00_r, u00_i); \n"
"    qcomplex u01 = (qcomplex)(u01_r, u01_i); \n"
"    qcomplex u10 = (qcomplex)(u10_r, u10_i); \n"
"    qcomplex u11 = (qcomplex)(u11_r, u11_i); \n"
"    \n"
"    state[idx0] = cadd(cmul(u00, a), cmul(u01, b)); \n"
"    state[idx1] = cadd(cmul(u10, a), cmul(u11, b)); \n"
"} \n"
"\n"
"__kernel void apply_controlled_gate(__global qcomplex *state, \n"
"                                    unsigned int control_qubit, \n"
"                                    unsigned int target_qubit, \n"
"                                    qreal u00_r, qreal u00_i, \n"
"                                    qreal u01_r, qreal u01_i, \n"
"                                    qreal u10_r, qreal u10_i, \n"
"                                    qreal u11_r, qreal u11_i) \n"
"{ \n"
"    size_t i = get_global_id(0); \n"
"    size_t target_bit = 1UL << target_qubit; \n"
"    size_t control_bit = 1UL << control_qubit; \n"
"    size_t mask_low = target_bit - 1; \n"
"    size_t mask_high = ~mask_low; \n"
"    \n"
"    size_t idx0 = (i & mask_low) | ((i & mask_high) << 1); \n"
"    \n"
"    if ((idx0 & control_bit) != 0) { \n"
"        size_t idx1 = idx0 | target_bit; \n"
"        qcomplex a = state[idx0]; \n"
"        qcomplex b = state[idx1]; \n"
"        \n"
"        qcomplex u00 = (qcomplex)(u00_r, u00_i); \n"
"        qcomplex u01 = (qcomplex)(u01_r, u01_i); \n"
"        qcomplex u10 = (qcomplex)(u10_r, u10_i); \n"
"        qcomplex u11 = (qcomplex)(u11_r, u11_i); \n"
"        \n"
"        state[idx0] = cadd(cmul(u00, a), cmul(u01, b)); \n"
"        state[idx1] = cadd(cmul(u10, a), cmul(u11, b)); \n"
"    } \n"
"} \n"
"\n"
"__kernel void apply_toffoli_gate(__global qcomplex *state, \n"
"                                 unsigned int c1, unsigned int c2, \n"
"                                 unsigned int target, \n"
"                                 qreal u00_r, qreal u00_i, \n"
"                                 qreal u01_r, qreal u01_i, \n"
"                                 qreal u10_r, qreal u10_i, \n"
"                                 qreal u11_r, qreal u11_i) \n"
"{ \n"
"    size_t i = get_global_id(0); \n"
"    size_t t_bit = 1UL << target; \n"
"    size_t c1_bit = 1UL << c1; \n"
"    size_t c2_bit = 1UL << c2; \n"
"    size_t mask_low = t_bit - 1; \n"
"    size_t mask_high = ~mask_low; \n"
"    \n"
"    size_t idx0 = (i & mask_low) | ((i & mask_high) << 1); \n"
"    \n"
"    if ((idx0 & c1_bit) && (idx0 & c2_bit)) { \n"
"        size_t idx1 = idx0 | t_bit; \n"
"        qcomplex a = state[idx0]; \n"
"        qcomplex b = state[idx1]; \n"
"        \n"
"        qcomplex u00 = (qcomplex)(u00_r, u00_i); \n"
"        qcomplex u01 = (qcomplex)(u01_r, u01_i); \n"
"        qcomplex u10 = (qcomplex)(u10_r, u10_i); \n"
"        qcomplex u11 = (qcomplex)(u11_r, u11_i); \n"
"        \n"
"        state[idx0] = cadd(cmul(u00, a), cmul(u01, b)); \n"
"        state[idx1] = cadd(cmul(u10, a), cmul(u11, b)); \n"
"    } \n"
"} \n";

static const char *CACHE_FILENAME = "qc2_kernel_cache.bin";

static int save_program_binary(cl_program prog, const char *filename) {
    cl_uint num_devices;
    clGetProgramInfo(prog, CL_PROGRAM_NUM_DEVICES, sizeof(cl_uint), &num_devices, NULL);
    if (num_devices == 0) return 0;

    // Get sizes
    size_t *sizes = malloc(sizeof(size_t) * num_devices);
    clGetProgramInfo(prog, CL_PROGRAM_BINARY_SIZES, sizeof(size_t) * num_devices, sizes, NULL);

    // Get binaries
    unsigned char **binaries = malloc(sizeof(unsigned char*) * num_devices);
    for (cl_uint i = 0; i < num_devices; i++) {
        binaries[i] = malloc(sizes[i]);
    }
    clGetProgramInfo(prog, CL_PROGRAM_BINARIES, sizeof(unsigned char*) * num_devices, binaries, NULL);

    // Save the first binary (assuming 1 device for now)
    FILE *fp = fopen(filename, "wb");
    if (fp) {
        fwrite(binaries[0], 1, sizes[0], fp);
        fclose(fp);
    }

    // Cleanup
    for (cl_uint i = 0; i < num_devices; i++) free(binaries[i]);
    free(binaries);
    free(sizes);
    return (fp != NULL);
}

static unsigned char* load_program_binary(const char *filename, size_t *size) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) return NULL;

    fseek(fp, 0, SEEK_END);
    long pos = ftell(fp);
    if (pos < 0) {
        fclose(fp);
        return NULL;
    }
    *size = (size_t)pos;
    rewind(fp);

    unsigned char *binary = malloc(*size);
    if (fread(binary, 1, *size, fp) != *size) {
        free(binary);
        fclose(fp);
        return NULL;
    }
    fclose(fp);
    return binary;
}

/**
 * @brief Initializes OpenCL context, queue, and program.
 *
 * Compiles the kernels.
 *
 * @return 1 on success.
 */
int init_opencl(void) {
    if (cl_ctx) return 1;

    cl_int err;
    cl_uint num_platforms;
    err = clGetPlatformIDs(0, NULL, &num_platforms);
    if (err != CL_SUCCESS || num_platforms == 0) return 0;

    cl_platform_id platform;
    err = clGetPlatformIDs(1, &platform, NULL);
    if (err != CL_SUCCESS) return 0;

    cl_device_id device;
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
    if (err != CL_SUCCESS) {
        err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &device, NULL);
        if (err != CL_SUCCESS) return 0;
    }

    cl_ctx = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
    if (!cl_ctx) return 0;

    #if defined(CL_TARGET_OPENCL_VERSION) && (CL_TARGET_OPENCL_VERSION >= 200)
        cl_queue_properties properties[] = {0};
        cl_queue = clCreateCommandQueueWithProperties(cl_ctx, device, properties, &err);
    #else
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        cl_queue = clCreateCommandQueue(cl_ctx, device, 0, &err);
        #pragma GCC diagnostic pop
    #endif

    if (!cl_queue) return 0;

    // Try loading binary
    size_t bin_size;
    unsigned char *bin = load_program_binary(CACHE_FILENAME, &bin_size);

    int built_from_source = 0;

    if (bin) {
        cl_int binary_status;
        cl_prog = clCreateProgramWithBinary(cl_ctx, 1, &device, &bin_size, (const unsigned char**)&bin, &binary_status, &err);
        free(bin);

        if (err == CL_SUCCESS && binary_status == CL_SUCCESS) {
            err = clBuildProgram(cl_prog, 1, &device, NULL, NULL, NULL);
            if (err != CL_SUCCESS) {
                 // Binary load failed (e.g. driver change), fallback to source
                clReleaseProgram(cl_prog);
                cl_prog = NULL;
            }
        } else {
            if (cl_prog) clReleaseProgram(cl_prog);
            cl_prog = NULL;
        }
    }

    if (!cl_prog) {
        // Build from source
        const char *header;
        #if defined(QC2_USE_DOUBLE)
            header = "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
                        "typedef double qreal;\n"
                        "typedef double2 qcomplex;\n";
        #elif defined(QC2_USE_FLOAT)
            header = "#pragma OPENCL EXTENSION cl_khr_fp32 : enable\n"
                        "typedef float qreal;\n"
                        "typedef float2 qcomplex;\n";
        #else
            #error "Precision not defined."
        #endif

        const char *sources[] = {header, CL_KERNEL_SOURCE};
        cl_prog = clCreateProgramWithSource(cl_ctx, 2, sources, NULL, &err);
        if (!cl_prog) return 0;

        err = clBuildProgram(cl_prog, 1, &device, NULL, NULL, NULL);
        if (err != CL_SUCCESS) {
            size_t len;
            clGetProgramBuildInfo(cl_prog, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
            char *log = malloc(len);
            clGetProgramBuildInfo(cl_prog, device, CL_PROGRAM_BUILD_LOG, len, log, NULL);
            fprintf(stderr, "OpenCL Build Error:\n%s\n", log);
            free(log);
            return 0;
        }
        built_from_source = 1;
    }

    if (built_from_source) {
        save_program_binary(cl_prog, CACHE_FILENAME);
    }

    k_apply_gate = clCreateKernel(cl_prog, "apply_gate", &err);
    if (!k_apply_gate) return 0;

    k_apply_controlled_gate = clCreateKernel(cl_prog, "apply_controlled_gate", &err);
    if (!k_apply_controlled_gate) return 0;

    k_apply_toffoli_gate = clCreateKernel(cl_prog, "apply_toffoli_gate", &err);
    if (!k_apply_toffoli_gate) return 0;

    return 1;
}

/**
 * @brief Releases OpenCL resources.
 */
void cleanup_opencl(void) {
    if (k_apply_gate) clReleaseKernel(k_apply_gate);
    if (k_apply_controlled_gate) clReleaseKernel(k_apply_controlled_gate);
    if (k_apply_toffoli_gate) clReleaseKernel(k_apply_toffoli_gate);
    if (cl_prog) clReleaseProgram(cl_prog);
    if (cl_queue) clReleaseCommandQueue(cl_queue);
    if (cl_ctx) clReleaseContext(cl_ctx);

    cl_ctx = NULL;
    cl_queue = NULL;
    cl_prog = NULL;
    k_apply_gate = NULL;
    k_apply_controlled_gate = NULL;
}

/**
 * @brief Computes device buffer from host amplitudes.
 *
 * @param reg QuantumRegister
 */
void sync_to_device(QuantumRegister *reg) {
    if (!reg->device_amplitudes) return;
    clEnqueueWriteBuffer(cl_queue, reg->device_amplitudes, CL_TRUE, 0, reg->dim * sizeof(cfloat), reg->amplitudes, 0, NULL, NULL);
    reg->device_dirty = 0;
}

/**
 * @brief Copies device amplitudes back to host.
 *
 * Only performs copy if device memory is marked dirty.
 *
 * @param reg QuantumRegister
 */
void sync_to_host(QuantumRegister *reg) {
    if (!reg->device_amplitudes) return;
    if (!reg->device_dirty) return;
    clEnqueueReadBuffer(cl_queue, reg->device_amplitudes, CL_TRUE, 0, reg->dim * sizeof(cfloat), reg->amplitudes, 0, NULL, NULL);
    reg->device_dirty = 0;
}

/**
 * @brief Allocates and initializes OpenCL buffer for a register.
 *
 * @param reg QuantumRegister
 * @return 1 on success
 */
int qc2_opencl_init_register(QuantumRegister *reg) {
    if (!cl_ctx) return 0;
    cl_int err;
    reg->device_amplitudes = clCreateBuffer(cl_ctx, CL_MEM_READ_WRITE, reg->dim * sizeof(cfloat), NULL, &err);
    if (err != CL_SUCCESS) return 0;

    reg->device_dirty = 0;

    // Initialize with zeros or default pattern
    cfloat zero_pattern = Q_0_0 + Q_0_0 * I;
    if (CL_TARGET_OPENCL_VERSION >= 120) {
        err = clEnqueueFillBuffer(cl_queue, reg->device_amplitudes, &zero_pattern, sizeof(cfloat), 0, reg->dim * sizeof(cfloat), 0, NULL, NULL);
    } else {
        sync_to_device(reg); // Fallback
        return 1;
    }

    if (err != CL_SUCCESS) {
        sync_to_device(reg);
    } else {
        // Set first element to 1.0 (since we assume reg->amplitudes[0] is 1.0)
        cfloat initial_state = reg->amplitudes[0];
        clEnqueueWriteBuffer(cl_queue, reg->device_amplitudes, CL_TRUE, 0, sizeof(cfloat), &initial_state, 0, NULL, NULL);
    }
    return 1;
}

/**
 * @brief Frees OpenCL buffer for a register.
 *
 * @param reg QuantumRegister
 */
void qc2_opencl_free_register(QuantumRegister *reg) {
    if (reg->device_amplitudes) {
        clReleaseMemObject(reg->device_amplitudes);
        reg->device_amplitudes = NULL;
    }
}

/**
 * @brief Resets OpenCL buffer to |0...0>.
 *
 * @param reg QuantumRegister
 */
void qc2_opencl_reset_register(QuantumRegister *reg) {
    if (!reg->device_amplitudes) return;
    cl_int err;
    cfloat zero_pattern = Q_0_0 + Q_0_0 * I;

    if (CL_TARGET_OPENCL_VERSION >= 120) {
        err = clEnqueueFillBuffer(cl_queue, reg->device_amplitudes, &zero_pattern, sizeof(cfloat), 0, reg->dim * sizeof(cfloat), 0, NULL, NULL);
    } else {
        sync_to_device(reg);
        return;
    }

    if (err != CL_SUCCESS) {
        sync_to_device(reg);
    } else {
        cfloat initial_state = reg->amplitudes[0];
        clEnqueueWriteBuffer(cl_queue, reg->device_amplitudes, CL_TRUE, 0, sizeof(cfloat), &initial_state, 0, NULL, NULL);
    }
}

/**
 * @brief Enqueues the single qubit gate kernel.
 *
 * @param reg QuantumRegister
 * @param target Target qubit
 * @param u00 Matrix element
 * @param u01 Matrix element
 * @param u10 Matrix element
 * @param u11 Matrix element
 * @return 1 on success
 */
int apply_gate_opencl(QuantumRegister *reg, int target,
                        cfloat u00, cfloat u01, cfloat u10, cfloat u11) {
    if (!reg->device_amplitudes) return 0;

    cl_int err;
    qfloat args[8] = {
        Q_CREAL(u00), Q_CIMAG(u00),
        Q_CREAL(u01), Q_CIMAG(u01),
        Q_CREAL(u10), Q_CIMAG(u10),
        Q_CREAL(u11), Q_CIMAG(u11)
    };

    unsigned int target_arg = (unsigned int)target;

    int arg_idx = 0;
    clSetKernelArg(k_apply_gate, (cl_uint)arg_idx++, sizeof(cl_mem), &reg->device_amplitudes);
    clSetKernelArg(k_apply_gate, (cl_uint)arg_idx++, sizeof(unsigned int), &target_arg);
    for(int i=0; i<8; i++) clSetKernelArg(k_apply_gate, (cl_uint)arg_idx++, sizeof(qfloat), &args[i]);

    size_t global_work_size = reg->dim / 2;
    err = clEnqueueNDRangeKernel(cl_queue, k_apply_gate, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
    if (err != CL_SUCCESS) return 0;

    reg->device_dirty = 1;
    return 1;
}

/**
 * @brief Enqueues the controlled gate kernel.
 *
 * @param reg QuantumRegister
 * @param control Control qubit
 * @param target Target qubit
 * @param u00 Matrix element
 * @param u01 Matrix element
 * @param u10 Matrix element
 * @param u11 Matrix element
 * @return 1 on success
 */
int apply_controlled_gate_opencl(QuantumRegister *reg, int control, int target,
                                cfloat u00, cfloat u01, cfloat u10, cfloat u11) {
    if (!reg->device_amplitudes) return 0;

    cl_int err;
    qfloat args[8] = {
        Q_CREAL(u00), Q_CIMAG(u00),
        Q_CREAL(u01), Q_CIMAG(u01),
        Q_CREAL(u10), Q_CIMAG(u10),
        Q_CREAL(u11), Q_CIMAG(u11)
    };

    unsigned int control_arg = (unsigned int)control;
    unsigned int target_arg = (unsigned int)target;

    int arg_idx = 0;
    clSetKernelArg(k_apply_controlled_gate, (cl_uint)arg_idx++, sizeof(cl_mem), &reg->device_amplitudes);
    clSetKernelArg(k_apply_controlled_gate, (cl_uint)arg_idx++, sizeof(unsigned int), &control_arg);
    clSetKernelArg(k_apply_controlled_gate, (cl_uint)arg_idx++, sizeof(unsigned int), &target_arg);
    for(int i=0; i<8; i++) clSetKernelArg(k_apply_controlled_gate, (cl_uint)arg_idx++, sizeof(qfloat), &args[i]);

    size_t global_work_size = reg->dim / 2;
    err = clEnqueueNDRangeKernel(cl_queue, k_apply_controlled_gate, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
    if (err != CL_SUCCESS) return 0;

    reg->device_dirty = 1;
    return 1;
}

int apply_toffoli_gate_opencl(QuantumRegister *reg, int control1, int control2, int target,
                                cfloat u00, cfloat u01, cfloat u10, cfloat u11) {
    if (!reg->device_amplitudes) return 0;

    cl_int err;
    qfloat args[8] = {
        Q_CREAL(u00), Q_CIMAG(u00),
        Q_CREAL(u01), Q_CIMAG(u01),
        Q_CREAL(u10), Q_CIMAG(u10),
        Q_CREAL(u11), Q_CIMAG(u11)
    };

    unsigned int c1_arg = (unsigned int)control1;
    unsigned int c2_arg = (unsigned int)control2;
    unsigned int t_arg = (unsigned int)target;

    int arg_idx = 0;
    clSetKernelArg(k_apply_toffoli_gate, (cl_uint)arg_idx++, sizeof(cl_mem), &reg->device_amplitudes);
    clSetKernelArg(k_apply_toffoli_gate, (cl_uint)arg_idx++, sizeof(unsigned int), &c1_arg);
    clSetKernelArg(k_apply_toffoli_gate, (cl_uint)arg_idx++, sizeof(unsigned int), &c2_arg);
    clSetKernelArg(k_apply_toffoli_gate, (cl_uint)arg_idx++, sizeof(unsigned int), &t_arg);
    for(int i=0; i<8; i++) clSetKernelArg(k_apply_toffoli_gate, (cl_uint)arg_idx++, sizeof(qfloat), &args[i]);

    size_t global_work_size = reg->dim / 2;
    err = clEnqueueNDRangeKernel(cl_queue, k_apply_toffoli_gate, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
    if (err != CL_SUCCESS) return 0;

    reg->device_dirty = 1;
    return 1;
}
#endif
