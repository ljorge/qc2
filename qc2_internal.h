// QC v2.0
// This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public License along with this library; if not, see <https://www.gnu.org/licenses/>.
#ifndef _QC2_INTERNAL_H_
#define _QC2_INTERNAL_H_

#include "qc2.h"

// -- Memory Abstraction --
void* qc2_malloc(size_t size);
void qc2_free(void *ptr);

// -- Math Helpers --
#if defined(QC2_USE_TRIG_TABLES)
#define Q_TRIG_TABLE_SIZE 360000

qfloat fast_sin(qfloat theta);
qfloat fast_cos(qfloat theta);
int init_trig_tables(void);
void cleanup_trig_tables(void);
#else
#define fast_sin Q_SIN
#define fast_cos Q_COS
#endif

// -- OpenCL Internals --
#if defined(QC2_USE_OPENCL)
int init_opencl(void);
void cleanup_opencl(void);
void sync_to_device(QuantumRegister *reg);
void sync_to_host(QuantumRegister *reg);
int apply_gate_opencl(QuantumRegister *reg, int target,
                    cfloat u00, cfloat u01, cfloat u10, cfloat u11);
int apply_controlled_gate_opencl(QuantumRegister *reg, int control, int target,
                                cfloat u00, cfloat u01, cfloat u10, cfloat u11);
int apply_toffoli_gate_opencl(QuantumRegister *reg, int control1, int control2, int target,
                                cfloat u00, cfloat u01, cfloat u10, cfloat u11);
int apply_fredkin_gate_opencl(QuantumRegister *reg, int control, int target1, int target2);
int apply_ccz_gate_opencl(QuantumRegister *reg, int control1, int control2, int target);
int qc2_opencl_init_register(QuantumRegister *reg);
void qc2_opencl_free_register(QuantumRegister *reg);
void qc2_opencl_reset_register(QuantumRegister *reg);
#endif

#endif // _QC2_INTERNAL_H_
