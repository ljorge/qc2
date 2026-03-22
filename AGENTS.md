# QC2 Library - Documentation

## Project Overview

**QC2** is a high-performance, precision-agnostic quantum computing simulation library written in C. It was developed to provide a robust foundation for simulating quantum circuits with an emphasis on performance, accuracy, and code clarity.

## Origin and Evolution

The project started as a basic quantum simulation tool and has evolved into a highly optimized library capable of handling complex multi-qubit operations and exhaustive state verification.

**Key Modifications:**

- Refactored from individual qubits to a unified **State Vector** (`QuantumRegister`) architecture to support entanglement.
- Abstracted floating-point precision using `qfloat` and `cfloat` types.
- Implemented a comprehensive test suite with parametric sweeps for universal gates.

## Architecture

### Initialization & Memory Management

**Important Lifecycle:**

1. **Initialize**: `qc_init()` must be called once (thread-safe). Returns `1` on success, `0` on failure (e.g., memory allocation).

2. **Usage**: create/manipulate registers.

3. **Cleanup**: `qc_destroy()` frees dynamic tables.

### Quantum Register

The core structure is `QuantumRegister`, representing a system of $N$ qubits.

- **State Vector:** Stores $2^N$ complex amplitudes (`cfloat`).
- **Memory Layout:** Contiguous array, 64-byte aligned using `posix_memalign` (Linux/Posix) or `_aligned_malloc` (Windows) to maximize CPU cache efficiency and AVX vectorization potential.

### Precision Abstraction

The library supports seamless switching between precisions via compile-time constants:

- **`QC_USE_FLOAT`**: Standard `float` (fastest, lower precision).
- **`QC_USE_DOUBLE`**: Standard `double` (balanced).
- **`QC_USE_LONG_DOUBLE`**: Extended precision `long double` (highest accuracy).

Defined via `qfloat`, `cfloat` types and `Q_*` macros (`Q_SIN`, `Q_PI`, `Q_1_0`, etc.).

## Optimizations

### 1. Algorithmic Bit-Insertion

Gate operations use a highly optimized single-loop approach that iterates over $N/2$ pairs of amplitudes. This avoids branching logic dependent on qubit indices.

- **Formulation:** Instead of iterating all $2^N$ states and checking bits, we iterate $0 \dots 2^{N-1}-1$ and "insert" the target bit via bitwise operations:

  ```c
  idx0 = (i & mask_low) | ((i & mask_high) << 1);
  idx1 = idx0 | target_bit;
  ```

- **Benefit:** Guarantees perfect partitionability for parallel execution.

### 2. OpenMP Parallelism

All core operations, including register initialization and gate application, are guarded with `#pragma omp parallel for`. This ensures that for large qubit counts:

- Initialization time is reduced.

- Gate application scales linearly with thread count.

### 3. Trigonometric Caching (Optimized)

Rotation gates (`RX`, `RY`, `RZ`, `U3`) utilize a pre-computed lookup table (`sin_table`) with high resolution (default 360,000 entries).

**Memory Optimization:** The library now uses a *single* sine table for all trigonometric functions. `cos(x)` is calculated as `sin(x + \pi/2)` by shifting the index, effectively halving the memory footprint for lookup tables.

### 4. Compiler Optimization

The build system enables aggressive optimizations:

- `-O3`: Max optimization level.
- `-march=native`: Instructions tuned for the host CPU (AVX2, AVX-512).
- `-flto`: Link Time Optimization.

## Platform and Build

**Platform:** Linux (Debian/Ubuntu optimized).

**Build System:** `Makefile` driven.

```bash
# Standard Test Run
make test

# High Precision with Custom Step
make test PRECISION=LONG_DOUBLE U3_TEST_STEP_DEG=1
```

### Configuration Variables

You can pass these variables to `make` to configure the build (e.g., `make test PRECISION=FLOAT`).

| Variable | Default | Description |
| :--- | :--- | :--- |
| `COMPILER` | `gcc` | C Compiler to use (`gcc` or `clang`). |
| `PRECISION` | `DOUBLE` | Floating point precision: `FLOAT`, `DOUBLE`, `LONG_DOUBLE`. |
| `OPENMP` | `0` | Enable OpenMP parallelization (`1` or `0`). |
| `OPENCL` | `1` | Enable OpenCL acceleration (`1` or `0`). |
| `OPENCL_VERSION` | `120` | OpenCL target version (e.g. `120` for 1.2). |
| `ALIGNED_MEMORY` | `1` | Enable aligned memory allocation (`1` or `0`). |
| `USE_TRIG_TABLES` | `0` | Use pre-computed lookup tables for trigonometry (`1` or `0`). |
| `DEBUG` | `0` | Enable debug symbols and disable optimizations (`1` or `0`). |
| `U3_TEST_STEP_DEG` | `15` | Step size (degrees) for U3 sweep test. |

> [!NOTE]
> **OpenCL Restrictions:**
>
> - `LONG_DOUBLE` precision is **NOT** supported when `OPENCL=1`.
> - The library automatically detects and optimizes for OpenCL versions both `< 2.0` and `>= 2.0`.

## Testing Strategy

Tests are located in `tests/test_suite.c`.

- **Unit Tests:** Verify individual gates (Pauli, Hadamard, Phase, CNOT, etc.) against expected probabilities.
- **Combinatorial Sweep:** The `U3` gate is tested against millions of parameter combinations ($\theta, \phi, \lambda$) to ensure numerical stability across the Bloch sphere.
- **Tolerance:** Comparisons use `Q_TEST_TOLERANCE` (0.01) to account for floating-point drift.

## Documentation

The project uses **Doxygen** to generate API documentation.

- **Process:**
  1. Header files (`.h`) document public APIs.
  2. Source files (`.c`) document internal implementations.
  3. `make doc` generates HTML documentation in `doc/html`.

## Code Structure (For AI Reference)

### File Organization

| File | Purpose |
|------|---------|
| `qc2.h` | Public API declarations, type definitions, constants |
| `qc2_internal.h` | Internal function declarations, OpenCL interfaces |
| `qc2_core.c` | Register management, initialization, measurement, utilities |
| `qc2_gates.c` | Gate implementations (single, multi-qubit) |
| `qc2_opencl.c` | OpenCL kernels and GPU acceleration |

### Implementation Patterns

- **Bit-insertion**: Used in `apply_1q_gate()` and `apply_controlled_gate()` for O(N/2) parallelism
- **OpenMP**: `#pragma omp parallel for` on all loops over state vector
- **Memory**: 64-byte aligned using `posix_memalign` (or `_aligned_malloc` on Windows)
- **Types**: `qfloat` (real), `cfloat` (complex) - precision abstracted via macros

### Adding New Gates

1. Declare in `qc2.h` with Doxygen comments including matrix
2. Implement in `qc2_gates.c` using `apply_1q_gate()` or `apply_controlled_gate()`
3. Add OpenCL kernel in `qc2_opencl.c` if GPU support needed
4. Add test in `tests/test_suite.c`

### Precision Constants

Use `Q_*` macros defined in `qc2.h`: `Q_0_0`, `Q_0_25`, `Q_0_5`, `Q_0_75`, `Q_1_0`, `Q_2_0`, `Q_4_0`, `Q_PI`, `Q_E`, etc.
