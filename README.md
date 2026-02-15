# QC2: Quantum Computing Simulator library for C

![QC2](image/qc2.png)

QC2 is a minimalistic C library for simulating quantum circuits and algorithms.

## Project Overview

**QC2** is a high-performance, precision-agnostic quantum computing simulation library written in C. It was developed to provide a robust foundation for simulating quantum circuits with an emphasis on performance, accuracy, and code clarity.

### Origin and Evolution

This project is a derivative work based on ![QC](image/qc.png)
[QC](https://github.com/Null-byte-00/qc).
It started as a basic quantum simulation tool and has evolved into a fully optimized library capable of handling complex multi-qubit operations and exhaustive state verification.

**License:**

This project is licensed under the **GNU LGPL v2.1** (Lesser General Public License), same as the original work.

**Development Note:**

This project is a derivative work created and modified using **Gemini AI** under human supervision. Community feedback, improvements, and bug fixes are highly encouraged.

**Key Modifications:**

- Refactored from individual qubits to a unified **State Vector** (`QuantumRegister`) architecture to support entanglement.
- Abstracted floating-point precision using `qfloat` and `cfloat` types.
- Implemented a comprehensive test suite with parametric sweeps for universal gates.

---

## Basics

## Initialization & Cleanup

**Crucial:** You must call `qc_init()` before using any other function. This sets up the random number generator and allocates high-performance lookup tables.

At the end of your program, call `qc_destroy()` to free allocated memory.

```c
int main() {
    if (!qc_init()) {
        fprintf(stderr, "Initialization failed\n");
        return 1;
    }

    // Optional: Reset seed for reproducibility (qc_init calls this with time(NULL) by default)
    qc_seed(12345);

    // ... code ...

    qc_destroy();
    return 0;
}
```

### Defining a Register for N Qubits

To simulate a quantum system, you first create a `QuantumRegister`. This structure holds the state vector for `N` qubits (dimension `2^N`).

```c
// Create a register with 2 qubits, initialized to state |00>
QuantumRegister *reg = create_register(2);
```

### Randomness

By default, the library seeds the random number generator using the current time. For reproducible results (debugging), you can manually set the seed:

```c
qc_seed(42); // Seed with a fixed value
QuantumRegister *reg = create_register(2);
```

---

## Gates

All gates now operate on the register and take the target qubit index(es) as arguments.

### Single Qubit Gates

**Pauli-X (NOT), Y, Z:**

```c
q_pauli_x(reg, 0); // Apply X to qubit 0
q_pauli_y(reg, 0);
q_pauli_z(reg, 0);
```

**Hadamard & Phase (S) & T:**

```c
q_hadamard(reg, 0);
q_phase(reg, 0);
q_t_gate(reg, 0);
```

**Rotation Gates (RX, RY, RZ):**

```c
q_rx(reg, 0, 3.14159); // Radians
q_ry(reg, 0, 3.14159/2.0);
q_rz(reg, 0, 3.14159);
```

**Generic Gate (U3):**

```c
// u3(theta, phi, lambda)
q_u3(reg, 0, 3.14, 0.0, 1.57);
```

### Multi-Qubit Gates

**Controlled NOT (CNOT) & Controlled Z (CZ):**

```c
q_cnot(reg, 0, 1); // Control: 0, Target: 1
q_cz(reg, 0, 1);
```

**Swap Gate:**

```c
q_swap(reg, 0, 1);
```

**Toffoli (CCNOT) Gate:**

```c
q_ccnot(reg, 0, 1, 2); // Controls: 0, 1; Target: 2
```

---

## Measuring

Measuring a qubit collapses the state vector of the entire register based on the probability of the outcome.

**Getting probabilities (without collapsing):**

```c
float p1 = get_prob_one(reg, 0);
```

**Measuring (collapses state):**

```c
int result = measure(reg, 0); // Returns 0 or 1
```

---

## Architecture details

### Precision Abstraction

The library supports seamless switching between precisions via compile-time constants:

- **`QC_USE_FLOAT`**: Standard `float` (fastest, lower precision).
- **`QC_USE_DOUBLE`**: Standard `double` (balanced).
- **`QC_USE_LONG_DOUBLE`**: Extended precision `long double` (highest accuracy).

### Optimizations

1. **Algorithmic Bit-Insertion**: Gate operations use a highly optimized single-loop approach that iterates over $N/2$ pairs of amplitudes, avoiding branching logic.

2. **OpenMP Parallelism**: All core functions (initialization, gates) are parallelized with `#pragma omp parallel for`, ensuring scalability.

3. **Trigonometric Caching (Optimized)**: Uses a single shared lookup table (`sin_table`) for both sine and cosine calculations (via $\pi/2$ phase shift), reducing memory usage by 50%.

4. **Compiler Optimization**: The build system enables aggressive optimizations (`-O3`, `-march=native`, `-flto`).

---

## Platform and Build

**Platform:**

- **Linux:** Native support (Debian/Ubuntu optimized).
- **Windows:** Supported via MinGW/WSL (automatically handles aligned memory).

**Build System:** `Makefile` driven.

### Documentation

To generate the API documentation (HTML):

```bash
make doc
# Open doc/html/index.html
```

### Configuration Variables

You can pass these variables to `make` to configure the build:

| Variable           | Default  | Description                                                   |
| :----------------- | :------- | :------------------------------------------------------------ |
| `COMPILER`         | `gcc`    | C Compiler to use (`gcc` or `clang`).                         |
| `PRECISION`        | `DOUBLE` | Floating point precision: `FLOAT`, `DOUBLE`, `LONG_DOUBLE`.   |
| `OPENMP`           | `0`      | Enable OpenMP parallelization (`1` or `0`).                   |
| `OPENCL`           | `1`      | Enable OpenCL acceleration (`1` or `0`).                      |
| `OPENCL_VERSION`   | `120`    | OpenCL target version (e.g. `120` for 1.2).                   |
| `ALIGNED_MEMORY`   | `1`      | Enable aligned memory allocation (`1` or `0`).                |
| `USE_TRIG_TABLES`  | `0`      | Use pre-computed lookup tables for trigonometry (`1` or `0`). |
| `DEBUG`            | `0`      | Enable debug symbols and disable optimizations (`1` or `0`).  |
| `U3_TEST_STEP_DEG` | `15`     | Step size (degrees) for U3 sweep test.                        |

> [!NOTE]
> **OpenCL Restrictions:**
>
> - `LONG_DOUBLE` precision is **NOT** supported when `OPENCL=1`.
> - The library automatically detects and optimizes for OpenCL versions both `< 2.0` and `>= 2.0`.

```bash
# Generate test datasets (Required for QML tests)
python3 scripts/preprocess_qml.py

# Standard Test Run
make test

# High Precision with Custom Step
make test PRECISION=LONG_DOUBLE U3_TEST_STEP_DEG=1
```

## Example

**Bell State (|00> + |11>):**

```c
#include <stdio.h>
#include "qc.h"

int main() {
    qc_init();

    // 1. Create Register (State |00>)
    QuantumRegister *reg = create_register(2);

    // 2. Apply Hadamard to Qubit 0 -> (|00> + |10>)
    q_hadamard(reg, 0);

    // 3. Apply CNOT (Control 0, Target 1) -> (|00> + |11>)
    q_cnot(reg, 0, 1);

    // 4. Measure
    int m0 = measure(reg, 0);
    int m1 = measure(reg, 1);

    printf("Result: %d%d\n", m0, m1); // Always 00 or 11

    destroy_register(reg);
    qc_destroy();
    return 0;
}
```
