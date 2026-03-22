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

This project is a derivative work created and modified using **Antigravity/Gemini AI** and **Opencode/BigPickle** under human supervision. Community feedback, improvements, and bug fixes are highly encouraged.

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

All gates operate on the register and take the target qubit index(es) as arguments.

### Single-Qubit Gates

| Function | Description | Example |
| ---------- | ----------- | ------- |
| `q_identity(reg, qubit)` | Identity (I) | `q_identity(reg, 0);` |
| `q_pauli_x(reg, qubit)` | Pauli-X (NOT, X) | `q_pauli_x(reg, 0);` |
| `q_pauli_y(reg, qubit)` | Pauli-Y (Y) | `q_pauli_y(reg, 0);` |
| `q_pauli_z(reg, qubit)` | Pauli-Z (Z) | `q_pauli_z(reg, 0);` |
| `q_hadamard(reg, qubit)` | Hadamard (H) | `q_hadamard(reg, 0);` |
| `q_phase(reg, qubit)` | Phase (S) | `q_phase(reg, 0);` |
| `q_t_gate(reg, qubit)` | T gate | `q_t_gate(reg, 0);` |
| `q_rx(reg, qubit, theta)` | Rotation X | `q_rx(reg, 0, 3.14159);` |
| `q_ry(reg, qubit, theta)` | Rotation Y | `q_ry(reg, 0, 3.14159);` |
| `q_rz(reg, qubit, theta)` | Rotation Z | `q_rz(reg, 0, 3.14159);` |
| `q_u3(reg, qubit, theta, phi, lambda)` | U3 gate | `q_u3(reg, 0, 1.0, 0.5, 0.25);` |
| `q_sqrt_x(reg, qubit)` | Square-root X (√X) | `q_sqrt_x(reg, 0);` |
| `q_sqrt_z(reg, qubit)` | Square-root Z (√Z) | `q_sqrt_z(reg, 0);` |
| `q_u1(reg, qubit, lambda)` | Phase rotation U1 | `q_u1(reg, 0, 3.14159);` |
| `q_u2(reg, qubit, phi, lambda)` | U2 gate | `q_u2(reg, 0, 0.0, 1.57);` |
| `q_p(reg, qubit)` | Phase gate (P) | `q_p(reg, 0);` |
| `q_p_dagger(reg, qubit)` | Phase dagger (P†) | `q_p_dagger(reg, 0);` |

### Two-Qubit Gates

| Function | Description | Example |
| ---------- | ----------- | ------- |
| `q_cnot(reg, control, target)` | Controlled-NOT (CNOT) | `q_cnot(reg, 0, 1);` |
| `q_cy(reg, control, target)` | Controlled-Y (CY) | `q_cy(reg, 0, 1);` |
| `q_cz(reg, control, target)` | Controlled-Z (CZ) | `q_cz(reg, 0, 1);` |
| `q_cp(reg, control, target, theta)` | Controlled-Phase (CP) | `q_cp(reg, 0, 1, 3.14159);` |
| `q_ch(reg, control, target)` | Controlled-Hadamard (CH) | `q_ch(reg, 0, 1);` |
| `q_cs(reg, control, target)` | Controlled-Phase S (CS) | `q_cs(reg, 0, 1);` |
| `q_ct(reg, control, target)` | Controlled-T (CT) | `q_ct(reg, 0, 1);` |
| `q_crx(reg, control, target, theta)` | Controlled-RX (CRX) | `q_crx(reg, 0, 1, 3.14159);` |
| `q_cry(reg, control, target, theta)` | Controlled-RY (CRY) | `q_cry(reg, 0, 1, 3.14159);` |
| `q_crz(reg, control, target, theta)` | Controlled-RZ (CRZ) | `q_crz(reg, 0, 1, 3.14159);` |
| `q_cu1(reg, control, target, lambda)` | Controlled-U1 (CU1) | `q_cu1(reg, 0, 1, 3.14159);` |
| `q_swap(reg, q1, q2)` | SWAP | `q_swap(reg, 0, 1);` |
| `q_iswap(reg, q1, q2)` | iSWAP | `q_iswap(reg, 0, 1);` |
| `q_sqrt_cnot(reg, control, target)` | Square-root CNOT (√CNOT) | `q_sqrt_cnot(reg, 0, 1);` |
| `q_rxx(reg, q1, q2, theta)` | RXX | `q_rxx(reg, 0, 1, 3.14159);` |
| `q_ryy(reg, q1, q2, theta)` | RYY | `q_ryy(reg, 0, 1, 3.14159);` |
| `q_rzz(reg, q1, q2, theta)` | RZZ | `q_rzz(reg, 0, 1, 3.14159);` |
| `q_ecr(reg, q1, q2)` | ECR | `q_ecr(reg, 0, 1);` |
| `q_cz_equiv(reg, control, target)` | CZ (alias) | `q_cz_equiv(reg, 0, 1);` |

### Three-Qubit Gates

| Function | Description | Example |
| ---------- | ----------- | ------- |
| `q_ccnot(reg, c1, c2, target)` | Toffoli (CCNOT) | `q_ccnot(reg, 0, 1, 2);` |
| `q_fredkin(reg, control, t1, t2)` | Fredkin (CSWAP) | `q_fredkin(reg, 0, 1, 2);` |
| `q_cch(reg, c1, c2, target)` | C-Hadamard (CCH) | `q_cch(reg, 0, 1, 2);` |
| `q_ccz(reg, c1, c2, target)` | CCZ | `q_ccz(reg, 0, 1, 2);` |
| `q_ccp(reg, c1, c2, target, theta)` | CCP | `q_ccp(reg, 0, 1, 2, 3.14159);` |
| `q_toffoli(reg, c1, c2, target)` | Toffoli | `q_toffoli(reg, 0, 1, 2);` |

### Advanced Measurements

| Function | Description | Example |
| ---------- | ----------- | ------- |
| `q_measure_multiple(reg, qubits, n)` | Measure multiple qubits | `int* res = q_measure_multiple(reg, qubits, 3);` |
| `q_measure_partial(reg, qubits, n)` | Get probabilities without collapsing | `qfloat* probs = q_measure_partial(reg, qubits, 3);` |
| `q_measure_basis(reg, qubit, basis)` | Measure in basis X/Y/Z | `int result = q_measure_basis(reg, 0, 2);` |

### Utilities

| Function | Description | Example |
| ---------- | ----------- | ------- |
| `q_normalize(reg)` | Normalize state vector to unit norm | `q_normalize(reg);` |
| `q_fidelity(reg1, reg2)` | Calculate fidelity between two states | `qfloat f = q_fidelity(reg1, reg2);` |
| `q_partial_trace(reg, trace_q, n, remain_q, m)` | Compute partial trace probabilities (caller must free) | `qfloat* probs = q_partial_trace(reg, tq, 1, rq, 2); free(probs);` |
| `q_copy(reg)` | Create a copy of a quantum register (caller must destroy_register) | `QuantumRegister* copy = q_copy(reg); destroy_register(copy);` |

### Aliases (Lowercase)

For compatibility, lowercase aliases are provided:

```c
#define q_i               q_identity
#define q_swap_not        q_sqrt_x
#define q_cswap           q_fredkin
#define q_sqrt_cx         q_sqrt_cnot
#define q_toffoli_gate    q_toffoli
#define q_p               q_phase
#define q_cz_equiv        q_cz
#define q_toffoli         q_ccnot
```

### Precision-Abstracted Constants

To ensure cross-precision compatibility, use the `Q_*` macros for all mathematical values and functions:

| Macro | Description |
| :--- | :--- |
| `Q_PI` | $\pi$ (Pi) |
| `Q_PI_HALF` | $\pi/2$ |
| `Q_PI_QUARTER` | $\pi/4$ |
| `Q_E` | $e$ (Euler's number) |
| `Q_I` | $i$ (Imaginary unit) |
| `Q_0_0` .. `Q_4_0` | Literal constants (0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 4.0) |
| `Q_SIN(x)` / `Q_COS(x)` | Precision-aware sine and cosine |
| `Q_SQRT(x)` | Precision-aware square root |
| `Q_CEXP(x)` | Complex exponential |
| `Q_FMT` | Format string for printf (e.g., `"f"`, `"lf"`, `"Lf"`) |
| `Q_TEST_TOLERANCE` | Standard 0.01 tolerance for unit tests |

---

## Gate Matrices

### Single-Qubit Gate Matrices

#### Identity (I)

```text
│ 1  0 │
│ 0  1 │
```

#### Pauli-X (X, NOT)

```text
│ 0  1 │
│ 1  0 │
```

#### Pauli-Y (Y)

```text
│ 0  -i │
│ i   0 │
```

#### Pauli-Z (Z)

```text
│ 1   0 │
│ 0  -1 │
```

#### Hadamard (H)

```text
│ 1/√2   1/√2 │
│ 1/√2  -1/√2 │
```

#### Phase (S)

```text
│ 1  0 │
│ 0  i │
```

#### T Gate

```text
│ 1      0    │
│ 0  e^(iπ/4) │
```

#### Rotation X (RX)

```text
│   cos(θ/2)   -i·sin(θ/2) │
│ -i·sin(θ/2)    cos(θ/2)  │
```

#### Rotation Y (RY)

```text
│ cos(θ/2)  -sin(θ/2) │
│ sin(θ/2)   cos(θ/2) │
```

#### Rotation Z (RZ)

```text
│ e^(-iθ/2)      0    │
│     0      e^(iθ/2) │
```

#### U3 Gate

```text
│    cos(θ/2)        -e^(iλ)·sin(θ/2)  │
│ e^(iφ)·sin(θ/2)  e^(i(φ+λ))·cos(θ/2) │
```

#### Square-root X (√X)

```text
│ (1+i)/2  (1-i)/2 │
│ (1-i)/2  (1+i)/2 │
```

#### Square-root Z (√Z)

```text
│ 1      0    │
│ 0  e^(iπ/4) │
```

#### U1 Gate

```text
│ 1     0   │
│ 0  e^(iλ) │
```

#### U2 Gate

```text
│   1/√2         -e^(iλ)/√2  │
│ e^(iφ)/√2    e^(i(φ+λ))/√2 │
```

#### Phase Gate (P)

```text
│ 1  0 │
│ 0  i │
```

#### Phase Dagger (P†)

```text
│ 1   0 │
│ 0  -i │
```

### Two-Qubit Gate Matrices

#### CNOT (Controlled-NOT)

```text
│ 1  0  0  0 │
│ 0  1  0  0 │
│ 0  0  0  1 │
│ 0  0  1  0 │
```

#### CY (Controlled-Y)

```text
│ 1  0  0   0 │
│ 0  1  0   0 │
│ 0  0  0  -i │
│ 0  0  i   0 │
```

#### CZ (Controlled-Z)

```text
│ 1  0  0  0 │
│ 0  1  0  0 │
│ 0  0  1  0 │
│ 0  0  0 -1 │
```

#### CP (Controlled-Phase)

```text
│ 1  0  0    0    │
│ 0  1  0    0    │
│ 0  0  1    0    │
│ 0  0  0  e^(iθ) │
```

#### CH (Controlled-Hadamard)

```text
│ 1    0    0      0  │
│ 0    1    0      0  │
│ 0    0  1/√2   1/√2 │
│ 0    0  1/√2  -1/√2 │
```

#### CS (Controlled-Phase S)

```text
│ 1  0  0  0 │
│ 0  1  0  0 │
│ 0  0  1  0 │
│ 0  0  0  i │
```

#### CT (Controlled-T)

```text
│ 1  0  0      0    │
│ 0  1  0      0    │
│ 0  0  1      0    │
│ 0  0  0  e^(iπ/4) │
```

#### CRX (Controlled-RX)

```text
│ 1  0       0           0      │
│ 0  1       0           0      │
│ 0  0    cos(θ/2)  -i·sin(θ/2) │
│ 0  0  -i·sin(θ/2)   cos(θ/2)  │
```

#### CRY (Controlled-RY)

```text
│ 1  0     0          0     │
│ 0  1     0          0     │
│ 0  0  cos(θ/2)  -sin(θ/2) │
│ 0  0  sin(θ/2)   cos(θ/2) │
```

#### CRZ (Controlled-RZ)

```text
│ 1  0      0          0    │
│ 0  1      0          0    │
│ 0  0  e^(-iθ/2)      0    │
│ 0  0      0      e^(iθ/2) │
```

#### CU1 (Controlled-U1)

```text
│ 1  0  0    0    │
│ 0  1  0    0    │
│ 0  0  1    0    │
│ 0  0  0  e^(iλ) │
```

#### SWAP

```text
│ 1  0  0  0 │
│ 0  0  1  0 │
│ 0  1  0  0 │
│ 0  0  0  1 │
```

#### iSWAP

```text
│ 1  0   0  0 │
│ 0  0   i  0 │
│ 0  i   0  0 │
│ 0  0   0  1 │
```

#### √CNOT (Square-root CNOT)

```text
│ 1  0     0        0    │
│ 0  1     0        0    │
│ 0  0  (1+i)/2  (1-i)/2 │
│ 0  0  (1-i)/2  (1+i)/2 │
```

#### RXX

```text
│   cos(θ/2)        0            0       -i·sin(θ/2) │
│      0         cos(θ/2)   -i·sin(θ/2)       0      │
│      0        -i·sin(θ/2)   cos(θ/2)        0      │
│ -i·sin(θ/2)       0            0         cos(θ/2)  │
```

#### RYY

```text
│  cos(θ/2)       0           0        i·sin(θ/2) │
│    0         cos(θ/2)   -i·sin(θ/2)       0     │
│    0        -i·sin(θ/2)   cos(θ/2)        0     │
│ i·sin(θ/2)      0           0         cos(θ/2)  │
```

#### RZZ

```text
│ e^(-iθ/2)    0        0         0     │
│    0      e^(iθ/2)    0         0     │
│    0         0     e^(iθ/2)     0     │
│    0         0        0     e^(-iθ/2) │
```

#### ECR

```text
│      0           0      1/√2   i/√2 │
│      0           0      i/√2   1/√2 │
│    1/√2      -i/√2         0      0 │
│   -i/√2       1/√2         0      0 │
```

### Three-Qubit Gate Matrices

#### CCNOT (Toffoli)

```text
│ 1  0  0  0  0  0  0  0 │
│ 0  1  0  0  0  0  0  0 │
│ 0  0  1  0  0  0  0  0 │
│ 0  0  0  1  0  0  0  0 │
│ 0  0  0  0  1  0  0  0 │
│ 0  0  0  0  0  1  0  0 │
│ 0  0  0  0  0  0  0  1 │
│ 0  0  0  0  0  0  1  0 │
```

Intercambia los estados |111⟩ y |110⟩

#### Fredkin (CSWAP)

```text
│ 1  0  0  0  0  0  0  0 │
│ 0  1  0  0  0  0  0  0 │
│ 0  0  1  0  0  0  0  0 │
│ 0  0  0  1  0  0  0  0 │
│ 0  0  0  0  1  0  0  0 │
│ 0  0  0  0  0  0  1  0 │
│ 0  0  0  0  0  1  0  0 │
│ 0  0  0  0  0  0  0  1 │
```

Intercambia los qubits objetivo cuando el control es 1

#### CCH (Controlled-Controlled-Hadamard)

```text
│ 1  0  0  0  0  0    0    0   │
│ 0  1  0  0  0  0    0    0   │
│ 0  0  1  0  0  0    0    0   │
│ 0  0  0  1  0  0    0    0   │
│ 0  0  0  0  1  0    0    0   │
│ 0  0  0  0  0  1    0    0   │
│ 0  0  0  0  0  0  1/√2  1/√2 │
│ 0  0  0  0  0  0  1/√2 -1/√2 │
```

#### CCZ

```text
│ 1  0  0  0  0  0  0  0 │
│ 0  1  0  0  0  0  0  0 │
│ 0  0  1  0  0  0  0  0 │
│ 0  0  0  1  0  0  0  0 │
│ 0  0  0  0  1  0  0  0 │
│ 0  0  0  0  0  1  0  0 │
│ 0  0  0  0  0  0  1  0 │
│ 0  0  0  0  0  0  0 -1 │
```

Aplica fase -1 al estado |111⟩

#### CCP (Controlled-Controlled-Phase)

```text
│ 1  0  0  0  0  0  0   0    │
│ 0  1  0  0  0  0  0   0    │
│ 0  0  1  0  0  0  0   0    │
│ 0  0  0  1  0  0  0   0    │
│ 0  0  0  0  1  0  0   0    │
│ 0  0  0  0  0  1  0   0    │
│ 0  0  0  0  0  0  1   0    │
│ 0  0  0  0  0  0  0 e^(iθ) │
```

Aplica fase e^(iθ) al estado |111⟩

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
