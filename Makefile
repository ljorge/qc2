# Compiler Selection gcc/clang (default: gcc)
COMPILER ?= gcc
CC = $(COMPILER)

# Command Definitions
PYTHON ?= python3
AR ?= ar
RM ?= rm -f
RMR ?= rm -rf
MKDIR ?= mkdir -p
ECHO ?= echo
LCOV ?= lcov
GENHTML ?= genhtml
DOXYGEN ?= doxygen
CPPCHECK ?= cppcheck

# Base CFLAGS
CFLAGS = -Wall -Wextra -pedantic -Wshadow -Wconversion -Wunreachable-code -I. -DQC2_COMPILER_NAME=\"$(CC)\"

# Base LDFLAGS
LDFLAGS = -lm

# Configuration Defaults
PRECISION ?= DOUBLE
OPENMP ?= 1
OPENCL ?= 0
OPENCL_VERSION ?= 200
ALIGNED_MEMORY ?= 1
USE_TRIG_TABLES ?= 1
U3_TEST_STEP_DEG ?= 10
DEBUG ?= 0

# Precision flags
ifeq ($(PRECISION),FLOAT)
    CFLAGS += -DQC2_USE_FLOAT
else ifeq ($(PRECISION),DOUBLE)
    CFLAGS += -DQC2_USE_DOUBLE
else ifeq ($(PRECISION),LONG_DOUBLE)
    CFLAGS += -DQC2_USE_LONG_DOUBLE
else
    $(error Invalid PRECISION. Use FLOAT, DOUBLE, or LONG_DOUBLE)
endif

# OpenMP flags
ifeq ($(OPENMP),1)
    CFLAGS += -fopenmp -DQC2_USE_OPENMP
    LDFLAGS += -fopenmp
endif

# OpenCL flags
ifeq ($(OPENCL), 1)
	CFLAGS += -DQC2_USE_OPENCL -DCL_TARGET_OPENCL_VERSION=$(OPENCL_VERSION)
	LDFLAGS += -lOpenCL
endif

# Aligned Memory
ifeq ($(ALIGNED_MEMORY), 0)
	CFLAGS += -DQC2_NO_ALIGNED_MEMORY
endif

# Trig Tables
ifeq ($(USE_TRIG_TABLES), 1)
	CFLAGS += -DQC2_USE_TRIG_TABLES
endif

# U3 Test Step
CFLAGS += -DQC2_U3_TEST_STEP_DEG=$(U3_TEST_STEP_DEG)

# Debug / Optimization Logic
ifeq ($(DEBUG),1)
    # Debug Mode: No optimizations, debug symbols
    CFLAGS += -O0 -g -DQC2_DEBUG
else
    # Release Mode: Max optimizations
    # -native: optimizes for local machine architecture
    # -flto: Link Time Optimization
    CFLAGS += -O3 -march=native -flto
endif

# Object Files
QC2_OBJS = qc2_core.o qc2_gates.o qc2_opencl.o

# Linking Variables
LIB_NAME = qc2
LIB_STATIC = lib$(LIB_NAME).a

# Targets
.PHONY: all clean test doc

all: $(LIB_STATIC)

# Static Library
$(LIB_STATIC): $(QC2_OBJS)
	$(AR) rcs $@ $^

# Build individual modules
qc2_core.o: qc2_core.c qc2.h qc2_internal.h
	$(CC) $(CFLAGS) -c qc2_core.c -o qc2_core.o

qc2_gates.o: qc2_gates.c qc2.h qc2_internal.h
	$(CC) $(CFLAGS) -c qc2_gates.c -o qc2_gates.o

qc2_opencl.o: qc2_opencl.c qc2.h qc2_internal.h
	$(CC) $(CFLAGS) -c qc2_opencl.c -o qc2_opencl.o

# Output test executable to tests/bin/
# Link against static library
tests/bin/test_suite: tests/test_suite.c $(LIB_STATIC)
	@$(MKDIR) tests/bin
	$(CC) $(CFLAGS) tests/test_suite.c -L. -l$(LIB_NAME) -o tests/bin/test_suite $(LDFLAGS)

tests/bin/bell_state: tests/bell_state.c $(LIB_STATIC)
	@$(MKDIR) tests/bin
	$(CC) $(CFLAGS) tests/bell_state.c -L. -l$(LIB_NAME) -o tests/bin/bell_state $(LDFLAGS)

tests/bin/benchmark: tests/benchmark.c $(LIB_STATIC)
	@$(MKDIR) tests/bin
	$(CC) $(CFLAGS) tests/benchmark.c -L. -l$(LIB_NAME) -o tests/bin/benchmark $(LDFLAGS)

tests/bin/ghz_state: tests/ghz_state.c $(LIB_STATIC)
	@$(MKDIR) tests/bin
	$(CC) $(CFLAGS) tests/ghz_state.c -L. -l$(LIB_NAME) -o tests/bin/ghz_state $(LDFLAGS)

tests/bin/bv_algorithm: tests/bv_algorithm.c $(LIB_STATIC)
	@$(MKDIR) tests/bin
	$(CC) $(CFLAGS) tests/bv_algorithm.c -L. -l$(LIB_NAME) -o tests/bin/bv_algorithm $(LDFLAGS)

tests/bin/dj_algorithm: tests/dj_algorithm.c $(LIB_STATIC)
	@$(MKDIR) tests/bin
	$(CC) $(CFLAGS) tests/dj_algorithm.c -L. -l$(LIB_NAME) -o tests/bin/dj_algorithm $(LDFLAGS)

tests/bin/qft_test: tests/qft_test.c $(LIB_STATIC)
	@$(MKDIR) tests/bin
	$(CC) $(CFLAGS) tests/qft_test.c -L. -l$(LIB_NAME) -o tests/bin/qft_test $(LDFLAGS)



test: tests/bin/test_suite tests/bin/bell_state tests/bin/ghz_state tests/bin/bv_algorithm tests/bin/dj_algorithm tests/bin/qft_test tests/bin/qml_real
	@$(ECHO) "***************** ./tests/bin/test_suite *****************"
	./tests/bin/test_suite
	@$(ECHO) "***************** ./tests/bin/bell_state *****************"
	./tests/bin/bell_state
	@$(ECHO) "***************** ./tests/bin/ghz_state (4 Qubits) *******"
	./tests/bin/ghz_state 4
	@$(ECHO) "***************** ./tests/bin/bv_algorithm (5 Qubits) ****"
	./tests/bin/bv_algorithm 5
	@$(ECHO) "***************** ./tests/bin/dj_algorithm (8 Qubits) ****"
	./tests/bin/dj_algorithm 8
	@$(ECHO) "***************** ./tests/bin/qft_test (6 Qubits) ********"
	./tests/bin/qft_test 6
	@$(ECHO) "***************** ./tests/bin/qml_real ***********************"
	./tests/bin/qml_real

tests/bin/qml_real: tests/qml_real.c $(LIB_STATIC) | train_samples.bin
	@$(MKDIR) tests/bin
	$(CC) $(CFLAGS) tests/qml_real.c -L. -l$(LIB_NAME) -o tests/bin/qml_real $(LDFLAGS)

train_samples.bin:
	@$(ECHO) "Downloading/Generating QML Datasets..."
	$(PYTHON) scripts/preprocess_qml.py

benchmark: tests/bin/benchmark
	@$(ECHO) "***************** ./tests/bin/benchmark *****************"
	./tests/bin/benchmark 30

doc:
	$(DOXYGEN) Doxyfile

serve: doc
	@$(ECHO) "Serving documentation at http://localhost:8000"
	$(PYTHON) -m http.server -d doc/html 8000

# Cleanup (including generated data)
clean:
	$(RM) *.o *.a tests/bin/*
	$(RMR) doc html
	$(RM) *.gcno *.gcda *.gcov coverage.info
	$(RMR) coverage_report qc2_kernel_cache.bin
	$(RM) centroids.bin train_samples.bin test_samples.bin curated.tar.gz*
	$(RMR) temp_dataset

# Code Quality
cppcheck:
	$(CPPCHECK) --force --enable=all --inconclusive --std=c99 --suppress=missingIncludeSystem $(filter -D% -I%,$(CFLAGS)) .

coverage: CFLAGS += -fprofile-arcs -ftest-coverage
coverage: LDFLAGS += --coverage
coverage: clean test
	$(LCOV) --capture --directory . --output-file coverage.info
	$(LCOV) --remove coverage.info '/usr/*' --output-file coverage.info
	$(GENHTML) coverage.info --output-directory coverage_report
	@$(ECHO) "Coverage report generated in coverage_report/index.html"
