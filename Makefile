# Fortran Compiler
FC = gfortran

FFLAGS = -Wall -g -fcheck=all -O3 -std=f2008 -cpp

# Optimzation Flags
# FFLAGS += -D_FMA -mfma

# Dirs
SRC_DIR = src
TEST_DIR = tests
DATA_DIR = test_data

# Library name
LIB = libeftdot.a
OBJS = eftdot.o

PYTHON = python3

all: $(LIB)

$(LIB): $(OBJS)
	ar rcs $@ $^

# Compile the library
eftdot.o: $(SRC_DIR)/eftdot.f90
	$(FC) $(FFLAGS) -c $< -o $@

test_accuracy: $(TEST_DIR)/test.f90 $(LIB)
	$(FC) $(FFLAGS) $< -I. -L. -leftdot -o $@

tests: test_accuracy
	@echo "[1/3] Generating ill-conditioned data..."
	$(PYTHON) $(TEST_DIR)/gen_dot_data.py
	@echo "[2/3] Running precision benchmarks..."
	./test_accuracy
	@echo "[3/3] Plotting results..."
	$(PYTHON) $(TEST_DIR)/plot_res.py
	@echo "Validation complete. See accuracy_plot.pdf"

clean:
	rm -f *.o *.mod $(LIB) test_accuracy *.dat *.pdf input_data.bin results.csv *.png $(SRC_DIR)/*.mod
	rm -rf $(DATA_DIR)

.PHONY: all clean 