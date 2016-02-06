# makefile: makes the adFE2D program


current_dir = $(shell pwd)
SRC=$(current_dir)/src
OBJ=$(current_dir)/obj
BIN=$(current_dir)/bin

# Compiler
FF = gfortran
#FFlags = -Wall -fbounds-check

# Extra object files required by main program
objects=$(OBJ)/input.o $(OBJ)/globals.o $(OBJ)/construct.o $(OBJ)/legendre.o $(OBJ)/linsolver.o


$(BIN)/adFE2D: $(OBJ)/adFE2D.o $(objects)
	$(FF) $(FFlags) -o $@ $+
$(OBJ)/adFE2D.o: $(SRC)/adFE2D.f90 $(objects)
	$(FF) $(FFlags) -I$(OBJ) -c -o $@ $<
$(OBJ)/input.o: $(SRC)/input.f90 $(OBJ)/globals.o
	$(FF) $(FFlags) -J$(OBJ) -c -o $@ $<
$(OBJ)/globals.o: $(SRC)/globals.f90
	$(FF) $(FFlags) -J$(OBJ) -c -o $@ $<
$(OBJ)/construct.o: $(SRC)/construct.f90 $(OBJ)/globals.o $(OBJ)/legendre.o
	$(FF) $(FFlags) -J$(OBJ) -c -o $@ $<
$(OBJ)/legendre.o: $(SRC)/legendre.f90 $(OBJ)/linsolver.o
	$(FF) $(FFlags) -J$(OBJ) -c -o $@ $<
$(OBJ)/linsolver.o: $(SRC)/linsolver.f90 $(OBJ)/globals.o
	$(FF) $(FFlags) -J$(OBJ) -c -o $@ $<
clean:
	rm $(OBJ)/*.o $(OBJ)/*.mod $(BIN)/adFE2D
