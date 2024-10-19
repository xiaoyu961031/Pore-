OBJS   = defaults.o \
         matrixops.o \
         matrix.o \
         nonorthoSA.o 

FOPS = -g 

nonorthoSA: $(OBJS) 
	@echo "Linking nonorthoSA.exe ..."
	@ifort $(FOPS) $(OBJS) -o nonorthoSA.exe
	@echo "done"

clean:;	@ rm -f *.o *.mod


%.o: %.F90
	ifort -c $(FOPS) $^
