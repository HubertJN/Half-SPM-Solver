########################################
#MAKEFILE for PX915 Assignment, Group A#
########################################

#Fortran compiler & flags 
#Uncomment -02 to add compiler optimisation 
compiler=gfortran #-02 -Wall -Wextra

flags=`nf-config --fflags` #mpif90 -I/warwick/desktop/2018/software/libpng/1.6.37-GCCcore-8.3.0/include/

libraries=`nf-config --flibs` -llapack -fopenmp

main=main.f90

exe=test.out

object=input_output_netcdf.o fd.o pde.o


#Compile line
test.out: $(object)
	$(compiler) $(flags) $(object) $(main) $(libraries) -o $(exe)
	@#automatically execute commands
	./test.out
	python3 plots.py
	@#clean output files after visualisation is produced 
	@#make clean
	
	
#Checking if object files created 
#ifeq ("$(wildcard $(./fd.o))","")
#	@echo "FD solver ready"
#else
#	@echo "FD solver failed"
#endif	



ifeq ("$(wildcard $(./pde.o))","")
	@echo "Crank-Nicolson solver ready"
else
	@echo "Crank-Nicolson solver failed"
endif


ifeq ("$(wildcard $(./input_output_netcdf.o))","")
	@echo "Input ready"
else
	@echo "Input failed"
endif	

	
#Checking the compilation status
ifeq ("$(wildcard $(./test.out))","")
	@echo "Compilation sucessful"
else
	@echo "Compilation failed"
endif


#Purge build and output files, give no errors if they did not previously exist 
clean:
	rm -f *.o *.mod $(exe)
	rm -f SP_output.nc
	rm -f SP_check.chp
	@echo "Files removed"
	
#Remove the build files but keep the output 
checkpoint:
	rm -f *.o *.mod $(exe)
	@echo "Files removed"
	
#Rules for building object files
%.o: %.f90
	$(compiler) $(flags) -c $< $(libraries) -o $@
	

	
