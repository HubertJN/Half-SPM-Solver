########################################
#MAKEFILE for PX915 Assignment, Group A#
########################################

#Fortran compiler & flags 
compiler=gfortran #-Wall -Wextra

flags=`nf-config --fflags` #mpif90 -I/warwick/desktop/2018/software/libpng/1.6.37-GCCcore-8.3.0/include/

libraries=`nf-config --flibs` -llapack

main=main.f90

exe=test.out

object=input_output_netcdf.o pde.o



#Compile line
test.out: $(object)
	$(compiler) $(flags) $(object) $(main) $(libraries) -o $(exe)
	
#Checking if object files created 
ifeq ("$(wildcard $(./pde.o))","")
	@echo "Solver ready"
else
	@echo "Solver failed"
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


#Purge files, give no errors if they did not previously exist 
clean:
	rm -f *.o *.mod $(exe)
	rm -f SP_output.nc
	rm -f SP_check.chp
	
	
#Rules for building object files
%.o: %.f90
	$(compiler) $(flags) -c $< $(libraries) -o $@
	

	
