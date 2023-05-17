########################################
#MAKEFILE for PX915 Assignment, Group A#
########################################

#Fortran compiler & flags 
#Uncomment -O3 to add compiler optimisation 
compiler=gfortran #-O3 #-Wall -Wextra

#set to 1 to run in serial
num_threads=1

flags=`nf-config --fflags` 

ifeq ($(num_threads),1)	
	libraries=`nf-config --flibs` -llapack 
else 
	libraries=`nf-config --flibs` -llapack -fopenmp
endif

main=main.f90

exe=test.out

object=input_output_netcdf.o fd.o pde.o

#Compile line
test.out: $(object)
		$(compiler) $(flags) $(object) $(main) $(libraries) -o $(exe)

ifeq ($(num_threads),1)		
	@#automatically execute commands
	./test.out
	python3 plots.py
	@#clean output files after visualisation is produced 
	@#make clean
	@#echo "num_threads=" $(num_threads)
	@echo "Serial code running"

else
	@#automatically execute commands
	OMP_NUM_THREADS=$(num_threads) ./test.out
	python3 plots.py
	@#clean output files after visualisation is produced 
	@#make clean
	@echo "Parallelising. Number of threads = " $(num_threads)

endif


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
	@echo "Input successful"
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

# Generate Doxygen documentation
.PHONY: docs
docs:
	doxygen ./doxygen/Doxyfile
	(cd ./doxygen/output/latex && make)
	cp ./doxygen/output/latex/refman.pdf docs.pdf
	ln -s ./doxygen/output/html/index.html docs.html
