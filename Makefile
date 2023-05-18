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


#Checking if object files created 
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


#Rules for building object files
%.o: %.f90
	$(compiler) $(flags) -c $< $(libraries) -o $@
		

#Serial or parallel run
ifeq ($(num_threads),1)	
exe: 		
		@./test.out
		@echo "Serial code produced"
else
exe:
		@OMP_NUM_THREADS=$(num_threads) ./test.out 
		@echo "Parallel code produced. Number of threads = " $(num_threads)
		
endif
	
#Generate visualisation
visual: 
	python3 plots.py	



#Purge build and output files, give no errors if they did not previously exist 
clean:
	rm -f *.o *.mod $(exe)
	@echo "Files removed"



#Generate Doxygen documentation
.PHONY: docs
docs:
	-rm docs.html docs.pdf doxy-warns.log
	doxygen ./doxygen/Doxyfile
	(cd ./doxygen/output/latex && make)
	cp ./doxygen/output/latex/refman.pdf docs.pdf
	ln -s ./doxygen/output/html/index.html docs.html
