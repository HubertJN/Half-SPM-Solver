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

object=input_output_netcdf.o pde.o

#Compile line
test.out: $(object)
		$(compiler) $(flags) $(object) $(main) $(libraries) -o $(exe)
		-rm -f *.o *.mod
		chmod +x uq_code/sens_ana.sh
		chmod +x uq_code/up_code.sh

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
.PHONY: exe
ifeq ($(num_threads),1)	
exe: 		
		@./test.out
		@echo "Running Serial code"
else
exe:
		@OMP_NUM_THREADS=$(num_threads) ./test.out 
		@echo "Running Parallel code. Number of threads = " $(num_threads)
endif


#Generate visualisation
.PHONY: visual
visual: 
	gnome-terminal --tab -- python3 plots.py	


#Purge build and output files, give no errors if they did not previously exist
.PHONY: clean
clean:
	-rm -f *.o *.mod *.nc *chp $(exe)
	-rm -f uq_code/*.nc uq_code/*.csv uq_code/*.chp
	-rm -r uq_code/data_store_sens
	-rm -r uq_code/data_store_up
	@echo "Files removed"



#Generate Doxygen documentation
.PHONY: docs
docs:
	-rm docs.html docs.pdf doxy-warns.log
	doxygen ./doxygen/Doxyfile
	(cd ./doxygen/output/latex && make)
	cp ./doxygen/output/latex/refman.pdf docs.pdf
	ln -s ./doxygen/output/html/index.html docs.html

.PHONY: sensitive
sensitive:
	(cd ./uq_code && ./sens_ana.sh False)

.PHONY: vis_sens
vis_sens:
	(cd ./uq_code ; gnome-terminal --tab -- python3 visual_uq_res.py)

.PHONY: uncertain
uncertain:
	(cd ./uq_code && ./up_code.sh False)

.PHONY: vis_uncer
vis_uncer:
	(cd ./uq_code ; gnome-terminal --tab -- python3 visual_up_data.py False)

.PHONY: sens_uncer
sens_uncer:
	(cd ./uq_code && ./sens_ana.sh False)
	(cd ./uq_code && ./up_code.sh True)

.PHONY: vis_sens_uncer
vis_sens_uncer:
	(cd ./uq_code ; gnome-terminal --tab -- python3 visual_up_data.py True)

.PHONY: sens_uncer_sep
sens_uncer_sep:
	(cp ./uq_code/data_store_up/SPM_input_ori.nc ./uq_code/)
	(cd ./uq_code && python3 generate_inp_params.py True)
	(cd ./uq_code ; gnome-terminal --tab -- python3 visual_up_data.py True)

.PHONY: uncer_from_sens
uncer_from_sens:
	(cd ./uq_code && ./sens_ana.sh True)

.PHONY: vis_uncer_from_sens
vis_uncer_from_sens:
	(cd ./uq_code ; gnome-terminal --tab -- python3 vis_uncer_sens.py)
	
	
.PHONY: benchmarking
benchmarking:
	(cd ./benchmarking ; gnome-terminal --tab -- python3 benchmarking.py)		
	
#Create a virtual environment 	
.PHONY: virtual
virtual: 
	chmod +x datafitting/compile.sh
	#chmod +x compile_1.sh
	sudo apt install python3.10-venv
	python3 -m venv venv

#Install modules for the virtual environment
#must be done after venv is activated: source venv/bin/activate
.PHONY: mods
mods:
	pip install -r requirements.txt 
	jupyter contrib nbextension install --user
	jupyter nbextension enable codefolding/main
