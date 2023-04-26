compiler=gfortran #-Wall -Wextra

flags=`nf-config --fflags`

libraries=`nf-config --flibs`

inp=prog.f90

exe=test.out

object=input_output_netcdf.o

test.out: $(object)
	$(compiler) $(flags) $(inp) $(libraries) -o $(exe) $(object)

clean:
	rm -f *.o *.mod $(exe)

input_output_netcdf.o:
	$(compiler) $(flags) $(libraries) input_output_netcdf.f90 -c -o $@ $<
