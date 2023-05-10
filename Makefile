compiler=gfortran #-Wall -Wextra

flags=`nf-config --fflags`

libraries=`nf-config --flibs` -llapack

main=main.f90

exe=test.out

object=input_output_netcdf.o fd.o pde.o

test.out: $(object)
	$(compiler) $(flags) $(object) $(main) $(libraries) -o $(exe)

clean:
	rm -f *.o *.mod $(exe)
	rm SP_output.nc
	rm SP_check.chp

%.o: %.f90
	$(compiler) $(flags) -c $< $(libraries) -o $@
