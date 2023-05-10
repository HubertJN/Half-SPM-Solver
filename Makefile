compiler=gfortran #-Wall -Wextra

flags=`nf-config --fflags`

libraries=`nf-config --flibs` -llapack

main=main.f90

exe=test.out

object=input_output_netcdf.o pde.o

test.out: $(object)
	$(compiler) $(flags) $(object) $(main) $(libraries) -o $(exe)

clean:
	rm -f *.o *.mod $(exe)
	rm test.nc
	rm test.chp

%.o: %.f90
	$(compiler) $(flags) -c $< $(libraries) -o $@
