PROGRAM main
  !> @file main.f90
  !! @brief Main fortran file which calls other functions from other modules.
  !! 
  !! @details This main program calls the functions and subroutines from the pde.f90 file to solve the diffusion equation, given an
  !! input netcdf file. It contains a checkpoint system as a failsafe during simulation runtime so that the entire simulation is not
  !! lost in the event of an unforeseen error. 
  USE ISO_FORTRAN_ENV
  USE pde_solver
  USE input_output_netcdf
  
  IMPLICIT NONE
  
  !Overload function for voltage calculation
  !> @brief Voltage calculator interface
  !!
  !! @details Interface that wraps both the scalar and
  !! array voltage caclulation subroutines from the module pde_solver
  INTERFACE volt_calc
    !> @brief Scalar voltage calculator
    !!
    !! @details Calculates scalar voltage when given 
    !! a scalar input of concentration.
    MODULE PROCEDURE volt_scalar
    !> @brief Array voltage calculator
    !!
    !! @details Calculates array of voltages when given 
    !! an array input of concentration.
    MODULE PROCEDURE volt_array
  END INTERFACE volt_calc
  
  !> @var character loadname 
  !!
  !! Name of the input file.
  !> @var character outname
  !!
  !! Name of the output file.
  !> @var character checkname
  !!
  !! Name of the checkpoint file.
  !> @var real64 A(:,:)
  !!
  !! Left hand side coefficient matrix of the system of equations, please refer to formulation section.
  !> @var real64 B(:,:)
  !!
  !! Right hand side coefficient matrix of the system of equations, please refer to formulation section.
  !> @var real64 c_tmp(:)
  !!
  !! Temporary array of concentration that is loaded in from a checkpoint file.
  !> @var real64 volt(:,:)
  !!
  !! Voltage matrix, saving voltage arrays for each timestep intended to be 1D, but kept as 2D for writing
  !! to a netcdf file.
  !> @var real64 conc(:,:)
  !!
  !! Concentration matrix, saving concentration vectors for each timestep intended to be 2D.
  !> @var int32 quo
  !!
  !! Ratio of number of timesteps to number of outsteps (number of steps until writing to an output file).
  !> @var int32 step_prog
  !!
  !! Counter for the total number of steps.
  !> @var logical print_timing
  !!
  !! Logical toggle to output the timings of the simulation.
  !> @var int32 rate 
  !!
  !! Variable to store the CPU clock rate.
  !> @var int32 ts
  !!
  !! Variable to store the initial number of cpu clocks in a timing.
  !> @var int32 tc
  !!
  !! Variable to store final number of cpu clocks in a timing.
  !> @var int32 conc_time
  !!
  !! The number of CPU cycles taken to do the concentration calculation.
  !> @var int32 volt_time
  !!
  !! The number of CPU cycles taken to do the voltage calculation.
  !> @var int32 write_time
  !!
  !! The number of CPU cycles taken to write to an output file.
  !> @var int32 total_time
  !!
  !! The total number of CPU cycles for the entire simulation (including input and output).

  CHARACTER(LEN=12)               :: loadname = "SPM_input.nc", outname = "SP_output.nc", checkname='SP_check.chp'
  REAL(REAL64),       ALLOCATABLE :: A(:,:), B(:,:)
  REAL(REAL64),       ALLOCATABLE :: c_tmp(:), volt(:,:), conc(:,:)
  INTEGER(kind=int32)             :: i, j, quo, step_prog

  logical                         :: print_timing = .False.
  integer(kind=int32)             :: rate, ts, tc, conc_time, volt_time, write_time, total_time

  call system_clock(total_time,rate)
  
  !Import parameters from input file
  CALL import_input(loadname)

  call system_clock(tc,rate)
  if (print_timing .eqv. .true.) then
     write(0,'("Import Time",18x,": ",F15.6," seconds")')real(tc-total_time,kind=dp)/real(rate,kind=dp)
  end if
  ts = total_time
  
  !allocate arrays
  IF (ALLOCATED(conc)) THEN
    DEALLOCATE(conc)
  END IF
  
  IF (ALLOCATED(volt)) THEN
    DEALLOCATE(volt)
  END IF
    
  IF (ALLOCATED(A)) THEN
    DEALLOCATE(A)
  END IF
    
  IF(ALLOCATED(B)) THEN
    DEALLOCATE(B)
  END IF
  
  ALLOCATE(conc(space_steps,out_steps))
  ALLOCATE(A(space_steps,space_steps))
  ALLOCATE(B(space_steps,space_steps))

  call setup_crank_nicholson(A, B)

  !Assign concentration to an impossible number to check for errors
  conc = -1.0_REAL64
  
  !Get the number of chunks you want to split the concentration matrix into
  quo = CEILING(real(sim_steps, kind=real64)/real(out_steps, kind=real64))

  call system_clock(tc,rate)
  if (print_timing .eqv. .true.) then
     write(0,'("Setup Time",19x,": ",F15.6," seconds")')real(tc-ts,kind=dp)/real(rate,kind=dp)
  end if
  ts = tc
  
  !If starting from a checkpoint then allocate a temporary vector to store the checkpoint concentration
  !Otherwise initialise an output file and a checkpoint file and produce the initial concentration vector
  if (checkpoint .eqv. .True.) then
     IF (ALLOCATED(c_tmp)) THEN
        DEALLOCATE(c_tmp)
     END IF

     ALLOCATE(c_tmp(space_steps))
     call load_checkp(checkname, outname, c_tmp, volt_do)
     
  else
     call initiate_file(outname)
     call initiate_checkp(checkname)

     conc(:,1) = init_c
     
  end if

  call system_clock(tc,rate)
  if (print_timing .eqv. .true.) then
     write(0,'("Initiate out or check",8x,": ",F15.6," seconds")')real(tc-ts,kind=dp)/real(rate,kind=dp)
  end if
  ts = tc
  
  !Calculate parameters to be used later

  !If doing a voltage calculation you need to create a voltage variable in the output file, calculate it, and write it the output file
  !These additional steps are added in if voltage_do=.true.
  IF (volt_do .eqv. .False.) THEN

     !If we are starting from a checkpoint, we need the first row of conc to be the next evolution in the output file, so we perform one time step to get this
     if (checkpoint .eqv. .True.) then
        conc(:,1) = crank_nicholson(A, B, c_tmp)
     end if
    
    !Then we loop over all sim_steps in chunks of out_steps to save the total concentration matrix in blocks
    DO i = 1, (quo-1)
      DO j = 1, (out_steps-1)
         conc(:,(j+1)) = crank_nicholson(A, B, conc(:,j))
      END DO
      
      !Now we save the block to the output file keeping track of the number of steps performed including if we start from a checkpoint
      !And we save the final concentration vector to the checkpoint file overwriting the existing values
      step_prog = (((i-1)*out_steps)+1) + tot_steps
      CALL assign_exp_real(conc, output_id, step_prog, var_id_in=conc_out_id)
      call update_checkp(conc(:,out_steps), step_prog)

      !We then perform one more step to get the first row of the matrix for the next block of steps
      conc(:,1) = crank_nicholson(A, B, conc(:,out_steps))
   END DO

   !Once the above is completed for quo-1 blocks, for the final block we dont have to perform the additional time evolution at the end, so this is removed
   DO j = 1, (out_steps-1)
      conc(:,(j+1)) = crank_nicholson(A, B, conc(:,j))
   END DO
      
   !Now we save the block to the output file keeping track of the number of steps performed including if we start from a checkpoint
   !And we save the final concentration vector to the checkpoint file overwriting the existing values
   step_prog = (((i-1)*out_steps)+1) + tot_steps
   CALL assign_exp_real(conc, output_id, step_prog, var_id_in=conc_out_id)
   call update_checkp(conc(:,out_steps), step_prog)
    
  ELSE
    !If we are doing a voltage calculation as well, we allocate the voltage vector
     ALLOCATE(volt(1, out_steps))

     !If we aren't starting from a checkpoint, we need to create the voltage variable in the output file
     if (checkpoint .eqv. .False.) then
        call create_exp_var('volt', nf90_double, 1, output_id, units='V', act='add', var_id_out=volt_out_id)

     !If we are starting from a checkpoint we assume the voltage variable exists
     !If we are starting from a checkpoint, we need the first row of conc to be the next evolution in the output file, so we perform one time step to get this
     else
        conc(:,1) = crank_nicholson(A, B, c_tmp)
        
     end if

     call system_clock(tc,rate)
     if (print_timing .eqv. .true.) then
        write(0,'("Prepare Loop",17x,": ",F15.6," seconds")')real(tc-ts,kind=dp)/real(rate,kind=dp)
     end if
     ts=tc
     
     conc_time = 0
     volt_time = 0
     write_time = 0
     !Then we loop over all sim_steps in chunks of out_steps to save the total concentration matrix in blocks
     DO i = 1, (quo-1)
        DO j = 1, (out_steps-1)
           call system_clock(ts,rate)
           conc(:,(j+1)) = crank_nicholson(A, B, conc(:,j))
           call system_clock(tc,rate)
           conc_time = conc_time + (tc-ts)
        END DO
        
        call system_clock(ts,rate)
        !Now we take the conctration values at the edge of the vectors and use them to calculate the voltages
        volt(1,:) = volt_calc(conc(space_steps,:))
        call system_clock(tc,rate)
        volt_time = volt_time + (tc-ts)
        !Now we save the conc block and the voltage vector to the output file keeping track of the number of steps performed including if we start from a checkpoint
        !And we save the final concentration vector to the checkpoint file overwriting the existing values
        call system_clock(ts,rate)
        step_prog = (((i-1)*out_steps)+1) + tot_steps
        CALL assign_exp_real(conc, output_id, step_prog, var_id_in=conc_out_id)
        CALL assign_exp_real(volt, output_id, step_prog, var_id_in=volt_out_id)
        call update_checkp(conc(:,out_steps), step_prog)
        call system_clock(tc,rate)
        write_time = write_time + (tc-ts)
        conc(:,1) = crank_nicholson(A, B, conc(:,out_steps))
     END DO

     call system_clock(ts,rate)
     !Once the above is completed for quo-1 blocks, for the final block we dont have to perform the additional time evolution at the end, so this is removed
     DO j = 1, (out_steps-1)
        conc(:,(j+1)) = crank_nicholson(A, B, conc(:,j))
     END DO

     !Calculate the voltage vector for the final conc block
     volt(1,:) = volt_calc(conc(space_steps,:))

     !Now we save the conc block and the voltage vector to the output file keeping track of the number of steps performed including if we start from a checkpoint
     !And we save the final concentration vector to the checkpoint file overwriting the existing values
     step_prog = (((quo-1)*out_steps)+1) + tot_steps
     CALL assign_exp_real(conc, output_id, step_prog, var_id_in=conc_out_id)
     CALL assign_exp_real(volt, output_id, step_prog, var_id_in=volt_out_id)
     call update_checkp(conc(:,out_steps), step_prog)
     call system_clock(tc,rate)

     if (print_timing .eqv. .True.) then
        write(0,'("One Cycle",20x,": ",F15.6," seconds")')real(tc-ts,kind=dp)/real(rate,kind=dp)
        write(0,'("conc time",20x,": ",F15.6," seconds")')real(conc_time,kind=dp)/real(rate,kind=dp)
        write(0,'("volt time",20x,": ",F15.6," seconds")')real(volt_time,kind=dp)/real(rate,kind=dp)
        write(0,'("Write time",19x,": ",F15.6," seconds")')real(write_time,kind=dp)/real(rate,kind=dp)
     end if
     
     ts = tc
  END IF

  !Then we close the open output and checkpoint files
  call fin_in_out()

  call system_clock(tc,rate)
  if (print_timing .eqv. .true.) then
     write(0,'("Total time",19x,": ",F15.6," seconds")')real(tc-total_time,kind=dp)/real(rate,kind=dp)
  end if

END PROGRAM main
