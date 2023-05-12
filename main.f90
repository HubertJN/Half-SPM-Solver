PROGRAM main

  USE ISO_FORTRAN_ENV
  USE pde_solver
  use fd_solver
  USE input_output_netcdf
  
  IMPLICIT NONE
  
  !Overload function for voltage calculation
  INTERFACE volt_calc
    MODULE PROCEDURE volt_scalar
    MODULE PROCEDURE volt_array
  END INTERFACE volt_calc
  
  !TODO: main program:
  !test functions
  !essential variables
  !Read input
  !initialise/load from checkpoint
  !NO NEED TO CALCULATE STABILITY CRITERION BECAUSE CRANK-NICHOLSON IS UNCONDITIONALLY STABLE
  !Mesh
  !Parallel toggle
  !Time loop (propagate+voltage+write)
  
  !temporary
  CHARACTER(LEN=12) :: loadname = "SPM_input.nc", outname = "SP_output.nc", checkname='SP_check.chp'
  REAL(REAL64) :: a, iapp, flux_param
  REAL(REAL64), ALLOCATABLE :: conc(:,:)
  REAL(REAL64), ALLOCATABLE :: c_tmp(:), volt(:,:)
  INTEGER :: i, j, quo, step_prog
  LOGICAL :: checkpoint = .false., volt_do= .True. !test - default or not; voltage calcs or not
  
  !Import parameters from input file
  CALL import_input(loadname)
  
  !allocate arrays
  IF (ALLOCATED(conc)) THEN
    DEALLOCATE(conc)
  END IF
  
  IF (ALLOCATED(volt)) THEN
    DEALLOCATE(volt)
  END IF
  
  ALLOCATE(conc(space_steps,out_steps))

  !Assign concentration to an impossible number to check for errors
  conc = -1.0_REAL64
  !Get the number of chunks you want to split the concentration matrix into
  quo = CEILING(real(sim_steps, kind=real64)/real(out_steps, kind=real64)) 
  
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

  !Calculate parameters to be used later
  a = 3.0_REAL64*vol_per/(100.0_REAL64*rad)
  iapp = c_rate*dt/area
  flux_param = iapp/(a*farad*thick)

  !If doing a voltage calculation you need to create a voltage variable in the output file, calculate it, and write it the output file
  !These additional steps are added in if voltage_do=.true.
  IF (volt_do .eqv. .False.) THEN

     !If we are starting from a checkpoint, we need the first row of conc to be the next evolution in the output file, so we perform one time step to get this
     if (checkpoint .eqv. .True.) then
        conc(:,1) = crank_nicholson(rad,dif_coef,flux_param,dt,c_tmp)
     end if
    
    !Then we loop over all sim_steps in chunks of out_steps to save the total concentration matrix in blocks
    DO i = 1, (quo-1)
      DO j = 1, (out_steps-1)
         conc(:,(j+1)) = crank_nicholson(rad,dif_coef,flux_param,dt,conc(:,j))
         !conc(:,(j+1)) = fd(rad,dif_coef,flux_param,dt,conc(:,j))
      END DO
      
      !Now we save the block to the output file keeping track of the number of steps performed including if we start from a checkpoint
      !And we save the final concentration vector to the checkpoint file overwriting the existing values
      step_prog = (((i-1)*out_steps)+1) + tot_steps
      CALL assign_exp_real(conc, output_id, step_prog, var_id_in=conc_out_id)
      call update_checkp(conc(:,out_steps), step_prog)

      !We then perform one more step to get the first row of the matrix for the next block of steps
      conc(:,1) = crank_nicholson(rad,dif_coef,flux_param,dt,conc(:, out_steps))
   END DO

   !Once the above is completed for quo-1 blocks, for the final block we dont have to perform the additional time evolution at the end, so this is removed
   DO j = 1, (out_steps-1)
      conc(:,(j+1)) = crank_nicholson(rad,dif_coef,flux_param,dt,conc(:,j))
      !conc(:,(j+1)) = fd(rad,dif_coef,flux_param,dt,conc(:,j))
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
        conc(:,1) = crank_nicholson(rad,dif_coef,flux_param,dt,c_tmp)
        
     end if

     !Then we loop over all sim_steps in chunks of out_steps to save the total concentration matrix in blocks
     DO i = 1, (quo-1)
        DO j = 1, (out_steps-1)
           conc(:,(j+1)) = crank_nicholson(rad,dif_coef,flux_param,dt,conc(:,j))
           !conc(:,(j+1)) = fd(rad,dif_coef,flux_param,dt,conc(:,j))
        END DO
      
        !Now we take the conctration values at the edge of the vectors and use them to calculate the voltages
        volt(1,:) = volt_calc(conc(space_steps,:), gas_con, temp, farad, iapp, a, thick, rr_coef, max_c)

        !Now we save the conc block and the voltage vector to the output file keeping track of the number of steps performed including if we start from a checkpoint
        !And we save the final concentration vector to the checkpoint file overwriting the existing values
        step_prog = (((i-1)*out_steps)+1) + tot_steps
        CALL assign_exp_real(conc, output_id, step_prog, var_id_in=conc_out_id)
        CALL assign_exp_real(volt, output_id, step_prog, var_id_in=volt_out_id)
        call update_checkp(conc(:,out_steps), step_prog)

        conc(:,1) = crank_nicholson(rad,dif_coef,flux_param,dt,conc(:,out_steps))
        !conc(:,(j+1)) = fd(rad,dif_coef,flux_param,dt,conc(:,j))
     END DO

     !Once the above is completed for quo-1 blocks, for the final block we dont have to perform the additional time evolution at the end, so this is removed
     DO j = 1, (out_steps-1)
        conc(:,(j+1)) = crank_nicholson(rad,dif_coef,flux_param,dt,conc(:,j))
        !conc(:,(j+1)) = fd(rad,dif_coef,flux_param,dt,conc(:,j))
     END DO

     !Calculate the voltage vector for the final conc block
     volt(1,:) = volt_calc(conc(space_steps,:), gas_con, temp, farad, iapp, a, thick, rr_coef, max_c)

     !Now we save the conc block and the voltage vector to the output file keeping track of the number of steps performed including if we start from a checkpoint
     !And we save the final concentration vector to the checkpoint file overwriting the existing values
     step_prog = (((quo-1)*out_steps)+1) + tot_steps
     CALL assign_exp_real(conc, output_id, step_prog, var_id_in=conc_out_id)
     CALL assign_exp_real(volt, output_id, step_prog, var_id_in=volt_out_id)
     call update_checkp(conc(:,out_steps), step_prog)
    
  END IF

  !Then we close the open output and checkpoint files
  call fin_in_out()


END PROGRAM main
