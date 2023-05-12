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
  
  !Import file and initialise mesh
  CALL import_input(loadname)
  
  !allocate
  IF (ALLOCATED(conc)) THEN
    DEALLOCATE(conc)
  END IF
  
  IF (ALLOCATED(volt)) THEN
    DEALLOCATE(volt)
  END IF
  
  ALLOCATE(conc(space_steps,out_steps))
  
  conc = -1.0_REAL64
  quo = CEILING(real(sim_steps, kind=real64)/real(out_steps, kind=real64)) 
  
  !initial conditions and derived parameters
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
  
  a = 3.0_REAL64*vol_per/(100.0_REAL64*rad)
  iapp = c_rate*dt/area
  flux_param = iapp/(a*farad*thick)
  
  IF (volt_do .eqv. .False.) THEN

     if (checkpoint .eqv. .True.) then
        conc(:,1) = crank_nicholson(rad,dif_coef,flux_param,dt,c_tmp)
     end if
    
    !loop to evolve
    DO i = 1, (quo-1)
      DO j = 1, (out_steps-1)
         conc(:,(j+1)) = crank_nicholson(rad,dif_coef,flux_param,dt,conc(:,j))
         !conc(:,(j+1)) = fd(rad,dif_coef,flux_param,dt,conc(:,j))
      END DO
      
      !save
      step_prog = (((i-1)*out_steps)+1) + tot_steps
      CALL assign_exp_real(conc, output_id, step_prog, var_id_in=conc_out_id)
      call update_checkp(conc(:,out_steps), step_prog)

      conc(:,1) = crank_nicholson(rad,dif_coef,flux_param,dt,conc(:, out_steps))
   END DO

   DO j = 1, (out_steps-1)
      conc(:,(j+1)) = crank_nicholson(rad,dif_coef,flux_param,dt,conc(:,j))
      !conc(:,(j+1)) = fd(rad,dif_coef,flux_param,dt,conc(:,j))
   END DO
      
   !save
   step_prog = (((i-1)*out_steps)+1) + tot_steps
   CALL assign_exp_real(conc, output_id, step_prog, var_id_in=conc_out_id)
   call update_checkp(conc(:,out_steps), step_prog)
    
  ELSE
    !include voltage calcs
     ALLOCATE(volt(1, out_steps))
     
     if (checkpoint .eqv. .False.) then
        call create_exp_var('volt', nf90_double, 1, output_id, units='V', act='add', var_id_out=volt_out_id)
        
     else
        conc(:,1) = crank_nicholson(rad,dif_coef,flux_param,dt,c_tmp)
        
     end if
    
    DO i = 1, (quo-1)
      DO j = 1, (out_steps-1)
         conc(:,(j+1)) = crank_nicholson(rad,dif_coef,flux_param,dt,conc(:,j))
         !conc(:,(j+1)) = fd(rad,dif_coef,flux_param,dt,conc(:,j))
      END DO
      
      !Actual calculation
      volt(1,:) = volt_calc(conc(space_steps,:), gas_con, temp, farad, iapp, a, thick, rr_coef, max_c)
      
      !save
      step_prog = (((i-1)*out_steps)+1) + tot_steps
      print*, step_prog
      CALL assign_exp_real(conc, output_id, step_prog, var_id_in=conc_out_id)
      CALL assign_exp_real(volt, output_id, step_prog, var_id_in=volt_out_id)
      call update_checkp(conc(:,out_steps), step_prog)

      conc(:,1) = crank_nicholson(rad,dif_coef,flux_param,dt,conc(:,out_steps))
      !conc(:,(j+1)) = fd(rad,dif_coef,flux_param,dt,conc(:,j))
   END DO

   DO j = 1, (out_steps-1)
      conc(:,(j+1)) = crank_nicholson(rad,dif_coef,flux_param,dt,conc(:,j))
      !conc(:,(j+1)) = fd(rad,dif_coef,flux_param,dt,conc(:,j))
   END DO
   
   volt(1,:) = volt_calc(conc(space_steps,:), gas_con, temp, farad, iapp, a, thick, rr_coef, max_c)

   step_prog = (((quo-1)*out_steps)+1) + tot_steps
   CALL assign_exp_real(conc, output_id, step_prog, var_id_in=conc_out_id)
   CALL assign_exp_real(volt, output_id, step_prog, var_id_in=volt_out_id)
   call update_checkp(conc(:,out_steps), step_prog)
    
END IF

    call fin_in_out()


END PROGRAM main
