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
  LOGICAL :: checkpoint = .False., volt_do= .False.                                        !test - default or not; voltage calcs or not
  
  !Import file and initialise mesh
  CALL import_input(loadname)
  
  !allocate
  IF (ALLOCATED(conc)) THEN
    DEALLOCATE(conc)
  END IF
  
  IF (ALLOCATED(c_tmp)) THEN
    DEALLOCATE(c_tmp)
  END IF
  
  IF (ALLOCATED(volt)) THEN
    DEALLOCATE(volt)
  END IF
  
  ALLOCATE(conc(space_steps,out_steps))
  ALLOCATE(c_tmp(space_steps))
  
  conc = -1.0_REAL64
  volt_do = .FALSE.
  quo = CEILING(real(sim_steps, kind=real64)/real(out_steps, kind=real64)) 
  
  !initial conditions and derived parameters
  if (checkpoint .eqv. .True.) then
     call load_checkp(checkname, c_tmp)
  else
     call initiate_file(outname)
     call initiate_checkp(checkname)
     
     !c_tmp = init_c
     c_tmp = 0.0_DP
     c_tmp(1) = init_c
     c_tmp(20) = init_c
  end if
  
  a = 3.0_REAL64*vol_per/(100.0_REAL64*rad)
  iapp = c_rate*dt/area
  flux_param = iapp/(a*farad*thick)
  
  IF (.NOT. volt_do) THEN
    
    !loop to evolve
    DO i = 1, quo
      DO j = 1, out_steps
         c_tmp = crank_nicholson(rad,dif_coef,flux_param,dt,c_tmp)
         !c_tmp = fd(rad,dif_coef,flux_param,dt,c_tmp)
         conc(:,j) = c_tmp
         print*, c_tmp(1), c_tmp(5), c_tmp(10), c_tmp(15), c_tmp(20)
      END DO
      
      !save
      step_prog = ((i-1)*out_steps+1)
      CALL save_exp_real('conc', conc, outname, step_prog)
      call update_checkp(checkname, conc(:,out_steps), step_prog)
    END DO
    
  ELSE
    !include voltage calcs
     ALLOCATE(volt(1, out_steps))
     
     if (checkpoint .eqv. .False.) then
        call create_exp_var('volt', nf90_double, 1, units='V', act='add', file_name=outname)
     end if
    
    DO i = 1, quo
      DO j = 1, out_steps
        c_tmp = crank_nicholson(rad,dif_coef,flux_param,dt,c_tmp)
        conc(:,j) = c_tmp
      END DO
      
      !Actual calculation
      volt(1,:) = volt_calc(conc(space_steps,:), gas_con, temp, farad, iapp, a, thick, rr_coef, max_c)
      
      !save
      step_prog = ((i-1)*out_steps+1)
      CALL save_exp_real('conc', conc, outname, step_prog)
      CALL save_exp_real('volt', volt, outname, step_prog)
      call update_checkp(checkname, conc(:,out_steps), step_prog)
    END DO
    
  END IF
    


END PROGRAM main
