PROGRAM main

  USE ISO_FORTRAN_ENV
  USE pde_solver
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
  CHARACTER(LEN=7) :: loadname = "test.nc", outname = "tout.nc"
  REAL(REAL64) :: a, iapp, flux_param
  REAL(REAL64), ALLOCATABLE :: conc(:,:)
  REAL(REAL64), ALLOCATABLE :: c_tmp(:), volt(:,:)
  INTEGER :: i, j, quo
  LOGICAL :: test, volt_do                                        !test - default or not; voltage calcs or not
  
  !Import file and initialise mesh
  CALL import_input(loadname,test)
  
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
  quo = CEILING(sim_steps/out_steps) 
  
  !initial conditions and derived parameters
  c_tmp = init_c
  a = 3.0_REAL64*vol_perc/(100.0_REAL64*rad)
  iapp = c_rate*dt/area
  flux_param = iapp/(a*farad*thick)

  call initiate_file(outname)
  
  IF (.NOT. volt_do) THEN
    
    !loop to evolve
    DO i = 1:quo
      DO j = 1:out_steps
        c_tmp = crank_nicholson(rad,dif_coef,flux_param,dt,c_tmp)
        conc(:,j) = c_tmp
      END DO
      
      !save
      CALL save_exp_real('conc',conc,outname,(i-1)*out_steps+1)
    END DO
    
  ELSE
    !include voltage calcs
     ALLOCATE(volt(1, out_steps))
     call create_exp_var('volt', nf90_double, 1, units='V', act='add', file_name=outname)
    
    DO i = 1:quo
      DO j = 1:out_steps
        c_tmp = crank_nicholson(rad,dif_coef,flux_param,dt,c_tmp)
        conc(:,j) = c_tmp
      END DO
      
      !Actual calculation
      volt(1,:) = volt_calc(conc(space_steps,:), gas_con, temp, farad, iapp, a, thick, rr_coef, max_c)
      
      !save
      CALL save_exp_real('conc',conc,outname,(i-1)*out_steps+1)
      CALL save_exp_real('volt',volt,outname,(i-1)*out_steps+1)
    END DO
    
  END IF
    


END PROGRAM main
