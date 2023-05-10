MODULE fd_solver

  USE ISO_FORTRAN_ENV
  
  IMPLICIT NONE
  
  CONTAINS
  
  FUNCTION fd(rad,dif_coef,flux_param,dt,c_cur)
  
    !solves the diffusion equation with constant diffusion coefficient
    !using the usual normal central finite difference algorithm
    !STABILITY IS IMPORTANT HERE: D*dt/(dr)^2 < 1/2
    !via the LAPACK library dgtsv function for tridiagonal matrices
    !evolves the given state by one timestep
    
    !rad - max radius of the geometry
    !dif_coef - the diffusion coefficient
    !flux_param - equals to iapp/(aFL), needs to be calculated beforehand!
    !dt - the timestep
    !c_cur - the current concentration vector in in out
    
    REAL(REAL64), INTENT(IN) :: rad, dif_coef, flux_param, dt
    REAL(REAL64), DIMENSION(:), INTENT(IN) :: c_cur
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: fd
    REAL(REAL64) :: dr, ri, ai                                                  !spacestep, current radius, parameter
    INTEGER :: i, n                                                             !loop variables and size of array
    
    !Get size of input array
    n = SIZE(c_cur)
    
    IF (n < 2) THEN
      PRINT *, "Warning: Invalid input - input array size is less than 2. Terminating." 
      STOP
    END IF
    
    !Check status of allocatables and deallocate as required
    IF (ALLOCATED(fd)) THEN
      DEALLOCATE(fd)
    END IF
    
    !Allocate arrays
    ALLOCATE(fd(n))
    
    !initialise arrays
    fd = 0.0_REAL64
    
    !calculate space-step and evolution parameter
    dr = rad/(REAL(n-1, KIND=REAL64))
    ai = dt*dif_coef/(dr*dr)
    
    !set new array values
    
    DO i=2,n-1
      !current radius and dimensionless parameter ai
      ri = REAL((i-1),KIND=REAL64)*dr
      fd(i) = ai*(1.0_REAL64 + dr/ri)*c_cur(i+1) + (1.0_REAL64 - 2.0_REAL64*ai)*c_cur(i) + ai*(1.0_REAL64 - dr/ri)*c_cur(i-1)
      
    END DO
    
    !smooth boundary conditions at r=0
    fd(1) = 2.0_REAL64*ai*c_cur(2) + (1.0_REAL64 - 2.0_REAL64*ai)*c_cur(1)
    
    !flux boundary conditions
    ri = REAL((n-1),KIND=REAL64)*dr
    fd(n) = 2.0_REAL64*ai*c_cur(n-1) + (1.0_REAL64 - 2.0_REAL64*ai)*c_cur(n) + 2.0_REAL64*dt*flux_param*(1.0_REAL64 + dr/ri)/dr
    
  END FUNCTION fd
    
END MODULE fd_solver
