MODULE pde_solver

  USE ISO_FORTRAN_ENV
  
  IMPLICIT NONE
  
  CONTAINS
  
  FUNCTION crank_nicholson(rad,dif_coef,flux_param,dt,c_cur)
  
    !solves the diffusion equation with constant diffusion coefficient
    !using the Crank-Nicholson algorithm
    !via the LAPACK library dgtsv function for tridiagonal matrices
    !evolves the given state by one timestep
    
    !rad - max radius of the geometry
    !dif_coef - the diffusion coefficient
    !flux_param - equals to iapp/(aFL), needs to be calculated beforehand!
    !dt - the timestep
    !c_cur - the current concentration vector 
    
    REAL(REAL64), INTENT(IN) :: rad, dif_coef, flux_param, dt
    REAL(REAL64), DIMENSION(:), INTENT(IN) :: c_cur
    REAL(REAL64), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: crank_nicholson
    REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: B                              !evolution matrix
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: AL, A, AU, rhs		        !lhs and rhs of equation
    REAL(REAL64) :: dr, ri, ai                                                  !spacestep, current radius, parameter
    INTEGER :: i                                                                !loop variables
    INTEGER :: n, info                                                          !size of array and dgtsv info
    
    !Get size of input array
    n = SIZE(c_cur)
    
    IF (n < 2) THEN
      PRINT *, "Warning: Invalid input - input array size is less than 2. Terminating." 
      STOP
    END IF
    
    !Check status of allocatables and deallocate as required
    IF (ALLOCATED(AL)) THEN
      DEALLOCATE(AL)
    END IF
    
    IF (ALLOCATED(A)) THEN
      DEALLOCATE(A)
    END IF
    
    IF (ALLOCATED(AU)) THEN
      DEALLOCATE(AU)
    END IF
    
    IF(ALLOCATED(B)) THEN
      DEALLOCATE(B)
    END IF
    
    IF(ALLOCATED(rhs)) THEN
      DEALLOCATE(rhs)
    END IF
    
    IF (ALLOCATED(crank_nicholson)) THEN
      DEALLOCATE(crank_nicholson)
    END IF
    
    !Allocate arrays
    ALLOCATE(AL(n-1))
    ALLOCATE(A(n))
    ALLOCATE(AU(n-1))
    ALLOCATE(B(n,n))
    ALLOCATE(rhs(n))
    ALLOCATE(crank_nicholson(n))
    
    !initialise arrays
    AL = 0.0_REAL64
    A = 0.0_REAL64
    AU = 0.0_REAL64
    B = 0.0_REAL64
    rhs = 0.0_REAL64
    crank_nicholson = 0.0_REAL64
    
    !calculate space-step and evolution parameter
    dr = rad/(REAL(n-1, KIND=REAL64))
    
    !set array values
    !so that A*crank = B*c_cur is the equation of evolution
    
    DO i=2,n-1
      !current radius and dimensionless parameter ai
      ri = REAL((i-1),KIND=REAL64)*dr
      ai = dt*dif_coef/(2.0_REAL64*dr*ri)
      
      !LHS
      AL(i-1) = ai*(1.0_REAL64 - ri/dr)
      A(i) = 1.0_REAL64 + ai*2.0_REAL64*ri/dr
      AU(i) = ai*(-1.0_REAL64 - ri/dr)
      
      !RHS
      B(i,i-1) = ai*(ri/dr - 1.0_REAL64)
      B(i,i) = 1.0_REAL64 - ai*2.0_REAL64*ri/dr
      B(i,i+1) = ai*(1.0_REAL64 + ri/dr)
      
    END DO
    
    !smooth boundary conditions at r=0
    ai = dt*dif_coef/(dr*dr)
    
    A(1) = 1.0_REAL64 + ai
    AU(1) = -1.0_REAL64*ai
    B(1,1) = 1.0_REAL64 - ai
    B(1,2) = ai
    
    !flux boundary conditions
    AL(n-1) = -1.0_REAL64*ai
    A(n) = 1.0_REAL64 + ai
    B(n,n-1) = ai
    B(n,n) = 1.0_REAL64 - ai
    
    !generate RHS
    rhs = MATMUL(B,c_cur)
    rhs(n) = rhs(n) - 2.0_REAL64*flux_param*(dt/rad + dt/dr)
    
    !dgtsv solver
    CALL dgtsv(n,1,AL,A,AU,rhs,n,info)
    
    !Error handling - crude at this stage
    IF (info < 0) THEN
      PRINT *, 'Argument number', -info, 'in dgtsv is illegal' 
      ERROR STOP
    ELSE IF (info > 0) THEN
      PRINT *, 'Determinant of matrix A is zero in Ax=b'
      ERROR STOP
    END IF
    
    !Assign solution to output matrix
    crank_nicholson = rhs
    
  END FUNCTION crank_nicholson
  
  
  FUNCTION U_scalar(x)
  
    !Accepts SCALAR x - stoichiometry
    !Outputs U(c)
    
    REAL(REAL64), INTENT(IN) :: x
    REAL(REAL64), INTENT(OUT) :: U_scalar
    
    !U+ implemented (U- can be implemented later if needed)
    U_scalar = -0.8090_REAL64*x + 4.4875_REAL64 - 0.0428_REAL64*TANH(18.5138_REAL64*(x-0.5542_REAL64)) &
               -17.7326_REAL64*TANH(15.7890_REAL64*(x-0.3117_REAL64)) & 
               +17.5842_REAL64*TANH(15.9308_REAL64*(x-0.3120_REAL64))
                 
  END FUNCTION U_scalar
  
  
  FUNCTION U_arr(x)
    
    !Accepts 1D ARRAY x - stoichiometry
    !Outputs U(c) as 1D array
    
    REAL(REAL64), DIMENSION(:), INTENT(IN) :: x
    REAL(REAL64), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: U_arr
    INTEGER :: size_arr
    
    IF(ALLOCATED(U_arr)) THEN
      DEALLOCATE(U_arr)
    END IF
    
    size_arr = SIZE(x)
    
    ALLOCATE(U_arr(size_arr))
    
    !U+ implemented (U- can be implemented later if needed)
    U_arr = -0.8090_REAL64*x + 4.4875_REAL64 - 0.0428_REAL64*TANH(18.5138_REAL64*(x-0.5542_REAL64)) &
            -17.7326_REAL64*TANH(15.7890_REAL64*(x-0.3117_REAL64)) & 
            +17.5842_REAL64*TANH(15.9308_REAL64*(x-0.3120_REAL64))
    
  END FUNCTION U_arr
  
  
  FUNCTION volt_scalar(cin, Rg, T, F, iapp, a, L, K, cmax)
  
    !Calculates scalar voltage when given 
    !a SCALAR INPUT of concentration
    !and ESSENTIAL PARAMETERS
    
    !Rg = 'gas_con'
    !T = 'temp'
    !F = 'farad'
    !L = 'thick'
    !K = 'rr_coef'
    !cmax = 'max_c'
    
    !TODO: Neaten up (prevent writing so many params, maybe incorporate into module itself?)
    REAL(REAL64), INTENT(IN) :: cin, Rg, T, F, iapp, a, L, K, cmax
    REAL(REAL64) :: arsinh
    REAL(REAL64), INTENT(OUT) :: volt_scalar 
    
    !Calculates arsinh part
    arsinh = F*K*SQRT((cin/cmax)*(1.0_REAL64 - cin/cmax))          !jc
    arsinh = iapp/(a*L*arsinh)                                     !argument of arsinh
    arsinh = LOG(arsinh + SQRT(arsinh**2.0_REAL64 + 1.0_REAL64))   !arsinh
    
    volt_scalar = U_scalar(cin/cmax) - (2.0_REAL64*Rg*T/F)*arsinh
    
  END FUNCTION volt_scalar
  
  
  FUNCTION volt_array(arrin, Rg, T, F, iapp, a, L, K, cmax)
    
    !Calculates an array of voltages when given 
    !an ARRAY INPUT of concentrations
    !and ESSENTIAL PARAMETERS
    
    !TODO: Neaten up (prevent writing so many params, maybe incorporate into module itself?)
    REAL(REAL64), DIMENSION(:), INTENT(IN) :: arrin
    REAL(REAL64), INTENT(IN) :: Rg, T, F, iapp, a, L, K, cmax
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: arsinh
    REAL(REAL64), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: volt_array 
    INTEGER :: size_arr
    
    IF (ALLOCATED(volt_array)) THEN
      DEALLOCATE(volt_array)
    END IF
    
    IF (ALLOCATED(arsinh)) THEN
      DEALLOCATE(arsinh)
    END IF
    
    size_arr = SIZE(arrin)
    
    ALLOCATE(volt_array(size_arr))
    ALLOCATE(arsinh(size_arr))
    
    !Calculates arsinh part
    arsinh = F*K*SQRT((arrin/cmax)*(1.0_REAL64 - arrin/cmax))       !jc
    arsinh = iapp/(a*L*arsinh)                                      !argument of arsinh
    arsinh = LOG(arsinh + SQRT(arsinh**2.0_REAL64 + 1.0_REAL64))    !arsinh
    
    volt_array = U_arr(arrin/cmax) - (2.0_REAL64*Rg*T/F)*arsinh
    
  END FUNCTION volt_array
    
END MODULE pde_solver
