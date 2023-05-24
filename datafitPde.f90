MODULE datafitPde
  !> @file datafitPde.f90
  !! @brief Program file for solving the diffusion equation and voltage calculation, called for data fitting. 
  !!
  !! @details This contains the subroutines necessary for evolving the state of the system under diffusion and 
  !! calculating the voltage for each simulation timestep. The program contains one singular function which outputs the voltage
  !! so that it can be called in python using f2py and f90wrap for optimisation purposes. 

  USE ISO_FORTRAN_ENV
  !$  use omp_lib
  
  IMPLICIT NONE
  
CONTAINS

  !> @brief Function which runs the Crank-Nicholson solver and calculates the voltage. 
  !!
  !! @details Solves the diffusion equation with a constant diffusion coefficient
  !! and a constant \f$ i_{app} \f$
  !! using the Crank-Nicholson algorithm.
  !! The function uses the LAPACK library, calling the dgesv function for solving systems
  !! of linear equations to evolve the given state by one timestep.
  !!
  !! @brief Function to calculate the positive electrode OCV curve, \f$ U(c) \f$.
  !!
  !! @details OCV stands for Open Circuit Voltage. Please refer to the paper, 
  !! Chang-Hui Chen et. al. 2020 J. Electrochem Soc. 167 080534
  !!
  !! @param[in] x                The stoichiometry, as an array
  !!
  !! @param[in] A             The left hand side coefficient matrix of the system 
  !! of equations (please refer to formulation section
  !! for specific details on the matrix elements)
  !! @param[in] B             The right hand side coefficient matrix of the system
  !! of equations (please refer to formulation section
  !! for specific details on the matrix elements)
  !! @param[in] c0         The current time concentration array inputed into the
  !! Crank-Nicholson solver. 
  !!
  !! @var real64 rhs_const
  !!
  !! Rescaled source at the flux boundary given by: \f$ \left(dt + \frac{dt}{dr} \right) \frac{200 R i_{app}}{3 \epsilon_{actk} F L} \f$
  !> @var real64 volt_con_ial
  !!
  !! This is a combination of parameters given by: \f$ \frac{100 i_{app} R}{3 \epsilon_{actk} L} \f$
  !> @var real64 volt_con_rtf
  !!
  !! This is a combination of parameters given by: \f$ \frac{2 R_g T}{F} \f$
  !> @var real64 mod_dif
  !!
  !! The rescaled diffusion coefficient given by: \f$ \frac{D}{R^2} \f$ 
  !! ipiv - pivotal arrays used by dgesv (link to dgesv lapack) 
  FUNCTION crank_nicholson(n, totalTime, D, R, volPer, iapp, F, L, Rg, T, K, maxCon, c0, dt) RESULT(voltArray)

    INTEGER(4), INTENT(IN)                                     :: n, totalTime
    REAL(8), DIMENSION(:,:), ALLOCATABLE                       :: A, B
    REAL(8),   DIMENSION(:,:), ALLOCATABLE                     :: A_mod
    REAL(8), DIMENSION(n), INTENT(IN)                          :: c0
    REAL(8), DIMENSION(:),   ALLOCATABLE                       :: c_cur, rhs
    REAL(8)                                                    :: U_scalar
    REAL(8), DIMENSION(totalTime)                              :: voltArray
    
    REAL(8), INTENT(IN)                                        :: D, R, volPer, iapp, F, L, Rg, T, K, maxCon, dt
    REAL(8)                                                    :: ai, ri, num, fluxParam, dr, &
                                                                       rhsConst, voltConIal, voltConRtf, modD 
    REAL(8)                                                    ::  arsinh, div_const !div_const is the stoichiometry.
    INTEGER(4)                                                 :: i, info, time
    INTEGER, DIMENSION(n)                                      :: ipiv
    
    print *, 'Flux b.c. rescaled: '
    print *, 'Diffusion coefficient is: ', D

    ALLOCATE(A(n,n))
    ALLOCATE(B(n,n))

    IF (n < 2) THEN
      PRINT *, "Warning: Invalid input - input array size is less than 2. Terminating." 
      STOP
    END IF
    
   !----------------------Setup matrices for the system of linear equations: ---------------------
    dr = 1.0d0/(REAL(n-1, 8))

    num = 3.0d0*volPer/(100.0d0*R)
    modD = D/(R**2)
    
    fluxParam = iapp/(num*F*L*R)
    voltConIal = iapp/(num*L)
    voltConRtf = (2.0d0*Rg*T)/F

    A = 0.0d0
    B = 0.0d0

    ai = dt*modD/(2.0d0*dr*dr)

    DO i=2,n-1
      !current radius and dimensionless parameter ai
      ri = REAL((i-1),KIND=REAL64)*dr
      
      !LHS
      A(i,i-1) = ai*(-1.0d0 + dr/ri)
      A(i,i) = 1.0d0 + 2.0d0*ai
      A(i,i+1) = ai*(-1.0d0 - dr/ri)
      
      !RHS
      B(i,i-1) = ai*(1.0d0 - dr/ri)
      B(i,i) = 1.0d0 - 2.0d0*ai
      B(i,i+1) = ai*(1.0d0 + dr/ri)      
   END DO

   !smooth boundary conditions at r=0
    
    A(1,1) = 1.0d0 + 2.0d0*ai
    A(1,2) = -2.0d0*ai
    B(1,1) = 1.0d0 - 2.0d0*ai
    B(1,2) = 2.0d0*ai
    
    !flux boundary conditions
    A(n,n-1) = -2.0d0*ai
    A(n,n) = 1.0d0 + 2.0d0*ai
    B(n,n-1) = 2.0d0*ai
    B(n,n) = 1.0d0 - 2.0d0*ai

    rhsConst = 2.0d0*fluxParam*(dt + dt/dr)
    !------------------------------------------------------------------------------------------
    !-------------------------------Run the solver over the time array: -----------------------
    ALLOCATE(c_cur(n))
    ALLOCATE(rhs(n))
    c_cur = 0.0d0
    rhs = 0.0d0

    DO time = 1, totalTime

      IF (ALLOCATED(A_mod)) THEN
        DEALLOCATE(A_mod)
      END IF
      
      ALLOCATE(A_mod(n,n))
      
      
      A_mod = A
      !Use initial concentration at first timestep, otherwise use previour timestep concentration. 
      IF(time == 1) THEN
        rhs = MATMUL(B, c0)
      ELSE
        rhs = MATMUL(B, c_cur)
      END IF
   
      rhs(n) = rhs(n) - rhsConst
      
      !dgesv solver
      CALL dgesv(n,1,A_mod,n,ipiv,rhs,n,info)
      
      !Error handling - crude at this stage
      IF (info < 0) THEN
        PRINT *, 'Argument number', -info, 'in dgtsv is illegal' 
        ERROR STOP
      ELSE IF (info > 0) THEN
        PRINT *, 'Determinant of matrix A is zero in Ax=b'
        ERROR STOP
      END IF
      
      !Assign solution to output matrix
      c_cur = rhs
   
      !--------------------------Calculate stoichiometry: -----------------------------------------------
     ! IF (c_cur(n) < 0.0) THEN
      !  PRINT *, 'Concentration at the boundary has reached minimum at iteration time ', time
     !   voltArray(time) = voltArray(time-1)
     ! ELSE
      div_const = c_cur(n)/maxCon
      !------------------------U+ implemented (U- can be implemented later if needed): -----------------
      U_scalar = -0.8090d0*div_const + 4.4875d0 - 0.0428d0*TANH(18.5138d0*(div_const-0.5542d0)) &
              -17.7326d0*TANH(15.7890d0*(div_const-0.3117d0)) & 
              +17.5842d0*TANH(15.9308d0*(div_const-0.3120d0))
      !----------------------------Calculate arsinh part: ----------------------------------------------
      arsinh = F*K*SQRT(div_const - (div_const**2))
      arsinh = voltConIal/arsinh
      arsinh = ASINH(arsinh)
      !----------------------------Calculate the voltage: ----------------------------------------------
      !Checking precisions:
      voltArray(time) = U_scalar - (voltConRtf*arsinh)
     ! END IF
    END DO

    !DEALLOCATE(c_cur); DEALLOCATE(A_mod); DEALLOCATE(rhs); DEALLOCATE(A);
    !DEALLOCATE(B); DEALLOCATE(voltArray);
    
  END FUNCTION crank_nicholson

END MODULE datafitPde
