MODULE pde_solver
  !> @file pde.f90
  !! @brief Module file for solving the diffusion equation and voltage calculation. 
  !!
  !! @details This contains the subroutines necessary for evolving the state of the system under diffusion and 
  !! calculating the voltage for each simulation timestep. 

  USE ISO_FORTRAN_ENV
  use input_output_netcdf
  !$  use omp_lib
  
  IMPLICIT NONE
  !> @var real64 rhs_const
  !!
  !! Rescaled source at the flux boundary given by: \f$ \left(dt + \frac{dt}{dr} \right) \frac{200 R i_{app}}{3 \epsilon_{actk} F L} \f$.
  !> @var real64 volt_con_ial
  !!
  !! This is a combination of parameters given by: \f$ \frac{100 i_{app} R}{3 \epsilon_{actk} L} \f$.
  !> @var real64 volt_con_rtf
  !!
  !! This is a combination of parameters given by: \f$ \frac{2 R_g T}{F} \f$.
  !> @var real64 mod_dif
  !!
  !! The rescaled diffusion coefficient given by: \f$ \frac{D}{R^2} \f$.
  real(kind=real64) :: rhs_const, volt_con_ial, volt_con_rtf, mod_dif
  
CONTAINS

  !> @brief Subroutine to setup the matrices of the system of equations. 
  !!
  !! @details This subroutine takes in the corresponding matrices and assigns their elements as the appropriate coefficients in the 
  !! discretised diffusion equation using the Crank-Nicholson scheme.
  !!
  !! @param[in] A             The left hand side coefficient matrix of the system 
  !! of equations (please refer to formulation section
  !! for specific details on the matrix elements)
  !! @param[in] B             The right hand side coefficient matrix of the system
  !! of equations (please refer to formulation section
  !! for specific details on the matrix elements)
  subroutine setup_crank_nicholson(A, B)

    REAL(REAL64), DIMENSION(space_steps,space_steps), intent(inout) :: A, B

    REAL(REAL64)                                                    :: ai, ri, num, flux_param, dr
    integer(int32)                                                  :: n, i

    n = space_steps
    IF (n < 2) THEN
      PRINT *, "Warning: Invalid input - input array size is less than 2. Terminating." 
      STOP
    END IF

    dr = 1.0_REAL64/(REAL(space_steps-1, KIND=REAL64))

    num = 3.0_REAL64*vol_per/(100.0_REAL64*rad)
    mod_dif = dif_coef/(rad**2)
    
    flux_param = iapp/(num*farad*thick*rad)
    volt_con_ial = iapp/(num*thick)
    volt_con_rtf = (2.0_REAL64*gas_con*temp)/farad

    A = 0.0_REAL64
    B = 0.0_REAL64

    ai = dt*mod_dif/(2.0_REAL64*dr*dr)

    DO i=2,n-1
      !current radius and dimensionless parameter ai
      ri = REAL((i-1),KIND=REAL64)*dr
      
      !LHS
      A(i,i-1) = ai*(-1.0_REAL64 + dr/ri)
      A(i,i) = 1.0_REAL64 + 2.0_REAL64*ai
      A(i,i+1) = ai*(-1.0_REAL64 - dr/ri)
      
      !RHS
      B(i,i-1) = ai*(1.0_REAL64 - dr/ri)
      B(i,i) = 1.0_REAL64 - 2.0_REAL64*ai
      B(i,i+1) = ai*(1.0_REAL64 + dr/ri)      
   END DO

   !smooth boundary conditions at r=0
    
    A(1,1) = 1.0_REAL64 + 2.0_REAL64*ai
    A(1,2) = -2.0_REAL64*ai
    B(1,1) = 1.0_REAL64 - 2.0_REAL64*ai
    B(1,2) = 2.0_REAL64*ai
    
    !flux boundary conditions
    A(n,n-1) = -2.0_REAL64*ai
    A(n,n) = 1.0_REAL64 + 2.0_REAL64*ai
    B(n,n-1) = 2.0_REAL64*ai
    B(n,n) = 1.0_REAL64 - 2.0_REAL64*ai

    rhs_const = 2.0_REAL64*flux_param*(dt + dt/dr)
    
  end subroutine setup_crank_nicholson
    
  !> @brief Crank-Nicholson function.
  !!
  !! @details Solves the diffusion equation with a constant diffusion coefficient
  !! and a constant \f$ i_{app} \f$
  !! using the Crank-Nicholson algorithm.
  !! The function uses the LAPACK library, calling the dgesv function for solving systems
  !! of linear equations to evolve the given state by one timestep.
  !!
  !! @param[in] A             The left hand side coefficient matrix of the system 
  !! of equations (please refer to formulation section
  !! for specific details on the matrix elements)
  !! @param[in] B             The right hand side coefficient matrix of the system
  !! of equations (please refer to formulation section
  !! for specific details on the matrix elements)
  !! @param[in] c_cur         The current time concentration array inputed into the
  !! Crank-Nicholson solver
  FUNCTION crank_nicholson(A, B, c_cur)
  
    
    REAL(REAL64),   DIMENSION(:),                INTENT(IN) :: c_cur
    REAL(REAL64),   DIMENSION(:,:), ALLOCATABLE, intent(in) :: A, B
     
    REAL(REAL64),   DIMENSION(:),   ALLOCATABLE             :: crank_nicholson
    REAL(REAL64),   DIMENSION(:,:), ALLOCATABLE             :: A_mod, B_mod
    REAL(REAL64),   DIMENSION(:),   ALLOCATABLE             :: rhs
    INTEGER, DIMENSION(space_steps)                         :: ipiv
    INTEGER(INT32)                                          :: n, info, i, j
    
    ! Get size of input array
    n = space_steps
    
    IF(ALLOCATED(rhs)) THEN
      DEALLOCATE(rhs)
    END IF
    
    IF (ALLOCATED(crank_nicholson)) THEN
      DEALLOCATE(crank_nicholson)
    END IF

    IF (ALLOCATED(A_mod)) THEN
      DEALLOCATE(A_mod)
    END IF
    
    ALLOCATE(rhs(n))
    ALLOCATE(crank_nicholson(n))
    ALLOCATE(A_mod(n,n))
    
    rhs = 0.0_REAL64
    crank_nicholson = 0.0_REAL64
    
    !B_mod = B
    A_mod = A
    
    !generate RHS
    !$OMP parallel do default (shared) private(i,j)
    Do i=1, space_steps
       Do j=1, space_steps
          rhs(i) = rhs(i) + (B(i, j) * c_cur(j))
       end do
    end do
    
    !rhs = MATMUL(B,c_cur)
   
    rhs(n) = rhs(n) - rhs_const
    
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
    crank_nicholson = rhs
    
  END FUNCTION crank_nicholson
  
  
   
  !> @brief Function to calculate the positive electrode OCV curve, \f$ U(c) \f$.
  !!
  !! @details OCV stands for Open Circuit Voltage. Please refer to the paper, 
  !! Chang-Hui Chen et. al. 2020 J. Electrochem Soc. 167 080534, https:/doi.org/10.1149/1945-7111/ab9050.
  !!
  !! @param[in] x                The stoichiometry

  FUNCTION U_scalar(x)
  
    !Accepts SCALAR x - stoichiometry
    !Outputs U(c)
    
    REAL(REAL64), INTENT(IN) :: x
    REAL(REAL64) :: U_scalar
    
    !U+ implemented (U- can be implemented later if needed)
    U_scalar = -0.8090_REAL64*x + 4.4875_REAL64 - 0.0428_REAL64*TANH(18.5138_REAL64*(x-0.5542_REAL64)) &
               -17.7326_REAL64*TANH(15.7890_REAL64*(x-0.3117_REAL64)) & 
               +17.5842_REAL64*TANH(15.9308_REAL64*(x-0.3120_REAL64))
    ! and ESSENTIAL PARAMETERS
                 
  END FUNCTION U_scalar
  
  
  !> @brief Array version for the U_scalar function.
  !!
  !! @param x       Array version of the stoichiometry

  FUNCTION U_arr(x)
    
    !Accepts 1D ARRAY x - stoichiometry
    !Outputs U(c) as 1D array
    
    REAL(REAL64), DIMENSION(:), INTENT(IN) :: x
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: U_arr
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
  
  !> @brief Scalar voltage calculator
  !!
  !! @details Calculates scalar voltage when given 
  !! a scalar input of concentration.
  !! 
  !! @param[in] cin           input concentration
  
  FUNCTION volt_scalar(cin)
    
    REAL(REAL64), INTENT(IN) :: cin
    REAL(REAL64) :: arsinh, div_const
    REAL(REAL64) :: volt_scalar
    
    !Calculates arsinh part
    div_const = cin/max_c
    
    arsinh = farad*rr_coef*SQRT(div_const - (div_const**2))
    arsinh = volt_con_ial/arsinh
    arsinh = ASINH(arsinh)
    
    
    
    volt_scalar = U_scalar(div_const) - (volt_con_rtf*arsinh)
    
  END FUNCTION volt_scalar
 
  !> @brief Array voltage calculator
  !!
  !! @details Calculates array of voltages when given 
  !! an array input of concentration.
  !! 
  !! @param[in] arrin           array input of concentrations

  FUNCTION volt_array(arrin)
    
    !Calculates an array of voltages when given 
    !an ARRAY INPUT of concentrations
    
   
    REAL(REAL64), DIMENSION(:), INTENT(IN) :: arrin
    
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: arsinh, div_const
    REAL(REAL64), DIMENSION(:), ALLOCATABLE :: volt_array
    INTEGER :: size_arr, i

    
    IF (ALLOCATED(volt_array)) THEN
      DEALLOCATE(volt_array)
    END IF
    
    IF (ALLOCATED(arsinh)) THEN
      DEALLOCATE(arsinh)
   END IF

   IF (ALLOCATED(div_const)) THEN
      DEALLOCATE(div_const)
   END IF
    
   size_arr = SIZE(arrin)
    
   ALLOCATE(volt_array(size_arr))
   ALLOCATE(arsinh(size_arr))
   ALLOCATE(div_const(size_arr))

   !Calculates arsinh part
   !div_const = arrin/max_c
    
   !arsinh = farad*rr_coef*SQRT(div_const - (div_const**2))
   !arsinh = volt_con_ial/arsinh
   !arsinh = ASINH(arsinh)
    
   !volt_array = U_arr(div_const) - (volt_con_rtf*arsinh)
   
   !$OMP parallel do default (shared) private(i, arsinh)
   do i = 1, size_arr
      volt_array(i) = volt_scalar(arrin(i))
   end do 
   
   
    
 END FUNCTION volt_array
    
END MODULE pde_solver
