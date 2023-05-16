MODULE pde_solver

  USE ISO_FORTRAN_ENV
  use input_output_netcdf
  !$  use omp_lib
  
  IMPLICIT NONE
  real(kind=real64) :: rhs_const, volt_con_ial, volt_con_rtf
  
CONTAINS

  subroutine setup_crank_nicholson(AL, A, AU, B)
    implicit none

    REAL(REAL64), DIMENSION(:),   intent(inout) :: AL, A, AU
    REAL(REAL64), DIMENSION(:,:), intent(inout) :: B

    REAL(REAL64)                                :: ai, ri, num, iapp, flux_param, div_const, dr
    integer(int32)                              :: n, i

    n = space_steps
    IF (n < 2) THEN
      PRINT *, "Warning: Invalid input - input array size is less than 2. Terminating." 
      STOP
    END IF

    dr = rad/(REAL(space_steps-1, KIND=REAL64))

    num = 3.0_REAL64*vol_per/(100.0_REAL64*rad)
    iapp = c_rate*dt/area
    
    flux_param = iapp/(num*farad*thick)
    volt_con_ial = iapp/(num*thick)
    volt_con_rtf = (2.0_REAL64*gas_con*temp)/farad

    AL = 0.0_REAL64
    A = 0.0_REAL64
    AU = 0.0_REAL64
    B = 0.0_REAL64

    DO i=2,n-1
      !current radius and dimensionless parameter ai
      ri = REAL((i-1),KIND=REAL64)*dr
      ai = dt*dif_coef/(2.0_REAL64*dr*ri)

      div_const = ai*(ri/dr)
      
      !LHS
      AL(i-1) = ai - div_const
      A(i) = 1.0_REAL64 + 2.0_REAL64*div_const
      AU(i) = -ai - div_const
      
      !RHS
      B(i,i-1) = div_const - ai
      B(i,i) = 1.0_REAL64 - 2.0_REAL64*div_const
      B(i,i+1) = ai + div_const
      
   END DO

   !smooth boundary conditions at r=0
    ai = dt*dif_coef/(dr*dr)
    
    A(1) = 1.0_REAL64 + ai
    AU(1) = -ai
    B(1,1) = 1.0_REAL64 - ai
    B(1,2) = ai
    
    !flux boundary conditions
    AL(n-1) = -ai
    A(n) = 1.0_REAL64 + ai
    B(n,n-1) = ai
    B(n,n) = 1.0_REAL64 - ai

    rhs_const = 2.0_REAL64*flux_param*(dt/rad + dt/dr)
    
  end subroutine setup_crank_nicholson
    
  FUNCTION crank_nicholson(AL, A, AU, B, c_cur)
  
    !solves the diffusion equation with constant diffusion coefficient
    !using the Crank-Nicholson algorithm
    !via the LAPACK library dgtsv function for tridiagonal matrices
    !evolves the given state by one timestep
    
    !rad - max radius of the geometry
    !dif_coef - the diffusion coefficient
    !dt - the timestep
    !c_cur - the current concentration vector in in out
    
    REAL(REAL64),   DIMENSION(:),                INTENT(IN) :: c_cur, AL, A, AU
    REAL(REAL64),   DIMENSION(:,:), ALLOCATABLE, intent(in) :: B
     
    REAL(REAL64),   DIMENSION(:),   ALLOCATABLE             :: crank_nicholson
    REAL(REAL64),   DIMENSION(:,:), ALLOCATABLE             :: B_mod
    REAL(REAL64),   DIMENSION(:),   ALLOCATABLE             :: AL_mod, A_mod, AU_mod, rhs
    INTEGER(INT32)                                          :: n, info, i, j
    
    !Get size of input array
    n = space_steps
    
    IF(ALLOCATED(rhs)) THEN
      DEALLOCATE(rhs)
    END IF
    
    IF (ALLOCATED(crank_nicholson)) THEN
      DEALLOCATE(crank_nicholson)
    END IF
    
    ALLOCATE(rhs(n))
    ALLOCATE(crank_nicholson(n))
    
    rhs = 0.0_REAL64
    crank_nicholson = 0.0_REAL64
    
    !call omp_set_num_threads(4)
    
    
    !generate RHS
    !$OMP parallel do default (shared) private(i,j)
    Do i=1, space_steps
       Do j=1, space_steps
          rhs(i) = rhs(i) + (B(j, i) * c_cur(j))
       end do
    end do
    
    !rhs = MATMUL(B,c_cur)
   
    rhs(n) = rhs(n) - rhs_const

    AL_mod = AL
    A_mod = A
    AU_mod = AU
    
    !dgtsv solver
    CALL dgtsv(n,1,AL_mod,A_mod,AU_mod,rhs,n,info)
    
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
    REAL(REAL64) :: U_scalar
    
    !U+ implemented (U- can be implemented later if needed)
    U_scalar = -0.8090_REAL64*x + 4.4875_REAL64 - 0.0428_REAL64*TANH(18.5138_REAL64*(x-0.5542_REAL64)) &
               -17.7326_REAL64*TANH(15.7890_REAL64*(x-0.3117_REAL64)) & 
               +17.5842_REAL64*TANH(15.9308_REAL64*(x-0.3120_REAL64))
                 
  END FUNCTION U_scalar
  
  
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
  
  
  FUNCTION volt_scalar(cin)!, Rg, T, F, iapp, a, L, K, cmax)
  
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
  
  
  FUNCTION volt_array(arrin)
    
    !Calculates an array of voltages when given 
    !an ARRAY INPUT of concentrations
    !and ESSENTIAL PARAMETERS
    
   
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
  ! div_const = arrin/max_c
    
   !arsinh = farad*rr_coef*SQRT(div_const - (div_const**2))
  ! arsinh = volt_con_ial/arsinh
   !arsinh = ASINH(arsinh)
    
   !volt_array = U_arr(div_const) - (volt_con_rtf*arsinh)
   
   !$OMP parallel do default (shared) private(i, arsinh)
   do i = 1, size_arr
      volt_array(i) = volt_scalar(arrin(i))
   end do 
   
   
    
 END FUNCTION volt_array
    
END MODULE pde_solver
