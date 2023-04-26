module input_output_netcdf
  
  use netcdf
  use iso_fortran_env
  
  implicit none

  integer,           parameter :: DP=kind(1.0D0)
  
  real(kind=Real64), parameter :: farad = 96485.3321233100184_DP !C/mol
  real(kind=Real64), parameter :: gas_con = 8.31446261815324_DP !J/(K.mol)
  
  real(kind=Real64)            :: temp, rad, thick, rr_coef, dif_coef
  real(kind=Real64)            :: init_c, max_c, c_rate
  integer(kind=int32)          :: sim_steps, out_steps
  
contains
  
  subroutine error_check(ierr)
    implicit none
    
    integer, intent(in) :: ierr
    
    if (ierr /= nf90_noerr) then
       print*, trim(nf90_strerror(ierr))
       stop
    end if
    
  end subroutine error_check

  subroutine assign_int(var_name, var, file_id)
    implicit none

    character(len=*),    intent(in)  :: var_name
    integer(kind=int32), intent(in)  :: file_id
    
    integer(kind=int32), intent(out) :: var

    integer(kind=int32)              :: var_id, ierr

    ierr = nf90_inq_varid(file_id, var_name, var_id)
    if (ierr /= 0) print*, var_name
    call error_check(ierr)

    ierr = nf90_get_var(file_id, var_id, var)
    call error_check(ierr)

  end subroutine assign_int

  subroutine assign_real(var_name, var, file_id)
    implicit none

    character(len=*),    intent(in)  :: var_name
    integer(kind=int32), intent(in)  :: file_id
    
    real(kind=real64),   intent(out) :: var

    integer(kind=int32)              :: var_id, ierr

    ierr = nf90_inq_varid(file_id, var_name, var_id)
    if (ierr /= 0) print*, var_name
    call error_check(ierr)

    ierr = nf90_get_var(file_id, var_id, var)
    call error_check(ierr)

  end subroutine assign_real
  
  subroutine import_input(file_name)
    implicit none

    character(len=*), intent(in) :: file_name
    
    integer(kind=int32)          :: ierr, file_id, var_id
    
    ierr = nf90_open(file_name, NF90_NOWRITE, file_id)
    call error_check(ierr)
    
    call assign_int('sim_steps', sim_steps, file_id)
    call assign_int('out_steps', out_steps, file_id)
    call assign_real('temp',          temp, file_id)
    call assign_real('rad',            rad, file_id)
    call assign_real('thick',        thick, file_id)
    call assign_real('rr_coef',    rr_coef, file_id)
    call assign_real('dif_coef',  dif_coef, file_id)
    call assign_real('init_c',      init_c, file_id)
    call assign_real('max_c',        max_c, file_id)
    call assign_real('c_rate',      c_rate, file_id)
    
    ierr = nf90_close(file_id)
    call error_check(ierr)
    
  end subroutine import_input
 
end module input_output_netcdf

!ierr = nf90_inquire_variable(file_id, var_id, dimids=dim_id)
    !call error_check(ierr)
    
    !ierr = nf90_inquire_dimension(file_id, dim_id(1), len=dim_len)
    !call error_check(ierr)

!if (dim_len < 4) ERROR STOP 'IMPROPER NUMBER OF VARIABLES!'

!integer, dimension(nf90_max_var_dims) :: dim_id = 0
