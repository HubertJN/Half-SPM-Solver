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

  !This subroutine takes in an integer error code from NetCDF (ierr), prints out the associated error, and stops the code.
  !If there is no error, nothing happens
  subroutine error_check(ierr)
    implicit none
    
    integer, intent(in) :: ierr
    
    if (ierr /= nf90_noerr) then
       print*, trim(nf90_strerror(ierr))
       stop
    end if
    
  end subroutine error_check

  !This subroutine reads (act=r) or writes (act=w) a single integer (var), from/to a variable with name (var_name) in a NetCDF file with id (file_id)
  subroutine assign_int(var_name, var, file_id, act)
    implicit none

    character(len=*),    intent(in)    :: var_name, act
    integer(kind=int32), intent(in)    :: file_id
    
    integer(kind=int32), intent(inout) :: var

    integer(kind=int32)                :: var_id, ierr

    ierr = nf90_inq_varid(file_id, var_name, var_id)
    if (ierr /= 0) print*, var_name
    call error_check(ierr)

    if (act == 'r') then
       ierr = nf90_get_var(file_id, var_id, var)
       call error_check(ierr)

    else if (act == 'w') then
       ierr = nf90_put_var(file_id, var_id, var)
       call error_check(ierr)

    end if
    
  end subroutine assign_int

  !This subroutine reads (act=r) or writes (act=w) a single real number (var), from/to a variable with name (var_name) in a NetCDF file with id (file_id)
  subroutine assign_real(var_name, var, file_id, act)
    implicit none

    character(len=*),    intent(in)    :: var_name, act
    integer(kind=int32), intent(in)    :: file_id
    
    real(kind=real64),   intent(inout) :: var

    integer(kind=int32)                :: var_id, ierr

    ierr = nf90_inq_varid(file_id, var_name, var_id)
    if (ierr /= 0) print*, var_name
    call error_check(ierr)

    if (act == 'r') then
       ierr = nf90_get_var(file_id, var_id, var)
       call error_check(ierr)

    else if (act == 'w') then
       ierr = nf90_put_var(file_id, var_id, var)
       call error_check(ierr)

    end if
  end subroutine assign_real

  !This subroutine opens the netcdf file with name (file_name), reads the input values, writes them to global variables, and closes the file
  subroutine import_input(file_name)
    implicit none

    character(len=*), intent(in) :: file_name
    
    integer(kind=int32)          :: ierr, file_id, var_id
    
    ierr = nf90_open(file_name, NF90_NOWRITE, file_id)
    call error_check(ierr)

    call assign_int('sim_steps', sim_steps, file_id, 'r')
    call assign_int('out_steps', out_steps, file_id, 'r')
    call assign_real('temp',          temp, file_id, 'r')
    call assign_real('rad',            rad, file_id, 'r')
    call assign_real('thick',        thick, file_id, 'r')
    call assign_real('rr_coef',    rr_coef, file_id, 'r')
    call assign_real('dif_coef',  dif_coef, file_id, 'r')
    call assign_real('init_c',      init_c, file_id, 'r')
    call assign_real('max_c',        max_c, file_id, 'r')
    call assign_real('c_rate',      c_rate, file_id, 'r')
    
    ierr = nf90_close(file_id)
    call error_check(ierr)
    
  end subroutine import_input

  !This subroutine creates a single value variable called (var_name), with data type (var_typ (f90_int or f90_double)), in a netcdf file with id (file_id)
  !You can optionally prescribe units to the variable (units)
  !If the file is NOT in definition mode, so already exists, you can use (act='add') to add the variable to an existing netcdf file with id (file_id)
  subroutine create_sing_var(var_name, var_typ, file_id, units, act)
    implicit none

    integer(kind=int32), intent(in)           :: file_id, var_typ
    character(len=*),    intent(in)           :: var_name
    character(len=*),    intent(in), optional :: act, units

    integer(kind=int32)                       :: ierr, dim_id, var_id

    if (act == 'add') then
       ierr = nf90_redef(file_id)
       call error_check(ierr)
    end if
    
    ierr = nf90_def_dim(file_id, var_name, 1, dim_id)
    call error_check(ierr)

    ierr = nf90_def_var(file_id, var_name, var_typ, dim_id, var_id)
    call error_check(ierr)

    if (present(units)) then
       ierr = nf90_put_att(file_id, var_id, 'Units', units)
       call error_check(ierr)
    end if

    if (act == 'add') then
       ierr = nf90_enddef(file_id)
       call error_check(ierr)
    end if
    
  end subroutine create_sing_var

  !This subroutine creates an expanding variable called (var_name), with dimension (var_len x undefined), and with data type (var_typ (f90_int or f90_double)), in a netcdf file with id (file_id)
  !You can optionally prescribe units to the variable (units)
  !If the file is NOT in definition mode, so already exists, you can use (act='add') to add the variable to an existing netcdf file with id (file_id)
  subroutine create_exp_var(var_name, var_typ, var_len, file_id, units, act)
    implicit none

    integer(kind=int32), intent(in)           :: file_id, var_typ, var_len
    character(len=*),    intent(in)           :: var_name
    character(len=*),    intent(in), optional :: act, units

    integer(kind=int32)                       :: ierr, dim_id(2), var_id

    if (act == 'add') then
       ierr = nf90_redef(file_id)
       call error_check(ierr)
    end if

    ierr = nf90_def_dim(file_id, var_name // '_x', var_len, dim_id(1))
    call error_check(ierr)

    ierr = nf90_def_dim(file_id, var_name // '_y', nf90_unlimited, dim_id(2))
    call error_check(ierr)

    ierr = nf90_def_var(file_id, var_name, var_typ, dim_id, var_id)
    call error_check(ierr)

    if (present(units)) then
       ierr = nf90_put_att(file_id, var_id, 'Units', units)
       call error_check(ierr)
    end if

    if (act == 'add') then
       ierr = nf90_enddef(file_id)
       call error_check(ierr)
    end if

  end subroutine create_exp_var

  !This subroutine writes the integer 2D array (var), to the variable named (var_name), in the netcdf file with id (file_id), at position (it)
  !This should be used to write integer arrays (var), to a variable (var_name), with an infinite dimension, (var_len x undefined)
  subroutine assign_exp_int(var_name, var, file_id, it)
    implicit none

    character(len=*),    intent(in)    :: var_name
    integer(kind=int32), intent(in)    :: file_id, it
    
    integer(kind=int32), intent(in)    :: var(:,:)

    integer(kind=int32)                :: var_id, ierr

    ierr = nf90_inq_varid(file_id, var_name, var_id)
    if (ierr /= 0) print*, var_name
    call error_check(ierr)

    ierr = nf90_put_var(file_id, var_id, var, (/1,it/))
    call error_check(ierr)

  end subroutine assign_exp_int

  !This subroutine writes the real 2D array (var), to the variable named (var_name), in the netcdf file with id (file_id), at position (it)
  !This should be used to write real arrays (var), to a variable (var_name), with an infinite dimension, (var_len x undefined)
  subroutine assign_exp_real(var_name, var, file_id, it)
    implicit none

    character(len=*),    intent(in)    :: var_name
    integer(kind=int32), intent(in)    :: file_id, it
    
    real(kind=real64), intent(in)      :: var(:,:)

    integer(kind=int32)                :: var_id, ierr

    ierr = nf90_inq_varid(file_id, var_name, var_id)
    if (ierr /= 0) print*, var_name
    call error_check(ierr)

    ierr = nf90_put_var(file_id, var_id, var, (/1,it/))
    call error_check(ierr)

  end subroutine assign_exp_real

  !This subroutine creates a netcdf file named (file_name), initiates input variables, saves input variables, and initiates variables to write output data to
  subroutine initiate_file(file_name)
    implicit none
    
    character(len=*), intent(in) :: file_name

    integer(kind=int32)          :: ierr, file_id
    
    ierr = nf90_create(file_name, nf90_clobber, file_id)
    call error_check(ierr)

    call create_sing_var('sim_steps', nf90_int,    file_id)
    call create_sing_var('out_steps', nf90_int,    file_id)
    call create_sing_var('temp',      nf90_double, file_id, 'K')
    call create_sing_var('rad',       nf90_double, file_id, 'm')
    call create_sing_var('thick',     nf90_double, file_id, 'm')
    call create_sing_var('rr_coef',   nf90_double, file_id)
    call create_sing_var('dif_coef',  nf90_double, file_id)
    call create_sing_var('init_c',    nf90_double, file_id)
    call create_sing_var('max_c',     nf90_double, file_id)
    call create_sing_var('c_rate',    nf90_double, file_id)

    ierr = nf90_enddef(file_id)
    call error_check(ierr)

    call assign_int('sim_steps', sim_steps, file_id, 'w')
    call assign_int('out_steps', out_steps, file_id, 'w')
    call assign_real('temp',          temp, file_id, 'w')
    call assign_real('rad',            rad, file_id, 'w')
    call assign_real('thick',        thick, file_id, 'w')
    call assign_real('rr_coef',    rr_coef, file_id, 'w')
    call assign_real('dif_coef',  dif_coef, file_id, 'w')
    call assign_real('init_c',      init_c, file_id, 'w')
    call assign_real('max_c',        max_c, file_id, 'w')
    call assign_real('c_rate',      c_rate, file_id, 'w')

    ierr = nf90_close(file_id)
    call error_check(ierr)

  end subroutine initiate_file

  !This subroutine opens an existing netcdf file named (file_name), and writes a single integer value (var) to an existing variable named (var_name)
  !if (act='new'), this assumes the variable doesn't already exist and will create a single integer variable called (var_name), with units (units), and write the single integer (var) to this variable
  subroutine save_int(var_name, var, file_name, units, act)
    implicit none

    character(len=*),    intent(in)           :: var_name, file_name
    character(len=*),    intent(in), optional :: act, units

    integer(kind=int32), intent(inout)        :: var

    integer(kind=int32)                       :: ierr, file_id

    ierr = nf90_open(file_name, NF90_WRITE, file_id)
    call error_check(ierr)

    if (act == 'new') then
       call create_sing_var(var_name, nf90_int, file_id, units, act='add')
    end if
    
    call assign_int(var_name, var, file_id, 'w')

    ierr = nf90_close(file_id)
    call error_check(ierr)
    
  end subroutine save_int

  !This subroutine opens an existing netcdf file named (file_name), and writes a single real value (var) to an existing variable named (var_name)
  !if (act='new'), this assumes the variable doesn't already exist and will create a single integer variable called (var_name), with units (units), and write the single integer (var) to this variable
  subroutine save_real(var_name, var, file_name, units, act)
    implicit none

    character(len=*),  intent(in)           :: var_name, file_name
    character(len=*),  intent(in), optional :: act, units

    real(kind=real64), intent(inout)        :: var

    integer(kind=int32)                     :: ierr, file_id

    ierr = nf90_open(file_name, NF90_WRITE, file_id)
    call error_check(ierr)

    if (act == 'new') then
       call create_sing_var(var_name, nf90_double, file_id, units, act='add')
    end if
    
    call assign_real(var_name, var, file_id, 'w')

    ierr = nf90_close(file_id)
    call error_check(ierr)
    
  end subroutine save_real

  !This subroutine opens an existing netcdf file named (file_name), and writes the integer 2D array (var), to the variable named (var_name) at position (it)
  !This should be used to write integer arrays (var), to a variable (var_name), with an infinite dimension, (var_len x undefined)
  !if (act='new'), this assumes the variable doesn't already exist and will create an integer 2D array variable called (var_name), with shape (var_len x undefined), with units (units),...
  !...and then write the 2D integer array (var) to this variable
  subroutine save_exp_int(var_name, var, file_name, it, units, act)
    implicit none

    character(len=*),    intent(in)           :: var_name, file_name
    integer(kind=int32), intent(in)           :: var(:,:), it
    character(len=*),    intent(in), optional :: units, act

    integer(kind=int32)                       :: ierr, file_id, var_len(2)

    ierr = nf90_open(file_name, NF90_WRITE, file_id)
    call error_check(ierr)

    if (act == 'new') then
       var_len = shape(var)
       
       call create_exp_var(var_name, nf90_int, var_len(1), file_id, units, act='add')
    end if

    call assign_exp_int(var_name, var, file_id, it)

    ierr = nf90_close(file_id)
    call error_check(ierr)
    
  end subroutine save_exp_int

  !This subroutine opens an existing netcdf file named (file_name), and writes the real 2D array (var), to the variable named (var_name) at position (it)
  !This should be used to write real arrays (var), to a variable (var_name), with an infinite dimension, (var_len x undefined)
  !if (act='new'), this assumes the variable doesn't already exist and will create a real 2D array variable called (var_name), with shape (var_len x undefined), with units (units),...
  !...and then write the 2D real array (var) to this variable
  subroutine save_exp_real(var_name, var, file_name, it, units, act)
    implicit none

    character(len=*),    intent(in)           :: var_name, file_name
    real(kind=real64),   intent(in)           :: var(:,:)
    integer(kind=int32), intent(in)           :: it
    character(len=*),    intent(in), optional :: units, act
    
    integer(kind=int32)                       :: ierr, file_id, var_len(2)

    ierr = nf90_open(file_name, NF90_WRITE, file_id)
    call error_check(ierr)

    if (act == 'new') then
       var_len = shape(var)
       
       call create_exp_var(var_name, nf90_double, var_len(1), file_id, units, act='add')
    end if

    call assign_exp_real(var_name, var, file_id, it)

    ierr = nf90_close(file_id)
    call error_check(ierr)
    
  end subroutine save_exp_real
  
end module input_output_netcdf

!Need to modify assign_exp... to read as well for the checkpoint system, otherwise complete.
