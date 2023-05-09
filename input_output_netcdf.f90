module input_output_netcdf
  
  use netcdf
  use iso_fortran_env
  
  implicit none

  integer,           parameter :: DP=kind(1.0D0)
  
  real(kind=Real64), parameter :: farad = 96485.3321233100184_DP !C/mol
  real(kind=Real64), parameter :: gas_con = 8.31446261815324_DP !J/(K.mol)
  
  real(kind=Real64)            :: temp, rad, thick, rr_coef, dif_coef, area
  real(kind=Real64)            :: init_c, max_c, c_rate, dt, vol_per, final_time
  integer(kind=int32)          :: sim_steps, out_steps, space_steps, tot_steps
  
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

  !This subroutine reads (act=r) or writes (act=w) a single integer (sing) or a vector (vect), from/to a variable with name (var_name) in a NetCDF file with id (file_id)
  subroutine assign_int(var_name, file_id, act, vect, sing)
    implicit none

    character(len=*),    intent(in)              :: var_name, act
    integer(kind=int32), intent(in)              :: file_id
    
    integer(kind=int32), intent(inout), optional :: vect(:)
    integer(kind=int32), intent(inout), optional :: sing

    integer(kind=int32)                          :: var_id, ierr

    ierr = nf90_inq_varid(file_id, var_name, var_id)
    if (ierr /= nf90_noerr) print*, var_name
    call error_check(ierr)

    if (act == 'r') then
       if (present(vect)) then
          ierr = nf90_get_var(file_id, var_id, vect)
          call error_check(ierr)
          
       else if (present(sing)) then
          ierr = nf90_get_var(file_id, var_id, sing)
          call error_check(ierr)

       else
          print*, 'Error reading from NetCDF file'
          Error Stop
          
       end if

    else if (act == 'w') then
       if (present(vect)) then
          ierr = nf90_put_var(file_id, var_id, vect)
          call error_check(ierr)
          
       else if (present(sing)) then
          ierr = nf90_put_var(file_id, var_id, sing)
          call error_check(ierr)

       else
          print*, 'Error writing to NetCDF file'
          Error Stop
          
       end if

    end if
    
  end subroutine assign_int

  !This subroutine reads (act=r) or writes (act=w) a single real number (sing) or vector (vect), from/to a variable with name (var_name) in a NetCDF file with id (file_id)
  subroutine assign_real(var_name, file_id, act, vect, sing)
    implicit none

    character(len=*),    intent(in)              :: var_name, act
    integer(kind=int32), intent(in)              :: file_id
    
    real(kind=real64),   intent(inout), optional :: vect(:)
    real(kind=real64),   intent(inout), optional :: sing

    integer(kind=int32)                          :: var_id, ierr

    ierr = nf90_inq_varid(file_id, var_name, var_id)
    if (ierr /= 0) print*, var_name
    call error_check(ierr)

    if (act == 'r') then
       if (present(vect)) then
          ierr = nf90_get_var(file_id, var_id, vect)
          call error_check(ierr)
          
       else if (present(sing)) then
          ierr = nf90_get_var(file_id, var_id, sing)
          call error_check(ierr)
          
       else
          print*, 'Error reading from NetCDF file'
          Error stop
       end if

       else if (act == 'w') then
          if (present(vect)) then
             ierr = nf90_put_var(file_id, var_id, vect)
             call error_check(ierr)
             
          else if (present(sing)) then
             ierr = nf90_put_var(file_id, var_id, sing)
             call error_check(ierr)
             
          else
             print*, 'Error writing to NetCDF file'
             Error stop
             
          end if
          
    end if
  end subroutine assign_real

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

  !This subroutine creates a variable called (var_name) with length (var_len) and data type (var_typ (f90_int or f90_double)), in a netcdf file with id (file_id)
  !You can optionally prescribe units to the variable (units)
  !If the file is NOT in definition mode, so already exists, you can use (act='add') to add the variable to an existing netcdf file with id (file_id)
  !Or if you want to create a single variable mid prgram, you can input the file_name and leave (file_id_in) blank
  subroutine create_sing_var(var_name, var_typ, var_len, file_id_in, units, act, file_name)
    implicit none

    integer(kind=int32), intent(in)              :: var_typ, var_len
    integer(kind=int32), intent(in),    optional :: file_id_in
    character(len=*),    intent(in)              :: var_name
    character(len=*),    intent(in),    optional :: act, units, file_name

    integer(kind=int32)                          :: ierr, dim_id, var_id, file_id

    if ((present(file_id_in) .eqv. .false.) .and. (present(file_name) .eqv. .false.)) then
       print*, 'Need file_id or file_name to create variable'
       error stop
    end if

    if (present(file_name)) then
       ierr = nf90_open(file_name, NF90_WRITE, file_id)
       call error_check(ierr)
    else
       file_id = file_id_in
    end if

    if (act == 'add') then
       ierr = nf90_redef(file_id)
       call error_check(ierr)
    end if
    
    ierr = nf90_def_dim(file_id, var_name, var_len, dim_id)
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

    if (present(file_name)) then
       ierr = nf90_close(file_id)
       call error_check(ierr)
    end if
    
  end subroutine create_sing_var

  !This subroutine creates an expanding variable called (var_name), with dimension (var_len x undefined), and with data type (var_typ (f90_int or f90_double)), in a netcdf file with id (file_id)
  !You can optionally prescribe units to the variable (units)
  !If the file is NOT in definition mode, so already exists, you can use (act='add') to add the variable to an existing netcdf file with id (file_id)
  !Or if you want to create a single variable mid prgram, you can input the file_name
  subroutine create_exp_var(var_name, var_typ, var_len, file_id_in, units, act, file_name)
    implicit none

    integer(kind=int32), intent(in)              :: var_typ, var_len
    integer(kind=int32), intent(in),    optional :: file_id_in
    character(len=*),    intent(in)              :: var_name
    character(len=*),    intent(in),    optional :: act, units, file_name

    integer(kind=int32)                          :: ierr, dim_id(2), var_id, file_id

    if ((present(file_id_in) .eqv. .false.) .and. (present(file_name) .eqv. .false.)) then
       print*, 'Need file_id or file_name to create variable'
       error stop
    end if

    if (present(file_name)) then
       ierr = nf90_open(file_name, NF90_WRITE, file_id)
       call error_check(ierr)
    else
       file_id = file_id_in
    end if
    
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

    if (present(file_name)) then
       ierr = nf90_close(file_id)
       call error_check(ierr)
    end if

  end subroutine create_exp_var

  !This subroutine opens an existing netcdf file named (file_name), and writes an integer vector (vect) or single number (sing) to an existing variable named (var_name)
  !if (act='new'), this assumes the variable doesn't already exist and will create an integer variable called (var_name), with length (var_len) and units (units)...
  !... and write the integer variable (vect or sing) to this variable
  subroutine save_int(var_name, file_name, vect, sing, units, act, var_len)
    implicit none

    character(len=*),    intent(in)           :: var_name, file_name
    character(len=*),    intent(in), optional :: act, units
    integer(kind=int32), intent(in), optional :: var_len
    integer(kind=int32), intent(inout), optional :: vect(:), sing

    integer(kind=int32)                       :: ierr, file_id

    ierr = nf90_open(file_name, NF90_WRITE, file_id)
    call error_check(ierr)

    if (act == 'new') then
       call create_sing_var(var_name, nf90_int, var_len, file_id, units, act='add')
    end if

    if (present(vect)) then
       call assign_int(var_name, file_id, 'w', vect=vect)
    else if (present(sing)) then
       call assign_int(var_name, file_id, 'w', sing=sing)
    end if
    
    ierr = nf90_close(file_id)
    call error_check(ierr)
    
  end subroutine save_int

  !This subroutine opens an existing netcdf file named (file_name), and writes a real vector (vect) or single number (sing) to an existing variable named (var_name)
  !if (act='new'), this assumes the variable doesn't already exist and will create a real variable called (var_name), with length (var_len) and units (units)...
  !... and write the real variable (vect or sing) to this variable
  subroutine save_real(var_name, file_name, vect, sing, units, act, var_len)
    implicit none

    character(len=*),    intent(in)             :: var_name, file_name
    character(len=*),    intent(in),   optional :: act, units
    integer(kind=int32), intent(in),   optional :: var_len
    real(kind=real64),   intent(inout), optional :: vect(:), sing

    integer(kind=int32)                         :: ierr, file_id

    ierr = nf90_open(file_name, NF90_WRITE, file_id)
    call error_check(ierr)

    if (act == 'new') then
       call create_sing_var(var_name, nf90_double, var_len, file_id, units, act='add')
    end if

    if (present(vect)) then
       call assign_real(var_name, file_id, 'w', vect=vect)
    else if (present(sing)) then
       call assign_real(var_name, file_id, 'w', sing=sing)
    end if

    ierr = nf90_close(file_id)
    call error_check(ierr)
    
  end subroutine save_real

  !This subroutine opens an existing netcdf file named (file_name), and writes the integer 2D array (var), to the variable named (var_name) at position (it)
  !This should be used to write integer arrays (var), to a variable (var_name), with an infinite dimension, (var_len x undefined)
  subroutine save_exp_int(var_name, var, file_name, it)
    implicit none

    character(len=*),    intent(in)           :: var_name, file_name
    integer(kind=int32), intent(in)           :: var(:,:), it

    integer(kind=int32)                       :: ierr, file_id, var_len(2)

    ierr = nf90_open(file_name, NF90_WRITE, file_id)
    call error_check(ierr)

    call assign_exp_int(var_name, var, file_id, it)

    ierr = nf90_close(file_id)
    call error_check(ierr)
    
  end subroutine save_exp_int

  !This subroutine opens an existing netcdf file named (file_name), and writes the real 2D array (var), to the variable named (var_name) at position (it)
  !This should be used to write real arrays (var), to a variable (var_name), with an infinite dimension, (var_len x undefined)
  !if (act='new'), this assumes the variable doesn't already exist and will create a real 2D array variable called (var_name), with shape (var_len x undefined), with units (units),...
  !...and then write the 2D real array (var) to this variable
  subroutine save_exp_real(var_name, var, file_name, it)
    implicit none

    character(len=*),    intent(in)           :: var_name, file_name
    real(kind=real64),   intent(in)           :: var(:,:)
    integer(kind=int32), intent(in)           :: it
    
    integer(kind=int32)                       :: ierr, file_id, var_len(2)

    ierr = nf90_open(file_name, NF90_WRITE, file_id)
    call error_check(ierr)

    call assign_exp_real(var_name, var, file_id, it)

    ierr = nf90_close(file_id)
    call error_check(ierr)
    
  end subroutine save_exp_real

  !This subroutine opens the netcdf file with name (file_name), reads the input values, writes them to global variables, and closes the file
  !As a test case, if no file_name is given a series of test values are prescribed instead
  subroutine import_input(file_name)
    implicit none

    character(len=*), intent(in), optional :: file_name
    
    integer(kind=int32)                    :: ierr, file_id, var_id

    if (present(file_name)) then
    
    ierr = nf90_open(file_name, NF90_NOWRITE, file_id)
    call error_check(ierr)

    call assign_int('sim_steps',   file_id, 'r', sing=sim_steps)
    call assign_int('out_steps',   file_id, 'r', sing=out_steps)
    call assign_int('space_steps', file_id, 'r', sing=space_steps)
    call assign_real('temp',       file_id, 'r', sing=temp)
    call assign_real('rad',        file_id, 'r', sing=rad)
    call assign_real('thick',      file_id, 'r', sing=thick)
    call assign_real('rr_coef',    file_id, 'r', sing=rr_coef)
    call assign_real('dif_coef',   file_id, 'r', sing=dif_coef)
    call assign_real('init_c',     file_id, 'r', sing=init_c)
    call assign_real('max_c',      file_id, 'r', sing=max_c)
    call assign_real('c_rate',     file_id, 'r', sing=c_rate)
    call assign_real('dt',         file_id, 'r', sing=dt)
    call assign_real('vol_per',    file_id, 'r', sing=vol_per)
    call assign_real('area',       file_id, 'r', sing=area)
    
    
    ierr = nf90_close(file_id)
    call error_check(ierr)

    else
       sim_steps = 1000
       out_steps = 5
       space_steps = 20 
       temp = 294.15_DP
       rad = 5.22e-6_DP
       thick = 75.6e-6_DP
       rr_coef = 3.42_DP
       dif_coef = 1.48e-15_DP
       init_c = 51765.0_DP
       max_c = 51765.0_DP
       c_rate = 0.001_DP
       dt = 2e-9_DP
       vol_per = 66.5_DP
       area = 0.1027_DP
    end if
    
  end subroutine import_input

  !This subroutine creates a netcdf file named (file_name), initiates input variables and saves variables that are available
  subroutine initiate_file(file_name)
    implicit none
    
    character(len=*), intent(in) :: file_name

    integer(kind=int32)          :: ierr, file_id
    
    ierr = nf90_create(file_name, nf90_netcdf4, file_id)
    call error_check(ierr)

    call create_sing_var('sim_steps',   nf90_int,    1, file_id)
    call create_sing_var('out_steps',   nf90_int,    1, file_id)
    call create_sing_var('space_steps', nf90_int,    1, file_id)
    call create_sing_var('temp',        nf90_double, 1, file_id, 'K')
    call create_sing_var('rad',         nf90_double, 1, file_id, 'm')
    call create_sing_var('thick',       nf90_double, 1, file_id, 'm')
    call create_sing_var('rr_coef',     nf90_double, 1, file_id, '$Am^{-2}(m^3mol^{-1})^{1.5}$')
    call create_sing_var('dif_coef',    nf90_double, 1, file_id, '$10^{-15} m^2 s^{-1}$')
    call create_sing_var('init_c',      nf90_double, 1, file_id, '$mol m^{-3}$')
    call create_sing_var('max_c',       nf90_double, 1, file_id, '$mol m^{-3}$')
    call create_sing_var('c_rate',      nf90_double, 1, file_id, 'A/s')
    call create_sing_var('dt',          nf90_double, 1, file_id, 's')
    call create_sing_var('vol_per',     nf90_double, 1, file_id, '%')
    call create_sing_var('area',        nf90_double, 1, file_id, '$m^2$')

    call create_exp_var('conc', nf90_real, space_steps, file_id, '$mol m^{-3}$')

    ierr = nf90_enddef(file_id)
    call error_check(ierr)

    call assign_int('sim_steps',   file_id, 'w', sing=sim_steps)
    call assign_int('out_steps',   file_id, 'w', sing=out_steps)
    call assign_int('space_steps', file_id, 'w', sing=space_steps)
    call assign_real('temp',       file_id, 'w', sing=temp)
    call assign_real('rad',        file_id, 'w', sing=rad)
    call assign_real('thick',      file_id, 'w', sing=thick)
    call assign_real('rr_coef',    file_id, 'w', sing=rr_coef)
    call assign_real('dif_coef',   file_id, 'w', sing=dif_coef)
    call assign_real('init_c',     file_id, 'w', sing=init_c)
    call assign_real('max_c',      file_id, 'w', sing=max_c)
    call assign_real('c_rate',     file_id, 'w', sing=c_rate)
    call assign_real('dt',         file_id, 'w', sing=dt)
    call assign_real('vol_per',    file_id, 'w', sing=vol_per)
    call assign_real('area',       file_id, 'w', sing=area)

    ierr = nf90_close(file_id)
    call error_check(ierr)

  end subroutine initiate_file
  
  !This subroutine creates a netcdf file named (file_name), initiates checkpoint specific variables and saves available variables
  subroutine initiate_checkp(file_name)
    implicit none
    
    character(len=*), intent(in) :: file_name

    integer(kind=int32)          :: ierr, file_id
    
    ierr = nf90_create(file_name, nf90_netcdf4, file_id)
    call error_check(ierr)

    call create_sing_var('tot_steps',   nf90_int,    1, file_id)
    call create_sing_var('space_steps', nf90_int,    1, file_id)
    !call create_sing_var('rad',         nf90_double, 1, file_id, 'm')
    !call create_sing_var('thick',       nf90_double, 1, file_id, 'm')
    !call create_sing_var('rr_coef',     nf90_double, 1, file_id, '$Am^{-2}(m^3mol^{-1})^{1.5}$')
    !call create_sing_var('dif_coef',    nf90_double, 1, file_id, '$10^{-15} m^2 s^{-1}$')
    !call create_sing_var('init_c',      nf90_double, 1, file_id, '$mol m^{-3}$')
    !call create_sing_var('max_c',       nf90_double, 1, file_id, '$mol m^{-3}$')
    !call create_sing_var('c_rate',      nf90_double, 1, file_id, 'A/s')
    !call create_sing_var('dt',          nf90_double, 1, file_id, 's')
    !call create_sing_var('vol_per',     nf90_double, 1, file_id, '%')
    !call create_sing_var('area',        nf90_double, 1, file_id, '$m^2$')

    call create_sing_var('time', nf90_real,    1, file_id)
    call create_exp_var('conc', nf90_real, space_steps, file_id, '$mol m^{-3}$')

    ierr = nf90_enddef(file_id)
    call error_check(ierr)

    call assign_int('space_steps', file_id, 'w', sing=space_steps)

    ierr = nf90_close(file_id)
    call error_check(ierr)
    
  end subroutine initiate_checkp

  !This subroutine reads from a checkpoint file named (file_name) and saves the concentration vector to (conc)
  !It also extracts other variables to keep track of the total simulation steps and total simulation time
  subroutine load_checkp(file_name, conc)
    implicit none
    
    character(len=*),    intent(in)    :: file_name
    real(kind=real64),   intent(inout) :: conc(:)
    
    integer(kind=int32)                :: ierr, file_id

    ierr = nf90_open(file_name, NF90_WRITE, file_id)
    call error_check(ierr)

    call assign_real('time', file_id, 'r', sing=final_time)
    call assign_real('conc', file_id, 'r', vect=conc)
    call assign_int('tot_steps', file_id, 'r', sing=tot_steps)

    ierr = nf90_close(file_id)
    call error_check(ierr)
    
  end subroutine load_checkp

  !This subroutine opens an already created chackpoint file named (file_name) and overwrites the concentration vector (conc) and number of time steps (step_num)
  subroutine update_checkp(file_name, conc, step_num)
    implicit none
    
    character(len=*),    intent(in)    :: file_name
    real(kind=real64),   intent(inout) :: conc(:)
    integer(kind=int32), intent(inout)    :: step_num

    integer(kind=int32)                :: ierr, file_id
    real(kind=real64)                  :: time

    ierr = nf90_open(file_name, NF90_WRITE, file_id)
    call error_check(ierr)

    time = (step_num * dt) + final_time
    call assign_real('time', file_id, 'w', sing=time)
    call assign_real('conc', file_id, 'w', vect=conc)
    call assign_int('tot_steps', file_id, 'w', sing=step_num)

    ierr = nf90_close(file_id)
    call error_check(ierr)
  end subroutine update_checkp
  
end module input_output_netcdf

!Need to modify assign_exp... to read as well for the checkpoint system, otherwise complete.
