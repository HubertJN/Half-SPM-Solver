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
  integer(kind=int32)          :: output_id, check_id
  integer(kind=int32)          :: volt_out_id, conc_out_id, conc_check_id, time_check_id, ts_check_id
  logical                      :: volt_do, checkpoint
  
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

  !This subroutine reads (act=r) or writes (act=w) a single integer (sing) or a vector (vect), from/to a variable with name (var_name) or id (var_id_in)...
  !...in a NetCDF file with id (file_id)
  !Optionally you can save the variable id by inputting a variable to store it (var_id_out)
  subroutine assign_int(file_id, act, vect, sing, var_name, var_id_in, var_id_out)
    implicit none

    character(len=*),    intent(in)              :: act
    integer(kind=int32), intent(in)              :: file_id
    character(len=*),    intent(in), optional    :: var_name
    integer(kind=int32), intent(in), optional    :: var_id_in
    integer(kind=int32), intent(out), optional   :: var_id_out
    
    integer(kind=int32), intent(inout), optional :: sing, vect(:)

    integer(kind=int32)                          :: var_id, ierr

    if (present(var_name)) then
       ierr = nf90_inq_varid(file_id, var_name, var_id)
       if (ierr /= nf90_noerr) print*, var_name
       call error_check(ierr)
       
    else
       var_id = var_id_in
       
    end if

    if (present(var_id_out)) then
       var_id_out = var_id
    end if

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

  !This subroutine reads (act=r) or writes (act=w) a single real number (sing) or a vector (vect), from/to a variable with name (var_name) or id (var_id_in)...
  !...in a NetCDF file with id (file_id)
  !Optionally you can save the variable id by inputting a variable to store it (var_id_out)
  subroutine assign_real(file_id, act, vect, sing, var_name, var_id_in, var_id_out)
    implicit none

    character(len=*),    intent(in)              :: act
    integer(kind=int32), intent(in)              :: file_id
    character(len=*),    intent(in),    optional :: var_name
    integer(kind=int32), intent(in),    optional :: var_id_in
    integer(kind=int32), intent(out),   optional :: var_id_out
    real(kind=real64),   intent(inout), optional :: sing, vect(:)

    integer(kind=int32)                          :: var_id, ierr

    if (present(var_name)) then
       ierr = nf90_inq_varid(file_id, var_name, var_id)
       if (ierr /= nf90_noerr) print*, var_name
       call error_check(ierr)

    else
       var_id = var_id_in
       
    end if

    if (present(var_id_out)) then
       var_id_out = var_id
    end if
    
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

  !This subroutine writes the integer 2D array (var), to the variable named (var_name) with variable id (var_id_in), in the netcdf file with id (file_id), at position (it)
  !This should be used to write integer arrays (var), to a variable (var_name), with an infinite dimension, (var_len x undefined)
  subroutine assign_exp_int(var, file_id, it, var_name, var_id_in)
    implicit none

    character(len=*),    intent(in), optional    :: var_name
    integer(kind=int32), intent(in), optional    :: var_id_in
    integer(kind=int32), intent(in)    :: file_id, it
    
    integer(kind=int32), intent(in)    :: var(:,:)

    integer(kind=int32)                :: var_id, ierr

    if (present(var_name)) then
       ierr = nf90_inq_varid(file_id, var_name, var_id)
       if (ierr /= nf90_noerr) print*, var_name
       call error_check(ierr)
       
    else
       var_id = var_id_in
       
    end if
    
    ierr = nf90_put_var(file_id, var_id, var, (/1,it/))
    call error_check(ierr)

  end subroutine assign_exp_int

  !This subroutine writes the real 2D array (var), to the variable named (var_name) with variable id (var_id_in), in the netcdf file with id (file_id), at position (it)
  !This should be used to write real arrays (var), to a variable (var_name), with an infinite dimension, (var_len x undefined)
  subroutine assign_exp_real(var, file_id, it, var_name, var_id_in)
    implicit none

    character(len=*),    intent(in), optional    :: var_name
    integer(kind=int32), intent(in), optional    :: var_id_in
    integer(kind=int32), intent(in)              :: file_id, it
    
    real(kind=real64), intent(in)      :: var(:,:)

    integer(kind=int32)                :: var_id, ierr

    if (present(var_name)) then
       ierr = nf90_inq_varid(file_id, var_name, var_id)
       if (ierr /= nf90_noerr) print*, var_name
       call error_check(ierr)
       
    else
       var_id = var_id_in
       
    end if
    
    ierr = nf90_put_var(file_id, var_id, var, (/1,it/))
    call error_check(ierr)

  end subroutine assign_exp_real

  !This subroutine creates a variable called (var_name) with length (var_len) and data type (var_typ (f90_int or f90_double)), in a netcdf file with id (file_id)
  !You can optionally prescribe units to the variable (units), and if you want to save the variable id you can input a variable to store it (var_id_out)
  !If the file is NOT in definition mode, so already exists, you can use (act='add') to add the variable to an existing netcdf file with id (file_id)
  subroutine create_sing_var(var_name, var_typ, var_len, file_id, units, act, var_id_out)
    implicit none

    integer(kind=int32), intent(in)              :: var_typ, var_len, file_id
    character(len=*),    intent(in)              :: var_name
    character(len=*),    intent(in),    optional :: act, units
    integer(kind=int32), intent(out),   optional :: var_id_out
    integer(kind=int32)                          :: ierr, dim_id, var_id

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

    if (present(var_id_out)) then
       var_id_out = var_id
    end if
    
    if (act == 'add') then
       ierr = nf90_enddef(file_id)
       call error_check(ierr)
    end if
    
  end subroutine create_sing_var

  !This subroutine creates an expanding variable called (var_name), with dimension (var_len x undefined), and with data type (var_typ (f90_int or f90_double)), in a netcdf file with id (file_id)
  !You can optionally prescribe units to the variable (units) and if you want to save the variable id you can input a variable to store it (var_id_out)
  !If the file is NOT in definition mode, so already exists, you can use (act='add') to add the variable to an existing netcdf file with id (file_id)
  subroutine create_exp_var(var_name, var_typ, var_len, file_id, units, act, var_id_out)
    implicit none

    integer(kind=int32), intent(in)              :: var_typ, var_len, file_id
    character(len=*),    intent(in)              :: var_name
    character(len=*),    intent(in),    optional :: act, units
    integer(kind=int32), intent(out),   optional :: var_id_out
    integer(kind=int32)                          :: ierr, dim_id(2), var_id
    
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

    if (present(var_id_out)) then
       var_id_out = var_id
    end if

    if (act == 'add') then
       ierr = nf90_enddef(file_id)
       call error_check(ierr)
    end if

  end subroutine create_exp_var

  !This subroutine writes an integer vector (vect) or single number (sing) to an existing variable named (var_name) in a netcdf file with id (file_id)
  !if (act='new'), this assumes the variable doesn't already exist and will create an integer variable called (var_name), with length (var_len) and units (units)...
  !... and write the integer variable (vect or sing) to this variable
  subroutine save_int(var_name, file_id, vect, sing, units, act, var_len)
    implicit none

    character(len=*),    intent(in)           :: var_name
    integer(kind=int32), intent(in)           :: file_id
    character(len=*),    intent(in), optional :: act, units
    integer(kind=int32), intent(in), optional :: var_len
    integer(kind=int32), intent(inout), optional :: vect(:), sing

    integer(kind=int32)                       :: ierr

    if (act == 'new') then
       call create_sing_var(var_name, nf90_int, var_len, file_id, units, act='add')
    end if

    if (present(vect)) then
       call assign_int(file_id, 'w', vect=vect, var_name=var_name)
    else if (present(sing)) then
       call assign_int(file_id, 'w', sing=sing, var_name=var_name)
    end if
    
  end subroutine save_int

  !This subroutine writes a real vector (vect) or single number (sing) to an existing variable named (var_name) in a netcdf file with id (file_id)
  !if (act='new'), this assumes the variable doesn't already exist and will create a real variable called (var_name), with length (var_len) and units (units)...
  !... and write the real variable (vect or sing) to this variable
  subroutine save_real(var_name, file_id, vect, sing, units, act, var_len)
    implicit none

    character(len=*),    intent(in)             :: var_name
    integer(kind=int32), intent(in)             :: file_id
    character(len=*),    intent(in),   optional :: act, units
    integer(kind=int32), intent(in),   optional :: var_len
    real(kind=real64),   intent(inout), optional :: vect(:), sing

    integer(kind=int32)                         :: ierr

    if (act == 'new') then
       call create_sing_var(var_name, nf90_double, var_len, file_id, units, act='add')
    end if

    if (present(vect)) then
       call assign_real(file_id, 'w', vect=vect, var_name=var_name)
    else if (present(sing)) then
       call assign_real(file_id, 'w', sing=sing, var_name=var_name)
    end if
    
  end subroutine save_real

  !This subroutine closes NetCDF output and checkpoint files that are open
  subroutine fin_in_out()
    implicit none

    integer(kind=int32) :: ierr

    ierr = nf90_close(output_id)
    call error_check(ierr)
    
    ierr = nf90_close(check_id)
    call error_check(ierr)
    
  end subroutine fin_in_out

  !This subroutine opens the netcdf file with name (file_name), reads the input values, writes them to global variables, and closes the file
  !As a test case, if no file_name is given a series of test values are prescribed instead
  subroutine import_input(file_name)
    implicit none

    character(len=*), intent(in), optional :: file_name
    
    integer(kind=int32)                    :: ierr, file_id, var_id, volt_do_int, checkpoint_int

    if (present(file_name)) then
    
    ierr = nf90_open(file_name, NF90_NOWRITE, file_id)
    call error_check(ierr)

    call assign_int(file_id, 'r', sing=sim_steps, var_name='sim_steps')
    call assign_int(file_id, 'r', sing=out_steps, var_name='out_steps')
    call assign_int(file_id, 'r', sing=space_steps, var_name='space_steps')
    call assign_real(file_id, 'r', sing=temp, var_name='temp')
    call assign_real(file_id, 'r', sing=rad, var_name='rad')
    call assign_real(file_id, 'r', sing=thick, var_name='thick')
    call assign_real(file_id, 'r', sing=rr_coef, var_name='rr_coef')
    call assign_real(file_id, 'r', sing=dif_coef, var_name='dif_coef')
    call assign_real(file_id, 'r', sing=init_c, var_name='init_c')
    call assign_real(file_id, 'r', sing=max_c, var_name='max_c')
    call assign_real(file_id, 'r', sing=c_rate, var_name='c_rate')
    call assign_real(file_id, 'r', sing=dt, var_name='dt')
    call assign_real(file_id, 'r', sing=vol_per, var_name='vol_per')
    call assign_real(file_id, 'r', sing=area, var_name='area')
    call assign_int(file_id, 'r', sing=volt_do_int, var_name='volt_do')
    call assign_int(file_id, 'r', sing=checkpoint_int, var_name='checkpoint')

    volt_do = volt_do_int
    checkpoint = checkpoint_int
    
    
    ierr = nf90_close(file_id)
    call error_check(ierr)
    print*, volt_do, checkpoint
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
       volt_do = .True.
       checkpoint = .False.
    end if
    
  end subroutine import_input

  !This subroutine creates a netcdf file named (file_name), initiates input variables and saves variables that are available
  subroutine initiate_file(file_name)
    implicit none
    
    character(len=*), intent(in) :: file_name

    integer(kind=int32)          :: ierr
    
    ierr = nf90_create(file_name, nf90_netcdf4, output_id)
    call error_check(ierr)

    call create_sing_var('sim_steps',   nf90_int,    1, output_id)
    call create_sing_var('out_steps',   nf90_int,    1, output_id)
    call create_sing_var('space_steps', nf90_int,    1, output_id)
    call create_sing_var('temp',        nf90_double, 1, output_id, 'K')
    call create_sing_var('rad',         nf90_double, 1, output_id, 'm')
    call create_sing_var('thick',       nf90_double, 1, output_id, 'm')
    call create_sing_var('rr_coef',     nf90_double, 1, output_id, '$Am^{-2}(m^3mol^{-1})^{1.5}$')
    call create_sing_var('dif_coef',    nf90_double, 1, output_id, '$10^{-15} m^2 s^{-1}$')
    call create_sing_var('init_c',      nf90_double, 1, output_id, '$mol m^{-3}$')
    call create_sing_var('max_c',       nf90_double, 1, output_id, '$mol m^{-3}$')
    call create_sing_var('c_rate',      nf90_double, 1, output_id, 'A/s')
    call create_sing_var('dt',          nf90_double, 1, output_id, 's')
    call create_sing_var('vol_per',     nf90_double, 1, output_id, '%')
    call create_sing_var('area',        nf90_double, 1, output_id, '$m^2$')

    call create_exp_var('conc', nf90_double, space_steps, output_id, '$mol m^{-3}$', var_id_out=conc_out_id)

    ierr = nf90_enddef(output_id)
    call error_check(ierr)

    call assign_int(output_id, 'w', sing=sim_steps, var_name='sim_steps')
    call assign_int(output_id, 'w', sing=out_steps, var_name='out_steps')
    call assign_int(output_id, 'w', sing=space_steps, var_name='space_steps')
    call assign_real(output_id, 'w', sing=temp, var_name='temp')
    call assign_real(output_id, 'w', sing=rad, var_name='rad')
    call assign_real(output_id, 'w', sing=thick, var_name='thick')
    call assign_real(output_id, 'w', sing=rr_coef, var_name='rr_coef')
    call assign_real(output_id, 'w', sing=dif_coef, var_name='dif_coef')
    call assign_real(output_id, 'w', sing=init_c, var_name='init_c')
    call assign_real(output_id, 'w', sing=max_c, var_name='max_c')
    call assign_real(output_id, 'w', sing=c_rate, var_name='c_rate')
    call assign_real(output_id, 'w', sing=dt, var_name='dt')
    call assign_real(output_id, 'w', sing=vol_per, var_name='vol_per')
    call assign_real(output_id, 'w', sing=area, var_name='area')

    tot_steps = 0
    final_time = 0.0

  end subroutine initiate_file
  
  !This subroutine creates a netcdf file named (file_name), initiates checkpoint specific variables and saves available variables
  subroutine initiate_checkp(file_name)
    implicit none
    
    character(len=*), intent(in) :: file_name

    integer(kind=int32)          :: ierr
    
    ierr = nf90_create(file_name, nf90_netcdf4, check_id)
    call error_check(ierr)

    call create_sing_var('tot_steps',   nf90_int,    1, check_id, var_id_out=ts_check_id)
    call create_sing_var('space_steps', nf90_int,    1, check_id)
    call create_sing_var('time',        nf90_double, 1, check_id, var_id_out=time_check_id)
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

    call create_exp_var('conc', nf90_double, space_steps, check_id, '$mol m^{-3}$', var_id_out=conc_check_id)

    ierr = nf90_enddef(check_id)
    call error_check(ierr)

    call assign_int(check_id, 'w', sing=space_steps, var_name='space_steps')
    
  end subroutine initiate_checkp

  !This subroutine reads from a checkpoint file named (file_name) and saves the concentration vector to (conc)
  !It extracts other variables to keep track of the total simulation steps and total simulation time
  !It also opens an old netcdf output file and gets variable ids for the variables it will write new data to
  subroutine load_checkp(check_file_name, out_file_name, conc, volt_do)
    implicit none
    
    character(len=*),    intent(in)    :: check_file_name, out_file_name
    real(kind=real64),   intent(inout) :: conc(:)
    logical,             intent(in)    :: volt_do
    integer(kind=int32)                :: ierr

    ierr = nf90_open(check_file_name, NF90_WRITE, check_id)
    call error_check(ierr)

    call assign_real(check_id, 'r', sing=final_time, var_name='time', var_id_out=time_check_id)
    call assign_real(check_id, 'r', vect=conc, var_name='conc', var_id_out=conc_check_id)
    call assign_int(check_id, 'r', sing=tot_steps, var_name='tot_steps', var_id_out=ts_check_id)
    print*, tot_steps

    ierr = nf90_open(out_file_name, NF90_WRITE, output_id)
    call error_check(ierr)

    ierr = nf90_inq_varid(output_id, 'conc', conc_out_id)
    if (ierr /= nf90_noerr) print*, 'conc'
    call error_check(ierr)

    if (volt_do .eqv. .True.) then
       ierr = nf90_inq_varid(output_id, 'volt', volt_out_id)
       if (ierr /= nf90_noerr) print*, 'volt'
       call error_check(ierr)
    end if
    
  end subroutine load_checkp

  !This subroutine overwrites the concentration vector (conc) and number of time steps (step_num) in the checkpoint NetCDF file
  subroutine update_checkp(conc, step_num)
    implicit none
    
    real(kind=real64),   intent(inout) :: conc(:)
    integer(kind=int32), intent(inout) :: step_num

    integer(kind=int32)                :: ierr, steps
    real(kind=real64)                  :: time

    time = (step_num * dt) + final_time
    steps = step_num + out_steps - 1
    call assign_real(check_id, 'w', sing=time, var_id_in=time_check_id)
    call assign_real(check_id, 'w', vect=conc, var_id_in=conc_check_id)
    call assign_int(check_id, 'w', sing=steps, var_id_in=ts_check_id)
    
  end subroutine update_checkp
  
end module input_output_netcdf
