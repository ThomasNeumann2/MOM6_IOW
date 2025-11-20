!> Interface for sediment properties
module MOM_sediment_interface

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING, callTree_enter, callTree_leave
use MOM_file_parser,   only : get_param, log_version, param_file_type, log_param
use MOM_grid,          only : ocean_grid_type
use MOM_hor_index,     only : hor_index_type
use MOM_dyn_horgrid,   only : dyn_horgrid_type
use MOM_io,            only : file_exists, get_var_sizes, read_variable
use MOM_io,            only : vardesc, var_desc, slasher, MOM_read_data
use MOM_time_manager,  only : time_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_restart,       only : register_restart_pair, register_restart_field, MOM_restart_CS

implicit none ; private

#include <MOM_memory.h>

public MOM_sediment_interface_init ! Public interface to fully initialize the sediment routines
public sediment_register_restarts  ! Public interface to register sediment restart fields
public sediment_end                ! clean up in driver
public update_sediment             ! updates sediment in set_MOM
!public initialize_sediment_from_file ! Public interface to initialize sediment properties field

!> Container for all sediment related parameters
type, public :: sediment_parameters_CS ; private
  logical, public :: UseSediment = .false.  !< Flag to enable sediment features
  real, allocatable, dimension(:,:), public :: prop ! sediment propertiy type
  real, allocatable, dimension(:,:), public :: tt   ! for testing
                                   ! timing of diagnostic output.
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
  !>@{ Diagnostic handles
  integer :: id_prop = -1, id_tt = -1
  !>@}
end type sediment_parameters_CS


contains

subroutine MOM_sediment_interface_init(Time, G, US, param_file, CS, diag)
  type(time_type), target, intent(in)    :: Time       !< Model time
  type(ocean_grid_type),   intent(inout) :: G          !< Grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Input parameter structure
  type(sediment_parameters_CS), pointer  :: CS         !< Wave parameter control structure
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostic Pointer

  ! Local variables
  character(len=40)  :: mdl = "MOM_sediment_interface" !< This module's name.
  logical :: use_sediment 
  character(len=200) :: sed_config
# include "version_variable.h"

!!deallocate(CS)
  ! Dummy Check
  if (.not. associated(CS)) then
    call MOM_error(FATAL, "sediment_interface_init called without an associated control structure.")
    return
  endif

  call get_param(param_file, mdl, "USE_SEDIMENT", use_sediment, &
                 "If true, enables surface wave modules.", default=.false.)

  CS%UseSediment = use_sediment
  if (.not. use_sediment) return

  call log_version(param_file, mdl, version)

  CS%diag => diag
  CS%Time => Time

  call get_param(param_file, mdl, "SED_CONFIG", sed_config, &
                "This specifies how sediment properties are specified: \n"//&
                " \t file - read sediment information from the file \n"//&
                " \t none - set to zero, unknown sediment.",&
                default="none")

  select case ( trim(sed_config) )
    case ("file");      call  initialize_sediment_from_file(CS%prop, G, param_file, US)
    case ("none");      CS%prop = 0.0
    case default ;      CS%prop = 0.0
  end select

! prop already allocated in register_restarts

  CS%id_prop = register_diag_field('ocean_model','prop', &
               CS%diag%axesT1,Time,'Sediment type properties', 'none')

end subroutine MOM_sediment_interface_init


!> Read gridded sediment from file
subroutine initialize_sediment_from_file(S, G, param_file, US)
  type(ocean_grid_type),            intent(in)  :: G          !< Grid structure
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(out) :: S !< Sediment properties
  type(param_file_type),            intent(in)  :: param_file !< Parameter file structure
  type(unit_scale_type),            intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables
  character(len=200) :: filename, sed_file, inputdir ! Strings for file/path
  character(len=200) :: sed_varname                  ! Variable name in file
  character(len=40)  :: mdl = "initialize_sediment_from_file" ! This subroutine's name.

  call callTree_enter(trim(mdl)//"(), MOM_sediment_interface.F90")

  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "SED_FILE", sed_file, &
                 "The file from which the sediment properties is read", &
                 default="sediment.nc")
  call get_param(param_file, mdl, "SED_VARNAME", sed_varname, &
                 "The name of the sediment variable in SED_FILE.", &
                 default="sedtype")

  filename = trim(inputdir)//trim(sed_file)
  call log_param(param_file, mdl, "INPUTDIR/SED_FILE", filename)

  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_sediment_from_file: Unable to open "//trim(filename))

  S(:,:) = 0.0  ! Initializing to zero, meaning unknown properties
  call MOM_read_data(filename, trim(sed_varname), S, G%Domain)

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_sediment_from_file

!> Register sediment restart fields. To be called before MOM_sediment_interface_init
!> need to be called before MOM_sediment_interface_init and allocates
!> fields for restart
!!!!! not  used jet, probably later if sediment properties are not sttic
!!!!! e.g. roughness length
subroutine sediment_register_restarts(CS, HI, US, param_file, restart_CSp)
  type(sediment_parameters_CS),  pointer  :: CS           !< Sediment parameter control structure
  type(hor_index_type),    intent(inout)  :: HI           !< Grid structure
  type(unit_scale_type),    intent(in)    :: US           !< A dimensional unit scaling type
  type(param_file_type),    intent(in)    :: param_file   !< Input parameter structure
  type(MOM_restart_CS),     pointer       :: restart_CSp  !< Restart structure, data intent(inout)

! local variables
logical   :: use_sed

character(len=40)  :: mdl = "MOM_sediment_interface"  !< This module's name.

if (associated(CS)) then
  call MOM_error(FATAL, "sediment_register_restarts: Called with initialized sediment control structure")
endif

allocate(CS)

!!!>>>>> FOR NOW WE RETURN HERE  <<<<<<<!!!!!
!!return

call get_param(param_file, mdl, "USE_SEDIMENT", use_sed, &
     "If true, enables sediment modules.", do_not_log=.true., default=.false.)

if (use_sed) then
   allocate(CS%prop(HI%isd:HI%ied,HI%jsd:HI%jed), source=0.0)

   call register_restart_field(CS%prop(:,:), "prop", .false., restart_CSp, &
        "sediment type", "none")
endif

end subroutine sediment_register_restarts

!> update sediment parameters
!>
subroutine update_sediment(CS, G)
    type(sediment_parameters_CS), pointer       :: CS    !< Sediment parameter Control structure
    type(ocean_grid_type), intent(inout) :: G   !< Grid structure

    if (CS%id_prop > 0) then
      call post_data(CS%id_prop, CS%prop, CS%diag)
    endif

end subroutine update_sediment

!> clean up
subroutine sediment_end(CS)
  type(sediment_parameters_CS),  pointer  :: CS           !< Sediment parameter control structure

  if (allocated(CS%prop))     deallocate(CS%prop)

  deallocate( CS )

end subroutine sediment_end



end module MOM_sediment_interface




