!> Configures the ISOMIP test case.
module ISOMIP_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_ALE_sponge, only : ALE_sponge_CS, set_up_ALE_sponge_field, initialize_ALE_sponge
use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, MOM_read_data, slasher
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA
use regrid_consts, only : REGRIDDING_SIGMA_SHELF_ZSTAR
implicit none ; private

#include <MOM_memory.h>

character(len=40) :: mdl = "ISOMIP_initialization" !< This module's name.

! The following routines are visible to the outside world
public ISOMIP_initialize_topography
public ISOMIP_initialize_thickness
public ISOMIP_initialize_temperature_salinity
public ISOMIP_initialize_sponges

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> Initialization of topography for the ISOMIP configuration
subroutine ISOMIP_initialize_topography(D, G, param_file, max_depth, US)
  type(dyn_horgrid_type),          intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                   intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth [Z ~> m]
  type(unit_scale_type),           intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables
  real :: min_depth       ! The minimum depth of the ocean [Z ~> m].
  ! The following variables are used to set up the bathymetry in the ISOMIP example.
  real :: bmax            ! maximum depth of bedrock topography [Z ~> m]
  real :: b0, b2, b4, b6  ! first, second, third and fourth bedrock topography coeffs [Z ~> m]
  real :: xbar            ! characteristic along-flow length scale of the bedrock [L ~> m]
  real :: dc              ! depth of the trough compared with side walls [Z ~> m].
  real :: fc              ! characteristic width of the side walls of the channel [L ~> m]
  real :: wc              ! half-width of the trough [L ~> m]
  real :: ly              ! domain width (across ice flow) [L ~> m]
  real :: bx, by          ! The x- and y- contributions to the bathymetric profiles at a point [Z ~> m]
  real :: xtil            ! x-positon normalized by the characteristic along-flow length scale [nondim]
  logical :: is_2D        ! If true, use a 2D setup
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "ISOMIP_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call MOM_mesg("  ISOMIP_initialization.F90, ISOMIP_initialize_topography: setting topography", 5)

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "ISOMIP_2D", is_2D, 'If true, use a 2D setup.', default=.false.)
  call get_param(param_file, mdl, "ISOMIP_MAX_BEDROCK", bmax, &
                 "Maximum depth of bedrock topography in the ISOMIP configuration.", &
                 units="m", default=720.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "ISOMIP_TROUGH_DEPTH", dc, &
                 "Depth of the trough compared with side walls in the ISOMIP configuration.", &
                 units="m", default=500.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "ISOMIP_BEDROCK_LENGTH", xbar, &
                 "Characteristic along-flow length scale of the bedrock in the ISOMIP configuration.", &
                 units="m", default=300.0e3, scale=US%m_to_L)
  call get_param(param_file, mdl, "ISOMIP_TROUGH_WIDTH", wc, &
                 "Half-width of the trough in the ISOMIP configuration.", &
                 units="m", default=24.0e3, scale=US%m_to_L)
  call get_param(param_file, mdl, "ISOMIP_DOMAIN_WIDTH", ly, &
                 "Domain width (across ice flow) in the ISOMIP configuration.", &
                 units="m", default=80.0e3, scale=US%m_to_L)
  call get_param(param_file, mdl, "ISOMIP_SIDE_WIDTH", fc, &
                 "Characteristic width of the side walls of the channel in the ISOMIP configuration.", &
                 units="m", default=4.0e3, scale=US%m_to_L)

  if (G%grid_unit_to_L <= 0.) call MOM_error(FATAL, "ISOMIP_initialization.F90: " //&
          "ISOMIP_initialize_topography is only set to work with Cartesian axis units.")

  ! The following variables should be transformed into runtime parameters.
  b0 = -150.0*US%m_to_Z ; b2 = -728.8*US%m_to_Z ; b4 = 343.91*US%m_to_Z ; b6 = -50.57*US%m_to_Z

  if (is_2D) then
    do j=js,je ; do i=is,ie
      ! For the 2D setup take a slice through the middle of the domain
      xtil = G%geoLonT(i,j)*G%grid_unit_to_L / xbar
      ! xtil = 450.0e3*US%m_to_L / xbar
      bx = b0 + b2*xtil**2 + b4*xtil**4 + b6*xtil**6

      by = 2.0 * dc / (1.0 + exp(2.0*wc / fc))

      D(i,j) = -max(bx+by, -bmax)
      if (D(i,j) > max_depth) D(i,j) = max_depth
      if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
    enddo ; enddo

  else
    do j=js,je ; do i=is,ie
      ! 3D setup
      ! ===== TEST =====
      !if (G%geoLonT(i,j)<500.) then
      !  xtil = 500.0e3*US%m_to_L / xbar
      !else
      !  xtil = G%geoLonT(i,j)*G%grid_unit_to_L / xbar
      !endif
      ! ===== TEST =====

      xtil = G%geoLonT(i,j)*G%grid_unit_to_L / xbar

      bx = b0 + b2*xtil**2 + b4*xtil**4 + b6*xtil**6
      by = (dc / (1.0 + exp(-2.*(G%geoLatT(i,j)*G%grid_unit_to_L - 0.5*ly - wc) / fc))) + &
           (dc / (1.0 + exp(2.*(G%geoLatT(i,j)*G%grid_unit_to_L - 0.5*ly + wc) / fc)))

      D(i,j) = -max(bx+by, -bmax)
      if (D(i,j) > max_depth) D(i,j) = max_depth
      if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
    enddo ; enddo
  endif

end subroutine ISOMIP_initialize_topography

!> Initialization of thicknesses
subroutine ISOMIP_initialize_thickness ( h, depth_tot, G, GV, US, param_file, tv, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure to parse for model parameter values
  type(thermo_var_ptrs),   intent(in)  :: tv          !< A structure containing pointers to any
                                                      !! available thermodynamic fields, including
                                                      !! the eqn. of state.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.
  ! Local variables
  real :: e0(SZK_(GV)+1)  ! The resting interface heights, in depth units [Z ~> m],
                          !  usually negative because it is positive upward.
  real :: eta1D(SZK_(GV)+1) ! Interface height relative to the sea surface
                            ! positive upward, in depth units [Z ~> m].
  integer :: i, j, k, is, ie, js, je, nz
  real    :: min_thickness     ! Minimum layer thicknesses [Z ~> m]
  real    :: S_sur, S_bot      ! Surface and bottom salinities [S ~> ppt]
  real    :: T_sur, T_bot      ! Surface and bottom temperatures [C ~> degC]
  real    :: rho_sur, rho_bot  ! Surface and bottom densities [R ~> kg m-3]
  real    :: rho_range    ! The range of densities [R ~> kg m-3]
  !character(len=256) :: mesg  ! The text of an error message
  character(len=40) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.just_read) &
    call MOM_mesg("ISOMIP_initialization.F90, ISOMIP_initialize_thickness: setting thickness")

  call get_param(param_file, mdl,"MIN_THICKNESS", min_thickness, &
                 'Minimum layer thickness', units='m', default=1.e-3, do_not_log=just_read, scale=US%m_to_Z)
  call get_param(param_file, mdl,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)

  select case ( coordinateMode(verticalCoordinate) )

  case ( REGRIDDING_LAYER, REGRIDDING_RHO ) ! Initial thicknesses for isopycnal coordinates
    call get_param(param_file, mdl, "ISOMIP_T_SUR", t_sur, &
                   "Temperature at the surface (interface)", &
                   units="degC", default=-1.9, scale=US%degC_to_C, do_not_log=just_read)
    call get_param(param_file, mdl, "ISOMIP_S_SUR", s_sur, &
                   "Salinity at the surface (interface)", &
                   units="ppt", default=33.8, scale=US%ppt_to_S, do_not_log=just_read)
    call get_param(param_file, mdl, "ISOMIP_T_BOT", t_bot, &
                   "Temperature at the bottom (interface)", &
                   units="degC", default=-1.9, scale=US%degC_to_C, do_not_log=just_read)
    call get_param(param_file, mdl, "ISOMIP_S_BOT", s_bot, &
                   "Salinity at the bottom (interface)", &
                   units="ppt", default=34.55, scale=US%ppt_to_S, do_not_log=just_read)

    if (just_read) return ! All run-time parameters have been read, so return.

    ! Compute min/max density using T_SUR/S_SUR and T_BOT/S_BOT
    call calculate_density(T_sur, S_sur, 0.0, rho_sur, tv%eqn_of_state)
    ! write(mesg,*) 'Surface density is:', rho_sur
    ! call MOM_mesg(mesg,5)
    call calculate_density(T_bot, S_bot, 0.0, rho_bot, tv%eqn_of_state)
    ! write(mesg,*) 'Bottom density is:', rho_bot
    ! call MOM_mesg(mesg,5)
    rho_range = rho_bot - rho_sur
    ! write(mesg,*) 'Density range is:', rho_range
    ! call MOM_mesg(mesg,5)

    ! Construct notional interface positions
    e0(1) = 0.
    do K=2,nz
      e0(k) = -G%max_depth * ( 0.5 * ( GV%Rlay(k-1) + GV%Rlay(k) ) - rho_sur ) / rho_range
      e0(k) = min( 0., e0(k) ) ! Bound by surface
      e0(k) = max( -G%max_depth, e0(k) ) ! Bound by possible deepest point in model
      ! write(mesg,*) 'G%max_depth,GV%Rlay(k-1),GV%Rlay(k),e0(k)', &
      !     G%max_depth,GV%Rlay(k-1),GV%Rlay(k),e0(k)
      ! call MOM_mesg(mesg,5)
    enddo
    e0(nz+1) = -G%max_depth

    ! Calculate thicknesses
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -depth_tot(i,j)
      do k=nz,1,-1
        eta1D(k) = e0(k)
        if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_Z)) then
          eta1D(k) = eta1D(k+1) + GV%Angstrom_Z
          h(i,j,k) = GV%Angstrom_Z
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
    enddo ; enddo

  case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA_SHELF_ZSTAR )   ! Initial thicknesses for z coordinates
    if (just_read) return ! All run-time parameters have been read, so return.
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -depth_tot(i,j)
      do k=nz,1,-1
        eta1D(k) =  -G%max_depth * real(k-1) / real(nz)
        if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
          eta1D(k) = eta1D(k+1) + min_thickness
          h(i,j,k) = min_thickness
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
    enddo ; enddo

  case ( REGRIDDING_SIGMA )             ! Initial thicknesses for sigma coordinates
    if (just_read) return ! All run-time parameters have been read, so return.
    do j=js,je ; do i=is,ie
      h(i,j,:) = depth_tot(i,j) / real(nz)
    enddo ; enddo

  case default
      call MOM_error(FATAL,"isomip_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

end subroutine ISOMIP_initialize_thickness

!> Initial values for temperature and salinity
subroutine ISOMIP_initialize_temperature_salinity ( T, S, h, depth_tot, G, GV, US, param_file, &
                                                    eqn_of_state, just_read)
  type(ocean_grid_type),                     intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< Vertical grid structure
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< Potential temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< Salinity [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h  !< Layer thickness [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)  :: depth_tot  !< The nominal total bottom-to-top
                                                               !! depth of the ocean [Z ~> m]
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  type(EOS_type),                            intent(in)  :: eqn_of_state !< Equation of state structure
  logical,                                   intent(in)  :: just_read !< If true, this call will
                                                      !! only read parameters without changing T & S.
  ! Local variables
  real      :: rho_sur, rho_bot  ! Surface and bottom densities [R ~> kg m-3]
  real      :: xi0, xi1 ! Heights in depth units [Z ~> m].
  real      :: S_sur, S_bot ! Salinity at the surface and bottom [S ~> ppt]
  real      :: T_sur, T_bot ! Temperature at the surface and bottom [C ~> degC]
  real      :: dT_dz  ! Vertical gradient of temperature [C Z-1 ~> degC m-1].
  real      :: dS_dz  ! Vertical gradient of salinity [S Z-1 ~> ppt m-1].
  real :: T0(SZK_(GV))       ! A profile of temperatures [C ~> degC]
  real :: S0(SZK_(GV))       ! A profile of salinities [S ~> ppt]
  real :: drho_dT(SZK_(GV))  ! Derivative of density with temperature [R C-1 ~> kg m-3 degC-1].
  real :: drho_dS(SZK_(GV))  ! Derivative of density with salinity [R S-1 ~> kg m-3 ppt-1].
  real :: rho_guess(SZK_(GV)) ! Potential density at T0 & S0 [R ~> kg m-3].
  real :: pres(SZK_(GV))     ! An array of the reference pressure [R L2 T-2 ~> Pa]. (zero here)
  real :: drho_dT1           ! A prescribed derivative of density with temperature [R C-1 ~> kg m-3 degC-1]
  real :: drho_dS1           ! A prescribed derivative of density with salinity [R S-1 ~> kg m-3 ppt-1].
  real :: T_ref              ! Default value for other temperatures [C ~> degC]
  real :: S_ref              ! Default value for other salinities [S ~> ppt]
  logical :: fit_salin       ! If true, accept the prescribed temperature and fit the salinity.
  !real :: rho_tmp    ! A temporary density used for debugging [R ~> kg m-3]
  !character(len=256) :: mesg ! The text of an error message
  character(len=40) :: verticalCoordinate
  integer   :: i, j, k, is, ie, js, je, nz, itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  pres(:) = 0.0

  call get_param(param_file, mdl, "REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)
  call get_param(param_file, mdl, "ISOMIP_T_SUR",t_sur, &
                 "Temperature at the surface (interface)", &
                 units="degC", default=-1.9, scale=US%degC_to_C, do_not_log=just_read)
  call get_param(param_file, mdl, "ISOMIP_S_SUR", s_sur, &
                 "Salinity at the surface (interface)", &
                 units="ppt", default=33.8, scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl, "ISOMIP_T_BOT", t_bot, &
                 "Temperature at the bottom (interface)", &
                 units="degC", default=-1.9, scale=US%degC_to_C, do_not_log=just_read)
  call get_param(param_file, mdl, "ISOMIP_S_BOT", s_bot, &
                 "Salinity at the bottom (interface)", &
                 units="ppt", default=34.55, scale=US%ppt_to_S, do_not_log=just_read)

  call calculate_density(T_sur, S_sur, 0.0, rho_sur, eqn_of_state)
  ! write(mesg,*) 'Density in the surface layer:', rho_sur
  ! call MOM_mesg(mesg,5)
  call calculate_density(T_bot, S_bot, 0.0, rho_bot, eqn_of_state)
  ! write(mesg,*) 'Density in the bottom layer::', rho_bot
  ! call MOM_mesg(mesg,5)

  select case ( coordinateMode(verticalCoordinate) )

    case (  REGRIDDING_RHO, REGRIDDING_ZSTAR, REGRIDDING_SIGMA_SHELF_ZSTAR, REGRIDDING_SIGMA )
      if (just_read) return ! All run-time parameters have been read, so return.

      dS_dz = (S_sur - S_bot) / G%max_depth
      dT_dz = (T_sur - T_bot) / G%max_depth
      do j=js,je ; do i=is,ie
        xi0 = -depth_tot(i,j)
        do k = nz,1,-1
          xi0 = xi0 + 0.5 * h(i,j,k)     ! Depth in middle of layer
          S(i,j,k) = S_sur + dS_dz * xi0
          T(i,j,k) = T_sur + dT_dz * xi0
          xi0 = xi0 + 0.5 * h(i,j,k)     ! Depth at top of layer
        enddo
      enddo ; enddo

    case ( REGRIDDING_LAYER )
      call get_param(param_file, mdl, "FIT_SALINITY", fit_salin, &
                  "If true, accept the prescribed temperature and fit the "//&
                  "salinity; otherwise take salinity and fit temperature.", &
                  default=.false., do_not_log=just_read)
      call get_param(param_file, mdl, "DRHO_DS", drho_dS1, &
                  "Partial derivative of density with salinity.", &
                  units="kg m-3 ppt-1", scale=US%kg_m3_to_R*US%S_to_ppt, &
                  fail_if_missing=.not.just_read, do_not_log=just_read)
      call get_param(param_file, mdl, "DRHO_DT", drho_dT1, &
                  "Partial derivative of density with temperature.", &
                  units="kg m-3 K-1", scale=US%kg_m3_to_R*US%C_to_degC, &
                  fail_if_missing=.not.just_read, do_not_log=just_read)
      call get_param(param_file, mdl, "T_REF", T_Ref, &
                  "A reference temperature used in initialization.", &
                  units="degC", scale=US%degC_to_C, fail_if_missing=.not.just_read, do_not_log=just_read)
      call get_param(param_file, mdl, "S_REF", S_Ref, &
                  "A reference salinity used in initialization.", &
                  units="ppt", default=35.0, scale=US%ppt_to_S, do_not_log=just_read)
      if (just_read) return ! All run-time parameters have been read, so return.

      ! write(mesg,*) 'read drho_dS, drho_dT', drho_dS1, drho_dT1
      ! call MOM_mesg(mesg,5)

      dS_dz = (S_sur - S_bot) / G%max_depth
      dT_dz = (T_sur - T_bot) / G%max_depth

      do j=js,je ; do i=is,ie
        xi0 = 0.0
        do k = 1,nz
          !T0(k) = T_Ref; S0(k) = S_Ref
          xi1 = xi0 + 0.5 * h(i,j,k)
          S0(k) = S_sur - dS_dz * xi1
          T0(k) = T_sur - dT_dz * xi1
          xi0 = xi0 + h(i,j,k)
          ! write(mesg,*) 'S,T,xi0,xi1,k',S0(k),T0(k),xi0,xi1,k
          ! call MOM_mesg(mesg,5)
        enddo

        call calculate_density_derivs(T0, S0, pres, drho_dT, drho_dS, eqn_of_state, (/1,1/) )
        ! write(mesg,*) 'computed drho_dS, drho_dT', drho_dS(1), drho_dT(1)
        ! call MOM_mesg(mesg,5)
        call calculate_density(T0(1), S0(1), pres(1), rho_guess(1), eqn_of_state)

        if (fit_salin) then
          ! A first guess of the layers' salinity.
          do k=nz,1,-1
            S0(k) = max(0.0, S0(1) + (GV%Rlay(k) - rho_guess(1)) / drho_dS1)
          enddo
          ! Refine the guesses for each layer.
          do itt=1,6
            call calculate_density(T0, S0, pres, rho_guess, eqn_of_state)
            call calculate_density_derivs(T0, S0, pres, drho_dT, drho_dS, eqn_of_state)
            do k=1,nz
              S0(k) = max(0.0, S0(k) + (GV%Rlay(k) - rho_guess(k)) / drho_dS1)
            enddo
          enddo

        else
          ! A first guess of the layers' temperatures.
          do k=nz,1,-1
            T0(k) = T0(1) + (GV%Rlay(k) - rho_guess(1)) / drho_dT1
          enddo

          do itt=1,6
            call calculate_density(T0, S0, pres, rho_guess, eqn_of_state)
            call calculate_density_derivs(T0, S0, pres, drho_dT, drho_dS, eqn_of_state)
            do k=1,nz
              T0(k) = T0(k) + (GV%Rlay(k) - rho_guess(k)) / drho_dT(k)
            enddo
          enddo
        endif

        do k=1,nz
          T(i,j,k) = T0(k) ; S(i,j,k) = S0(k)
        enddo

      enddo ; enddo

    case default
      call MOM_error(FATAL,"isomip_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

  ! for debugging
  !i=G%iec; j=G%jec
  !do k = 1,nz
  !  call calculate_density(T(i,j,k), S(i,j,k),0.0,rho_tmp,eqn_of_state, scale=US%kg_m3_to_R)
  !  write(mesg,*) 'k,h,T,S,rho,Rlay',k,US%Z_to_m*h(i,j,k),US%C_to_degC*T(i,j,k),US%S_to_ppt*S(i,j,k),rho_tmp,GV%Rlay(k)
  !  call MOM_mesg(mesg,5)
  !enddo

end subroutine ISOMIP_initialize_temperature_salinity

!> Sets up the the inverse restoration time (Idamp), and
! the values towards which the interface heights and an arbitrary
! number of tracers should be restored within each sponge.
subroutine ISOMIP_initialize_sponges(G, GV, US, tv, depth_tot, PF, use_ALE, CSp, ACSp)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure containing pointers to any available
                                              !! thermodynamic fields, potential temperature and
                                              !! salinity or mixed layer density.
                                              !! Absent fields have NULL ptrs.
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: depth_tot !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in) :: PF   !< A structure to parse for model parameter values
  logical,                 intent(in) :: use_ALE !< If true, indicates model is in ALE mode
  type(sponge_CS),         pointer    :: CSp  !< Layer-mode sponge structure
  type(ALE_sponge_CS),     pointer    :: ACSp !< ALE-mode sponge structure
  ! Local variables
  real :: T(SZI_(G),SZJ_(G),SZK_(GV)) ! A temporary array for temp [C ~> degC]
  real :: S(SZI_(G),SZJ_(G),SZK_(GV)) ! A temporary array for salt [S ~> ppt]
  ! real :: RHO(SZI_(G),SZJ_(G),SZK_(GV)) ! A temporary array for RHO [R ~> kg m-3]
  real :: dz(SZI_(G),SZJ_(G),SZK_(GV)) ! Sponge layer thicknesses in height units [Z ~> m]
  real :: Idamp(SZI_(G),SZJ_(G))    ! The sponge damping rate [T-1 ~> s-1]
  real :: TNUDG                     ! Nudging time scale [T ~> s]
  real :: S_sur, S_bot              ! Surface and bottom salinities in the sponge region [S ~> ppt]
  real :: T_sur, T_bot              ! Surface and bottom temperatures in the sponge region [C ~> degC]
  real :: T_ref                     ! Default value for other temperatures [C ~> degC]
  real :: S_ref                     ! Default value for other salinities [S ~> ppt]
  real :: rho_sur, rho_bot          ! Surface and bottom densities [R ~> kg m-3]
  real :: rho_range                 ! The range of densities [R ~> kg m-3]
  real :: dT_dz                     ! Vertical gradient of temperature [C Z-1 ~> degC m-1]
  real :: dS_dz                     ! Vertical gradient of salinity [S Z-1 ~> ppt m-1]

  real :: e0(SZK_(GV)+1)            ! The resting interface heights [Z ~> m], usually
                                    ! negative because it is positive upward.
  real :: eta1D(SZK_(GV)+1)         ! Interface height relative to the sea surface, positive upward [Z ~> m].
  real :: eta(SZI_(G),SZJ_(G),SZK_(GV)+1) ! A temporary array for interface heights [Z ~> m].
  real :: min_depth                 ! The minimum depth of the ocean [Z ~> m]
  real :: min_thickness             ! The minimum layer thickness [Z ~> m]
  real :: xi0                       ! Interface heights in depth units [Z ~> m], usually negative.
  !real :: rho_tmp                   ! A temporary density used for debugging [R ~> kg m-3]
  character(len=40) :: verticalCoordinate, filename, state_file
  character(len=40) :: temp_var, salt_var, eta_var, inputdir

  character(len=40)  :: mdl = "ISOMIP_initialize_sponges" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call get_param(PF, mdl, "MIN_THICKNESS", min_thickness, "Minimum layer thickness", &
                 units="m", default=1.e-3, scale=US%m_to_Z)

  call get_param(PF, mdl, "REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
            default=DEFAULT_COORDINATE_MODE)

  call get_param(PF, mdl, "ISOMIP_TNUDG", TNUDG, "Nudging time scale for sponge layers", &
                 units="days", default=0.0, scale=86400.0*US%s_to_T)

  call get_param(PF, mdl, "T_REF", T_ref, "Reference temperature", &
                 units="degC", default=10.0, scale=US%degC_to_C, do_not_log=.true.)

  call get_param(PF, mdl, "S_REF", s_ref, "Reference salinity", &
                 units="ppt", default=35.0, scale=US%ppt_to_S, do_not_log=.true.)

  call get_param(PF, mdl, "ISOMIP_S_SUR_SPONGE", s_sur, &
                 "Surface salinity in sponge layer.", &
                 units="ppt", default=US%S_to_ppt*S_ref, scale=US%ppt_to_S)

  call get_param(PF, mdl, "ISOMIP_S_BOT_SPONGE", s_bot, &
                 "Bottom salinity in sponge layer.", &
                 units="ppt", default=US%S_to_ppt*S_ref, scale=US%ppt_to_S)

  call get_param(PF, mdl, "ISOMIP_T_SUR_SPONGE", t_sur, &
                 "Surface temperature in sponge layer.", &
                 units="degC", default=US%C_to_degC*T_ref, scale=US%degC_to_C)

  call get_param(PF, mdl, "ISOMIP_T_BOT_SPONGE", t_bot, &
                 "Bottom temperature in sponge layer.", &
                 units="degC", default=US%C_to_degC*T_ref, scale=US%degC_to_C)

  T(:,:,:) = 0.0 ; S(:,:,:) = 0.0 ; Idamp(:,:) = 0.0 !; RHO(:,:,:) = 0.0

!   Set up sponges for ISOMIP configuration
  call get_param(PF, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0, scale=US%m_to_Z)

  if (associated(CSp)) call MOM_error(FATAL, &
        "ISOMIP_initialize_sponges called with an associated control structure.")
  if (associated(ACSp)) call MOM_error(FATAL, &
        "ISOMIP_initialize_sponges called with an associated ALE-sponge control structure.")

  !  Here the inverse damping time [T-1 ~> s-1], is set. Set Idamp to 0
  !  wherever there is no sponge, and the subroutines that are called
  !  will automatically set up the sponges only where Idamp is positive
  !  and mask2dT is 1.

  do j=js,je ; do i=is,ie
    if (depth_tot(i,j) <= min_depth) then
      Idamp(i,j) = 0.0
    elseif (G%geoLonT(i,j) >= 790.0 .AND. G%geoLonT(i,j) <= 800.0) then
      Idamp(i,j) = (1.0/TNUDG) * max(0.0, (G%geoLonT(i,j)-790.0) / (800.0-790.0))
    else
      Idamp(i,j) = 0.0
    endif

  enddo ; enddo

  ! Compute min/max density using T_SUR/S_SUR and T_BOT/S_BOT
  call calculate_density(T_sur, S_sur, 0.0, rho_sur, tv%eqn_of_state)
  !write (mesg,*) 'Surface density in sponge:', rho_sur
  ! call MOM_mesg(mesg,5)
  call calculate_density(T_bot, S_bot, 0.0, rho_bot, tv%eqn_of_state)
  !write (mesg,*) 'Bottom density in sponge:', rho_bot
  ! call MOM_mesg(mesg,5)
  rho_range = rho_bot - rho_sur
  !write (mesg,*) 'Density range in sponge:', rho_range
  ! call MOM_mesg(mesg,5)

  if (use_ALE) then

    select case ( coordinateMode(verticalCoordinate) )

      case ( REGRIDDING_RHO )
        ! Construct notional interface positions
        e0(1) = 0.
        do K=2,nz
          e0(k) = -G%max_depth * ( 0.5 * ( GV%Rlay(k-1) + GV%Rlay(k) ) - rho_sur ) / rho_range
          e0(k) = min( 0., e0(k) ) ! Bound by surface
          e0(k) = max( -G%max_depth, e0(k) ) ! Bound by possible deepest point in model
          ! write(mesg,*) 'G%max_depth,GV%Rlay(k-1),GV%Rlay(k),e0(k)',&
          !       G%max_depth,GV%Rlay(k-1),GV%Rlay(k),e0(k)
          ! call MOM_mesg(mesg,5)
        enddo
        e0(nz+1) = -G%max_depth

        ! Calculate thicknesses
        do j=js,je ; do i=is,ie
          eta1D(nz+1) = -depth_tot(i,j)
          do k=nz,1,-1
            eta1D(k) = e0(k)
            if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_Z)) then
              eta1D(k) = eta1D(k+1) + GV%Angstrom_Z
              dz(i,j,k) = GV%Angstrom_Z
            else
              dz(i,j,k) = eta1D(k) - eta1D(k+1)
            endif
          enddo
        enddo ; enddo

      case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA_SHELF_ZSTAR )  ! Initial thicknesses for z coordinates
        do j=js,je ; do i=is,ie
          eta1D(nz+1) = -depth_tot(i,j)
          do k=nz,1,-1
            eta1D(k) =  -G%max_depth * real(k-1) / real(nz)
            if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
              eta1D(k) = eta1D(k+1) + min_thickness
              dz(i,j,k) = min_thickness
            else
              dz(i,j,k) = eta1D(k) - eta1D(k+1)
            endif
          enddo
        enddo ; enddo

      case ( REGRIDDING_SIGMA )             ! Initial thicknesses for sigma coordinates
        do j=js,je ; do i=is,ie
          dz(i,j,:) = depth_tot(i,j) / real(nz)
        enddo ; enddo

      case default
         call MOM_error(FATAL,"ISOMIP_initialize_sponges: "// &
         "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

    end select

    dS_dz = (S_sur - S_bot) / G%max_depth
    dT_dz = (T_sur - T_bot) / G%max_depth
    do j=js,je ; do i=is,ie
      xi0 = -depth_tot(i,j)
      do k = nz,1,-1
        xi0 = xi0 + 0.5 * dz(i,j,k)  ! Depth in middle of layer
        S(i,j,k) = S_sur + dS_dz * xi0
        T(i,j,k) = T_sur + dT_dz * xi0
        xi0 = xi0 + 0.5 * dz(i,j,k)  ! Depth at top of layer
      enddo
    enddo ; enddo

    ! for debugging
    !i=G%iec; j=G%jec
    !do k = 1,nz
    !  call calculate_density(T(i,j,k), S(i,j,k), 0.0, rho_tmp, tv%eqn_of_state, scale=US%kg_m3_to_R)
    !  write(mesg,*) 'Sponge - k,h,T,S,rho,Rlay',k,h(i,j,k),T(i,j,k),S(i,j,k),rho_tmp,GV%Rlay(k)
    !  call MOM_mesg(mesg,5)
    !enddo

    ! This call sets up the damping rates and interface heights in the sponges.
    call initialize_ALE_sponge(Idamp, G, GV, PF, ACSp, dz, nz, data_h_is_Z=.true.)

    !   Now register all of the fields which are damped in the sponge.   !
    ! By default, momentum is advected vertically within the sponge, but !
    ! momentum is typically not damped within the sponge.                !

    !  The remaining calls to set_up_sponge_field can be in any order. !
    if ( associated(tv%T) ) call set_up_ALE_sponge_field(T, G, GV, tv%T, ACSp, 'temp', &
        sp_long_name='temperature', sp_unit='degC s-1')
    if ( associated(tv%S) ) call set_up_ALE_sponge_field(S, G, GV, tv%S, ACSp, 'salt', &
        sp_long_name='salinity', sp_unit='g kg-1 s-1')


  else ! layer mode
    ! 1) Read eta, salt and temp from IC file
    call get_param(PF, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    ! GM: get two different files, one with temp and one with salt values
    ! this is work around to avoid having wrong values near the surface
    ! because of the FIT_SALINITY option. To get salt values right in the
    ! sponge, FIT_SALINITY=False. The oposite is true for temp. One can
    ! combined the *correct* temp and salt values in one file instead.
    call get_param(PF, mdl, "ISOMIP_SPONGE_FILE", state_file, &
              "The name of the file with temps., salts. and interfaces to "//&
              "damp toward.", fail_if_missing=.true.)
    call get_param(PF, mdl, "SPONGE_PTEMP_VAR", temp_var, &
              "The name of the potential temperature variable in "//&
              "SPONGE_STATE_FILE.", default="Temp")
    call get_param(PF, mdl, "SPONGE_SALT_VAR", salt_var, &
              "The name of the salinity variable in "//&
              "SPONGE_STATE_FILE.", default="Salt")
    call get_param(PF, mdl, "SPONGE_ETA_VAR", eta_var, &
              "The name of the interface height variable in "//&
              "SPONGE_STATE_FILE.", default="eta")

    !read temp and eta
    filename = trim(inputdir)//trim(state_file)
    if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
          "ISOMIP_initialize_sponges: Unable to open "//trim(filename))
    call MOM_read_data(filename, eta_var, eta(:,:,:), G%Domain, scale=US%m_to_Z)
    call MOM_read_data(filename, temp_var, T(:,:,:), G%Domain, scale=US%degC_to_C)
    call MOM_read_data(filename, salt_var, S(:,:,:), G%Domain, scale=US%ppt_to_S)

    ! for debugging
    !i=G%iec; j=G%jec
    !do k = 1,nz
    !  call calculate_density(T(i,j,k), S(i,j,k), 0.0, rho_tmp, tv%eqn_of_state, scale=US%kg_m3_to_R)
    !  write(mesg,*) 'Sponge - k,eta,T,S,rho,Rlay',k,eta(i,j,k),T(i,j,k),&
    !              S(i,j,k),rho_tmp,GV%Rlay(k)
    !  call MOM_mesg(mesg,5)
    !enddo

    ! Set the sponge damping rates so that the model will know where to
    ! apply the sponges, along with the interface heights.
    call initialize_sponge(Idamp, eta, G, PF, CSp, GV)
    ! Apply sponge in tracer fields
    call set_up_sponge_field(T, tv%T, G, GV, nz, CSp)
    call set_up_sponge_field(S, tv%S, G, GV, nz, CSp)

  endif

end subroutine ISOMIP_initialize_sponges

!> \namespace isomip_initialization
!!
!! See this paper for details: http://www.geosci-model-dev-discuss.net/8/9859/2015/gmdd-8-9859-2015.pdf
end module ISOMIP_initialization
