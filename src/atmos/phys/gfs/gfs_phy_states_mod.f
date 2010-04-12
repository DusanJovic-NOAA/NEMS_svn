       MODULE GFS_Phy_States_Mod

!BOP

! !MODULE: GFS_Phy_States_Mod --- Define Physics Import/Export states

! 
!  june 2005 		weiyu yang             initial code.
!  february 2007 	hann-ming henry juang  
!			for gfs physics and gaussian grid DATA.
!  March 2009           Weiyu Yang, modified for the ensemble NEMS run.
!  2009/10/05           Sarah Lu, grid_gr unfolded to 3D
!
!!USEs:
!
      USE ESMF_Mod                 ! the esmf library.

! the derived TYPE of the internal state.
!----------------------------------------
      USE gfs_physics_internal_state_mod    ! Physics internal state
      USE gfs_physics_namelist_mod          ! Physics configuration

! routines which can be used to add a fortran array to 
! an esmf state and to get a fortran array from an 
! esmf state
!-----------------------------------------------------
      USE gfs_physics_grid_create_mod,    ONLY: mgrid
      USE gfs_physics_add_get_state_mod
      USE gfs_physics_err_msg_mod
      USE tracer_const, ONLY: cpi, ri

      TYPE(GFS_Phy_State_Namelist) :: cf

      IMPLICIT none

! !REVISION HISTORY:
!  Sarah Lu  2009-08-04  First version
!  Sarah Lu  2009-10-12  Port to the latest trunk
!  Sarah Lu  2009-10-16  Tracer bundle added; (shum, oz, cld) removed
!  Sarah Lu  2009-11-13  2D diag fields added to export state (for GOCART)
!  Sarah Lu  2009-12-08  3D diag fields added to export state (for GOCART)
!  Sarah Lu  2009-12-15  DQDT added to export state (for GOCART)
!  Sarah Lu  2010-02-09  Get ri, cpi attribute from import tracer bundle
!  Sarah Lu  2010-02-11  add set_phy_attribute to insert tracer_config to export state;
!                        add get_phy_attribute to retreive ri, cpi from import state
!  Sarah Lu  2010-04-09  restore fArr2D allowing the linkage to GOCART
!
      CONTAINS

      SUBROUTINE gfs_physics_import2internal(imp_gfs_phy, int_state, rc)

! this subroutine can be used to update the initial condition 
! fields of the internal state from the esmf inport state.
!------------------------------------------------------------

! every possible import state has its own turn-on/turn-off switch flag
! which can be used to fit different interface requirement from different
! outside grid component systems.
!------------------------------------------------------------------------

      TYPE(esmf_state)                                          :: imp_gfs_phy  
      TYPE(gfs_physics_internal_state), POINTER,  INTENT(inout) :: int_state 
      INTEGER, OPTIONAL,                          INTENT(out)   :: rc     

      TYPE(ESMF_Field)         :: Field
      TYPE(ESMF_FieldBundle)   :: Bundle

      INTEGER                  :: rc1, rcfinal
      INTEGER                  :: i, status

      REAL, DIMENSION(:, :),    POINTER :: FArr2D
      REAL, DIMENSION(:, :, :), POINTER :: FArr3D

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

      CALL esmf_logwrite(						&
           " update internal state with the esmf import state", 	&
            esmf_log_info, rc = rc1)

      cf = int_state%esmf_sta_list

! get the surface orography array from the esmf import state.
!------------------------------------------------------------
      IF(cf%z_import == 1) THEN
          CALL getf90arrayfromstate(imp_gfs_phy, 'hs', FArr2D, 0, rc = rc1)
          IF(int_state%grid_aldata) THEN
              int_state%grid_fld%z =  FArr2D
          ELSE
              int_state%grid_fld%z => FArr2D
          END IF
          CALL gfs_physics_err_msg(rc1, 'retrieve Farray from field -z', rcfinal)
      END IF

! get the surface pressure array from the esmf import state.
! for the detailed comments for every computational steps
! please refer to the surface orography array code.
!-----------------------------------------------------------
      IF(cf%ps_import == 1) THEN
          CALL getf90arrayfromstate(imp_gfs_phy, 'ps', FArr2D, 0, rc = rc1)
          IF(int_state%grid_aldata) THEN
              int_state%grid_fld%ps =  FArr2D
          ELSE
              int_state%grid_fld%ps => FArr2D
          END IF
          CALL gfs_physics_err_msg(rc1, "retrieve Farray from field -ps", rcfinal)
      END IF

! get the temperature array from the esmf import state.
!------------------------------------------------------
      IF(cf%temp_import == 1) THEN
          CALL getf90arrayfromstate(imp_gfs_phy, 't', FArr3D, 0, rc = rc1)
          IF(int_state%grid_aldata) THEN
              int_state%grid_fld%t =  FArr3D
          ELSE
              int_state%grid_fld%t => FArr3D
          END IF
          CALL gfs_physics_err_msg(rc1, "retrieve Farray from field -t", rcfinal)
      END IF

! get the zonal-wind array from the esmf import state.
! for detailed line by line comments please refer to 
! the temperature code.
!-----------------------------------------------------
      IF(cf%u_import == 1) THEN
          CALL getf90arrayfromstate(imp_gfs_phy, 'u', FArr3D, 0, rc = rc1)
          IF(int_state%grid_aldata) THEN
              int_state%grid_fld%u =  FArr3D
          ELSE
              int_state%grid_fld%u => FArr3D
          END IF
          CALL gfs_physics_err_msg(rc1, "retrieve Farray from field -u", rcfinal)
      END IF

! get the meridian-wind array from the esmf import state.
!-----------------------------------------------------
      IF(cf%v_import == 1) THEN
          CALL getf90arrayfromstate(imp_gfs_phy, 'v', FArr3D, 0, rc = rc1)
          IF(int_state%grid_aldata) THEN
              int_state%grid_fld%v =  FArr3D
          ELSE
              int_state%grid_fld%v => FArr3D
          END IF
          CALL gfs_physics_err_msg(rc1, "retrieve Farray from field -v", rcfinal)
      END IF

! get the tracer array from the esmf import state.
!-------------------------------------------------
      IF(cf%tracer_import == 1) THEN
          CALL ESMF_StateGet(imp_gfs_phy, 'tracers', Bundle, rc = rc1 )
          CALL gfs_physics_err_msg(rc1, 'retrieve Ebundle from state', rcfinal)
          DO i = 1, int_state%ntrac
              IF(ASSOCIATED(FArr3D)) NULLIFY(FArr3D)
              CALL ESMF_FieldBundleGet(Bundle,                            &
                                       trim(int_state%gfs_phy_tracer%vname(i)), &
                                       Field,                             &
                                       rc = rc1)
              CALL gfs_physics_err_msg(rc1, 'retrieve Efield from bundle', rcfinal)
              CALL ESMF_FieldGet(Field, FArray = FArr3D, localDE = 0, rc = rc1)
              CALL gfs_physics_err_msg(rc1, 'retrieve Farray from field', rcfinal)
              IF(int_state%grid_aldata) THEN
                  int_state%grid_fld%tracers(i)%flds =  FArr3D
              ELSE
                  int_state%grid_fld%tracers(i)%flds => FArr3D
              END IF
          END DO

!  --- Retrieve ri/cpi from import state; fill in tracer_const local arrays
          if ( int_state%start_step ) then
            call get_phy_attribute
!  ---      inputs:  (in scope variables)
!  ---      outputs: (in scope variables)
          endif

      END IF

! get the pressure array from the esmf import state.
!-------------------------------------------------------------
      IF(cf%p_import == 1) THEN
          CALL getf90arrayfromstate(imp_gfs_phy, 'p', FArr3D, 0, rc = rc1)
          IF(int_state%grid_aldata) THEN
              int_state%grid_fld%p =  FArr3D
          ELSE
              int_state%grid_fld%p => FArr3D
          END IF
          CALL gfs_physics_err_msg(rc1, "retrieve Farray from field -p", rcfinal)
      END IF

!

! get the pressure layer depth (dp) array from the esmf import state.
!-------------------------------------------------------------
      IF(cf%dp_import == 1) THEN
          CALL getf90arrayfromstate(imp_gfs_phy, 'dp', FArr3D, 0, rc = rc1)
          IF(int_state%grid_aldata) THEN
              int_state%grid_fld%dp =  FArr3D
          ELSE
              int_state%grid_fld%dp => FArr3D
          END IF
          CALL gfs_physics_err_msg(rc1, "retrieve Farray from field -dp", rcfinal)
      END IF

!
! get the omega (dpdt) array from the esmf import state.
!-------------------------------------------------------------
      IF(cf%dpdt_import == 1) THEN
          CALL getf90arrayfromstate(imp_gfs_phy, 'dpdt', FArr3D, 0, rc = rc1)
          IF(int_state%grid_aldata) THEN
              int_state%grid_fld%dpdt =  FArr3D
          ELSE
              int_state%grid_fld%dpdt => FArr3D
          END IF
          CALL gfs_physics_err_msg(rc1, "retrieve Farray from field -dpdt", rcfinal)
      END IF

!
!
! PRINT out the final error signal message and put it to rc.
!-----------------------------------------------------------
      IF(PRESENT(rc)) CALL gfs_physics_err_msg_final(rcfinal, &
          "gfs_physics_import2internal", rc)

     contains
! =================

!-----------------------------
      subroutine get_phy_attribute

!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

! ---  Retrieve ri and cpi
       CALL ESMF_AttributeGet(Bundle                         &  !<-- Tracer bundle
               ,name ='cpi_dryair'                           &  !<-- Name of the attribute to retrieve
               ,value = int_state%gfs_phy_tracer%cpi(0)      &  !<-- Value of the attribute
               ,rc   =RC)
       call gfs_physics_err_msg(rc,'retrieve cpi(0) attribute from phy_imp',rcfinal)

       CALL ESMF_AttributeGet(Bundle                         &  !<-- Tracer bundle
               ,name ='ri_dryair'                            &  !<-- Name of the attribute to retrieve
               ,value = int_state%gfs_phy_tracer%ri(0)       &  !<-- Value of the attribute
               ,rc   =RC)
       call gfs_physics_err_msg(rc,'retrieve ri(0) attribute from phy_imp',rcfinal)

       CALL ESMF_AttributeGet(Bundle                         &  !<-- Tracer bundle
               ,name ='cpi'                                  &  !<-- Name of the attribute to retrieve
               ,count= int_state%ntrac                       &  !<-- Number of values in the attribute
               ,valueList = int_state%gfs_phy_tracer%cpi(1:int_state%ntrac)  &!<-- Value of the attribute
               ,rc   =RC)
       call gfs_physics_err_msg(rc,'retrieve cpi(:) attribute from phy_imp',rcfinal)

       CALL ESMF_AttributeGet(Bundle                         &  !<-- Tracer bundle
               ,name ='ri'                                   &  !<-- Name of the attribute to retrieve
               ,count= int_state%ntrac                       &  !<-- Number of values in the attribute
               ,valueList = int_state%gfs_phy_tracer%ri(1:int_state%ntrac)  &!<-- Value of the attribute
               ,rc   =RC)
       call gfs_physics_err_msg(rc,'retrieve ri(:) attribute from phy_imp',rcfinal)

! ---  Fill in ri/cpi local array
       allocate(ri(0:int_state%ntrac), stat=status)
         if( status .ne. 0 ) print *, 'ERROR: Fail to allocate ri'
       allocate(cpi(0:int_state%ntrac), stat=status)
         if( status .ne. 0 ) print *, 'ERROR: Fail to allocate cpi'
       cpi(0:int_state%ntrac) =                                       &
              int_state%gfs_phy_tracer%cpi(0:int_state%ntrac)
       ri(0:int_state%ntrac) =                                        &
              int_state%gfs_phy_tracer%ri(0:int_state%ntrac)

      return
      end subroutine get_phy_attribute

      END SUBROUTINE gfs_physics_import2internal


! =========================================================================


      SUBROUTINE gfs_physics_internal2export(int_state, exp_gfs_phy, rc)

! 
! this subroutine will be changed that all export
! esmf states will get data directly from the gfs internal structure data arrays
! 

!
!!uses:
!
!
! !input/output variables and parameters:
!----------------------------------------
      TYPE(esmf_state),                           INTENT(inout) :: exp_gfs_phy 
      TYPE(gfs_physics_internal_state), POINTER,  INTENT(inout) :: int_state 
      INTEGER, OPTIONAL,                          INTENT(out)   :: rc     

! local array size parameter of the esmf export state arrays.
!------------------------------------------------------------
      INTEGER               :: rc1     ! error signal work variable.
      INTEGER               :: rcfinal ! the final error signal variable.
      INTEGER               :: k

      TYPE(ESMF_Field)                              :: Field
      TYPE(ESMF_FieldBundle), SAVE                  :: Bundle

      REAL(KIND = kind_evod), DIMENSION(:, :),    POINTER ::    &
                   hold_z,    hold_ps

      REAL(KIND = kind_evod), DIMENSION(:, :, :), POINTER ::    &
                   hold_temp, hold_u,      hold_v,              &
                   hold_p,      hold_dp,   hold_dpdt

      REAL, DIMENSION(:, :, :), POINTER :: FArr3D
      REAL, DIMENSION(:, :),    POINTER :: FArr2D
      integer                  :: i, j

!
!*  2-D/3-D quantities for GOCART gridded component
!*  nfld_2d and nfld_3d define the number of 2D and 3D fields
!*  the field names are specified in vname_2d and vname_3d
!
! TROPP: tropopause_pressure_based_on_blended_estimate --
! LWI  : land-ocean-ice_mask             ==> slmsk, fice
! ZPBL : Planetary boundary layer height ==> hpbl
! FRLAKE: fraction_of_lake
! FRACI: ice_covered_fraction_of_tile
! WET1 : surface_soil_wetness            ==> smc(1), stype
! LAI  : leaf_area_index                 ==> vtype
! GRN  : greeness_fraction               ==> vfrac
! CN_PRCP: Surface Conv. rain flux needed by land ==> rainc
! NCN_PRCP: Non-convective precipitation ==> rain
! PS   : surface_pressure                ==>  ps
! SH   : sensible_heat_flux_from_turbulence ==> dtsfci
! TA   : surface_temperature_from_surface ==> tsea
! TSOIL1: soil_temperatures_layer_1       ==> stc(1)
! U10M : 10-meter_eastward_wind           ==> u10m
! V10M : 10-meter_northward_wind          ==> v10m
! USTAR: surface_velocity_scale           ==> uustar
! Z0H: surface_roughness_for_heat         ==> zorl
!
      integer, parameter      :: nfld_2d = 16
      integer, parameter      :: nfld_3d = 2
      integer(ESMF_KIND_I4), allocatable, save :: lonsperlar_r(:)
      integer, save           :: lonr, lats_node_r, lats_node_r_max
      integer                 :: ilat
      character*50            :: msg
      character*8             :: vname, vname_2d(nfld_2d), &
                                 vname_3d(nfld_3d)

      data vname_2d /'slmsk', 'fice', 'hpbl', 'smc1',     &
                     'stype', 'vtype', 'vfrac', 'rainc', &
                     'rain', 'dtsfci', 'tsea', 'stc1',   &
                     'u10m', 'v10m',  'ustar','zorl'/

      data vname_3d /'fcld','dqdt'/

      LOGICAL, SAVE :: first
      DATA first/.true./
      SAVE hold_z,    hold_ps,     hold_temp, hold_u, hold_v, &
           hold_p,    hold_dp,     hold_dpdt

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

      IF(.NOT. first) RETURN

      cf = int_state%esmf_sta_list

! PRINT out the information.
!-----------------------------------------------
      PRINT*, 'do int_state to exp_gfs_phy'

      CALL esmf_logwrite("begining to put the esmf export state.", &
                    esmf_log_info, rc = rc1)

! orography field. gaussian grid
!----------------------------------------------------------------

      IF(cf%z_export == 1) THEN
          hold_z => int_state%grid_fld%z
          CALL addf90arraytostate(exp_gfs_phy, mgrid, 'hs',             & 
                                  hold_z, rc = rc1)
          CALL gfs_physics_err_msg(rc1, "done z_export.", rcfinal)
      END IF

! surface pressure
!----------------------------------------------------------------

      IF(cf%ps_export == 1) THEN
          hold_ps => int_state%grid_fld%ps
          CALL addf90arraytostate(exp_gfs_phy, mgrid, 'ps',             & 
                                  hold_ps, rc = rc1)
          CALL gfs_physics_err_msg(rc1, "put to esmf state - ps_ex", rcfinal)
      END IF

! add the temperature fileds into the esmf export state.
!-------------------------------------------------------

      IF(cf%temp_export == 1) THEN
          hold_temp => int_state%grid_fld%t
          CALL addf90arraytostate(exp_gfs_phy, mgrid, 't',             & 
                                  hold_temp, rc = rc1)
          CALL gfs_physics_err_msg(rc1, "put to esmf state - t_ex", rcfinal)
      END IF

! to add the u field into the esmf export state.  for the detailed
! description comments please refer to the temperature field.
!--------------------------------------------------------------------------

      IF(cf%u_export == 1) THEN
          hold_u => int_state%grid_fld%u
          CALL addf90arraytostate(exp_gfs_phy, mgrid, 'u',             & 
                                  hold_u, rc = rc1)
          CALL gfs_physics_err_msg(rc1, "put to esmf state - u_ex", rcfinal)
      END IF

! add the v field into the esmf export state.
!----------------------------------------------------

      IF(cf%v_export == 1) THEN
          hold_v => int_state%grid_fld%v
          CALL addf90arraytostate(exp_gfs_phy, mgrid, 'v',             &
                                  hold_v, rc = rc1)
          CALL gfs_physics_err_msg(rc1, "put to esmf state - v_ex", rcfinal)
      END IF

! add the tracer fields into the esmf export state.
!---------------------------------------------------
      IF(cf%tracer_export == 1) THEN
          Bundle = ESMF_FieldBundleCreate(name = 'tracers', grid = mgrid, rc = rc1)
          CALL gfs_physics_err_msg(rc1, "create empty fieldbundle", rcfinal)
          DO k = 1, int_state%ntrac
              NULLIFY(fArr3D)
              fArr3D => int_state%grid_fld%tracers(k)%flds
              Field  = ESMF_FieldCreate(mgrid, fArr3D,                    &
                  name = trim(int_state%gfs_phy_tracer%vname(k)), rc = rc1)
              CALL ESMF_FieldBundleAdd(Bundle, Field, rc = rc1)
              CALL gfs_physics_err_msg(rc1, "add field to bundle, vname=1", rcfinal)
          END DO

!  --- Insert tracer_config to tracer bundle (to be exported to other component)
       call set_phy_attribute
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

          CALL ESMF_StateAdd(exp_gfs_phy, Bundle, rc = rc1)
          CALL gfs_physics_err_msg(rc1, "add to esmf state - tracer", rcfinal)
      END IF
! add layer pressure (pp) fileds into the esmf export state.
!-------------------------------------------------------

      IF(cf%p_export == 1) THEN
          hold_p => int_state%grid_fld%p
          CALL addf90arraytostate(exp_gfs_phy, mgrid, 'p',       & 
                                  hold_p, rc = rc1)
          CALL gfs_physics_err_msg(rc1,                          &
                  "put to esmf state - p_ex",rcfinal)
      END IF

! add pressure depth (dp) fileds into the esmf export state.
!-------------------------------------------------------

      IF(cf%dp_export == 1) THEN
          hold_dp => int_state%grid_fld%dp
          CALL addf90arraytostate(exp_gfs_phy, mgrid, 'dp',         & 
                                  hold_dp, rc = rc1)
          CALL gfs_physics_err_msg(rc1,                             &
                  "put to esmf state - dp_ex",rcfinal)
      END IF


! add omega (dpdt) fileds into the esmf export state.
!-------------------------------------------------------

      IF(cf%dpdt_export == 1) THEN
          hold_dpdt => int_state%grid_fld%dpdt
          CALL addf90arraytostate(exp_gfs_phy, mgrid, 'dpdt',        & 
                                  hold_dpdt, rc = rc1)
          CALL gfs_physics_err_msg(rc1,                             &
                  "put to esmf state - dpdt_ex",rcfinal)
      END IF

!!****************************************************************************
!!****************************************************************************
!! add 2D diag fields (for GOCART) to esmf export state
!! optioal -- check lgocart flag
!!------------------------------------------------------------

      lab_if_gocart :  if (int_state%lgocart ) then         !!! <======= optional

      lonr            = int_state%lonr
      lats_node_r     = int_state%lats_node_r
      lats_node_r_max = int_state%lats_node_r_max

      if ( .not. allocated (lonsperlar_r)) then
        allocate ( lonsperlar_r(lats_node_r_max))
        lonsperlar_r(:) = 0
        do i = 1, lats_node_r
          ilat = int_state%global_lats_r(int_state%ipt_lats_node_r-1+i)
          lonsperlar_r(i) =  int_state%lonsperlar(ilat)
        enddo
      endif

      CALL ESMF_AttributeSet(state=exp_gfs_phy             &  !<-- The physics export state
                            ,name ='lats_node_r'           &  !<-- Name of the attribute to insert
                            ,value= lats_node_r            &  !<-- Value of the attribute
                            ,rc   =RC)

      CALL ESMF_AttributeSet(state=exp_gfs_phy             &
                            ,name ='lats_node_r_max'       &
                            ,value= lats_node_r_max        &
                            ,rc   =RC)

      CALL ESMF_AttributeSet(state=exp_gfs_phy             &
                            ,name ='lonr'                  &
                            ,value= lonr                   &
                            ,rc   =RC)

      CALL ESMF_AttributeSet(state=exp_gfs_phy             &  !<-- The physics export state
                            ,name ='lonsperlar_r'          &  !<-- Name of the attribute to insert
                            ,count= lats_node_r_max        &  !<-- Number of values in the attribute
                            ,valueList =lonsperlar_r       &  !<-- Value of the attribute
                            ,rc   =RC)

! loop through the 2D diag fields

      lab_do_2d : DO i = 1, nfld_2d

        vname  = trim(vname_2d(i))
        if(associated(fArr2D)) nullify(fArr2D)

        SELECT CASE (vname)

!!        LWI: land-ocean-ice_mask
          CASE ('slmsk')   ! Land-sea mask (1=land; 0=sea)
            fArr2D => int_state%sfc_fld%slmsk

          CASE ('fice')    ! Ice concentration (ice>0; no ice=0)
            fArr2D => int_state%sfc_fld%fice

!!        Planetary boundary layer height (m)
          CASE ('hpbl')    ! Boundary layer height (m)
            fArr2D => int_state%flx_fld%hpbl

!!        surface_soil_wetness
          CASE ('smc1')
            fArr2D => int_state%sfc_fld%smc(1,:,:)
          CASE ('stype')
            fArr2D => int_state%sfc_fld%stype

!!        LAI: leaf_area_index
          CASE ('vtype')
            fArr2D => int_state%sfc_fld%vtype

!!        GRN: greeness_fraction
          CASE ('vfrac')
            fArr2D => int_state%sfc_fld%vfrac

!!        CN_PRCP: Surface Conv. rain flux (kg/m^2/s)
          CASE ('rainc')
            fArr2D => int_state%flx_fld%rainc

!!        NCN_PRCP: Non-convective precipitation rate (kg/m^2/s)
          CASE ('rain')
            fArr2D => int_state%flx_fld%rain

!!        SHFX: sensible_heat_flux_from_turbulence (W m-2)
          CASE ('dtsfci')
            fArr2D => int_state%flx_fld%dtsfci

!!        TA: Surface Air Temperature (K)
          CASE ('tsea')
            fArr2D => int_state%sfc_fld%tsea

!!        TSOIL1: soil_temperatures_layer_1 (k)
          CASE ('stc1')
            fArr2D => int_state%sfc_fld%stc(1,:,:)

!!        U10M: 10-meter_eastward_wind (m s-1)
          CASE ('u10m')
            fArr2D => int_state%flx_fld%u10m

!!        V10M: 10-meter_northward_wind (m s-1)
          CASE ('v10m')
            fArr2D => int_state%flx_fld%v10m

!!        USTAR: surface_velocity_scale (m s-1)
          CASE ('ustar')
            fArr2D => int_state%sfc_fld%uustar

!!        Z0H: surface_roughness_for_heat (m)
          CASE ('zorl')
            fArr2D => int_state%sfc_fld%zorl

        END SELECT

        msg   = "Create ESMF Field from "//vname
        field = ESMF_FieldCreate(name=vname, grid=mgrid,       &
                fArray=fArr2D, gridToFieldMap=(/1,2,0/),       &
                indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
        call gfs_physics_err_msg(rc,msg,rcfinal)

        msg   = "Add to Physics Export State"
        call ESMF_StateAdd(exp_gfs_phy, field, rc=rc)
        call gfs_physics_err_msg(rc, msg ,rcfinal)

      ENDDO   lab_do_2D

! loop through the 3D diag fields

      lab_do_3d : DO i = 1, nfld_3d

        vname  = trim(vname_3d(i))
        if(associated(fArr3D)) nullify(fArr3D)

        SELECT CASE (vname)

          CASE ('fcld')   ! cloud cover
            fArr3D => int_state%g3d_fld%fcld

          CASE ('dqdt')   ! total moisture tendency
            fArr3D => int_state%g3d_fld%dqdt

        END SELECT

        msg   = "Create ESMF Field from "//vname
        field = ESMF_FieldCreate(name=vname, grid=mgrid,       &
                fArray=fArr3D,indexFlag=ESMF_INDEX_DELOCAL,rc=rc)
        call gfs_physics_err_msg(rc,msg,rcfinal)

        msg   = "Add to Physics Export State"
        call ESMF_StateAdd(exp_gfs_phy, field, rc=rc)
        call gfs_physics_err_msg(rc, msg ,rcfinal)

      ENDDO   lab_do_3D

      endif lab_if_gocart

!
!! print out the final error signal information and put it to the rc.
!!-------------------------------------------------------------------
      IF(PRESENT(rc)) CALL gfs_physics_err_msg_final(rcfinal, &
          "gfs_physics_internal2export", rc)

      first = .false.

      contains
! =================

!-----------------------------
      subroutine set_phy_attribute
!.............................
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

!  ---  locals:
       TYPE(ESMF_Logical)      :: doing_SU, doing_SS, &
                                  doing_DU, doing_OC, doing_BC

! ---  Set attributes for ri and cpi
       msg   = "insert cpi(0) attribute to phy_exp"
       CALL ESMF_AttributeSet(Bundle                         &  !<-- Phy export state tracer bundle
               ,name ='cpi_dryair'                           &  !<-- Name of the attribute to insert
               ,value = int_state%gfs_phy_tracer%cpi(0)      &  !<-- Value of the attribute
               ,rc   =RC)
       call gfs_physics_err_msg(rc,msg,rcfinal)

       msg   = "inset ri(0) attribute to phy_exp"
       CALL ESMF_AttributeSet(Bundle                         &  !<-- Phy export state tracer bundle
               ,name ='ri_dryair'                            &  !<-- Name of the attribute to insert
               ,value = int_state%gfs_phy_tracer%ri(0)       &  !<-- Value of the attribute
               ,rc   =RC)
       call gfs_physics_err_msg(rc,msg,rcfinal)

       msg   = "insert cpi(:) attribute to phy_exp"
       CALL ESMF_AttributeSet(Bundle                         &  !<-- Phy export state tracer bundle
               ,name ='cpi'                                  &  !<-- Name of the attribute array
               ,count= int_state%ntrac                       &  !<-- Length of array being inserted
               ,valueList = int_state%gfs_phy_tracer%cpi(1:int_state%ntrac)  &!<-- The array being inserted
               ,rc   =RC)
       call gfs_physics_err_msg(rc,msg,rcfinal)

       msg   = "insert ri(:) attribute to phy_exp"
       CALL ESMF_AttributeSet(Bundle                         &  !<-- Phy export state tracer bundle
               ,name ='ri'                                   &  !<-- Name of the attribute array
               ,count= int_state%ntrac                       &  !<-- Length of array being inserted
               ,valueList = int_state%gfs_phy_tracer%ri(1:int_state%ntrac)  &!<-- The array being inserted
               ,rc   =RC)
       call gfs_physics_err_msg(rc,msg,rcfinal)

! ---  Set attribute for ntrac
       msg   = "insert ntrac attribute to phy_exp"
       CALL ESMF_AttributeSet(Bundle                         &  !<-- Phy export state tracer bundle
               ,name ='ntrac'                                &  !<-- Name of the attribute to insert
               ,value = int_state%gfs_phy_tracer%ntrac        &  !<-- Value of the attribute
               ,rc   =RC)
       call gfs_physics_err_msg(rc,msg,rcfinal)

! ---  Set attributes for doing_[DU,SS,SU,OC,BC]
       doing_DU = ESMF_False
       if ( int_state%gfs_phy_tracer%doing_DU ) doing_DU = ESMF_True
       msg   = "insert logical doing_DU attribute to phy_exp"
       CALL ESMF_AttributeSet(Bundle                          &  !<-- Phy export state tracer bundle
               ,name  ='doing_DU'                             &  !<-- Name of the logical
               ,value = doing_DU                              &  !<-- The logical being inserted into the Bundle
               ,rc   =RC)
       call gfs_physics_err_msg(rc,msg,rcfinal)

       doing_SU = ESMF_False
       if ( int_state%gfs_phy_tracer%doing_SU ) doing_SU = ESMF_True
       msg   = "insert logical doing_SU attribute to phy_exp"
       CALL ESMF_AttributeSet(Bundle                          &  !<-- Phy export state tracer bundle
               ,name  ='doing_SU'                             &  !<-- Name of the logical
               ,value = doing_SU                              &  !<-- The logical being inserted into the Bundle
               ,rc   =RC)
       call gfs_physics_err_msg(rc,msg,rcfinal)

       doing_SS = ESMF_False
       if ( int_state%gfs_phy_tracer%doing_SS ) doing_SS = ESMF_True
       msg   = "insert logical doing_SS attribute to phy_exp"
       CALL ESMF_AttributeSet(Bundle                          &  !<-- Phy export state tracer bundle
               ,name  ='doing_SS'                             &  !<-- Name of the logical
               ,value = doing_SS                              &  !<-- The logical being inserted into the Bundle
               ,rc   =RC)
       call gfs_physics_err_msg(rc,msg,rcfinal)

       doing_OC = ESMF_False
       if ( int_state%gfs_phy_tracer%doing_OC ) doing_OC = ESMF_True
       msg   = "insert logical doing_OC attribute to phy_exp"
       CALL ESMF_AttributeSet(Bundle                          &  !<-- Phy export state tracer bundle
               ,name  ='doing_OC'                             &  !<-- Name of the logical
               ,value = doing_OC                              &  !<-- The logical being inserted into the Bundle
               ,rc   =RC)
       call gfs_physics_err_msg(rc,msg,rcfinal)

       doing_BC = ESMF_False
       if ( int_state%gfs_phy_tracer%doing_BC ) doing_BC = ESMF_True
       msg   = "insert logical doing_BC attribute to phy_exp"
       CALL ESMF_AttributeSet(Bundle                          &  !<-- Phy export state tracer bundle
               ,name  ='doing_BC'                             &  !<-- Name of the logical
               ,value = doing_BC                              &  !<-- The logical being inserted into the Bundle
               ,rc   =RC)
       call gfs_physics_err_msg(rc,msg,rcfinal)

      RETURN
      end subroutine set_phy_attribute

      END SUBROUTINE gfs_physics_internal2export

! ==========================================================================

      END MODULE GFS_Phy_States_Mod
