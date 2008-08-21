!
! !module: gfs_physics_run_mod --- run module of the grided
!                              component of the gfs physics.
!
! !description: gfs run module.
!
! !revision history:
!
!  november 2004      weiyu yang initial code.
!  may      2005      weiyu yang, updated to the new version of gfs.
!  janusry  2007      hann-ming henry juang for gfs dynamics only
!  july     2007      shrinivas moorthi for gfs physics only
!  november 2007      hann-ming henry juang continue for gfs physics
!
!
! !interface:
!
      module gfs_physics_run_mod
!
!!uses:
!
      use gfs_physics_internal_state_mod

      implicit none

      contains

      subroutine gfs_physics_run(gis_phy, rc)

      type(gfs_physics_internal_state), pointer, intent(inout) :: gis_phy
      integer, optional,                          intent(out)   :: rc

      real , save	:: timestep=0.0
      integer		rc1, k 
!***********************************************************************
!
!     lsfwd      logical true during a first forward step
!     lssav      logical true during a step for which
!                diagnostics are accumulated
!     lscca      logical true during a step for which convective clouds
!                are calculated from convective precipitation rates
!     phour      real forecast hour at the end of the previous timestep
!
!     lsout controls freq. of output
!     nsout controls freq. of output in time steps
!     fhout controls freq. of output in hours
!     nszer time steps between zeroing fluxes
!
!***********************************************************************

       print *,' enter gfs_physics_run '
!
       rc1 = 0
!
       call common_to_physics_vars(gis_phy%grid_gr(1,gis_phy%g_ps),      &
                                   gis_phy%grid_gr(1,gis_phy%g_t ),      &
                                   gis_phy%grid_gr(1,gis_phy%g_q ),      &
                                   gis_phy%grid_gr(1,gis_phy%g_u ),      &
                                   gis_phy%grid_gr(1,gis_phy%g_v ),      &
                                   gis_phy%grid_gr(1,gis_phy%g_p ),      &
                                   gis_phy%grid_gr(1,gis_phy%g_dp ),     &
                                   gis_phy%grid_gr(1,gis_phy%g_dpdt ),   &
                                   gis_phy%global_lats_r,                &
                                   gis_phy%lonsperlar)

       print *,' end of common_to_physics_vars '

!
! ---------------------------------------------------------------------
! ======================================================================
!                     do one physics time step
! ---------------------------------------------------------------------
        if (gis_phy%lsoutd) then
!         gis_phy%LSOUT = MOD(gis_phy%kdt ,NSOUT).EQ.0
          gis_phy%LSOUT = MOD(gis_phy%kdt ,NSOUT).EQ.0 .or. gis_phy%kdt == 1
        else
          gis_phy%LSOUT = .false.
        endif
!     print *,' kdt in run =',gis_phy%kdt,' lsout=',gis_phy%LSOUT &
!    ,' lsoutd=',gis_phy%lsoutd,' lssavd=',gis_phy%lssavd,' lssav=',lssav
!
! ======================================================================
        call do_physics_one_step(                                         &
                 gis_phy%deltim,  gis_phy%kdt,     gis_phy%phour,         &
                 gis_phy%grid_gr, gis_phy%sfc_fld, gis_phy%flx_fld,       &
                 gis_phy%lats_nodes_r,   gis_phy%global_lats_r,           &
                 gis_phy%lonsperlar,                                      &
                 gis_phy%XLON,    gis_phy%XLAT,    gis_phy%COSZDG,        &
                 gis_phy%HPRIME,  gis_phy%SWH,     gis_phy%HLW,           &
                 gis_phy%FLUXR,   gis_phy%SFALB,                          &
                 gis_phy%SLAG,    gis_phy%SDEC,    gis_phy%CDEC,          &
                 gis_phy%OZPLIN,  gis_phy%JINDX1,  gis_phy%JINDX2,        &
                 gis_phy%DDY,                                             &
                 gis_phy%phy_f3d, gis_phy%phy_f2d, gis_phy%NBLCK,         &
                 gis_phy%ZHOUR,   gis_phy%N3,      gis_phy%N4,            &
                 gis_phy%lssavd,                                          &
                 gis_phy%LSOUT,   gis_phy%COLAT1,  gis_phy%CFHOUR1)
                        

       print *,' end of do_physics_one_step '

! =======================================================================

! update hour
!
      gis_phy%phour      = fhour

! --------------------------------------------------------------------------
       call physics_to_common_vars(gis_phy%grid_gr(1,gis_phy%g_ps),      &
                                   gis_phy%grid_gr(1,gis_phy%g_t ),      &
                                   gis_phy%grid_gr(1,gis_phy%g_q ),      &
                                   gis_phy%grid_gr(1,gis_phy%g_u ),      &
                                   gis_phy%grid_gr(1,gis_phy%g_v ),      &
                                   gis_phy%grid_gr(1,gis_phy%g_p ),      &
                                   gis_phy%grid_gr(1,gis_phy%g_dp ),     &
                                   gis_phy%grid_gr(1,gis_phy%g_dpdt ),   &
                                   gis_phy%global_lats_r,                &
                                   gis_phy%lonsperlar)
! --------------------------------------------------------------------------
       print *,' end of physics_to_common_vars '
!
      if(present(rc)) then
          rc = rc1
      end if

      print *,' end of gfs_physics_run '

      end subroutine gfs_physics_run

      end module gfs_physics_run_mod
