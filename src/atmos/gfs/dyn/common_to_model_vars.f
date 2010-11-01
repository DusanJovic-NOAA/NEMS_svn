      subroutine common_to_model_vars (psg,ttg,rqg,uug,vvg,
     &                                 global_lats_a,lonsperlat)
!!
!! hmhj - this routine change variables from common usage to model
!!        common usage are t=dry temperature (k), p is pascal, real winds
!!        model  usage are t=virtal temperature (k) or enthalpy, 
!!                         p is centibar, mapping winds
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_coordinate_def 
      use gfs_dyn_tracer_const
      use gfs_dyn_physcons, fv => con_fvirt, rerth => con_rerth,
     &              grav => con_g,  cp => con_cp , rd => con_rd
      implicit none
!!
!
      integer              global_lats_a(latg)
      integer                 lonsperlat(latg)
!
      REAL(KIND=KIND_EVOD) psg    (lonf,lats_node_a_max)
      REAL(KIND=KIND_EVOD) ttg    (lonf,lats_node_a_max,levs)
      REAL(KIND=KIND_EVOD) uug    (lonf,lats_node_a_max,levs)
      REAL(KIND=KIND_EVOD) vvg    (lonf,lats_node_a_max,levs)
      REAL(KIND=KIND_EVOD) rqg    (lonf,lats_node_a_max,levh)
!
      real(kind=kind_evod)   tfac(lonf,levs), sumq(lonf,levs), tkrt0
!
      integer              i,j,k,kk, nn, nnl
      integer              l,lan,lat
      integer              lons_lat
!
      real(kind=kind_evod), parameter :: qmin=1.0e-10, rkappa = cp / rd
      real(kind=kind_evod), parameter :: one=1.0, pa2cb=0.001
!
!--------------------------------------------------------------------
!
      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
!
! get factor for t
        if (thermodyn_id == 3) then
          do k=1,levs
            do i=1,lons_lat
              tfac(i,k) = 0.0
              sumq(i,k) = 0.0
            enddo
          enddo
          do nn=1,ntrac
            nnl = (nn-1)*levs
            if (cpi(nn) .ne. 0.0) then
              do k=1,levs
                do i=1,lons_lat
                  sumq(i,k) = sumq(i,k) + rqg(i,lan,nnl+k)
                  tfac(i,k) = tfac(i,k) + 
     &                       cpi(nn)*rqg(i,lan,nnl+k)
                enddo
              enddo
            endif
          enddo
          do k=1,levs
            do i=1,lons_lat
              tfac(i,k) = (one-sumq(i,k))*cpi(0) + tfac(i,k)
            enddo
          enddo
        else
          do k=1,levs
            do i=1,lons_lat
              tfac(i,k) = one + fv*max(rqg(i,lan,k),qmin) 
            enddo
          enddo
        endif
!
! get model t (virtual temperature or enthalpy)
        do k=1,levs
          do i=1,lons_lat
            ttg(i,lan,k) = ttg(i,lan,k) * tfac(i,k)
            uug(i,lan,k) = uug(i,lan,k) * coslat_a(lat)
            vvg(i,lan,k) = vvg(i,lan,k) * coslat_a(lat)
          enddo
        enddo
!
! get model ps (log surface pressure or surface pressure)
        if (gen_coord_hybrid) then   ! Ps is the prognostic variable
          do i=1,lons_lat
            psg (i,lan) =  psg(i,lan) * pa2cb 
          enddo
        else                         ! ln(Ps) is the prognostic variable
          do i=1,lons_lat
            psg(i,lan) = log( psg(i,lan) * pa2cb )
          enddo
        endif
!       call mymaxmin(psg(1,lan),lons_lat,lonf,1,' psg in com to mdl')
!
!       call mymaxmin(rqg(1,lan,1),lons_lat,lonf,1,' rqg in com to mdl')
      enddo
!
!      print *,' exit common_to_model_vars '
!!
      return
      end
