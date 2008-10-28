      SUBROUTINE LONLAT_PARA(global_lats_r,XLON,XLAT,lonsperlat)
!
c***********************************************************************
!
      USE MACHINE , ONLY : kind_grid

      use resol_def
      use layout1
      use gg_def
      use physcons, pi => con_pi
      implicit none
      integer i,j,lat
      integer                 lonsperlat(latg)
      real (kind=kind_grid) tpi,hpi,bphi
      PARAMETER (TPI=2.E0*PI,HPI=0.5E0*PI)
      integer              global_lats_r(latg)
      real (kind=kind_grid) XLON(lonf,lats_node_r)
      real (kind=kind_grid) XLAT(lonf,lats_node_r)
!
      xlon=0.
      xlat=0.
 
      DO j=1,lats_node_r
        lat = global_lats_r(ipt_lats_node_r-1+j)
        BPHI = TPI/lonsperlat(lat)
        if (lat.le.latg2) then
          DO i=1,lonsperlat(lat)
            XLON(I,J) = (i-1) * BPHI
            XLAT(I,J) = HPI - colrad_r(lat)
          ENDDO
        else
          DO i=1,lonsperlat(lat)
            XLON(I,J) =  (i-1) * BPHI
            XLAT(I,J) = colrad_r(lat)-HPI
          ENDDO
        endif
      ENDDO
 
      RETURN
      END
 