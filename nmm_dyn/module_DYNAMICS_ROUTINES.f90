!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        module module_dynamics_routines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
use module_control,only:klog,kint,kfpt,kdbl &
                       ,adv_upstream,adv_standard &
                       ,e_bdy,n_bdy,s_bdy,w_bdy &
                       ,read_global_sums,write_global_sums &
                       ,timef
!
use module_clocktimes
!
use module_dm_parallel,only : ids,ide,jds,jde &
                             ,ims,ime,jms,jme &
                             ,its,ite,jts,jte &
                             ,its_b1,ite_b1,ite_b2 &
                             ,its_b1_h1,ite_b1_h1,ite_b1_h2 &
                             ,its_b1_h2 &
                             ,its_h1,ite_h1,its_h2,ite_h2 &
                             ,jts_b1,jte_b1,jte_b2 &
                             ,jts_b1_h1,jte_b1_h1,jte_b1_h2 &
                             ,jts_b1_h2 &
                             ,jts_h1,jte_h1,jts_h2,jte_h2 &
                             ,ihalo,jhalo &
                             ,mpi_comm_comp,mype_share 
!
use module_exchange
use module_fltbnds
use module_constants
!
integer(kind=kint),save :: &
 iunit_advec_sums 
!
integer(kind=kint),private :: &
 mype
!
real(kind=kfpt),private :: &
 btim
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine pgforce &
(ntsd,first &
,lm &
,dt,rdyv &
,dsg2,pdsg1 &
,rdxv,wpdar &
,fis,pd,pdo &
,t,q,cw,dwdt &
,pint &
,rtop &
!---temporary arguments-------------------------------------------------
,div &
,pcne,pcnw &
,pcx,pcy &
,tcu,tcv)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),parameter:: &
 jw=01                       ! rows w/o correction next to poles

real(kind=kfpt),parameter:: &
 cfc=1.533 &                 ! adams-bashforth positioning in time
,bfc=1.-cfc                  ! adams bashforth positioning in time
!-----------------------------------------------------------------------
logical(kind=klog),intent(in):: &
 first                       ! first pass

integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,ntsd                        ! timestep

real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,rdyv                        ! 1/deltay

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures

real(kind=kfpt),dimension(jds:jde),intent(in):: &
 rdxv &                      ! 1/deltax
,wpdar                       ! divergence correction weight

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 fis &                       ! surface geopotential
,pd &                        ! sigma range pressure difference
,pdo                         ! old sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 t &                         ! temperature
,q &                         ! specific humidity
,cw &                        ! condensate
,dwdt                        ! nonhydrostatic correction factor

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(in):: &
 pint                        ! pressure at interfaces

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out):: &
 rtop                        ! RT/p

!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out) :: &
 div &                       ! horizontal mass divergence
,pcne &                      ! second term of pgf, ne direction
,pcnw &                      ! second term of pgf, nw direction
,pcx &                       ! second term of pgf, x direction
,pcy                         ! second term of pgf, y direction

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout) :: &
 tcu &                       ! time change of u
,tcv                         ! time change of v
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,jcl &                       ! lower bound for no divergence correction
,jch &                       ! upper bound for no divergence correction
,l                           ! index in p direction
                                                                                                                                              
real(kind=kfpt):: &
 apd &                       ! hydrostatic pressure difference at the point
,apelp &                     ! pressure at the point
,dfip &                      ! delta phi
,dfdp &                      ! dfi/dp
,fiup &                      ! geopotential at the upper interface
,pdnep &                     ! hydrostatic pressure difference at the point
,pdnwp &                     ! hydrostatic pressure difference at the point
,pdxp &                      ! hydrostatic pressure difference at the point
,pdyp &                      ! hydrostatic pressure difference at the point
,ppne &                      ! first term of pgf, ne direction
,ppnw &                      ! first term of pgf, nw direction
,ppx &                       ! first term of pgf, x direction
,ppy &                       ! first term of pgf, y direction
,rdu &                       !
,rdv &                       !
,rpdp &                      !
,wprp                        ! divergence modification weight at the point
                                                                                                                                              
real(kind=kfpt),dimension(its_h1:ite_h1,jts_h1:jte_h1):: &
 apel &                      ! scratch, pressure in the middle of the layer
,dfi &                       ! scratch, delta phi
,filo &                      ! scratch, geopotential at lower interface
,fim                         ! scratch, geopotential in the middle of the layer
real(kind=kfpt),dimension(its_h1:ite_h1,jts_h1:jte_h1,1:lm):: &
 apel_3d &                      ! scratch, 3d copy of pressure in the middle of the layer
,fim_3d  &                      ! scratch, 3d copy of geopotential in the middle of the layer
,dfi_3d                         ! scratch, 3d copy of delta phi

 
real(kind=kfpt),dimension(its_b1:ite_h1,jts_b1:jte_h1):: &
 pdne &                      ! hydrostatic pressure difference at the point
,pdnw                        ! hydrostatic pressure difference at the point
real(kind=kfpt),dimension(its_b1:ite_h1,jts:jte_h1):: &
 pdx                         ! hydrostatic pressure difference at the point
 
real(kind=kfpt),dimension(its:ite_h1,jts_b1:jte_h1):: &
 pdy                         ! hydrostatic pressure difference at the point
 
real(kind=kfpt),dimension(its_b1:ite_h1,jts_b1:jte_h1):: &
 pgne &                      ! scratch, pgf, ne direction
,pgnw                        ! scratch, pgf, nw direction
 
real(kind=kfpt),dimension(its_b1:ite_h1,jts:jte_h1):: &
 pgx                         ! scratch, pgf, x direction
 
real(kind=kfpt),dimension(its:ite_h1,jts_b1:jte_h1):: &
 pgy                         ! scratch, pgf, y direction
!-----------------
!-----------------
integer(kind=kint) :: &
 jstart &
,jstop &
,nth &
,omp_get_num_threads &
,omp_get_thread_num &
,tid 
!-----------------
!-----------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!do l=1,lm
!do j=1,jm
!  print*,j,l,u(ide-2,j,l),u(ide-1,j,l),u(ide,j,l),u(1,j,l),u(2,j,l)
!  print*,j,l,v(ide-2,j,l),v(ide-1,j,l),v(ide,j,l),v(1,j,l),v(2,j,l)
!enddo
!enddo
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      do j=jts_h1,jte_h1
        do i=its_h1,ite_h1
          filo(i,j)=fis(i,j)
        enddo
      enddo
!-----------------------------------------------------------------------
      do j=jts,jte_h1
        do i=its_b1,ite_h1
          pdx (i,j)=((pd (i-1,j)+pd (i,j))*cfc &
                    +(pdo(i-1,j)+pdo(i,j))*bfc)*0.5
        enddo
      enddo
!
      do j=jts_b1,jte_h1
        do i=its,ite_h1
          pdy (i,j)=((pd (i,j-1)+pd (i,j))*cfc &
                    +(pdo(i,j-1)+pdo(i,j))*bfc)*0.5
        enddo
      enddo
!
      do j=jts_b1,jte_h1
        do i=its_b1,ite_h1
          pdne(i,j)=((pd (i-1,j-1)+pd (i,j))*cfc &
                    +(pdo(i-1,j-1)+pdo(i,j))*bfc)*0.5
          pdnw(i,j)=((pd (i,j-1)+pd (i-1,j))*cfc &
                    +(pdo(i,j-1)+pdo(i-1,j))*bfc)*0.5
        enddo
      enddo
!-----------------------------------------------------------------------
!---vertical grand loop-------------------------------------------------
!-----------------------------------------------------------------------
!-----------------
!-----------------
!.......................................................................
!$omp parallel private(nth,tid,i,j,l,jstart,jstop,apelp,fiup,dfdp,dfip)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth, jts_h1, jte_h1, jstart, jstop)
!-----------------
!-----------------
!-----------------------------------------------------------------------
      do l=lm,1,-1
!-----------------------------------------------------------------------
        do j=jstart,jstop
          do i=its_h1,ite_h1
            apelp=(pint(i,j,l)+pint(i,j,l+1))*0.5
            apel_3d(i,j,l)=apelp
            dfdp=(q(i,j,l)*0.608+(1.-cw(i,j,l)))*t(i,j,l)*r/apelp
            dfip=dfdp*(dsg2(l)*pd(i,j)+pdsg1(l))
            rtop(i,j,l)=dfdp
            fiup=filo(i,j)+dfip
            dfi_3d(i,j,l)=dfip*0.5
            fim_3d(i,j,l)=(filo(i,j)+fiup)*0.5
            filo(i,j)=fiup
          enddo
        enddo
      enddo
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------
!.......................................................................
!$omp parallel do &
!$omp private (l,j,i,fim,apel,dfi,pdxp,ppx,pgx,pdyp,ppy,pgy, &
!$omp          pdnep,pdnwp,ppne,ppnw,pgne,pgnw,jcl,jch,wprp,rdv,rdu,apd,&
!$omp          rpdp)
!.......................................................................
!-----------------------------------------------------------------------
!---vertical grand loop-------------------------------------------------
!-----------------------------------------------------------------------
      vertical_loop: do l=lm,1,-1
        do j=jts_h1,jte_h1
          do i=its_h1,ite_h1
           apel(i,j)= apel_3d(i,j,l)
           fim(i,j) = fim_3d(i,j,l)
           dfi(i,j) = dfi_3d(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---pressure gradient force components, behind h points-----------------
!-----------------------------------------------------------------------
        do j=jts,jte_h1
          do i=its_b1,ite_h1
            pdxp=(dsg2(l)*pdx(i,j)+pdsg1(l))*0.5
            ppx=(fim(i,j)-fim(i-1,j))  &
               *(dwdt(i-1,j,l)+dwdt(i,j,l))*pdxp
            pcx(i,j,l)=(dfi (i-1,j  )+dfi (i  ,j  )) &
                      *(apel(i  ,j  )-apel(i-1,j  ))
            pgx(i,j)=ppx+pcx(i,j,l)
          enddo
        enddo
!
        do j=jts_b1,jte_h1
          do i=its,ite_h1
            pdyp=(dsg2(l)*pdy(i,j)+pdsg1(l))*0.5
            ppy=(fim(i,j)-fim(i,j-1))  &
               *(dwdt(i,j-1,l)+dwdt(i,j,l))*pdyp
            pcy(i,j,l)=(dfi (i  ,j-1)+dfi (i  ,j  )) &
                      *(apel(i  ,j  )-apel(i  ,j-1))
            pgy(i,j)=ppy+pcy(i,j,l)
          enddo
        enddo
!
        do j=jts_b1,jte_h1
          do i=its_b1,ite_h1
            pdnep=(dsg2(l)*pdne(i,j)+pdsg1(l))*0.5
            pdnwp=(dsg2(l)*pdnw(i,j)+pdsg1(l))*0.5
            ppne=(fim(i,j)-fim(i-1,j-1))  &
                *(dwdt(i-1,j-1,l)+dwdt(i,j,l))*pdnep
            ppnw=(fim(i-1,j)-fim(i,j-1))  &
                *(dwdt(i,j-1,l)+dwdt(i-1,j,l))*pdnwp
            pcne(i,j,l)=(dfi (i-1,j-1)+dfi (i  ,j  )) & 
                       *(apel(i  ,j  )-apel(i-1,j-1))
            pcnw(i,j,l)=(dfi (i  ,j-1)+dfi (i-1,j  )) & 
                       *(apel(i-1,j  )-apel(i  ,j-1))
            pgne(i,j)=ppne+pcne(i,j,l)
            pgnw(i,j)=ppnw+pcnw(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---divergence correction, janjic 1974 mwr, 1979 beitrage---------------
!-----------------------------------------------------------------------
        jcl=jds
        jch=jds-1+jw
!
        if(jts<=jch)then
          do j=max(jts,jcl),min(jte,jch)
            do i=its,ite
              div(i,j,l)=0.
            enddo
          enddo
        endif
!
        jcl=jde-jw+1
        jch=jde
!
        if(jte>=jcl)then
          do j=max(jts,jcl),min(jte,jch)
            do i=its,ite
              div(i,j,l)=0.
            enddo
          enddo
        endif
!
        jcl=jds+jw
        if(jts<=jcl.and.jte>=jcl)then
          wprp=wpdar(jcl)
          do i=its_b1,ite_b1
            div(i,jcl,l)=(pgx(i+1,jcl)-pgx(i,jcl)  &
                          +pgy(i  ,jcl+1)  &
                          -( pgne(i+1,jcl+1)  &
                            +pgnw(i  ,jcl+1))*0.5)*wprp
          enddo
        endif
!
        jch=jde-jw
        if(jts<=jch.and.jte>=jch)then
          wprp=wpdar(jch)
          do i=its_b1,ite_b1
            div(i,jch,l)=(pgx(i+1,jch)-pgx(i,jch)  &
                           -pgy(i  ,jch)  &
                           -(-pgne(i  ,jch)  &
                             -pgnw(i+1,jch))*0.5)*wprp
          enddo
        endif
!
        jcl=jds+jw+1
        jch=jde-jw-1
        if(jte>=jcl.and.jts<=jch)then
          do j=max(jts,jcl),min(jte,jch)
            wprp=wpdar(j)
            do i=its_b1,ite_b1
              div(i,j,l)=(pgx(i+1,j)-pgx(i,j)  &
                         +pgy(i,j+1)-pgy(i,j)  &
                        -(pgne(i+1,j+1)-pgne(i,j)  &
                         +pgnw(i,j+1)-pgnw(i+1,j))*0.5)*wprp
            enddo
          enddo
        endif
!-----------------------------------------------------------------------
!---first pass switch---------------------------------------------------
!-----------------------------------------------------------------------
        if(.not.first) then
!-----------------------------------------------------------------------
!---updating u and v due to pgf force, end of time-step for u and v-----
!-----------------------------------------------------------------------
          rdv=rdyv*dt
          do j=jts_b1,jte_b2
            rdu=rdxv(j)*dt
            do i=its_b1,ite_b2
              apd=(pd(i,j)+pd(i+1,j)+pd(i,j+1)+pd(i+1,j+1))*0.25
              rpdp=0.3333333333/(dsg2(l)*apd+pdsg1(l))
!
              tcu(i,j,l)=-((pgx(i+1,j)+pgx(i+1,j+1)) &
                          +(pgne(i+1,j+1)-pgnw(i+1,j+1))*0.5)*rdu*rpdp &
                         +tcu(i,j,l)
              tcv(i,j,l)=-((pgy(i,j+1)+pgy(i+1,j+1)) &
                          +(pgne(i+1,j+1)+pgnw(i+1,j+1))*0.5)*rdv*rpdp &
                         +tcv(i,j,l)
            enddo
          enddo
!-----------------------------------------------------------------------
        else
!-----------------------------------------------------------------------
          do j=jts_b1,jte_b2
            do i=its_b1,ite_b2
              tcu(i,j,l)=0.
              tcv(i,j,l)=0.
            enddo
          enddo
!-----------------------------------------------------------------------
        endif
!-----------------------------------------------------------------------
!
      enddo vertical_loop
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
                        endsubroutine pgforce
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine dht &
(lm &
,dyv &
,dsg2,pdsg1 &
,dxv,fcp,fdiv &
,pd,pdo &
,u,v &
,omgalf &
!---temporary arguments-------------------------------------------------
,pcne,pcnw,pcx,pcy,pfne,pfnw,pfx,pfy,div,tdiv)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 cfc=1.533 &                 ! adams-bashforth positioning in time
,bfc=1.-cfc                  ! adams bashforth positioning in time
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm                          ! total # of levels
 
real(kind=kfpt),intent(in):: &
 dyv                         ! deltay, v point
 
real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures
 
real(kind=kfpt),dimension(jds:jde),intent(in):: &
 dxv &                       ! deltax, v point
,fcp &                       !
,fdiv                        !
 
real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd &                        ! sigma range pressure difference
,pdo                         ! sigma range pressure difference
 
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 u &                         ! u wind component
,v                           ! v wind component
 
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out):: &
 omgalf                      ! omega-alfa (horizontal)
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in) :: &
 pcne &                      ! second term of pgf, ne direction
,pcnw &                      ! second term of pgf, nw direction
,pcx &                       ! second term of pgf, x direction
,pcy                         ! second term of pgf, y direction
 
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out) :: &
 pfne &                      ! mass flux, ne direction
,pfnw &                      ! mass flux, nw direction
,pfx &                       ! mass flux, x direction
,pfy &                       ! mass flux, y direction
,tdiv                        ! integrated horizontal mass divergence

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout) :: &
 div                         ! horizontal mass divergence
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction
 
real(kind=kfpt):: &
 dyp1 &                      !
,dyp &                       !
,dxp1 &                      !
,dxp &                       !
,fcpp &                      !
,fdp &                       !
,pdp &                       ! hydrostatic pressure difference at the point
,pdxp &                      ! hydrostatic pressure at the point
,pdyp &                      ! hydrostatic pressure at the point
,pdnep &                     ! hydrostatic pressure at the point
,pdnwp &                     ! hydrostatic pressure at the point
,udy &                       !
,vdx &                       !
,udy1 &                      !
,vdx1                        !
 
real(kind=kfpt),dimension(its_b1:ite_h2,jts_b1:jte_h2):: &
 pdne &                      ! hydrostatic pressure difference at the point
,pdnw &                      ! hydrostatic pressure difference at the point
,pdx &                       ! hydrostatic pressure difference at the point
,pdy &                       ! hydrostatic pressure difference at the point
,tne &                       ! temperature flux, ne direction
,tnw &                       ! temperature flux, nw direction
,tx &                        ! temperature flux, x direction
,ty                          ! temperature flux, y direction
 
!-----------------
!-----------------
integer(kind=kint) :: &
 jstart &
,jstop &
,nth &
,omp_get_num_threads &
,omp_get_thread_num &
,tid
!-----------------
!-----------------
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      dyp1=dyv
      dyp=dyp1*0.5
!
!-----------------
!-----------------
!.......................................................................
!$omp parallel &
!$omp private(nth,tid,i,j,jstart,jstop)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth,jts_b1,jte_h2,jstart,jstop)
!-----------------
!-----------------
!
      do j=jstart,jstop
        do i=its_b1,ite_h2
          pdx (i,j)=((pd (i-1,j)+pd (i,j))*cfc &
                    +(pdo(i-1,j)+pdo(i,j))*bfc)*0.5
        enddo
      enddo
!
      do j=jstart,jstop
        do i=its_b1,ite_h2
          pdy (i,j)=((pd (i,j-1)+pd (i,j))*cfc &
                    +(pdo(i,j-1)+pdo(i,j))*bfc)*0.5
        enddo
      enddo
!
      do j=jstart,jstop
        do i=its_b1,ite_h2
          pdne(i,j)=((pd (i-1,j-1)+pd (i,j))*cfc &
                    +(pdo(i-1,j-1)+pdo(i,j))*bfc)*0.5
          pdnw(i,j)=((pd (i,j-1)+pd (i-1,j))*cfc &
                    +(pdo(i,j-1)+pdo(i-1,j))*bfc)*0.5
        enddo
      enddo
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------
!-----------------------------------------------------------------------
!---vertical grand loop-------------------------------------------------
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp private (l,j,i,pdp,dxp1,dxp,pdxp,pdyp,pdnep,pdnwp,udy,vdx,udy1,vdx1,&
!$omp          tx,ty,tne,tnw,fdp,fcpp)
!.......................................................................
!-----------------------------------------------------------------------
      vertical_loop: do l=1,lm
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!---mass and pressure fluxes, on h points-------------------------------
!-----------------------------------------------------------------------
        pdp=pdsg1(l)
        do j=jts_b1,jte_h2
          dxp1=dxv(j-1)
          dxp=dxp1*0.5
          do i=its_b1,ite_h2
!-----------------------------------------------------------------------
            pdxp=dsg2(l)*pdx(i,j)+pdp
            pdyp=dsg2(l)*pdy(i,j)+pdp
            pdnep=dsg2(l)*pdne(i,j)+pdp
            pdnwp=dsg2(l)*pdnw(i,j)+pdp
!
            udy=(u(i-1,j-1,l)+u(i-1,j,l))*dyp
            vdx=(v(i-1,j-1,l)+v(i,j-1,l))*dxp
!
            pfx(i,j,l)=udy*pdxp
            pfy(i,j,l)=vdx*pdyp
!
            udy1=u(i-1,j-1,l)*dyp1
            vdx1=v(i-1,j-1,l)*dxp1
!
            pfne(i,j,l)=( udy1+vdx1)*pdnep
            pfnw(i,j,l)=(-udy1+vdx1)*pdnwp
!
            tx(i,j)=pcx(i,j,l)*udy
            ty(i,j)=pcy(i,j,l)*vdx
!
            tne(i,j)=pcne(i,j,l)*( udy1+vdx1)
            tnw(i,j)=pcnw(i,j,l)*(-udy1+vdx1)
          enddo
        enddo
!-----------------------------------------------------------------------
!---divergence and hor. pressure advection in t eq.---------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1_h1
          fdp=fdiv(j)
          fcpp=fcp(j)
          do i=its_b1,ite_b1_h1
            div(i,j,l)=((pfx (i+1,j  ,l)-pfx (i  ,j  ,l) &
                        +pfy (i  ,j+1,l)-pfy (i  ,j  ,l)) &
                       +(pfne(i+1,j+1,l)-pfne(i  ,j  ,l) &
                        +pfnw(i  ,j+1,l)-pfnw(i+1,j  ,l))*0.25)*fdp  &
                      +div(i,j,l)
            tdiv(i,j,l)=div(i,j,l)
            omgalf(i,j,l)=((tx (i  ,j  )+tx (i+1,j  ) &
                           +ty (i  ,j  )+ty (i  ,j+1)) &
                          +(tne(i+1,j+1)+tne(i  ,j  ) &
                           +tnw(i  ,j+1)+tnw(i+1,j  ))*0.25) &
                         *fcpp/(dsg2(l)*pd(i,j)+pdsg1(l))
          enddo
        enddo
!-----------------------------------------------------------------------
      enddo vertical_loop
!-----------------------------------------------------------------------
!.......................................................................
!$omp end parallel do
!.......................................................................
!-----------------
!-----------------
!.......................................................................
!$omp parallel &
!$omp private(nth,tid,i,j,jstart,jstop)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth,jts_b1,jte_b1,jstart,jstop)
!-----------------
!-----------------
      do l=2,lm
        do j=jstart,jstop
          do i=its_b1,ite_b1
            tdiv(i,j,l)=tdiv(i,j,l-1)+tdiv(i,j,l)
          enddo
        enddo
      enddo
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------
     
!-----------------------------------------------------------------------
!
                        endsubroutine dht
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine ddamp &
(lm &
,ddmpv,pdtop &
,dsg2,pdsg1 &
,ddmpu &
,pd,pdo &
,u,v &
!---temporary arguments-------------------------------------------------
,div)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 cfc=1.533 &                 ! adams-bashforth positioning in time
,bfc=1.-cfc                  ! adams bashforth positioning in time
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm                          ! total # of levels
 
real(kind=kfpt),intent(in):: &
 ddmpv &                     ! divergence damping, v component
,pdtop                       ! pressure coordinate depth
 
real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures
 
real(kind=kfpt),dimension(jds:jde),intent(in):: &
 ddmpu                       ! divergence damping, u direction
 
real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd &                        ! sigma range pressure difference
,pdo                         ! sigma range pressure difference
 
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 u &                         ! u wind component
,v                           ! v wind component
 
!---temporary arguments-------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 div                         ! horizontal mass divergence
 
!--local variables------------------------------------------------------
logical(kind=klog):: &
 extmod                      ! external mode divergence damping
 
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction
 
real(kind=kfpt):: &
 fcim &                      !
,fcxm &                      ! blow up factor for external mode damping
,rdpd &                      !
,dud &                       !
,dvd                         !
 
real(kind=kfpt),dimension(ims:ime,jms:jme):: &
 apd &                       !
,dive &                      !
,rddu &                      !
,rddv                        !
!-----------------
!-----------------
integer(kind=kint) :: &
 jstart &
,jstop &
,nth &
,omp_get_num_threads &
,omp_get_thread_num &
,tid
!-----------------
!-----------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      dvd=ddmpv
      do j=jts_b1,jte_b2
        do i=its_b1,ite_b2
          apd(i,j)=((pd (i  ,j  )+pd (i+1,j  ) &
                    +pd (i  ,j+1)+pd (i+1,j+1))*cfc &
                   +(pdo(i  ,j  )+pdo(i+1,j  ) &
                    +pdo(i  ,j+1)+pdo(i+1,j+1))*bfc)*0.25
        enddo
      enddo
!
!-----------------------------------------------------------------------
!---external mode-------------------------------------------------------
!-----------------------------------------------------------------------
!
      extmod=.false.
      extmod=.true.
      if(extmod) then
        fcxm=1.
        fcim=1.
!
!-----------------
!-----------------
!.......................................................................
!$omp parallel private(nth,tid, i,j,l,jstart,jstop)
!.......................................................................
        nth = omp_get_num_threads()
        tid = omp_get_thread_num()
        call looplimits(tid, nth,jts_b1,jte_b1_h2,jstart,jstop)
!-----------------
!-----------------
        do j=jstart,jstop
          do i=its_b1,ite_b1_h2
            dive(i,j)=div(i,j,1)
          enddo
        enddo
        do l=2,lm
        do j=jstart,jstop
          do i=its_b1,ite_b1_h2
              dive(i,j)=div(i,j,l)+dive(i,j)
            enddo
          enddo
        enddo
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------
!
        do j=jts_b1,jte_b2
          dud=ddmpu(j)
          do i=its_b1,ite_b2
            rdpd=fcxm/(apd(i,j)+pdtop)
            rddu(i,j)=(dive(i+1,j)+dive(i+1,j+1) &
                      -dive(i  ,j)-dive(i  ,j+1))*dud*rdpd
            rddv(i,j)=(dive(i,j+1)+dive(i+1,j+1) &
                      -dive(i,j  )-dive(i+1,j  ))*dvd*rdpd
          enddo
        enddo
!
        jstart = jts_b1
        jstop = jte_b2
!.......................................................................
!$omp parallel do private(l,j,i,dud,rdpd)
!.......................................................................
        do l=1,lm
          do j=jstart,jstop
            dud=ddmpu(j)
            do i=its_b1,ite_b2
              rdpd=fcim/(dsg2(l)*apd(i,j)+pdsg1(l))
              u(i,j,l)=((div(i+1,j,l)+div(i+1,j+1,l) &
                        -div(i  ,j,l)-div(i  ,j+1,l))*dud)*rdpd &
                      +rddu(i,j)+u(i,j,l)
              v(i,j,l)=((div(i,j+1,l)+div(i+1,j+1,l) &
                        -div(i,j  ,l)-div(i+1,j  ,l))*dvd)*rdpd &
                      +rddv(i,j)+v(i,j,l)
            enddo
          enddo
        enddo      ! end of vertical loop
!.......................................................................
!$omp end parallel do
!.......................................................................

!-----------------------------------------------------------------------
!
      else
!
!-----------------------------------------------------------------------
!---divergence damping--------------------------------------------------
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do private(l,j,i,dvd,dud,rdpd)
!.......................................................................
        do l=1,lm
          dvd=ddmpv
          do j=jts_b1,jte_b2
            dud=ddmpu(j)
            do i=its_b1,ite_b2
              rdpd=1./(dsg2(l)*apd(i,j)+pdsg1(l))
              u(i,j,l)=(div(i+1,j,l)+div(i+1,j+1,l) &
                       -div(i  ,j,l)-div(i  ,j+1,l)) &
                      *dud*rdpd+u(i,j,l)
              v(i,j,l)=(div(i,j+1,l)+div(i+1,j+1,l) &
                       -div(i,j  ,l)-div(i+1,j  ,l)) &
                      *dvd*rdpd+v(i,j,l)
            enddo
          enddo
        enddo      ! end of vertical loop
!.......................................................................
!$omp end parallel do
!.......................................................................
!-----------------------------------------------------------------------
!
      endif
!
!-----------------------------------------------------------------------
!
                        endsubroutine ddamp
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine pdtsdt &
(lm &
,dt &
,sg2 &
,pd &
,pdo,psdt &
,psgdt &
!---temporary arguments-------------------------------------------------
,div,tdiv)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm                          ! total # of levels

real(kind=kfpt),intent(in):: &
 dt                          ! time step

real(kind=kfpt),dimension(1:lm+1),intent(in):: &
 sg2                        !

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout):: &
 pd                          ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(out):: &
 pdo &                       ! sigma range pressure difference
,psdt                        ! hydrostatic pressure tendency

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(out):: &
 psgdt                       ! vertical mass flux
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 div &                       ! horizontal mass divergence
,tdiv                        ! integrated unfiltered mass divergence
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction
!-----------------
!-----------------
integer(kind=kint) :: &
 jstart &
,jstop &
,nth &
,omp_get_num_threads &
,omp_get_thread_num &
,tid
!-----------------
!-----------------
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!-----------------
!-----------------
!.......................................................................
!$omp parallel &
!$omp private(nth,tid,i,j,jstart,jstop)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth,jts_b1,jte_b1,jstart,jstop)
!-----------------
!-----------------
      do l=2,lm
        do j=jstart,jstop
          do i=its_b1,ite_b1
            div(i,j,l)=div(i,j,l-1)+div(i,j,l)
          enddo
        enddo
      enddo
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------

!-----------------------------------------------------------------------
      do j=jts_h2,jte_h2
        do i=its_h2,ite_h2
          psdt(i,j)=0.
          pdo(i,j)=pd(i,j)
        enddo
      enddo
!
!-----------------
!-----------------
!.......................................................................
!$omp parallel &
!$omp private(nth,tid,i,j,jstart,jstop)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth,jts_b1,jte_b1,jstart,jstop)
!-----------------
!-----------------
      do j=jstart,jstop
        do i=its_b1,ite_b1
          psdt(i,j)=-div(i,j,lm)
          pd(i,j)=psdt(i,j)*dt+pd(i,j)
        enddo
      enddo
!-----------------------------------------------------------------------
!---boundary conditions-------------------------------------------------
!-----------------------------------------------------------------------
      do l=1,lm-1
        do j=jstart,jstop
          do i=its_b1,ite_b1
            psgdt(i,j,l)=-(-tdiv(i,j,lm)*sg2(l+1)+tdiv(i,j,l))
          enddo
        enddo
      enddo
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------
!-----------------------------------------------------------------------
!
                        endsubroutine pdtsdt
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine adv1 &
(global,secadv &
,lm,lnsad,inpes,jnpes &
,dt,dyv,rdyh,rdyv &
,dsg2,pdsg1 &
,curv,dxv,fad,fah,rdxh,rdxv &
,f,pd,pdo &
,omgalf,psgdt &
,t,u,v &
,tp,up,vp &
!---temporary arguments-------------------------------------------------
,pfne,pfnw,pfx,pfy,tct,tcu,tcv)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 cfc=1.533 &                 ! adams-bashforth positioning in time
,bfc=1.-cfc &                ! adams bashforth positioning in time
,epscm=2.e-6 &               ! a floor value (not used)
,pfc=1.+4./6. &              ! 4th order momentum advection
,sfc=-1./6. &                ! 4th order momentum advection
,w1=1.0 &                    ! crank-nicholson uncentering
!,w1=0.0 &                    ! crank-nicholson uncentering
,w2=2.-w1                    ! crank-nicholson uncentering
 
logical(kind=klog),intent(in):: &
 global &                    ! global or regional
,secadv                      ! second order momentum advection
 
integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,lnsad &                     ! # of boundary lines w. upstream advection
,inpes &                     ! domain decomposition parameter
,jnpes                       ! domain decomposition parameter
 
real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,dyv &                       ! deltay
,rdyh &                      ! 1/deltay
,rdyv                        ! 1/deltay
 
real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures
 
real(kind=kfpt),dimension(jds:jde),intent(in):: &
 curv &                      ! curvature
,dxv &                       ! dxv
,fad &                       ! grid factor
,fah &                       ! grid factor
,rdxh &                      ! 1/deltax
,rdxv                        ! 1/deltax
 
real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 f &                         ! coriolis parameter
,pd &                        ! sigma range pressure difference
,pdo                         ! sigma range pressure difference
 
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(in):: &
 psgdt                       ! vertical mass flux
 
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 omgalf &                    !
,t &                         ! temperature
,u &                         ! u wind component
,v                           ! v wind component
 
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 tp &                        ! old temperature
,up &                        ! old u
,vp                          ! old v
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in) :: &
 pfne &                      ! mass flux, ne direction
,pfnw &                      ! mass flux, nw direction
,pfx &                       ! mass flux, x direction
,pfy                         ! mass flux, y direction
 
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out) :: &
 tct &                       ! time change of temperature
,tcu &                       ! time change of u
,tcv                         ! time change of v
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,iap &                       ! offset in x direction
,ibeg &                      ! starting i in some horiz advec loops
,iend &                      ! ending i in some horiz advec loops
,j &                         ! index in y direction
,jap &                       ! offset in y direction
,jbeg &                      ! starting j in some horiz advec loops
,jend &                      ! ending j in some horiz advec loops
,l                           ! index in p direction
 
real(kind=kfpt):: &
 cf &                        ! temporary
,cmt &                       ! temporary
,cmw &                       ! temporary
,crv &                       ! curvature
,dtq &                       ! dt/4
,dts &                       ! dt/16
,dux1 &                      ! momentum advection component
,duy1 &                      ! momentum advection component
,dvx1 &                      ! momentum advection component
,dvy1 &                      ! momentum advection component
,dxody &                     !
,dyodx &                     !
,emhp &                      ! scaling in x direction
,enhp &                      ! scaling in y direction
,emvp &                      ! scaling in x direction
,envp &                      ! scaling in y direction
,fadp &                      ! grid factor at the point
,fahp &                      ! temporary grid factor
,fdpp &                      ! temporary grid factor
,fp &                        ! coriolis parameter factor
,fpp &                       ! coriolis with curvature
,pp &                        ! scaled trajectory, x direction
,qq &                        ! scaled trajectory, y direction
,rdp &                       ! 1/deltap
,rdxp &                      !
,rdyp &                      !
,vvlo &                      ! vertical velocity, lower interface
,vvup &                      ! vertical velocity, upper interface
,pvvup                       ! vertical mass flux, upper interface
 
real(kind=kfpt),dimension(its:ite,jts:jte):: &
 pdop &                      ! hydrostatic pressure difference at v points
,pvvlo                       ! vertical mass flux, lower interface
 
real(kind=kfpt),dimension(its_b1:ite_b1_h1,jts_b1:jte_b1_h1):: &
 pfnex1 &                    ! average mass flux for momentum advection
,pfney1 &                    ! average mass flux for momentum advection
,pfnwx1 &                    ! average mass flux for momentum advection
,pfnwy1 &                    ! average mass flux for momentum advection
,pfxx1 &                     ! average mass flux for momentum advection
,pfxy1 &                     ! average mass flux for momentum advection
,pfyx1 &                     ! average mass flux for momentum advection
,pfyy1 &                     ! average mass flux for momentum advection
,ufnex1 &                    ! average mass flux for momentum advection
,ufney1 &                    ! average mass flux for momentum advection
,ufnwx1 &                    ! average mass flux for momentum advection
,ufnwy1 &                    ! average mass flux for momentum advection
,ufxx1 &                     ! average mass flux for momentum advection
,ufxy1 &                     ! average mass flux for momentum advection
,ufyx1 &                     ! average mass flux for momentum advection
,ufyy1 &                     ! average mass flux for momentum advection
,vfnex1 &                    ! average mass flux for momentum advection
,vfney1 &                    ! average mass flux for momentum advection
,vfnwx1 &                    ! average mass flux for momentum advection
,vfnwy1 &                    ! average mass flux for momentum advection
,vfxx1 &                     ! average mass flux for momentum advection
,vfxy1 &                     ! average mass flux for momentum advection
,vfyx1 &                     ! average mass flux for momentum advection
,vfyy1                       ! average mass flux for momentum advection
 
real(kind=kfpt),dimension(its_b1:ite_h1,jts_b1:jte_h1):: &
 tne &                       ! temperature flux, ne direction
,tnw &                       ! temperature flux, nw direction
,tx &                        ! temperature flux, x direction
,ty                          ! temperature flux, y direction
 
real(kind=kfpt),dimension(its_h1:ite_h1,jts_h1:jte_h1):: &
 t1                          ! extrapolated temperature between time levels
 
real(kind=kfpt),dimension(ims:ime,jms:jme):: &
 u2d &                       ! 4th order diagonal u between time levels
,v2d &                       ! 4th order diagonal v between time levels
 
!real(kind=kfpt),dimension(its_h1:ite_h1,jts_h1:jte_h1):: &
,u1d &                       ! extrapolated diagonal u between time levels
,v1d                         ! extrapolated diagonal v between time levels
 
real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 crt &                       ! vertical advection temporary
,rcmt &                      ! vertical advection temporary
,rstt                        ! vertical advection temporary
 
real(kind=kfpt),dimension(its_b1:ite_b2,jts_b1:jte_b2,1:lm):: &
 crw &                       ! vertical advection temporary
,rcmw &                      ! vertical advection temporary
,rstu &                      ! vertical advection temporary
,rstv                        ! vertical advection temporary
!-----------------
!-----------------
integer(kind=kint) :: &
 jstart &
,jstop &
,nth &
,omp_get_num_threads &
,omp_get_thread_num &
,tid
!-----------------
!-----------------
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      do j=jts,jte
        do i=its,ite
          pdop(i,j)=(pd(i,j)+pdo(i,j))*0.5
        enddo
      enddo
!-----------------------------------------------------------------------
!---crank-nicholson vertical advection----------------------------------
!-----------------------------------------------------------------------
!
      dtq=dt*0.25
!
!-----------------
!-----------------
!.......................................................................
!$omp parallel private(nth, tid, i,j,l,jstart, jstop,rdp,pvvup,vvup,vvlo,cf,cmt)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth, jts_b1, jte_b1, jstart, jstop)
!-----------------
!-----------------
       do j=jstart, jstop
        do i=its_b1,ite_b1
          pvvlo(i,j)=psgdt(i,j,1)*dtq
          vvlo=pvvlo(i,j)/(dsg2(1)*pdop(i,j)+pdsg1(1))
          cmt=-vvlo*w2+1.
!          if(abs(cmt).lt.epscm) cmt=epscm
          rcmt(i,j,1)=1./cmt
          crt(i,j,1)=vvlo*w2
          rstt(i,j,1)=(-vvlo*w1)*(t(i,j,2)-t(i,j,1))+t(i,j,1)
        enddo
      enddo
!
      do l=2,lm-1
       do j=jstart, jstop
          do i=its_b1,ite_b1
            rdp=1./(dsg2(l)*pdop(i,j)+pdsg1(l))
            pvvup=pvvlo(i,j)
            pvvlo(i,j)=psgdt(i,j,l)*dtq
!
            vvup=pvvup*rdp
            vvlo=pvvlo(i,j)*rdp
!
            cf=-vvup*w2*rcmt(i,j,l-1)
            cmt=-crt(i,j,l-1)*cf+((vvup-vvlo)*w2+1.)
!            if(abs(cmt).lt.epscm) cmt=epscm
            rcmt(i,j,l)=1./cmt
            crt(i,j,l)=vvlo*w2
            rstt(i,j,l)=-rstt(i,j,l-1)*cf+t(i,j,l) &
                        -((t(i,j,l)-t(i,j,l-1))*vvup &
                         +(t(i,j,l+1)-t(i,j,l))*vvlo)*w1
          enddo
        enddo
      enddo
!
       do j=jstart, jstop
        do i=its_b1,ite_b1
          pvvup=pvvlo(i,j)
          vvup=pvvup/(dsg2(lm)*pdop(i,j)+pdsg1(lm))
!
          cf=-vvup*w2*rcmt(i,j,lm-1)
          cmt=-crt(i,j,lm-1)*cf+(vvup*w2+1.)
!          if(abs(cmt).lt.epscm) cmt=epscm
          rcmt(i,j,lm)=1./cmt
          crt(i,j,lm)=0.
          rstt(i,j,lm)=-rstt(i,j,lm-1)*cf+t(i,j,lm) &
     &                 -(t(i,j,lm)-t(i,j,lm-1))*vvup*w1
!
          tct(i,j,lm)=rstt(i,j,lm)*rcmt(i,j,lm)-t(i,j,lm)
        enddo
      enddo
!
      do l=lm-1,1,-1
       do j=jstart, jstop
          do i=its_b1,ite_b1
            tct(i,j,l)=(-crt(i,j,l)*(t(i,j,l+1)+tct(i,j,l+1)) &
                        +rstt(i,j,l)) &
                      *rcmt(i,j,l)-t(i,j,l)
          enddo
        enddo
      enddo
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------

!.......................................................................
!$omp parallel do &
!$omp private (t1,tx,ty,tne,tnw,ibeg,iend,jbeg,jend,fahp,enhp,pp,qq, &
!$omp         iap,jap,emhp,j,i)
!.......................................................................
!-----------------------------------------------------------------------
      vertical_loop1: do l=1,lm
!-----------------------------------------------------------------------
        do j=jts_h1,jte_h1
          do i=its_h1,ite_h1
            t1(i,j)=t(i,j,l)*cfc+tp(i,j,l)*bfc
            tp(i,j,l)=t(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---temperature fluxes, on h points-------------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_h1
          do i=its_b1,ite_h1
            tx(i,j)=(t1(i,j)-t1(i-1,j))*pfx(i,j,l)
            ty(i,j)=(t1(i,j)-t1(i,j-1))*pfy(i,j,l)
            tne(i,j)=(t1(i,j)-t1(i-1,j-1))*pfne(i,j,l)
            tnw(i,j)=(t1(i-1,j)-t1(i,j-1))*pfnw(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---advection of temperature--------------------------------------------
!-----------------------------------------------------------------------
        if(adv_standard)then
          ibeg=max(its,ids+1+lnsad)
          iend=min(ite,ide-1-lnsad)
          jbeg=max(jts,jds+1+lnsad)
          jend=min(jte,jde-1-lnsad)
!
          do j=jbeg,jend
            fahp=fah(j)
            do i=ibeg,iend
            tct(i,j,l)=(((tx(i,j)+tx(i+1,j)+ty(i,j)+ty(i,j+1)) &
                        +(tne(i+1,j+1)+tne(i,j) &
                         +tnw(i,j+1)+tnw(i+1,j))*0.25)*fahp) &
                      /(dsg2(l)*pdop(i,j)+pdsg1(l)) &
                      +omgalf(i,j,l)+tct(i,j,l)
            enddo
          enddo
        endif
!-----------------------------------------------------------------------
!---regional branch-----------------------------------------------------
!-----------------------------------------------------------------------
        if(.not.global.and.adv_upstream) then
!-----------------------------------------------------------------------
          enhp=-dt*rdyh*0.25
!
!***  Upstream advection along southern rows
!
          do j=jts_b1,min(jte,jds+lnsad)
            emhp=-dt*rdxh(j)*0.25
            do i=its_b1,ite_b1
              pp=(u(i-1,j-1,l)+u(i,j-1,l)+u(i-1,j,l)+u(i,j,l))*emhp
              qq=(v(i-1,j-1,l)+v(i,j-1,l)+v(i-1,j,l)+v(i,j,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tct(i,j,l)=(t(i+iap,j,l)-t(i,j,l))*pp &
                        +(t(i,j+jap,l)-t(i,j,l))*qq &
                        +(t(i,j,l)-t(i+iap,j,l) &
                         -t(i,j+jap,l)+t(i+iap,j+jap,l))*pp*qq &
                        +omgalf(i,j,l)+tct(i,j,l)
            enddo
          enddo
!
!***  Upstream advection along northern rows
!
          do j=max(jts,jde-lnsad),jte_b1
            emhp=-dt*rdxh(j)*0.25
            do i=its_b1,ite_b1
              pp=(u(i-1,j-1,l)+u(i,j-1,l)+u(i-1,j,l)+u(i,j,l))*emhp
              qq=(v(i-1,j-1,l)+v(i,j-1,l)+v(i-1,j,l)+v(i,j,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tct(i,j,l)=(t(i+iap,j,l)-t(i,j,l))*pp &
                        +(t(i,j+jap,l)-t(i,j,l))*qq &
                        +(t(i,j,l)-t(i+iap,j,l) &
                         -t(i,j+jap,l)+t(i+iap,j+jap,l))*pp*qq &
                        +omgalf(i,j,l)+tct(i,j,l)
            enddo
          enddo
!
!***  Upstream advection along western rows
!
          do j=max(jts,jds+1+lnsad),min(jte,jde-1-lnsad)
            emhp=-dt*rdxh(j)*0.25
            do i=its_b1,min(ite,ids+lnsad)
              pp=(u(i-1,j-1,l)+u(i,j-1,l)+u(i-1,j,l)+u(i,j,l))*emhp
              qq=(v(i-1,j-1,l)+v(i,j-1,l)+v(i-1,j,l)+v(i,j,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tct(i,j,l)=(t(i+iap,j,l)-t(i,j,l))*pp &
                        +(t(i,j+jap,l)-t(i,j,l))*qq &
                        +(t(i,j,l)-t(i+iap,j,l) &
                         -t(i,j+jap,l)+t(i+iap,j+jap,l))*pp*qq &
                        +omgalf(i,j,l)+tct(i,j,l)
            enddo
          enddo
!
!***  Upstream advection along eastern rows
!
          do j=max(jts,jds+1+lnsad),min(jte,jde-1-lnsad)
            emhp=-dt*rdxh(j)*0.25
            do i=max(its,ide-lnsad),ite_b1
              pp=(u(i-1,j-1,l)+u(i,j-1,l)+u(i-1,j,l)+u(i,j,l))*emhp
              qq=(v(i-1,j-1,l)+v(i,j-1,l)+v(i-1,j,l)+v(i,j,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tct(i,j,l)=(t(i+iap,j,l)-t(i,j,l))*pp &
                        +(t(i,j+jap,l)-t(i,j,l))*qq &
                        +(t(i,j,l)-t(i+iap,j,l) &
                         -t(i,j+jap,l)+t(i+iap,j+jap,l))*pp*qq &
                        +omgalf(i,j,l)+tct(i,j,l)
            enddo
          enddo
!-----------------------------------------------------------------------
        endif ! regional lateral boundaries
!-----------------------------------------------------------------------
!
      enddo vertical_loop1
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!---crank-nicholson vertical advection----------------------------------
!-----------------------------------------------------------------------
!
      dts=dt*(0.25*0.25)
!
!-----------------
!-----------------
!.......................................................................
!$omp parallel private(nth,tid,i,j,l,jstart,jstop,vvlo,cmw,pvvup,vvup,cf,rdp)
!.......................................................................
       nth = omp_get_num_threads()
       tid = omp_get_thread_num()
       call looplimits(tid, nth, jts_b1, jte_b2, jstart, jstop)
!-----------------
!-----------------
      do j=jstart, jstop
        do i=its_b1,ite_b2
          pdop(i,j)=(pd (i,j)+pd (i+1,j)+pd (i,j+1)+pd (i+1,j+1) &
                    +pdo(i,j)+pdo(i+1,j)+pdo(i,j+1)+pdo(i+1,j+1))*0.125
          pvvlo(i,j)=(psgdt(i,j  ,1)+psgdt(i+1,j  ,1) &
                     +psgdt(i,j+1,1)+psgdt(i+1,j+1,1))*dts
          vvlo=pvvlo(i,j)/(dsg2(1)*pdop(i,j)+pdsg1(1))
          cmw=-vvlo*w2+1.
!          if(abs(cmw).lt.epscm) cmw=epscm
          rcmw(i,j,1)=1./cmw
          crw(i,j,1)=vvlo*w2
          rstu(i,j,1)=(-vvlo*w1)*(u(i,j,2)-u(i,j,1))+u(i,j,1)
          rstv(i,j,1)=(-vvlo*w1)*(v(i,j,2)-v(i,j,1))+v(i,j,1)
        enddo
      enddo
!
      do l=2,lm-1
        do j=jstart,jstop
          do i=its_b1,ite_b2
            rdp=1./(dsg2(l)*pdop(i,j)+pdsg1(l))
            pvvup=pvvlo(i,j)
            pvvlo(i,j)=(psgdt(i,j,l)+psgdt(i+1,j,l) &
                       +psgdt(i,j+1,l)+psgdt(i+1,j+1,l))*dts
            vvup=pvvup*rdp
            vvlo=pvvlo(i,j)*rdp
            cf=-vvup*w2*rcmw(i,j,l-1)
            cmw=-crw(i,j,l-1)*cf+((vvup-vvlo)*w2+1.)
!            if(abs(cmw).lt.epscm) cmw=epscm
            rcmw(i,j,l)=1./cmw
            crw(i,j,l)=vvlo*w2
            rstu(i,j,l)=-rstu(i,j,l-1)*cf+u(i,j,l) &
                        -((u(i,j,l)-u(i,j,l-1))*vvup &
                         +(u(i,j,l+1)-u(i,j,l))*vvlo)*w1
            rstv(i,j,l)=-rstv(i,j,l-1)*cf+v(i,j,l) &
                        -((v(i,j,l)-v(i,j,l-1))*vvup &
                         +(v(i,j,l+1)-v(i,j,l))*vvlo)*w1
          enddo
        enddo
      enddo
!
      do j=jstart,jstop
        do i=its_b1,ite_b2
          pvvup=pvvlo(i,j)
          vvup=pvvup/(dsg2(lm)*pdop(i,j)+pdsg1(lm))
          cf=-vvup*w2*rcmw(i,j,lm-1)
          cmw=-crw(i,j,lm-1)*cf+(vvup*w2+1.)
!          if(abs(cmw).lt.epscm) cmw=epscm
          rcmw(i,j,lm)=1./cmw
          crw(i,j,lm)=vvlo
          rstu(i,j,lm)=-rstu(i,j,lm-1)*cf+u(i,j,lm) &
                       -(u(i,j,lm)-u(i,j,lm-1))*vvup*w1
          rstv(i,j,lm)=-rstv(i,j,lm-1)*cf+v(i,j,lm) &
                       -(v(i,j,lm)-v(i,j,lm-1))*vvup*w1
          tcu(i,j,lm)=rstu(i,j,lm)*rcmw(i,j,lm)-u(i,j,lm)
          tcv(i,j,lm)=rstv(i,j,lm)*rcmw(i,j,lm)-v(i,j,lm)
        enddo
      enddo
!
      do l=lm-1,1,-1
        do j=jstart,jstop
          do i=its_b1,ite_b2
            tcu(i,j,l)=(-crw(i,j,l)*(u(i,j,l+1)+tcu(i,j,l+1)) &
                        +rstu(i,j,l)) &
                      *rcmw(i,j,l)-u(i,j,l)
            tcv(i,j,l)=(-crw(i,j,l)*(v(i,j,l+1)+tcv(i,j,l+1)) &
                        +rstv(i,j,l)) &
                      *rcmw(i,j,l)-v(i,j,l)
          enddo
        enddo
      enddo
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------

!-----------------------------------------------------------------------
!---grand vertical loop-------------------------------------------------
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do &
!$omp private (u1d,v1d,crv,fp,fpp,pfxx1,pfyx1,pfnex1,pfnwx1,pfxy1,pfyy1, &
!$omp          pfney1,pfnwy1,u2d,v2d,ufxx1,ufyx1,ufnex1,ufnwx1,          &
!$omp          ufxy1,ufyy1,ufney1,ufnwy1,vfxx1,vfyx1,vfnex1,vfnwx1,vfxy1, &
!$omp          vfyy1,vfney1,vfnwy1,rdyp,ibeg,iend,jbeg,jend,fadp,rdxp,    &                     
!$omp          fdpp,dux1,dvx1,duy1,dvy1,envp,emvp, pp,qq,iap,jap,j,i)
!.......................................................................
!-----------------------------------------------------------------------
!
      vertical_loop2: do l=1,lm
!
!-----------------------------------------------------------------------
!
        do j=jts_h2,jte_h2
          do i=its_h2,ite_h2
            u1d(i,j)=u(i,j,l)*cfc+up(i,j,l)*bfc
            v1d(i,j)=v(i,j,l)*cfc+vp(i,j,l)*bfc
!
            up(i,j,l)=u(i,j,l)
            vp(i,j,l)=v(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---coriolis force tendency---------------------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b2
          crv=curv(j)*dt
          do i=its_b1,ite_b2
            fp=f(i,j)*dt
            fpp=u1d(i,j)*crv+fp
!
            tcu(i,j,l)= fpp*v1d(i,j)+tcu(i,j,l)
            tcv(i,j,l)=-fpp*u1d(i,j)+tcv(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---diagonally averaged fluxes on v points------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1_h1
          do i=its_b1,ite_b1_h1
            pfxx1 (i,j)=pfx (i  ,j  ,l)+pfx (i+1,j+1,l)
            pfyx1 (i,j)=pfy (i  ,j  ,l)+pfy (i+1,j+1,l)
            pfnex1(i,j)=pfne(i  ,j  ,l)+pfne(i+1,j+1,l)
            pfnwx1(i,j)=pfnw(i  ,j  ,l)+pfnw(i+1,j+1,l)
!
            pfxy1 (i,j)=pfx (i+1,j  ,l)+pfx (i  ,j+1,l)
            pfyy1 (i,j)=pfy (i+1,j  ,l)+pfy (i  ,j+1,l)
            pfney1(i,j)=pfne(i+1,j  ,l)+pfne(i  ,j+1,l)
            pfnwy1(i,j)=pfnw(i+1,j  ,l)+pfnw(i  ,j+1,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---4th order momentum advection----------------------------------------
!-----------------------------------------------------------------------
        if(.not.secadv) then
          do j=jts_b1_h1,jte_b1_h1
            do i=its_b1_h1,ite_b1_h1
              u2d(i,j)=(u1d(i,j-1)+u1d(i-1,j) &
                       +u1d(i+1,j)+u1d(i,j+1))*sfc &
                      +(u1d(i,j)*pfc)
              v2d(i,j)=(v1d(i,j-1)+v1d(i-1,j) &
                       +v1d(i+1,j)+v1d(i,j+1))*sfc &
                      +(v1d(i,j)*pfc)
            enddo
          enddo
          if(global) then
            btim=timef()
            call swapwn(u2d,ims,ime,jms,jme,1,inpes)
            call swapwn(v2d,ims,ime,jms,jme,1,inpes)
            swapwn_tim=swapwn_tim+timef()-btim
!
            btim=timef()
            call polewn(u2d,v2d,ims,ime,jms,jme,1,inpes,jnpes)
            polewn_tim=polewn_tim+timef()-btim
          else
            if(s_bdy)then
              do i=ims,ime
                u2d(i,jds)=u1d(i,jds)
                v2d(i,jds)=v1d(i,jds)
              enddo
            endif
            if(n_bdy)then
              do i=ims,ime
                u2d(i,jde-1)=u1d(i,jde-1)
                v2d(i,jde-1)=v1d(i,jde-1)
                u2d(i,jde)=u1d(i,jde-1)
                v2d(i,jde)=v1d(i,jde-1)
              enddo
            endif
            if(w_bdy)then
              do j=jms,jme
                u2d(ids,j)=u1d(ids,j)
                v2d(ids,j)=v1d(ids,j)
              enddo
            endif
            if(e_bdy)then
              do j=jms,jme
                u2d(ide-1,j)=u1d(ide-1,j)
                v2d(ide-1,j)=v1d(ide-1,j)
                u2d(ide,j)=u1d(ide-1,j)
                v2d(ide,j)=v1d(ide-1,j)
              enddo
            endif
          endif
!
          do j=jts_h1,jte_h1
            do i=its_h1,ite_h1
              u1d(i,j)=u2d(i,j)
              v1d(i,j)=v2d(i,j)
            enddo
          enddo
        endif
!-----------------------------------------------------------------------
!---horizontal fluxes of momentum components on v points----------------
!-----------------------------------------------------------------------
        if(global) then
          if(s_bdy) then
            do i=its_b1,ite_b1_h1
              pfyx1 (i,jts_b1)=0.
              pfyy1 (i,jts_b1)=0.
              pfnex1(i,jts_b1)=0.
              pfnwx1(i,jts_b1)=0.
              pfney1(i,jts_b1)=0.
              pfnwy1(i,jts_b1)=0.
            enddo
          endif
!
          if(n_bdy) then
            do i=its_b1,ite_b1_h1
              pfyx1 (i,jte_b1_h1)=0.
              pfyy1 (i,jte_b1_h1)=0.
              pfnex1(i,jte_b1_h1)=0.
              pfnwx1(i,jte_b1_h1)=0.
              pfney1(i,jte_b1_h1)=0.
              pfnwy1(i,jte_b1_h1)=0.
            enddo
          endif
        endif
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1_h1
          do i=its_b1,ite_b1_h1
            ufxx1 (i,j)=(u1d(i  ,j  )-u1d(i-1,j  ))*pfxx1 (i,j)
            ufyx1 (i,j)=(u1d(i  ,j  )-u1d(i  ,j-1))*pfyx1 (i,j)
            ufnex1(i,j)=(u1d(i  ,j  )-u1d(i-1,j-1))*pfnex1(i,j)
            ufnwx1(i,j)=(u1d(i-1,j  )-u1d(i  ,j-1))*pfnwx1(i,j)
!
            ufxy1 (i,j)=(u1d(i  ,j  )-u1d(i-1,j  ))*pfxy1 (i,j)
            ufyy1 (i,j)=(u1d(i  ,j  )-u1d(i  ,j-1))*pfyy1 (i,j)
            ufney1(i,j)=(u1d(i  ,j  )-u1d(i-1,j-1))*pfney1(i,j)
            ufnwy1(i,j)=(u1d(i-1,j  )-u1d(i  ,j-1))*pfnwy1(i,j)
!
            vfxx1 (i,j)=(v1d(i  ,j  )-v1d(i-1,j  ))*pfxx1 (i,j)
            vfyx1 (i,j)=(v1d(i  ,j  )-v1d(i  ,j-1))*pfyx1 (i,j)
            vfnex1(i,j)=(v1d(i  ,j  )-v1d(i-1,j-1))*pfnex1(i,j)
            vfnwx1(i,j)=(v1d(i-1,j  )-v1d(i  ,j-1))*pfnwx1(i,j)
!
            vfxy1 (i,j)=(v1d(i  ,j  )-v1d(i-1,j  ))*pfxy1 (i,j)
            vfyy1 (i,j)=(v1d(i  ,j  )-v1d(i  ,j-1))*pfyy1 (i,j)
            vfney1(i,j)=(v1d(i  ,j  )-v1d(i-1,j-1))*pfney1(i,j)
            vfnwy1(i,j)=(v1d(i-1,j  )-v1d(i  ,j-1))*pfnwy1(i,j)
          enddo
        enddo
!-----------------------------------------------------------------------
!---advection of u1d and v1d--------------------------------------------
!-----------------------------------------------------------------------
        if(adv_standard) then
          ibeg=max(its,ids+1+lnsad)
          iend=min(ite,ide-2-lnsad)
          jbeg=max(jts,jds+1+lnsad)
          jend=min(jte,jde-2-lnsad)
!
          do j=jbeg,jend
            dxody=dxv(j)*rdyv
            dyodx=dyv   *rdxv(j)
            fadp=fad(j)
            do i=ibeg,iend
!
              fdpp=fadp/(dsg2(l)*pdop(i,j)+pdsg1(l))
!
              dux1= ufxx1 (i  ,j  )+ufxx1 (i+1,j  ) &
                   +ufyx1 (i  ,j  )+ufyx1 (i  ,j+1) &
                  +(ufnex1(i+1,j+1)+ufnex1(i  ,j  ) &
                   +ufnwx1(i  ,j+1)+ufnwx1(i+1,j  ))*0.25
!
              dvx1= vfxx1 (i  ,j  )+vfxx1 (i+1,j  ) &
                   +vfyx1 (i  ,j  )+vfyx1 (i  ,j+1) &
                  +(vfnex1(i+1,j+1)+vfnex1(i  ,j  ) &
                   +vfnwx1(i  ,j+1)+vfnwx1(i+1,j  ))*0.25
!
              duy1= ufxy1 (i  ,j  )+ufxy1 (i+1,j  ) &
                   +ufyy1 (i  ,j  )+ufyy1 (i  ,j+1) &
                  +(ufney1(i+1,j+1)+ufney1(i  ,j  ) &
                   +ufnwy1(i  ,j+1)+ufnwy1(i+1,j  ))*0.25
!
              dvy1= vfxy1 (i  ,j  )+vfxy1 (i+1,j  ) &
                   +vfyy1 (i  ,j  )+vfyy1 (i  ,j+1) &
                  +(vfney1(i+1,j+1)+vfney1(i  ,j  ) &
                   +vfnwy1(i  ,j+1)+vfnwy1(i+1,j  ))*0.25
!
              tcu(i,j,l)=((dvx1-dvy1)*dxody+(dux1+duy1))*fdpp &
                        +tcu(i,j,l)
              tcv(i,j,l)=((dux1-duy1)*dyodx+(dvx1+dvy1))*fdpp &
                        +tcv(i,j,l)
            enddo
          enddo
        endif
!-----------------------------------------------------------------------
!---regional branch-----------------------------------------------------
!-----------------------------------------------------------------------
!
        if(.not.global.and.adv_upstream) then
!
!-----------------------------------------------------------------------
!---upstream advection of momentum along lateral boundaries-------------
!-----------------------------------------------------------------------
!
          envp=-dt*rdyv
!
!***  Upstream advection along southern rows.
!
          jbeg=max(jts,jds+1)
          jend=min(jte,jds+lnsad)
!
          do j=jbeg,jend
            emvp=-dt*rdxv(j)
            do i=its_b1,ite_b2
              pp=u(i,j,l)*emvp
              qq=v(i,j,l)*envp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tcu(i,j,l)=(u(i+iap,j,l)-u(i,j,l))*pp &
                        +(u(i,j+jap,l)-u(i,j,l))*qq &
                        +(u(i,j,l)-u(i+iap,j,l) &
                         -u(i,j+jap,l)+u(i+iap,j+jap,l))*pp*qq &
                        +tcu(i,j,l)
              tcv(i,j,l)=(v(i+iap,j,l)-v(i,j,l))*pp &
                        +(v(i,j+jap,l)-v(i,j,l))*qq &
                        +(v(i,j,l)-v(i+iap,j,l) &
                         -v(i,j+jap,l)+v(i+iap,j+jap,l))*pp*qq &
                        +tcv(i,j,l)
            enddo
          enddo
!
!***  Upstream advection along northern rows.
!
          do j=max(jts,jde-1-lnsad),jte_b2
            emvp=-dt*rdxv(j)
            do i=its_b1,ite_b2
              pp=u(i,j,l)*emvp
              qq=v(i,j,l)*envp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tcu(i,j,l)=(u(i+iap,j,l)-u(i,j,l))*pp &
                        +(u(i,j+jap,l)-u(i,j,l))*qq &
                        +(u(i,j,l)-u(i+iap,j,l) &
                         -u(i,j+jap,l)+u(i+iap,j+jap,l))*pp*qq &
                        +tcu(i,j,l)
              tcv(i,j,l)=(v(i+iap,j,l)-v(i,j,l))*pp &
                        +(v(i,j+jap,l)-v(i,j,l))*qq &
                        +(v(i,j,l)-v(i+iap,j,l) &
                         -v(i,j+jap,l)+v(i+iap,j+jap,l))*pp*qq &
                        +tcv(i,j,l)
            enddo
          enddo
!
!***  Upstream advection along western rows.
!
          do j=max(jts,jds+1+lnsad),min(jte,jde-2-lnsad)
            emvp=-dt*rdxv(j)
            do i=its_b1,ids+lnsad
              pp=u(i,j,l)*emvp
              qq=v(i,j,l)*envp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tcu(i,j,l)=(u(i+iap,j,l)-u(i,j,l))*pp &
                        +(u(i,j+jap,l)-u(i,j,l))*qq &
                        +(u(i,j,l)-u(i+iap,j,l) &
                         -u(i,j+jap,l)+u(i+iap,j+jap,l))*pp*qq &
                        +tcu(i,j,l)
              tcv(i,j,l)=(v(i+iap,j,l)-v(i,j,l))*pp &
                        +(v(i,j+jap,l)-v(i,j,l))*qq &
                        +(v(i,j,l)-v(i+iap,j,l) &
                         -v(i,j+jap,l)+v(i+iap,j+jap,l))*pp*qq &
                        +tcv(i,j,l)
            enddo
          enddo
!
!***  Upstream advection along eastern rows.
!
          do j=max(jts,jds+1+lnsad),min(jte,jde-2-lnsad)
            emvp=-dt*rdxv(j)
            do i=max(its,ide-1-lnsad),ite_b2
              pp=u(i,j,l)*emvp
              qq=v(i,j,l)*envp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tcu(i,j,l)=(u(i+iap,j,l)-u(i,j,l))*pp &
                        +(u(i,j+jap,l)-u(i,j,l))*qq &
                        +(u(i,j,l)-u(i+iap,j,l) &
                         -u(i,j+jap,l)+u(i+iap,j+jap,l))*pp*qq &
                        +tcu(i,j,l)
              tcv(i,j,l)=(v(i+iap,j,l)-v(i,j,l))*pp &
                        +(v(i,j+jap,l)-v(i,j,l))*qq &
                        +(v(i,j,l)-v(i+iap,j,l) &
                         -v(i,j+jap,l)+v(i+iap,j+jap,l))*pp*qq &
                        +tcv(i,j,l)
            enddo
          enddo
!-----------------------------------------------------------------------
        endif ! regional lateral boundaries
!-----------------------------------------------------------------------
!
      enddo vertical_loop2
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
                        endsubroutine adv1
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine vtoa &
(lm &
,dt,ef4t,pt &
,sg2 &
,psdt &
,dwdt,rtop &
,omgalf &
,pint &
!---temporary arguments-------------------------------------------------
,tdiv,tct)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm                          ! total # of levels

real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,ef4t &                      ! vertical grid parameter
,pt                          ! pressure at the top of the model's atmosphere

real(kind=kfpt),dimension(1:lm+1),intent(in):: &
 sg2                         ! delta sigmas

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 psdt                        ! surface pressure tendency

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 dwdt &                      ! nonhydrostatic correction factor
,rtop                        ! RT/p

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 omgalf                      !

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(inout):: &
 pint                        ! pressure at interfaces
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 tdiv                        ! integrated horizontal mass divergence

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 tct                         ! time change of temperature
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction

real(kind=kfpt):: &
 dwdtp &                     ! nonhydrostatic correction factor at the point
,toa &                       ! omega-alpha temporary
,tpmp                        ! pressure temporary at the point

real(kind=kfpt),dimension(its:ite,jts:jte):: &
 tpm                         ! pressure temporary
!-----------------
!-----------------
integer(kind=kint) :: &
 jstart &
,jstop &
,nth &
,omp_get_num_threads &
,omp_get_thread_num &
,tid
!-----------------
!-----------------
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!      do l=1,lm !!!debugggggggg!!!!!!!
!        do j=jds,jde
!          do i=ids,ide
!            tct(i,j,l)=0.
!          enddo
!        enddo
!      enddo
!-----------------------------------------------------------------------
!-----------------
!-----------------
!.......................................................................
!$omp parallel private(nth, tid, i,j,l,jstart,jstop,dwdtp,toa)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth, jts_b1,jte_b1,jstart,jstop)
!-----------------
!-----------------
      do j=jstart,jstop
        do i=its_b1,ite_b1
          pint(i,j,1)=pt
          tpm(i,j)=pt+pint(i,j,2)
!
          dwdtp=dwdt(i,j,1)
!
          tpmp=pint(i,j,2)+pint(i,j,3)
          toa=-tdiv(i,j,1)*rtop(i,j,1)*dwdtp*ef4t
!
          omgalf(i,j,1)=omgalf(i,j,1)+toa
          tct(i,j,1)=tct(i,j,1)+toa
!
          pint(i,j,2)=psdt(i,j)*(sg2(1)+sg2(2))*dwdtp*dt &
                     +tpm(i,j)-pint(i,j,1)
!
          tpm(i,j)=tpmp
        enddo
      enddo
!
      do l=2,lm-1
        do j=jstart,jstop
          do i=its_b1,ite_b1
            dwdtp=dwdt(i,j,l)
!
            tpmp=pint(i,j,l+1)+pint(i,j,l+2)
            toa=-(tdiv(i,j,l-1)+tdiv(i,j,l))*rtop(i,j,l)*dwdtp*ef4t
!
            omgalf(i,j,l)=omgalf(i,j,l)+toa 
            tct(i,j,l)=tct(i,j,l)+toa
!
            pint(i,j,l+1)=psdt(i,j)*(sg2(l)+sg2(l+1))*dwdtp*dt &
                         +tpm(i,j)-pint(i,j,l)
!
            tpm(i,j)=tpmp
          enddo
        enddo
      enddo
!
      do j=jstart,jstop
        do i=its_b1,ite_b1
          dwdtp=dwdt(i,j,lm)
!
          toa=-(tdiv(i,j,lm-1)+tdiv(i,j,lm))*rtop(i,j,lm)*dwdtp*ef4t
!
          omgalf(i,j,lm)=omgalf(i,j,lm)+toa
          tct(i,j,lm)=tct(i,j,lm)+toa
!
          pint(i,j,lm+1)=psdt(i,j)*(sg2(lm)+sg2(lm+1))*dwdtp*dt &
                        +tpm(i,j)-pint(i,j,lm)
        enddo
      enddo
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------

!-----------------------------------------------------------------------
!
                        endsubroutine vtoa
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine updates &
(lm,s &
!---temporary arguments-------------------------------------------------
,tcs)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm                          ! total # of levels

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 s                           ! tracer
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 tcs                         ! tracer time change
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do private(l,j,i)
!.......................................................................
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            s(i,j,l)=s(i,j,l)+tcs(i,j,l)
!
            tcs(i,j,l)=0.
          enddo
        enddo
      enddo
!.......................................................................
!$omp end parallel do
!.......................................................................
!-----------------------------------------------------------------------
!
                        endsubroutine updates
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine updatet &
(lm,t &
!---temporary arguments-------------------------------------------------
,tct)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm                          ! total # of levels

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 t                           ! temperature
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 tct                         ! temperature time change
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do private(l,j,i)
!.......................................................................
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            t(i,j,l)=t(i,j,l)+tct(i,j,l)
!
            tct(i,j,l)=0.
          enddo
        enddo
      enddo
!.......................................................................
!$omp end parallel do
!.......................................................................
!-----------------------------------------------------------------------
!
                        endsubroutine updatet
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                            
!-----------------------------------------------------------------------
                        subroutine updateuv &
(lm,u,v &
!---temporary arguments-------------------------------------------------
,tcu,tcv)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 cflfc=1./(160.*160.)        ! cfl limit
           
integer(kind=kint),intent(in):: &
 lm                          ! total # of levels
           
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 u &                         ! u
,v                           ! v
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 tcu &                       ! u time change
,tcv                         ! v time change
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction
           
real(kind=kfpt):: &
 cflc &                      !
,rcflc                       !
           
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do private (l,j,i,cflc,rcflc)
!.......................................................................
      do l=1,lm
        do j=jts_b1,jte_b2
          do i=its_b1,ite_b2
            u(i,j,l)=u(i,j,l)+tcu(i,j,l)
            v(i,j,l)=v(i,j,l)+tcv(i,j,l)
!
            tcu(i,j,l)=0.
            tcv(i,j,l)=0.
!
            cflc=(u(i,j,l)*u(i,j,l)+v(i,j,l)*v(i,j,l))*cflfc
            if(cflc.gt.1.) then
              rcflc=sqrt(1./cflc)
              u(i,j,l)=u(i,j,l)*rcflc
              v(i,j,l)=v(i,j,l)*rcflc
            endif
          enddo
        enddo
      enddo
!.......................................................................
!$omp end parallel do
!.......................................................................
!-----------------------------------------------------------------------
!
                        endsubroutine updateuv
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                            
!-----------------------------------------------------------------------
                        subroutine hdiff &
(global,hydro,secdif &
,inpes,jnpes,lm,lpt2 &
,dyh,rdyh &
,dxv,rare,rdxh &
,sice,sm &
,hdacx,hdacy,hdacvx,hdacvy &
,w,z &
,cw,q,q2,t,u,v)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 scq2=50. &                  ! 2tke weighting factor
,defm=1.35e-3/2. &           ! deformation cap, /4. for smag2=0.4
,epsq=1.e-20 &               ! floor value for specific humidity
,epsq2=0.02 &                ! floor value for 2tke
,slopec=.05                  ! critical slope
           
!-----------------------------------------------------------------------
logical(kind=klog),intent(in):: &
 global &                    ! global or regional
,hydro &                     ! logical switch for nonhydrostatic dynamics
,secdif                      ! 2nd order diffusion
           
integer(kind=kint),intent(in):: &
 inpes &                     ! w-e # of subdomains
,jnpes &                     ! n-s # of subdomains
,lm &                        ! total # of levels
,lpt2                        ! # of levels in the pressure range
           
real(kind=kfpt),intent(in):: &
 dyh &                       !
,rdyh                        ! 1/deltay
 
real(kind=kfpt),dimension(jds:jde),intent(in):: &
 dxv &                       !
,rare &                      !
,rdxh                        ! 1/deltax
           
real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 hdacx &                     ! exchange coefficient for mass points
,hdacy &                     ! exchange coefficient for mass points
,hdacvx &                    ! exchange coefficient for velocity points
,hdacvy &                    ! exchange coefficient for velocity points
,sice &                      ! sea-ice mask
,sm                          ! sea mask
           
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 w &                         ! w wind component
,z                           ! height at mass points
           
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 cw &                        ! condensate
,q &                         ! specific humidity
,q2 &                        ! 2tke
,t &                         ! temperature
,u &                         ! u wind component
,v                           ! v wind component
           
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
logical(kind=klog):: &
 cilinx &                    ! coast/ice line
,ciliny                      ! coast/ice line
 
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction
           
real(kind=kfpt):: &
 defc &                      ! deformation floor
,def1 &                      ! component of deformation
,def2 &                      ! component of deformation
,def3 &                      ! component of deformation
,def4 &                      ! component of deformation
,defp &                      ! deformation at the point
,defs &                      ! component of deformation
,deft &                      ! component of deformation
!,defz &                      ! rotational component of deformation
,hkfx &                      ! def with slope factor
,hkfy &                      ! def with slope factor
,q2trm &                     !
,slopx &                     ! x slope
,slopy                       ! y slope

real(kind=kfpt),dimension(ims:ime,jms:jme,lm):: def3d          
real(kind=kfpt),dimension(ims:ime,jms:jme):: &
 cdif &                      ! condensate 2nd order diffusion
,cx &                        ! condensate difference, x direction
,cy &                        ! condensate difference, y direction
,def &                       ! deformation (in halo exchange)
,fmlx &                      ! x slope mask
,fmly &                      ! y slope mask
,hkx &                       ! deformation sum, x direction
,hky &                       ! deformation sum, y direction
,qdif &                      ! specific humidity 2nd order diffusion
,qx &                        ! specific humidity difference, x direction
,qy &                        ! specific humidity difference, y direction
,q2dif &                     ! 2tke 2nd order diffusion
,q2x &                       ! 2tke difference, x direction
,q2y &                       ! 2tke difference, y direction
,tdif &                      ! temperature 2nd order diffusion
,tx &                        ! temperature difference, x direction
,ty &                        ! temperature difference, y direction
,udif &                      ! u wind component 2nd order diffusion
,ux &                        ! u wind component difference, x direction
,uy &                        ! u wind component difference, y direction
,vdif &                      ! v wind component 2nd order diffusion
,vx &                        ! v wind component difference, x direction
,vy                          ! v wind component difference, y direction
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
     do l=1,lm
      if(s_bdy)then
        do i=ims,ime
          def3d(i,jds,l)=0.
        enddo
      endif
!
      if(n_bdy)then
        do i=ims,ime
          def3d(i,jde,l)=0.
        enddo
      endif
!
      if(w_bdy)then
        do j=jms,jme
          def3d(ids,j,l)=0.
        enddo
      endif
!
      if(e_bdy)then
        do j=jms,jme
          def3d(ide,j,l)=0.
        enddo
      endif
    enddo

!
!-----------------------------------------------------------------------
!---grand vertical loop-------------------------------------------------
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do &
!$omp private(l,defc,j,i,deft,defs,def1,def2,def3,def4,q2trm,defp)
!.......................................................................
!-----------------------------------------------------------------------
!
      vertical_loop_1: do l=1,lm
!
!-----------------------------------------------------------------------
        if(l.gt.lpt2/3) then
          defc=0.
        else
          defc=defm*0.01
        endif
!-----------------------------------------------------------------------
!        do j=jts_h1,jte_h2
!          do i=its_h1,ite_h2
!            q2(i,j,l)=max(q2(i,j,l),epsq2)
!          enddo
!        enddo
!-----------------------------------------------------------------------
        do j=jts_b1_h1,jte_b1_h2
          do i=its_b1_h1,ite_b1_h2
            deft=((u(i  ,j-1,l)+u(i  ,j  ,l) &
                  -u(i-1,j-1,l)-u(i-1,j  ,l))*dyh &
                 -(v(i-1,j  ,l)+v(i  ,j  ,l))*dxv(j  ) &
                 +(v(i-1,j-1,l)+v(i  ,j-1,l))*dxv(j-1))*rare(j)
            defs=((u(i-1,j  ,l)+u(i  ,j  ,l))*dxv(j  ) &
                 -(u(i-1,j-1,l)+u(i  ,j-1,l))*dxv(j-1) &
                 +(v(i  ,j-1,l)+v(i  ,j  ,l) &
                  -v(i-1,j-1,l)-v(i-1,j  ,l))*dyh     )*rare(j)
!            defz=(-(u(i-1,j  ,l)+u(i  ,j  ,l))*dxv(j  ) &
!                  +(u(i-1,j-1,l)+u(i  ,j-1,l))*dxv(j-1) &
!                  +(v(i  ,j-1,l)+v(i  ,j  ,l) &
!                   -v(i-1,j-1,l)-v(i-1,j  ,l))*dyh     )*rare(j)  *10.
!            if(defz.gt.0.) defz=0.
!
            def1=(w(i,j,l)-w(i,j-1,l))*rdyh
            def2=(w(i,j,l)-w(i-1,j,l))*rdxh(j)
            def3=(w(i+1,j,l)-w(i,j,l))*rdxh(j)
            def4=(w(i,j+1,l)-w(i,j,l))*rdyh
!
            if(q2(i,j,l).gt.epsq2) then
              q2trm=scq2*q2(i,j,l)*rare(j)
            else
              q2trm=0.
            endif
!
            defp=deft*deft+defs*defs &
!                +defz*defz &
                +def1*def1+def2*def2+def3*def3+def4*def4 &
                +q2trm
            defp=sqrt(defp+defp)
            defp=max (defp,defc)
            defp=min (defp,defm)
!
            def3d(i,j,l)=defp
          enddo
        enddo
!-----------------------------------------------------------------------
!
      enddo vertical_loop_1
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
      if(global) then
        btim=timef()
        call swaphn(def3d,ims,ime,jms,jme,lm,inpes)
        swaphn_tim=swaphn_tim+timef()-btim
!
        btim=timef()
        call polehn(def3d,ims,ime,jms,jme,lm,inpes,jnpes)
        polehn_tim=polehn_tim+timef()-btim
      endif
!
!     call halo_exch(def,1,1,1)
!-----------------------------------------------------------------------
!
      if(secdif) then ! 2nd order diffusion
!
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do &
!$omp private(l,j,i,cilinx,ciliny,slopx,slopy,fmlx,fmly,hkx,hky,hkfx,hkfy, &
!$omp         tx,qx,cx,q2x,ty,qy,cy,q2y,ux,vx,uy,vy,tdif,qdif,cdif,q2dif,  &
!$omp         udif,vdif)
!.......................................................................
!-----------------------------------------------------------------------
      vertical_loop_3: do l=1,lm
!-----------------------------------------------------------------------
        if(l.gt.lpt2.and..not.hydro) then
          do j=jts_b1,jte_h2
            do i=its_b1,ite_h2
              cilinx=sice(i-1,j).ne.sice(i,j) &
                 .or.sm  (i-1,j).ne.sm  (i,j)
              ciliny=sice(i,j-1).ne.sice(i,j) &
                 .or.sm  (i,j-1).ne.sm  (i,j)
              slopx=abs((z(i,j,l)-z(i-1,j,l))*rdxh(j))
              slopy=abs((z(i,j,l)-z(i,j-1,l))*rdyh   )
!
              if(slopx.le.slopec.or.cilinx) then
                fmlx(i,j)=1.
              else
                fmlx(i,j)=0.
              endif
!
              if(slopy.le.slopec.or.ciliny) then
                fmly(i,j)=1.
              else
                fmly(i,j)=0.
              endif
!
            enddo
          enddo
        else
          do j=jts_b1,jte_h2
            do i=its_b1,ite_h2
              fmlx(i,j)=1.
              fmly(i,j)=1.
            enddo
          enddo
        endif
!-----------------------------------------------------------------------
!---contributions behind mass points------------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_h2
          do i=its_b1,ite_h2
            hkx(i,j)=(def3d(i-1,j,l)+def3d(i,j,l))
            hky(i,j)=(def3d(i,j-1,l)+def3d(i,j,l))
            hkfx=hkx(i,j)*fmlx(i,j)
            hkfy=hky(i,j)*fmly(i,j)
!
            tx (i,j)=(t (i,j,l)-t (i-1,j,l))*hkfx
            qx (i,j)=(q (i,j,l)-q (i-1,j,l))*hkfx
            cx (i,j)=(cw(i,j,l)-cw(i-1,j,l))*hkfx
            q2x(i,j)=(q2(i,j,l)-q2(i-1,j,l))*hkfx
!
            ty (i,j)=(t (i,j,l)-t (i,j-1,l))*hkfy
            qy (i,j)=(q (i,j,l)-q (i,j-1,l))*hkfy
            cy (i,j)=(cw(i,j,l)-cw(i,j-1,l))*hkfy
            q2y(i,j)=(q2(i,j,l)-q2(i,j-1,l))*hkfy
          enddo
        enddo
!-----------------------------------------------------------------------
!---u,v, contributions, behind v points---------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1_h1
          do i=its_b1,ite_b1_h1
            ux(i,j)=(u(i,j,l)-u(i-1,j,l))*hky(i,j+1)
            vx(i,j)=(v(i,j,l)-v(i-1,j,l))*hky(i,j+1)
            uy(i,j)=(u(i,j,l)-u(i,j-1,l))*hkx(i+1,j)
            vy(i,j)=(v(i,j,l)-v(i,j-1,l))*hkx(i+1,j)
          enddo
        enddo
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            tdif (i,j)=(tx (i+1,j)-tx (i,j))*hdacx(i,j) &
                      +(ty (i,j+1)-ty (i,j))*hdacy(i,j)
            qdif (i,j)=(qx (i+1,j)-qx (i,j))*hdacx(i,j) &
                      +(qy (i,j+1)-qy (i,j))*hdacy(i,j)
            cdif (i,j)=(cx (i+1,j)-cx (i,j))*hdacx(i,j) &
                      +(cy (i,j+1)-cy (i,j))*hdacy(i,j)
            q2dif(i,j)=(q2x(i+1,j)-q2x(i,j))*hdacx(i,j) &
                      +(q2y(i,j+1)-q2y(i,j))*hdacy(i,j)
          enddo
        enddo
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b2
          do i=its_b1,ite_b2
            udif(i,j)=(ux(i+1,j)-ux(i,j))*hdacvx(i,j) &
                     +(uy(i,j+1)-uy(i,j))*hdacvy(i,j)
            vdif(i,j)=(vx(i+1,j)-vx(i,j))*hdacvx(i,j) &
                     +(vy(i,j+1)-vy(i,j))*hdacvy(i,j)
          enddo
        enddo
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-------------2-nd order diffusion--------------------------------------
!-----------------------------------------------------------------------
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              t (i,j,l)=t (i,j,l)+tdif (i,j)
              q (i,j,l)=q (i,j,l)+qdif (i,j)
              cw(i,j,l)=cw(i,j,l)+cdif (i,j)
              q2(i,j,l)=q2(i,j,l)+q2dif(i,j)
            enddo
          enddo
!-----------------------------------------------------------------------
!-------------2-nd order diffusion--------------------------------------
!-----------------------------------------------------------------------
          do j=jts_b1,jte_b2
            do i=its_b1,ite_b2
              u(i,j,l)=u(i,j,l)+udif(i,j)
              v(i,j,l)=v(i,j,l)+vdif(i,j)
            enddo
          enddo
!-----------------------------------------------------------------------
!
        enddo vertical_loop_3
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
        else ! 4th order diffusion
!-----------------------------------------------------------------------
!-------------4-th order diffusion--------------------------------------
!-----------------------------------------------------------------------
      vertical_loop_4: do l=1,lm
!-----------------------------------------------------------------------
        if(l.gt.lpt2.and..not.hydro) then
          do j=jts_b1,jte_h2
            do i=its_b1,ite_h2
              cilinx=sice(i-1,j).ne.sice(i,j) &
                 .or.sm  (i-1,j).ne.sm  (i,j)
              ciliny=sice(i,j-1).ne.sice(i,j) &
                 .or.sm  (i,j-1).ne.sm  (i,j)
              slopx=abs((z(i,j,l)-z(i-1,j,l))*rdxh(j))
              slopy=abs((z(i,j,l)-z(i,j-1,l))*rdyh   )
!
              if(slopx.le.slopec.or.cilinx) then
                fmlx(i,j)=1.
              else
                fmlx(i,j)=0.
              endif
!
              if(slopy.le.slopec.or.ciliny) then
                fmly(i,j)=1.
              else
                fmly(i,j)=0.
              endif
!
            enddo
          enddo
        else
          do j=jts_b1,jte_h2
            do i=its_b1,ite_h2
              fmlx(i,j)=1.
              fmly(i,j)=1.
            enddo
          enddo
        endif
!-----------------------------------------------------------------------
!---contributions behind mass points------------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_h2
          do i=its_b1,ite_h2
            hkx(i,j)=(def3d(i-1,j,l)+def3d(i,j,l))
            hky(i,j)=(def3d(i,j-1,l)+def3d(i,j,l))
            hkfx=hkx(i,j)*fmlx(i,j)
            hkfy=hky(i,j)*fmly(i,j)
!
            tx (i,j)=(t (i,j,l)-t (i-1,j,l))*hkfx
            qx (i,j)=(q (i,j,l)-q (i-1,j,l))*hkfx
            cx (i,j)=(cw(i,j,l)-cw(i-1,j,l))*hkfx
            q2x(i,j)=(q2(i,j,l)-q2(i-1,j,l))*hkfx
!
            ty (i,j)=(t (i,j,l)-t (i,j-1,l))*hkfy
            qy (i,j)=(q (i,j,l)-q (i,j-1,l))*hkfy
            cy (i,j)=(cw(i,j,l)-cw(i,j-1,l))*hkfy
            q2y(i,j)=(q2(i,j,l)-q2(i,j-1,l))*hkfy
          enddo
        enddo
!-----------------------------------------------------------------------
!---u,v, contributions, behind v points---------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1_h1
          do i=its_b1,ite_b1_h1
            ux(i,j)=(u(i,j,l)-u(i-1,j,l))*hky(i,j+1)
            vx(i,j)=(v(i,j,l)-v(i-1,j,l))*hky(i,j+1)
            uy(i,j)=(u(i,j,l)-u(i,j-1,l))*hkx(i+1,j)
            vy(i,j)=(v(i,j,l)-v(i,j-1,l))*hkx(i+1,j)
          enddo
        enddo
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            tdif (i,j)=(tx (i+1,j)-tx (i,j))*hdacx(i,j) &
                      +(ty (i,j+1)-ty (i,j))*hdacy(i,j)
            qdif (i,j)=(qx (i+1,j)-qx (i,j))*hdacx(i,j) &
                      +(qy (i,j+1)-qy (i,j))*hdacy(i,j)
            cdif (i,j)=(cx (i+1,j)-cx (i,j))*hdacx(i,j) &
                      +(cy (i,j+1)-cy (i,j))*hdacy(i,j)
            q2dif(i,j)=(q2x(i+1,j)-q2x(i,j))*hdacx(i,j) &
                      +(q2y(i,j+1)-q2y(i,j))*hdacy(i,j)
          enddo
        enddo
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b2
          do i=its_b1,ite_b2
            udif(i,j)=(ux(i+1,j)-ux(i,j))*hdacvx(i,j) &
                     +(uy(i,j+1)-uy(i,j))*hdacvy(i,j)
            vdif(i,j)=(vx(i+1,j)-vx(i,j))*hdacvx(i,j) &
                     +(vy(i,j+1)-vy(i,j))*hdacvy(i,j)
          enddo
        enddo
          if(global) then
            btim=timef()
            call swaphn(tdif,ims,ime,jms,jme,1,inpes)
            call swaphn(qdif,ims,ime,jms,jme,1,inpes)
            call swaphn(cdif,ims,ime,jms,jme,1,inpes)
            call swaphn(q2dif,ims,ime,jms,jme,1,inpes)
            swaphn_tim=swaphn_tim+timef()-btim
!
            btim=timef()
            call polehn(tdif,ims,ime,jms,jme,1,inpes,jnpes)
            call polehn(qdif,ims,ime,jms,jme,1,inpes,jnpes)
            call polehn(cdif,ims,ime,jms,jme,1,inpes,jnpes)
            call polehn(q2dif,ims,ime,jms,jme,1,inpes,jnpes)
            polehn_tim=polehn_tim+timef()-btim
!
            btim=timef()
            call swapwn(udif,ims,ime,jms,jme,1,inpes)
            call swapwn(vdif,ims,ime,jms,jme,1,inpes)
            swapwn_tim=swapwn_tim+timef()-btim
!
            btim=timef()
            call polewn(udif,vdif,ims,ime,jms,jme,1,inpes,jnpes)
            polewn_tim=polewn_tim+timef()-btim
          else
            if(s_bdy)then
              do i=ims,ime
                tdif(i,jds)=0.
                qdif(i,jds)=0.
                cdif(i,jds)=0.
                q2dif(i,jds)=0.
                udif(i,jds)=0.
                vdif(i,jds)=0.
              enddo
            endif
            if(n_bdy)then
              do i=ims,ime
                tdif(i,jde)=0.
                qdif(i,jde)=0.
                cdif(i,jde)=0.
                q2dif(i,jde)=0.
                udif(i,jde-1)=0.
                vdif(i,jde-1)=0.
                udif(i,jde)=0.
                vdif(i,jde)=0.
              enddo
            endif
            if(w_bdy)then
              do j=jms,jme
                tdif(ids,j)=0.
                qdif(ids,j)=0.
                cdif(ids,j)=0.
                q2dif(ids,j)=0.
                udif(ids,j)=0.
                vdif(ids,j)=0.
              enddo
            endif
            if(e_bdy)then
              do j=jms,jme
                tdif(ide,j)=0.
                qdif(ide,j)=0.
                cdif(ide,j)=0.
                q2dif(ide,j)=0.
                udif(ide-1,j)=0.
                vdif(ide-1,j)=0.
                udif(ide,j)=0.
                vdif(ide,j)=0.
              enddo
            endif
          endif
!
          btim=timef()
          call halo_exch( tdif,1,2,2)
          call halo_exch( qdif,1,2,2)
          call halo_exch( cdif,1,2,2)
          call halo_exch(q2dif,1,2,2)
          call halo_exch( udif,1,2,2)
          call halo_exch( vdif,1,2,2)
          exch_dyn_tim=exch_dyn_tim+timef()-btim
!---contributions behind mass points------------------------------------
          do j=jts_b1,jte_h1
            do i=its_b1,ite_h1
              hkfx=hkx(i,j)*fmlx(i,j)
              hkfy=hky(i,j)*fmly(i,j)
!
              tx (i,j)=(tdif (i,j)-tdif (i-1,j))*hkfx
              qx (i,j)=(qdif (i,j)-qdif (i-1,j))*hkfx
              cx (i,j)=(cdif (i,j)-cdif (i-1,j))*hkfx
              q2x(i,j)=(q2dif(i,j)-q2dif(i-1,j))*hkfx
!
              ty (i,j)=(tdif (i,j)-tdif (i,j-1))*hkfy
              qy (i,j)=(qdif (i,j)-qdif (i,j-1))*hkfy
              cy (i,j)=(cdif (i,j)-cdif (i,j-1))*hkfy
              q2y(i,j)=(q2dif(i,j)-q2dif(i,j-1))*hkfy
            enddo
          enddo
!---u,v, contributions, behind v points---------------------------------
          do j=jts_b1,jte_b1_h1
            do i=its_b1,ite_b1_h1
              ux(i,j)=(udif(i,j)-udif(i-1,j))*hky(i,j+1)
              vx(i,j)=(vdif(i,j)-vdif(i-1,j))*hky(i,j+1)
              uy(i,j)=(udif(i,j)-udif(i,j-1))*hkx(i+1,j)
              vy(i,j)=(vdif(i,j)-vdif(i,j-1))*hkx(i+1,j)
            enddo
          enddo
!-----------------------------------------------------------------------
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              t (i,j,l)=-(tx (i+1,j)-tx (i,j))*hdacx(i,j) &
                        -(ty (i,j+1)-ty (i,j))*hdacy(i,j) &
                       +t (i,j,l)
              q (i,j,l)=-(qx (i+1,j)-qx (i,j))*hdacx(i,j) &
                        -(qy (i,j+1)-qy (i,j))*hdacy(i,j) &
                       +q (i,j,l)
              cw(i,j,l)=-(cx (i+1,j)-cx (i,j))*hdacx(i,j) &
                        -(cy (i,j+1)-cy (i,j))*hdacy(i,j) &
                       +cw(i,j,l)
              q2(i,j,l)=-(q2x(i+1,j)-q2x(i,j))*hdacx(i,j) &
                        -(q2y(i,j+1)-q2y(i,j))*hdacy(i,j) &
                       +q2(i,j,l)
            enddo
          enddo
!-----------------------------------------------------------------------
          do j=jts_b1,jte_b2
            do i=its_b1,ite_b2
              u(i,j,l)=-(ux(i+1,j)-ux(i,j))*hdacx(i,j) &
                       -(uy(i,j+1)-uy(i,j))*hdacy(i,j) &
                      +u(i,j,l)
              v(i,j,l)=-(vx(i+1,j)-vx(i,j))*hdacx(i,j) &
                       -(vy(i,j+1)-vy(i,j))*hdacy(i,j) &
                      +v(i,j,l)
            enddo
          enddo
!-----------------------------------------------------------------------
      enddo vertical_loop_4
!-----------------------------------------------------------------------
     endif
!-----------------------------------------------------------------------
!
                        endsubroutine hdiff
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine cdzdt &
(global,hydro &
,lm &
,dt &
,dsg2,pdsg1 &
,fah &
,fis,pd,pdo &
,psgdt &
,cw,q,rtop,t &
,pint &
,dwdt,pdwdt,w &
,z &
!---temporary arguments-------------------------------------------------
,pfne,pfnw,pfx,pfy)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 cfc=1.533 &                 ! adams-bashforth positioning in time
,bfc=1.-cfc                  ! adams bashforth positioning in time
!-----------------------------------------------------------------------
logical(kind=klog),intent(in):: &
 global &                    ! global or regional
,hydro                       ! hydrostatic or nonhydrostatic
           
integer(kind=kint),intent(in):: &
 lm                          ! total # of levels
           
real(kind=kfpt),intent(in):: &
 dt                          ! dynamics time step
           
real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures
           
real(kind=kfpt),dimension(jds:jde),intent(in):: &
 fah                         !
           
real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 fis &                       ! surface geopotential
,pd &                        ! sigma range pressure difference
,pdo                         ! sigma range pressure difference
      
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(in):: &
 psgdt                       ! vertical mass flux
           
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 cw &                        ! condensate
,q &                         ! specific humidity
,rtop &                      ! rt/p
,t                           ! temperature
           
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(in):: &
 pint                        ! pressure at interfaces
           
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 dwdt &                      ! nonhydrostatic correction factor
,pdwdt &                     ! previous nonhydrostatic correction factor
,w                           ! w wind component
           
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out):: &
 z                           ! heights in the middle of the layers
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in) :: &
 pfne &                      ! mass flux, ne direction
,pfnw &                      ! mass flux, nw direction
,pfx &                       ! mass flux, x direction
,pfy                         ! mass flux, y direction
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction
           
real(kind=kfpt):: &
 dpup &                      !
,dw &                        ! w wind component increment
,dz &                        ! layer thickness
,fahp &                      !
,rg &                        !
,rdt &                       !
,trog &                      !
,wup &                       ! w wind component at upper interface
,zup                         ! height of upper interface
           
real(kind=kfpt),dimension(its_b1:ite_h1,jts_b1:jte_h1):: &
 zne &                       ! height flux, ne direction
,znw &                       ! height flux, nw direction
,zx &                        ! height flux, x direction
,zy                          ! height flux, y direction
           
real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1):: &
 tta &                       ! advection through upper interface
,ttb                         ! advection through lower interface
           
real(kind=kfpt),dimension(its_h1:ite_h1,jts_h1:jte_h1):: &
 wlo &                       ! w wind component at lower interface
,zlo                         ! height at lower interface
!-----------------
!-----------------
integer(kind=kint) :: &
 jstart &
,jstop &
,nth &
,omp_get_num_threads &
,omp_get_thread_num &
,tid
!-----------------
!-----------------
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rg=1./g
      rdt=1./dt
      trog=2.*r/g
!
!-----------------------------------------------------------------------
!-----------------
!-----------------
!.......................................................................
!$omp parallel private(nth, tid, i,j,l,jstart,jstop,dz,dw,zup,wup)
        nth = omp_get_num_threads()
        tid = omp_get_thread_num()
        call looplimits(tid, nth, jts_h1,jte_h1,jstart,jstop)
!.......................................................................
!-----------------
!-----------------
      do j=jstart,jstop
        do i=its_h1,ite_h1
          wlo(i,j)=0.
          zlo(i,j)=fis(i,j)*rg
        enddo
      enddo
!-----------------------------------------------------------------------
!---nonhydrostatic equation---------------------------------------------
!-----------------------------------------------------------------------
      do l=lm,1,-1
        do j=jstart,jstop
          do i=its_h1,ite_h1
            pdwdt(i,j,l)=dwdt(i,j,l)
            dwdt(i,j,l)=w(i,j,l)
!
            dz=(q(i,j,l)*0.608+(1.-cw(i,j,l)))*t(i,j,l)*trog &
              *(dsg2(l)*pd(i,j)+pdsg1(l))/(pint(i,j,l)+pint(i,j,l+1))
            dw=(dz-rtop(i,j,l)*(dsg2(l)*pdo(i,j)+pdsg1(l))*rg)*rdt
!
            zup=zlo(i,j)+dz
            wup=wlo(i,j)+dw
            z(i,j,l)=dz*0.5+zlo(i,j)
            w(i,j,l)=dw*0.5+wlo(i,j)
!
            zlo(i,j)=zup
            wlo(i,j)=wup
          enddo
        enddo
      enddo
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------
!-----------------------------------------------------------------------
!
      if(hydro) return
!
!-----------------------------------------------------------------------
!-----------------
!-----------------
!.......................................................................
!$omp parallel private(nth, tid, i,j,l,jstart,jstop)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth, jts_b1,jte_b1,jstart,jstop)
!-----------------
!-----------------
!
      do j=jstart,jstop
        do i=its_b1,ite_b1
          ttb(i,j)=0.
        enddo
      enddo
!
      do l=1,lm-1
        do j=jstart,jstop
          do i=its_b1,ite_b1
            tta(i,j)=(z(i,j,l+1)-z(i,j,l))*psgdt(i,j,l)*0.5
            w(i,j,l)=(tta(i,j)+ttb(i,j))/(dsg2(l)*pdo(i,j)+pdsg1(l)) &
                    +w(i,j,l)
            ttb(i,j)=tta(i,j)
          enddo
        enddo
      enddo
!
      do j=jstart,jstop
        do i=its_b1,ite_b1
          w(i,j,lm)=ttb(i,j)/(dsg2(lm)*pdo(i,j)+pdsg1(lm))+w(i,j,lm)
        enddo
      enddo
!
!-----------------------------------------------------------------------
!---grand horizontal loop-----------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------
!-----------------
!.......................................................................
!$omp do private(l,j,i,dpup,zx,zy,zne,znw,fahp)
!.......................................................................
!-----------------
!-----------------
!-----------------------------------------------------------------------
!
      vertical_loop: do l=1,lm
!
!-----------------------------------------------------------------------
!
        dpup=pdsg1(l)
!
!-----------------------------------------------------------------------
!-------------mass fluxes, on h points----------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_h1
          do i=its_b1,ite_h1
            zx(i,j)=(z(i,j,l)-z(i-1,j,l))*pfx(i,j,l)
            zy(i,j)=(z(i,j,l)-z(i,j-1,l))*pfy(i,j,l)
            zne(i,j)=(z(i,j,l)-z(i-1,j-1,l))*pfne(i,j,l)
            znw(i,j)=(z(i-1,j,l)-z(i,j-1,l))*pfnw(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---advection of height-------------------------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1
          fahp=-fah(j)/dt
          do i=its_b1,ite_b1
            w(i,j,l)=((zx(i,j)+zx(i+1,j)+zy(i,j)+zy(i,j+1)) &
                    +(zne(i+1,j+1)+zne(i,j) &
                     +znw(i,j+1)+znw(i+1,j))*0.25)*fahp &
                    /(dsg2(l)*pdo(i,j)+pdsg1(l)) &
                    +w(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---setting w on h points to 0. along boundaries------------------------
!-----------------------------------------------------------------------
        if(.not.global) then
!-----------------------------------------------------------------------
          if(s_bdy)then
            do j=jts,jts+1
              do i=its,ite
                w(i,j,l)=0.
              enddo
            enddo
          endif
!
          if(n_bdy)then
            do j=jte-1,jte
              do i=its,ite
                w(i,j,l)=0.
              enddo
            enddo
          endif
!
          if(w_bdy)then
            do j=jts,jte
              do i=its,its+1
                w(i,j,l)=0.
              enddo
            enddo
          endif
!
          if(e_bdy)then
            do j=jts,jte
              do i=ite-1,ite
                w(i,j,l)=0.
              enddo
            enddo
          endif
!-----------------------------------------------------------------------
        endif
!-----------------------------------------------------------------------
      enddo vertical_loop
!------------------
!------------------
!.......................................................................
!$omp end do
!$omp end parallel
!.......................................................................
!------------------
!------------------
!-----------------------------------------------------------------------
!
                        endsubroutine cdzdt
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine cdwdt &
(global,hydro,inpes,jnpes &
,lm,ntsd &
,dt,g &
,dsg2,pdsg1 &
,fah &
,pd,pdo &
,psgdt &
,dwdt,pdwdt,w &
,pint &
!---temporary arguments-------------------------------------------------
,pfx,pfy,pfne,pfnw)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),parameter:: &
 nsmud=00 &                  ! number of smoothing iterations
,lnsnh=02                    ! # of rows with smoothing along boundaries

real(kind=kfpt),parameter:: &
 epsfc=9.81 &                ! limiter value
,epsn=-epsfc &               ! floor value for vertical acceleration
,epsp=epsfc &                ! upper limit for vertical acceleration
,wa=0.125 &                  ! weighting factor
,wb=0.5 &                    ! weighting factor
!,wad=0.125 &                 ! lateral smoothing weight 
!,wad=0.0625 &                ! lateral smoothing weight
!,wad=0.075 &                 ! lateral smoothing weight
!,wad=0.050 &                 ! lateral smoothing weight
,wad=0. &                    ! lateral smoothing weight
!,wp=0.075 &                  ! time smoothing weight
,wp=0.                       ! time smoothing weight
!-----------------------------------------------------------------------
logical(kind=klog),intent(in):: &
 global &                    ! global or regional
,hydro                       ! hydrostatic or nonhydrostatic

integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,ntsd &                      ! time step
,inpes &                     ! tasks in x direction
,jnpes                       ! tasks in y direction

real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,g                           ! gravity

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures

real(kind=kfpt),dimension(jds:jde),intent(in):: &
 fah                         ! delta sigmas

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd &                        ! sigma range pressure difference
,pdo                         ! old sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(in):: &
 psgdt                       ! vertical mass flux

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 dwdt &                      ! nonhydrostatic correction factor
,pdwdt &                     ! previous nonhydrostatic correction factor
,w                           ! w wind component

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(inout):: &
 pint                        ! pressure
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in) :: &
 pfx &                       ! mass flux, x direction
,pfy &                       ! mass flux, y direction
,pfne &                      ! mass flux, ne direction
,pfnw                        ! mass flux, nw direction
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,imn &                       !
,imx &                       !
,j &                         ! index in y direction
,jmn &                       !
,jmx &                       !
,l &                         ! index in p direction
,lmn &                       !
,lmx &                       !
,kn &                        ! counter
,kp &                        ! counter
,ks                        ! smoothing counter

real(kind=kfpt):: &
 dwdtmn &                    ! minimum value of dwdt
,dwdtmx &                    ! maximum value of dwdt
,dwdtp &                     ! nonhydrostatic correction factor at the point
,fahp &                      ! grid factor
,rdt &                       ! 1/dt
,rg                          ! 1/g

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1):: &
 tta                         ! advection through upper interface
real(kind=kfpt),dimension(its:ite,jts:jte):: &
 ttb                         ! advection through lower interface
real(kind=kfpt),dimension(its_b1:ite_h1,jts_b1:jte_h1):: &
 wne &                       ! height flux, ne direction
,wnw &                       ! height flux, nw direction
,ww &                        ! temporary for lateral smoothing
,wx &                        ! height flux, x direction
,wy                          ! height flux, y direction
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
      if(hydro.or.ntsd.lt.2) then
!-----------------------------------------------------------------------
        do l=1,lm
          do j=jts,jte
            do i=its,ite
              dwdt(i,j,l)=1.
              pdwdt(i,j,l)=1.
              pint(i,j,l+1)=dsg2(l)*pd(i,j)+pdsg1(l)+pint(i,j,l)
            enddo
          enddo
        enddo
!-----------------------------------------------------------------------
        return
!-----------------------------------------------------------------------
      endif
!-----------------------------------------------------------------------
      if(.not.global) then 
!---smoothing of w on h points along boundaries-------------------------
        if(nsmud.gt.0) then
!-----------------------------------------------------------------------
          do ks=1,nsmud
!-----------------------------------------------------------------------
            do l=1,lm
!-----------------------------------------------------------------------
              if(s_bdy)then
                do j=jts+1,lnsnh
                  do i=its_b1,ite_b1
                    ww(i,j)=(w(i,j-1,l)+w(i-1,j,l) &
                            +w(i+1,j,l)+w(i,j+1,l))*wa &
                            +w(i,j,l)*wb       
                  enddo
                enddo
              endif
!
              if(n_bdy)then
                do j=jte-lnsnh+1,jte-1
                  do i=its_b1,ite_b1
                    ww(i,j)=(w(i,j-1,l)+w(i-1,j,l) &
                            +w(i+1,j,l)+w(i,j+1,l))*wa &
                            +w(i,j,l)*wb       
                  enddo
                enddo
              endif
!
              if(w_bdy)then
                do j=max(jts,jds+lnsnh),min(jte,jde-lnsnh)
                  do i=its+1,its-1+lnsnh
                    ww(i,j)=(w(i,j-1,l)+w(i-1,j,l) &
                            +w(i+1,j,l)+w(i,j+1,l))*wa &
                            +w(i,j,l)*wb       
                  enddo
                enddo
              endif
!
              if(e_bdy)then
                do j=max(jts,jds+lnsnh),min(jte,jde-lnsnh)
                  do i=ite-lnsnh+1,ite-1
                    ww(i,j)=(w(i,j-1,l)+w(i-1,j,l) &
                            +w(i+1,j,l)+w(i,j+1,l))*wa &
                            +w(i,j,l)*wb       
                  enddo
                enddo
              endif
!
              if(s_bdy)then
                do j=its+1,its-1+lnsnh
                  do i=its_b1,ite_b1
                    w(i,j,l)=ww(i,j)
                  enddo
                enddo
              endif
!
              if(n_bdy)then
                do j=jte-lnsnh+1,jte-1
                  do i=its_b1,ite_b1
                    w(i,j,l)=ww(i,j)
                  enddo
                enddo
              endif
!
              if(w_bdy)then
                do j=max(jts,jds+lnsnh),min(jte,jde-lnsnh)
                  do i=its+1,lnsnh
                    w(i,j,l)=ww(i,j)
                  enddo
                enddo
              endif
!
              if(e_bdy)then
                do j=max(jts,jds+lnsnh),min(jte,jde-lnsnh)
                  do i=ite-lnsnh+1,ite-1
                    w(i,j,l)=ww(i,j)
                  enddo
                enddo
              endif
!-----------------------------------------------------------------------
            enddo
!-----------------------------------------------------------------------
          enddo
!-----------------------------------------------------------------------
        endif ! end of smoothing of w along lateral boundaries
!-----------------------------------------------------------------------
!
      endif
!
!-----------------------------------------------------------------------
!---local derivative of w on h points-----------------------------------
!-----------------------------------------------------------------------
!
      rdt=1./dt
!
!.......................................................................
!$omp parallel do private(l,j,i)
!.......................................................................
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            dwdt(i,j,l)=(w(i,j,l)-dwdt(i,j,l))*rdt
          enddo
        enddo
      enddo
!.......................................................................
!$omp end parallel do 
!.......................................................................
!-----------------------------------------------------------------------
!---vertical advection of w---------------------------------------------
!-----------------------------------------------------------------------
      do j=jts,jte
        do i=its,ite
          ttb(i,j)=0.
        enddo
      enddo
!
      do l=1,lm-1
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            tta(i,j)=(w(i,j,l+1)-w(i,j,l))*psgdt(i,j,l)*0.5
            dwdt(i,j,l)=(tta(i,j)+ttb(i,j))/(dsg2(l)*pdo(i,j)+pdsg1(l)) &
                       +dwdt(i,j,l)
            ttb(i,j)=tta(i,j)
          enddo
        enddo
      enddo
!
      do j=jts,jte
        do i=its,ite
          dwdt(i,j,lm)=ttb(i,j)/(dsg2(lm)*pdo(i,j)+pdsg1(lm)) &
                      +dwdt(i,j,lm)
        enddo
      enddo
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private (l,j,i,wx,wy,wne,wnw,fahp)
!.......................................................................
!-----------------------------------------------------------------------
      do l=1,lm
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!---w fluxes, on h points-----------------------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_h1
          do i=its_b1,ite_h1
            wx(i,j)=(w(i,j,l)-w(i-1,j,l))*pfx(i,j,l)
            wy(i,j)=(w(i,j,l)-w(i,j-1,l))*pfy(i,j,l)
            wne(i,j)=(w(i,j,l)-w(i-1,j-1,l))*pfne(i,j,l)
            wnw(i,j)=(w(i-1,j,l)-w(i,j-1,l))*pfnw(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---advection of w------------------------------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1
          fahp=-fah(j)/dt
          do i=its_b1,ite_b1
            dwdt(i,j,l)=((wx(i,j)+wx(i+1,j)+wy(i,j)+wy(i,j+1)) &
                        +(wne(i+1,j+1)+wne(i,j) &
                         +wnw(i,j+1)+wnw(i+1,j))*0.25)*fahp &
                       /(dsg2(l)*pdo(i,j)+pdsg1(l)) &
                       +dwdt(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!
      enddo
!.......................................................................
!$omp end parallel do 
!.......................................................................
!
!-----------------------------------------------------------------------
      if(global) then
        btim=timef()
        call swaphn(dwdt,ims,ime,jms,jme,lm,inpes)
        swaphn_tim=swaphn_tim+timef()-btim
!
        btim=timef()
        call polehn(dwdt,ims,ime,jms,jme,lm,inpes,jnpes)
        polehn_tim=polehn_tim+timef()-btim
      endif
!
!-----------------------------------------------------------------------
!------------spatial filtering of dwdt----------------------------------
!-----------------------------------------------------------------------
!
      if(wad.gt.0.) then
!
!.......................................................................
!$omp parallel do private (l,j,i,ww)
!.......................................................................
        do l=1,lm
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              ww(i,j)=(dwdt(i,j-1,l)+dwdt(i-1,j,l) &
                      +dwdt(i+1,j,l)+dwdt(i,j+1,l)-dwdt(i,j,l)*4.0)
            enddo
          enddo
!
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              dwdt(i,j,l)=ww(i,j)*wad+dwdt(i,j,l)
            enddo
          enddo
        enddo
!.......................................................................
!$omp end parallel do 
!.......................................................................
      endif
!
!-----------------------------------------------------------------------
!
      if(global) then
        btim=timef()
        call swaphn(dwdt,ims,ime,jms,jme,lm,inpes)
        swaphn_tim=swaphn_tim+timef()-btim
!!!   if(mype==23)write(0,*)' in cdwdt swaphn_tim=',swaphn_tim
!
        btim=timef()
        call polehn(dwdt,ims,ime,jms,jme,lm,inpes,jnpes)
        polehn_tim=polehn_tim+timef()-btim
      endif
!-----------------------------------------------------------------------
      rg=1./g
!
      kn=0
      kp=0
      imn=0
      jmn=0
      lmn=0
      imx=0
      jmx=0
      lmx=0
      dwdtmx=0.
      dwdtmn=0.
!
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            dwdtp=dwdt(i,j,l)
            if(dwdtp.gt.dwdtmx) then
              dwdtmx=dwdtp
              imx=i
              jmx=j
              lmx=l
            endif
            if(dwdtp.lt.dwdtmn) then
              dwdtmn=dwdtp
              imn=i
              jmn=j
              lmn=l
            endif
            if(dwdtp.lt.epsn) then
              dwdtp=epsn
              kn=kn+1
            endif
            if(dwdtp.gt.epsp) then
              dwdtp=epsp
              kp=kp+1
            endif
            dwdt(i,j,l)=(dwdtp*rg+1.)*(1.-wp)+pdwdt(i,j,l)*wp
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!
      if(.not.global) then 
!
!-----------------------------------------------------------------------
!---setting dwdt on h points to 1. along boundaries---------------------
!-----------------------------------------------------------------------
        do l=1,lm
          if(s_bdy)then
            do j=jts,jts+1
              do i=its,ite
                dwdt(i,j,l)=1.
              enddo
            enddo
          endif
!
          if(n_bdy)then
            do j=jte-1,jte
              do i=its,ite
                dwdt(i,j,l)=1.
              enddo
            enddo
          endif
!
          if(w_bdy)then
            do j=jts,jte
              do i=its,its+1
                dwdt(i,j,l)=1.
              enddo
            enddo
          endif
!
          if(e_bdy)then
            do j=jts,jte
              do i=ite-1,ite
                dwdt(i,j,l)=1.
              enddo
            enddo
          endif
        enddo
!-----------------------------------------------------------------------
      endif ! regional 
!-----------------------------------------------------------------------
!
                        endsubroutine cdwdt
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine vsound &
(global,hydro &
,lm,ntsd &
,cp,dt,pt &
,dsg2,pdsg1 &
,pd &
,cw,q,rtop &
,dwdt,t,w &
,pint)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),parameter:: &
 nsmud=0                     ! number of smoothing iterations

real(kind=kfpt),parameter:: &
 wght=0.35                   ! first guess weight
!-----------------------------------------------------------------------
logical(kind=klog),intent(in):: &
 global &                    ! global or regional
,hydro                       ! hydrostatic or nonhydrostatic

integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,ntsd                        ! time step

real(kind=kfpt),intent(in):: &
 cp &                        ! cp
,dt &                        ! dynamics time step
,pt                          ! pressure at the top of the model's atmosphere

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd                          ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 cw &                        ! condensate
,q &                         ! specific humidity
,rtop                        ! rt/p

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 dwdt &                      ! nonhydrostatic correction factor
,t &                         ! previous nonhydrostatic correction factor
,w                           ! w wind component

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(inout):: &
 pint
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction

real(kind=kfpt):: &
 cappa &                     ! R/cp
,cofl &                      !
,dp &                        !
,delp &                      !
,dppl &                      !
,dpstr &                     !
,dptl &                      !
,fcc &                       ! 
,ffc &                       ! 
,gdt &                       ! g*dt
,gdt2 &                      ! gdt**2
,pp1 &                       !
,rcph &                      ! 
,rdpdn &                     !
,rdpup &                     !
,tfc &                       !
,tmp                         !

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1):: &
 dptu                        !

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 b1 &                        !
,b2 &                        !
,b3 &                        !
,c0 &                        !
,rdpp                        !

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm+1):: &
 chi &                       !
,coff &                      !
,dfrh &                      !
,pnp1 &                      !
,pone &                      !
,pstr                        !
!-----------------
!-----------------
integer(kind=kint) :: &
 jstart &
,jstop &
,nth &
,omp_get_num_threads &
,omp_get_thread_num &
,tid
!-----------------
!-----------------
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      if(hydro.or.ntsd.lt.2) return
!
!-----------------------------------------------------------------------
      cappa=r/cp
      gdt=g*dt
      gdt2=gdt*gdt
      ffc=-r*0.25/gdt2
      rcph=0.5/cp
!-----------------------------------------------------------------------
!-----------------
!-----------------
!.......................................................................
!$omp parallel &
!$omp private(nth,tid,i,j,l,jstart,jstop,dppl,dpstr,pp1,dp,tfc,fcc,&
!$omp         cofl,rdpdn,rdpup,tmp,dptl,delp )
!.......................................................................
         nth = omp_get_num_threads()
         tid = omp_get_thread_num()
         call looplimits(tid, nth, jts_b1, jte_b1, jstart, jstop)
!-----------------
!-----------------
 
       do j=jstart,jstop
        do i=its_b1,ite_b1
          pone(i,j,1)=pt
          pstr(i,j,1)=pt
          pnp1(i,j,1)=pt
        enddo
      enddo
!
      do l=2,lm+1
        do j=jstart,jstop
          do i=its_b1,ite_b1
            dppl=dsg2(l-1)*pd(i,j)+pdsg1(l-1)
            rdpp(i,j,l-1)=1./dppl
            dpstr=dwdt(i,j,l-1)*dppl
            pstr(i,j,l)=pstr(i,j,l-1)+dpstr
            pp1=pnp1(i,j,l-1)+dpstr
            pone(i,j,l)=pint(i,j,l)
            dp=(pp1-pone(i,j,l))*wght
            pnp1(i,j,l)=pone(i,j,l)+dp
            tfc=q(i,j,l-1)*0.608-cw(i,j,l-1)+1.
            fcc=(1.-cappa*tfc)*tfc*ffc
            cofl=t(i,j,l-1)*dppl*fcc &
                /((pnp1(i,j,l-1)+pnp1(i,j,l))*0.5)**2
            coff(i,j,l-1)=cofl
            dfrh(i,j,l)=(pstr(i,j,l-1)+pstr(i,j,l) &
                        -pone(i,j,l-1)-pone(i,j,l))*cofl
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      do l=2,lm
        do j=jstart,jstop
          do i=its_b1,ite_b1
            rdpdn=rdpp(i,j,l)
            rdpup=rdpp(i,j,l-1)
            b1(i,j,l)=coff(i,j,l-1)+rdpup
            b2(i,j,l)=coff(i,j,l-1)+coff(i,j,l)-rdpup-rdpdn
            b3(i,j,l)=coff(i,j,l)+rdpdn
            c0(i,j,l)=-dfrh(i,j,l)-dfrh(i,j,l+1)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
       do j=jstart,jstop
        do i=its_b1,ite_b1
          b2(i,j,lm)=b2(i,j,lm)+b3(i,j,lm)
        enddo
      enddo
!
      do l=3,lm
        do j=jstart,jstop
          do i=its_b1,ite_b1
            tmp=-b1(i,j,l)/b2(i,j,l-1)
            b2(i,j,l)=b3(i,j,l-1)*tmp+b2(i,j,l)
            c0(i,j,l)=c0(i,j,l-1)*tmp+c0(i,j,l)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
       do j=jstart,jstop
        do i=its_b1,ite_b1
          chi(i,j,1)=0.
          chi(i,j,lm)=c0(i,j,lm)/b2(i,j,lm)
          chi(i,j,lm+1)=chi(i,j,lm)
        enddo
      enddo
!
      do l=lm-1,2,-1
        do j=jstart,jstop
          do i=its_b1,ite_b1
            chi(i,j,l)=(-b3(i,j,l)*chi(i,j,l+1)+c0(i,j,l))/b2(i,j,l)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      do l=1,lm+1
       do j=jstart,jstop
           do i=its_b1,ite_b1
            pnp1(i,j,l)=chi(i,j,l)+pstr(i,j,l)
            pint(i,j,l)=pnp1(i,j,l)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      do j=jstart,jstop
        do i=its_b1,ite_b1
          dptu(i,j)=0.
        enddo
      enddo
!
      do l=1,lm
        do j=jstart,jstop
          do i=its_b1,ite_b1
            dptl=pnp1(i,j,l+1)-pone(i,j,l+1)
            t(i,j,l)=(dptu(i,j)+dptl)*rtop(i,j,l)*rcph+t(i,j,l)
            delp=(pnp1(i,j,l+1)-pnp1(i,j,l))*rdpp(i,j,l)
            w(i,j,l)=(delp-dwdt(i,j,l))*gdt+w(i,j,l)
            dwdt(i,j,l)=delp
            dptu(i,j)=dptl
          enddo
        enddo
      enddo
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------
!-----------------------------------------------------------------------
!
      if(.not.global) then
!
!-----------------------------------------------------------------------
!---setting dwdt on h points to 1. along boundaries---------------------
!-----------------------------------------------------------------------
        if(s_bdy)then
          do l=1,lm
            do j=jts,jts+1
              do i=its,ite
                dwdt(i,j,l)=1.
              enddo
            enddo
          enddo
        endif
!
        if(n_bdy)then
          do l=1,lm
            do j=jte-1,jte
              do i=its,ite
                dwdt(i,j,l)=1.
              enddo
            enddo
          enddo
        endif
!
        if(w_bdy)then
          do l=1,lm
            do j=jts,jte
              do i=its,its+1
                dwdt(i,j,l)=1.
              enddo
            enddo
          enddo
        endif
!
        if(e_bdy)then
          do l=1,lm
            do j=jts,jte
              do i=ite-1,ite
                dwdt(i,j,l)=1.
              enddo
            enddo
          enddo
        endif
!-----------------------------------------------------------------------
      endif ! regional
!-----------------------------------------------------------------------
!
                        endsubroutine vsound
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        subroutine adv2 &
(global &
,idtadt,kss,kse,lm,lnsad &
,dt,rdyh &
,dsg2,pdsg1 &
,fah,rdxh &
,pd,pdo &
,psgdt &
,up,vp &
,q2,indx_q2 &
,s,sp &
!---temporary arguments-------------------------------------------------
,pfne,pfnw,pfx,pfy,s1,tcs)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

include 'kind.inc'
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 cfc=1.533 &                 ! adams-bashforth positioning in time
,bfc=1.-cfc &                ! adams bashforth positioning in time
,epsq=1.e-20 &               ! floor value for specific humidity
,epsq2=0.02 &                ! floor value for 2tke
,pfc=1.+4./6. &              ! 4th order momentum advection
,sfc=-1./6. &                ! 4th order momentum advection
,epscm=2.e-6 &               ! a floor value (not used)
,w1=1.0 &                    ! crank-nicholson uncentering
!,w1=0.80 &                   ! crank-nicholson uncentering
,w2=2.-w1                    ! crank-nicholson uncentering

logical(kind=klog),intent(in):: &
 global

integer(kind=kint),intent(in):: &
 idtadt &                    !
,kse &                       ! terminal species index
,kss &                       ! initial species index
,lm &                        ! total # of levels
,lnsad &                     !
,indx_q2                     ! location of q2 in tracer arrays

real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,rdyh                        !

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures

real(kind=kfpt),dimension(jds:jde),intent(in):: &
 fah &                       ! grid factor
,rdxh                        !

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd &                        ! sigma range pressure difference
,pdo                         ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(in):: &
 psgdt                       ! vertical mass flux

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 up &                        !
,vp                          !

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
  q2                          ! 2.*TKE

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,kss:kse),intent(inout):: &
 s &                         ! tracers
,sp                          ! s at previous time level

!---temporary arguments-------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 pfne &                      ! mass flux, ne direction
,pfnw &                      ! mass flux, nw direction
,pfx &                       ! mass flux, x direction
,pfy                         ! mass flux, y direction

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,kss:kse),intent(inout):: &
 s1 &                        ! intermediate value of sqrt(s)
,tcs                         ! timechange of s

!--local variables------------------------------------------------------
integer(kind=kint):: &
 i &                         !
,iap &                       !
,ibeg &                      !
,iend &                      !
,j &                         !
,jap &                       !
,jbeg &                      !
,jend &                      !
,ks &                        !
,l                           !

real(kind=kfpt):: &
 cf &                        ! temporary
,cms &                       ! temporary
,dtq &                       ! dt/4
,emhp &                      !
,enhp &                      !
,fahp &                      ! temporary grid factor
,pp &                        !
,qq &                        !
,rdp &                       ! 1/deltap
,vvlo &                      ! vertical velocity, lower interface
,vvup &                      ! vertical velocity, upper interface
,pvvup                       ! vertical mass flux, upper interface

real(kind=kfpt),dimension(ims:ime,jms:jme):: &
 pdop &                      ! hydrostatic pressure difference at v points
,pvvlo &                     ! vertical mass flux, lower interface
,ss1 &                       ! extrapolated species between time levels 
,ssne &                      ! flux, ne direction
,ssnw &                      ! flux, nw direction
,ssx &                       ! flux, x direction
,ssy                         ! flux, y direction

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 crs &                       ! vertical advection temporary
,rcms                        ! vertical advection temporary

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,kss:kse):: &
 rsts                        ! vertical advection temporary

!-----------------
!-----------------
integer(kind=kint) :: &
 jstart &
,jstop &
,nth &
,omp_get_num_threads &
,omp_get_thread_num &
,tid
!-----------------
!-----------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
!-----------------
!-----------------
!.......................................................................
!$omp parallel private(nth,tid,i,j,l,ks,jstart,jstop)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth, jts_h1, jte_h1, jstart, jstop)
!-----------------
!-----------------
      do ks=kss,kse
        do l=1,lm
          do j=jstart,jstop
            do i=its_h1,ite_h1
              s(i,j,l,ks)=max(s(i,j,l,ks),epsq)
            enddo
          enddo
        enddo
      enddo
!
      do j=jstart,jstop
        do i=its_h1,ite_h1
          s(i,j,1,indx_q2)=max((q2 (i,j,1)+epsq2)*0.5,epsq2)
        enddo
      enddo
      do l=2,lm
        do j=jstart,jstop
          do i=its_h1,ite_h1
            s(i,j,l,indx_q2)=max((q2 (i,j,l)+q2 (i,j,l-1))*0.5,epsq2)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      do ks=kss,kse ! loop by species
!-----------------------------------------------------------------------
        do l=1,lm
          do j=jstart,jstop
            do i=its_h1,ite_h1
              s1(i,j,l,ks)=sqrt(s(i,j,l,ks))
            enddo
          enddo
        enddo
!-----------------------------------------------------------------------
!
      enddo ! end of the loop by species
!
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------
!-----------------------------------------------------------------------
!-----------------
!-----------------
!.......................................................................
!$omp parallel private(nth,tid,i,j,l,jstart,jstop,dtq,vvlo,cms,rdp,&
!$omp                  vvup,cf,pvvup)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth, jts_b1, jte_b1, jstart, jstop)
!-----------------
!-----------------
      do j=jstart,jstop
        do i=its_b1,ite_b1
          pdop(i,j)=(pd(i,j)+pdo(i,j))*0.5
        enddo
      enddo
!-----------------------------------------------------------------------
!---crank-nicholson vertical advection----------------------------------
!-----------------------------------------------------------------------
      dtq=dt*0.25*idtadt
      do j=jstart,jstop
        do i=its_b1,ite_b1
          pvvlo(i,j)=psgdt(i,j,1)*dtq
          vvlo=pvvlo(i,j)/(dsg2(1)*pdop(i,j)+pdsg1(1))
          cms=-vvlo*w2+1.
          rcms(i,j,1)=1./cms
          crs(i,j,1)=vvlo*w2
!
          do ks=kss,kse
            rsts(i,j,1,ks)=(-vvlo*w1) &
                          *(s1(i,j,2,ks)-s1(i,j,1,ks)) &
                          +s1(i,j,1,ks)
          enddo
        enddo
      enddo
      do l=2,lm-1
        do j=jstart,jstop
          do i=its_b1,ite_b1
            rdp=1./(dsg2(l)*pdop(i,j)+pdsg1(l))
            pvvup=pvvlo(i,j)
            pvvlo(i,j)=psgdt(i,j,l)*dtq
!
            vvup=pvvup*rdp
            vvlo=pvvlo(i,j)*rdp
!
            cf=-vvup*w2*rcms(i,j,l-1)
            cms=-crs(i,j,l-1)*cf+((vvup-vvlo)*w2+1.)
            rcms(i,j,l)=1./cms
            crs(i,j,l)=vvlo*w2
!
            do ks=kss,kse
              rsts(i,j,l,ks)=-rsts(i,j,l-1,ks)*cf+s1(i,j,l,ks) &
                             -(s1(i,j,l  ,ks)-s1(i,j,l-1,ks))*vvup*w1 &
                             -(s1(i,j,l+1,ks)-s1(i,j,l  ,ks))*vvlo*w1
            enddo
          enddo
        enddo
      enddo
      do j=jstart,jstop
        do i=its_b1,ite_b1
          pvvup=pvvlo(i,j)
          vvup=pvvup/(dsg2(lm)*pdop(i,j)+pdsg1(lm))
!
          cf=-vvup*w2*rcms(i,j,lm-1)
          cms=-crs(i,j,lm-1)*cf+(vvup*w2+1.)
          rcms(i,j,lm)=1./cms
          crs(i,j,lm)=0.
!
          do ks=kss,kse
            rsts(i,j,lm,ks)=-rsts(i,j,lm-1,ks)*cf+s1(i,j,lm,ks) &
                           -(s1(i,j,lm,ks)-s1(i,j,lm-1,ks))*vvup*w1
!
            tcs(i,j,lm,ks)=rsts(i,j,lm,ks)*rcms(i,j,lm)-s1(i,j,lm,ks)
          enddo
        enddo
      enddo
      do ks=kss,kse
        do l=lm-1,1,-1
          do j=jstart,jstop
            do i=its_b1,ite_b1
              tcs(i,j,l,ks)=(-crs(i,j,l)*(s1(i,j,l+1,ks)+tcs(i,j,l+1,ks)) &
                             +rsts(i,j,l,ks)) &
                           *rcms(i,j,l)-s1(i,j,l,ks)
            enddo
          enddo
        enddo
      enddo
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do &
!$omp private(ks,l,j,i,ss1,ssx,ssy,ssne,ssnw,ibeg,iend, &
!$omp         jbeg,jend,fahp)   
!.......................................................................
!-----------------------------------------------------------------------
!
      do ks=kss,kse ! loop by species
!
!-----------------------------------------------------------------------
        do l=1,lm
          do j=jts_h1,jte_h1
            do i=its_h1,ite_h1
              ss1(i,j)=s1(i,j,l,ks)*cfc+sp(i,j,l,ks)*bfc
              sp(i,j,l,ks)=s1(i,j,l,ks)
            enddo
          enddo
!---temperature fluxes, on h points-------------------------------------
          do j=jts_b1,jte_h1
            do i=its_b1,ite_h1
              ssx(i,j)=(ss1(i,j)-ss1(i-1,j))*pfx(i,j,l)
              ssy(i,j)=(ss1(i,j)-ss1(i,j-1))*pfy(i,j,l)
!
              ssne(i,j)=(ss1(i,j)-ss1(i-1,j-1))*pfne(i,j,l)
              ssnw(i,j)=(ss1(i-1,j)-ss1(i,j-1))*pfnw(i,j,l)
            enddo
          enddo
!---advection of species------------------------------------------------
          if(adv_standard)then
            ibeg=max(its,ids+1+lnsad)
            iend=min(ite,ide-1-lnsad)
            jbeg=max(jts,jds+1+lnsad)
            jend=min(jte,jde-1-lnsad)
!
            do j=jbeg,jend
              fahp=fah(j)*idtadt
              do i=ibeg,iend
                tcs(i,j,l,ks)=(((ssx (i  ,j  )+ssx (i+1,j  ) &
                                +ssy (i  ,j  )+ssy (i  ,j+1)) &
                               +(ssne(i+1,j+1)+ssne(i  ,j  ) &
                                +ssnw(i  ,j+1)+ssnw(i+1,j  ))*0.25) &
                               *fahp) &
                             /(dsg2(l)*pdop(i,j)+pdsg1(l)) &
                             +tcs(i,j,l,ks)
              enddo
            enddo
          endif
        enddo
!-----------------------------------------------------------------------
!
      enddo ! end of the loop by the species
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!---regional branch-----------------------------------------------------
!-----------------------------------------------------------------------
      if(.not.global.and.adv_upstream) then
!-----------------------------------------------------------------------
        enhp=-dt*rdyh*0.25*idtadt
!-----------------------------------------------------------------------
        do l=1,lm
!-----------------------------------------------------------------------
!
!***  Upstream advection along southern rows
!
          do j=jts_b1,min(jte,jds+lnsad)
            emhp=-dt*rdxh(j)*0.25*idtadt
            do i=its_b1,ite_b1
              pp=(up(i-1,j-1,l)+up(i  ,j-1,l) &
                 +up(i-1,j  ,l)+up(i  ,j  ,l))*emhp
              qq=(vp(i-1,j-1,l)+vp(i  ,j-1,l) &
                 +vp(i-1,j  ,l)+vp(i  ,j  ,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              do ks=kss,kse
                tcs(i,j,l,ks)=(s1(i+iap,j,l,ks)-s1(i,j,l,ks))*pp &
                             +(s1(i,j+jap,l,ks)-s1(i,j,l,ks))*qq &
                             +(s1(i,j,l,ks)-s1(i+iap,j,l,ks) &
                              -s1(i,j+jap,l,ks)+s1(i+iap,j+jap,l,ks)) &
                             *pp*qq &
                             +tcs(i,j,l,ks)
              enddo
            enddo
          enddo
!
!***  Upstream advection along northern rows
!
          do j=max(jts,jde-lnsad),jte_b1
            emhp=-dt*rdxh(j)*0.25*idtadt
            do i=its_b1,ite_b1
              pp=(up(i-1,j-1,l)+up(i  ,j-1,l) &
                 +up(i-1,j  ,l)+up(i  ,j  ,l))*emhp
              qq=(vp(i-1,j-1,l)+vp(i  ,j-1,l) &
                 +vp(i-1,j  ,l)+vp(i  ,j  ,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              do ks=kss,kse
                tcs(i,j,l,ks)=(s1(i+iap,j,l,ks)-s1(i,j,l,ks))*pp &
                             +(s1(i,j+jap,l,ks)-s1(i,j,l,ks))*qq &
                             +(s1(i,j,l,ks)-s1(i+iap,j,l,ks) &
                              -s1(i,j+jap,l,ks)+s1(i+iap,j+jap,l,ks)) &
                             *pp*qq &
                             +tcs(i,j,l,ks)
              enddo
            enddo
          enddo
!
!***  Upstream advection along western rows
!
          do j=max(jts,jds+1+lnsad),min(jte,jde-1-lnsad)
            emhp=-dt*rdxh(j)*0.25*idtadt
            do i=its_b1,min(ite,ids+lnsad)
              pp=(up(i-1,j-1,l)+up(i  ,j-1,l) &
                 +up(i-1,j  ,l)+up(i  ,j  ,l))*emhp
              qq=(vp(i-1,j-1,l)+vp(i  ,j-1,l) &
                 +vp(i-1,j  ,l)+vp(i  ,j  ,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              do ks=kss,kse
                tcs(i,j,l,ks)=(s1(i+iap,j,l,ks)-s1(i,j,l,ks))*pp &
                             +(s1(i,j+jap,l,ks)-s1(i,j,l,ks))*qq &
                             +(s1(i,j,l,ks)-s1(i+iap,j,l,ks) &
                              -s1(i,j+jap,l,ks)+s1(i+iap,j+jap,l,ks)) &
                             *pp*qq &
                             +tcs(i,j,l,ks)
              enddo
            enddo
          enddo
!
!***  Upstream advection along eastern rows
!
          do j=max(jts,jds+1+lnsad),min(jte,jde-1-lnsad)
            emhp=-dt*rdxh(j)*0.25*idtadt
            do i=max(its,ide-lnsad),ite_b1
              pp=(up(i-1,j-1,l)+up(i  ,j-1,l) &
                 +up(i-1,j  ,l)+up(i  ,j  ,l))*emhp
              qq=(vp(i-1,j-1,l)+vp(i  ,j-1,l) &
                 +vp(i-1,j  ,l)+vp(i  ,j  ,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              do ks=kss,kse
                tcs(i,j,l,ks)=(s1(i+iap,j,l,ks)-s1(i,j,l,ks))*pp &
                             +(s1(i,j+jap,l,ks)-s1(i,j,l,ks))*qq &
                             +(s1(i,j,l,ks)-s1(i+iap,j,l,ks) &
                              -s1(i,j+jap,l,ks)+s1(i+iap,j+jap,l,ks)) &
                             *pp*qq &
                             +tcs(i,j,l,ks)
              enddo
            enddo
          enddo
!-----------------------------------------------------------------------
        enddo
!-----------------------------------------------------------------------
      endif ! regional lateral boundaries
!-----------------------------------------------------------------------
!
                        endsubroutine adv2
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine mono &
(idtadt,kss,kse,lm &
,dsg2,pdsg1 &
,dare &
,pd &
,indx_q2 &
,s &
!---temporary arguments-------------------------------------------------
,s1,tcs)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
implicit none

include 'kind.inc'
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 epsq=1.e-20 &               ! floor value for specific humidity
,epsq2=0.02                  ! floor value for 2tke

integer(kind=kint),intent(in):: &
 idtadt &                    !
,kse &                       ! terminal species index
,kss &                       ! initial species index
,lm &                        ! total # of levels
,indx_q2                     ! location of q2 in 4-d tracer arrays

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures

real(kind=kfpt),dimension(jds:jde),intent(in):: &
 dare                        ! grid box area

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd                          ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,kss:kse),intent(inout):: &
 s                           ! s at previous time level

!---temporary arguments-------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,kss:kse),intent(inout):: &
 s1 &                        ! intermediate value of s
,tcs                         ! timechange of s
!--local variables------------------------------------------------------
integer(kind=kint):: &
 i &                         !
,ierr &                      !
,irecv &                     !
,j &                         !
,ks &                        !
,l &                         !
,lngth                       !

real(kind=kfpt):: &
 s1p &                       !
,smax &                      ! local maximum
,smin &                      ! local minimum
,smaxh &                     ! horizontal local maximum
,sminh &                     ! horizontal local minimum
,smaxv &                     ! vertical local maximum
,sminv &                     ! vertical local minimum
,sn &                        !
,steep                       !

real(kind=kdbl):: &
 dsp &                       !
,rfacs &                     !
,sfacs &                     !
,sumns &                     !
,sumps                       !

real(kind=kdbl):: &
 gsump &                     !
,xsump 

real(kind=kdbl),dimension(2*kss-1:2*kse):: &
 vgsums                      !

real(kind=kdbl),dimension(2*kss-1:2*kse,1:lm):: &
 gsums &                     ! sum of neg/pos changes all global fields
,xsums                       ! sum of neg/pos changes all local fields

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 dvol &                      ! grid box volume
,rdvol                       ! 1./grid box volume
!-----------------------------------------------------------------------
integer(kind=kint) :: &
 istat

logical(kind=klog) :: &
 opened

logical(kind=klog),save :: &
 sum_file_is_open=.false.

character(10) :: &
 fstatus
!-----------------------------------------------------------------------
real(kind=kdbl),save :: sumdrrw=0.
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      mype=mype_share
!-----------------------------------------------------------------------
!
      steep=1.-0.040*idtadt
!.......................................................................
!$omp parallel 
!$omp do private(l,j,i)
!.......................................................................
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            dvol (i,j,l)=(dsg2(l)*pd(i,j)+pdsg1(l))*dare(j)
            rdvol(i,j,l)=1./dvol(i,j,l)
          enddo
        enddo
      enddo
!.......................................................................
!$omp end do
!.......................................................................
!
!-----------------------------------------------------------------------
!---monotonization------------------------------------------------------
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp do private (ks,l,j,i,s1p,sminh,smaxh,sminv,smaxv,smin, &
!$omp             smax,sn,dsp)
!.......................................................................
!-----------------------------------------------------------------------
      do ks=kss,kse ! loop by species
!-----------------------------------------------------------------------
        do l=1,lm
          xsums(2*ks-1,l)=0.
          xsums(2*ks  ,l)=0.
          gsums(2*ks-1,l)=0.
          gsums(2*ks  ,l)=0.
!
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              s1p=(s1(i,j,l,ks)+tcs(i,j,l,ks))**2
              tcs(i,j,l,ks)=s1p-s(i,j,l,ks)
!
              sminh=min(s(i-1,j-1,l,ks) &
                       ,s(i  ,j-1,l,ks) &
                       ,s(i+1,j-1,l,ks) &
                       ,s(i-1,j  ,l,ks) &
                       ,s(i  ,j  ,l,ks) &
                       ,s(i+1,j  ,l,ks) &
                       ,s(i-1,j+1,l,ks) &
                       ,s(i  ,j+1,l,ks) &
                       ,s(i+1,j+1,l,ks))
              smaxh=max(s(i-1,j-1,l,ks) &
                       ,s(i  ,j-1,l,ks) &
                       ,s(i+1,j-1,l,ks) &
                       ,s(i-1,j  ,l,ks) &
                       ,s(i  ,j  ,l,ks) &
                       ,s(i+1,j  ,l,ks) &
                       ,s(i-1,j+1,l,ks) &
                       ,s(i  ,j+1,l,ks) &
                       ,s(i+1,j+1,l,ks))
!
              if(l.gt.1.and.l.lt.lm) then
                sminv=min(s(i,j,l-1,ks),s(i,j,l  ,ks),s(i,j,l+1,ks))
                smaxv=max(s(i,j,l-1,ks),s(i,j,l  ,ks),s(i,j,l+1,ks))
              elseif(l.eq.1) then
                sminv=min(s(i,j,l  ,ks),s(i,j,l+1,ks))
                smaxv=max(s(i,j,l  ,ks),s(i,j,l+1,ks))
              elseif(l.eq.lm) then
                sminv=min(s(i,j,l-1,ks),s(i,j,l  ,ks))
                smaxv=max(s(i,j,l-1,ks),s(i,j,l  ,ks))
              endif
!
              smin=min(sminh,sminv)
              smax=max(smaxh,smaxv)
!
              sn=s1p
              if(sn.gt.steep*smax) sn=smax
              if(sn.lt.     smin) sn=smin
!
              dsp=(sn-s1p)*dvol(i,j,l)
              s1(i,j,l,ks)=dsp
!
              if(dsp.gt.0.) then
                xsums(2*ks-1,l)=xsums(2*ks-1,l)+dsp
              else
                xsums(2*ks  ,l)=xsums(2*ks  ,l)+dsp
              endif
!
            enddo
          enddo
        enddo
!-----------------------------------------------------------------------
!
      enddo ! end of the loop by species
!.......................................................................
!$omp end do
!$omp end parallel
!.......................................................................
!
!-----------------------------------------------------------------------
!***  Global reductions
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Skip computing the global reduction if they are to be read in
!***  from another run to check bit reproducibility.
!-----------------------------------------------------------------------
      if(.not.read_global_sums)then
        lngth=(2*kse-2*kss+2)*lm
        call mpi_allreduce(xsums,gsums,lngth &
                          ,mpi_double_precision &
                          ,mpi_sum,mpi_comm_comp,irecv)
      endif
!-----------------------------------------------------------------------
!***  For bit reproducibility, read/write global sums.
!-----------------------------------------------------------------------
      bitsaf: if(read_global_sums.or.write_global_sums)then   !<--- NEVER SET BOTH READ AND WRITE TO .TRUE.
!
        if(.not.sum_file_is_open.and.mype==0)then
          open_unit_ad: do l=51,59
            inquire(l,opened=opened)
            if(.not.opened)then
              iunit_advec_sums=l
              if(read_global_sums)fstatus='OLD'
              if(write_global_sums)fstatus='REPLACE'
              open(unit=iunit_advec_sums,file='global_sums',status=fstatus &
                  ,form='UNFORMATTED',iostat=istat)
              sum_file_is_open=.true.
              exit open_unit_ad
            endif
          enddo open_unit_ad
          write(0,*)' mono opened iunit_advec_sums=',iunit_advec_sums
        endif
!
        if(write_global_sums.and.mype==0)then
          do ks=kss,kse
            do l=1,lm
              write(iunit_advec_sums) gsums(2*ks-1,l) &
                                     ,gsums(2*ks  ,l)
            enddo
          enddo
        endif
!
        if(read_global_sums)then
          if(mype==0)then
            do ks=kss,kse
              do l=1,lm
                read (iunit_advec_sums) gsums(2*ks-1,l) &
                                       ,gsums(2*ks  ,l) 
              enddo
            enddo
          endif
!
          call mpi_bcast(gsums,(kse-kss+1)*2*lm &
                        ,mpi_real,0,mpi_comm_comp,ierr)
!
        endif
!
      endif bitsaf
!----------------------------------------------------------------------
!---forced conservation after monotonization----------------------------
!----------------------------------------------------------------------
      do ks=kss,kse
        vgsums(2*ks-1)=0.
        vgsums(2*ks  )=0.
        do l=1,lm
          vgsums(2*ks-1)=gsums(2*ks-1,l)+vgsums(2*ks-1)
          vgsums(2*ks  )=gsums(2*ks  ,l)+vgsums(2*ks  )
        enddo
      enddo
!----------------------------------------------------------------------
!.......................................................................
!$omp parallel do private (ks,l,j,i,sumps,sumns,sfacs,rfacs,dsp)
!.......................................................................
      do ks=kss,kse
        sumps=vgsums(2*ks-1)
        sumns=vgsums(2*ks  )

!jaaif(mype.eq.0) write(0,*) 'sumps,sumns ',sumps,sumns

!
        if(sumps*(-sumns).gt.1.) then
          sfacs=-sumns/sumps
          rfacs=1./sfacs
        else
          sfacs=0.
          rfacs=0.
        endif
!
        do l=1,lm
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              dsp=s1(i,j,l,ks)*rdvol(i,j,l)
              if(sfacs.lt.1.) then
                if(dsp.gt.0.) dsp=dsp*sfacs
              endif
              tcs(i,j,l,ks)=tcs(i,j,l,ks)+dsp
            enddo
          enddo
        enddo
!-----------------------------------------------------------------------
!
      enddo ! end of the loop by species
!.......................................................................
!$omp end parallel do 
!.......................................................................
!
!-----------------------------------------------------------------------
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            tcs(i,j,l,indx_q2)=(dsg2(l)*pd(i,j)+pdsg1(l))*tcs(i,j,l,indx_q2)
          enddo
        enddo
      enddo
      do l=1,lm-1
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            tcs(i,j,l,indx_q2)=(tcs(i,j,l,indx_q2)+tcs(i,j,l+1,indx_q2)) &
                        /((dsg2(l  )*pd(i,j)+pdsg1(l  )) &
                         +(dsg2(l+1)*pd(i,j)+pdsg1(l+1)))
          enddo
        enddo
      enddo
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
          tcs(i,j,lm,indx_q2)=0.
        enddo
      enddo
!-----------------------------------------------------------------------
!zjwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
!-----------------------------------------------------------------------
      do ks=kss,kse ! loop by species
        do l=1,lm
          xsums(2*ks-1,l)=0.
          xsums(2*ks  ,l)=0.
          gsums(2*ks-1,l)=0.
          gsums(2*ks  ,l)=0.
!
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              xsums(2*ks-1,l)=xsums(2*ks-1,l)+s  (i,j,l,ks)*dvol(i,j,l)
              xsums(2*ks  ,l)=xsums(2*ks  ,l)+tcs(i,j,l,ks)*dvol(i,j,l)
            enddo
          enddo
        enddo
      enddo
!
      xsump=0.
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
          xsump=pd(i,j)*dare(j)+xsump
        enddo
      enddo
!-----------------------------------------------------------------------
!***  GLOBAL REDUCTION
!-----------------------------------------------------------------------
      lngth=1
      call mpi_allreduce(xsump,gsump,lngth &
                        ,mpi_double_precision &
                        ,mpi_sum,mpi_comm_comp,irecv)
!-----------------------------------------------------------------------
!***  END OF GLOBAL REDUCTION
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  GLOBAL REDUCTION
!-----------------------------------------------------------------------
      lngth=(2*kse-2*kss+2)*lm
      call mpi_allreduce(xsums,gsums,lngth &
                        ,mpi_double_precision &
                        ,mpi_sum,mpi_comm_comp,irecv)
!-----------------------------------------------------------------------
!***  END OF GLOBAL REDUCTION
!-----------------------------------------------------------------------
      do ks=kss,kse
        vgsums(2*ks-1)=0.
        vgsums(2*ks  )=0.
      enddo
!
      do ks=kss,kse
        do l=1,lm
          vgsums(2*ks-1)=gsums(2*ks-1,l)+vgsums(2*ks-1)
          vgsums(2*ks  )=gsums(2*ks  ,l)+vgsums(2*ks  )
        enddo
      enddo
!
      sumdrrw=vgsums(6)+sumdrrw
!
      if(mype.eq.0) then
        write(0,1000) (vgsums(ks),ks=2*kss-1,2*kse) &
                      ,gsump,sumdrrw
      endif
 1000 format('global vol sums ',10d13.5)
!zjmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
!-----------------------------------------------------------------------
!
                        endsubroutine mono
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine vadv2 &
(lm,idtad &
,dt &
,dsg2,pdsg1,psgml1,sgml2 &
,pd &
,psgdt &
,cw,q,q2,rrw &
!temporary argument passing
,e2)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
logical(kind=klog),parameter:: &
 traditional=.false.         !

integer(kind=kint),parameter:: &
 nsmud=0                     ! number of smoothing iterations

real(kind=kfpt),parameter:: &
 conserve_max=1.5 &          ! max limit on conservation ratios
,conserve_min=0.5 &          ! min limit on conservation ratios
,epsq=1.e-20 &               ! floor value for specific humidity
,epsq2=0.02 &                ! floor value for 2tke
,ff1=0.52500 &               ! antifiltering weighting factor
!ff1=0.50000 &               ! antifiltering weighting factor
,ff2=-0.64813 &              ! antifiltering weighting factor
,ff3=0.24520 &               ! antifiltering weighting factor
,ff4=-0.12189                ! antifiltering weighting factor
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,idtad                       ! timestep factor

real(kind=kfpt),intent(in):: &
 dt                          ! dynamics time step

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigma
,pdsg1 &                     ! delta pressure
,psgml1 &                    ! pressure at midlevels
,sgml2                       ! sigma at midlevels

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd                          ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(in):: &
 psgdt                       ! vertical mass flux

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 cw &                        ! condensate
,q &                         ! specific humidity
,q2 &                        ! 2tke
,rrw                         ! rt/p
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
!
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 e2                          ! 2TKE in the layers
!
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
!
!*** This next group of four arrays is dimensioned with haloes
!*** since they are used as primary arrays.
!
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm) :: &
 e1 &                        ! scratch, 2tke
,g1 &                        ! scratch, rrw
,q1 &                        ! scratch, specific humidity
,w1                          ! scratch, condensate
!
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l &                         ! index in p direction
,lap &                       ! l increment next to departure point
,llap                        ! vertical index next to departure point

real(kind=kfpt):: &
 addt &                      ! dt*idtad, time step
,afrp &                      !
,q1p &                       !
,w1p &                       !
,g1p &                       !
,e1p &                       !
,dqp &                       !
,dwp &                       !
,dgp &                       !
,dep &                       !
,dpdn &                      !
,dpup &                      !
,rdpdn &                     !
,rdpup &                     !
,d2pqq &                     !
,d2pqw &                     !
,d2pqg &                     !
,d2pqe &                     !
,ep &                        !
,e00 &                       !
,ep0 &                       !
,gp &                        !
,g00 &                       !
,gp0 &                       !
,qp &                        !
,q00 &                       !
,qp0 &                       !
,wp &                        !
,w00 &                       !
,wp0 &                       !
,pdsg &                      !
,psgdtp &                    ! vertical mass flux
,rfc &                       !
,rr                          !

logical(kind=klog),dimension(its_b1:ite_b1,jts_b1:jte_b1):: &
 bot                         !

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1):: &
 sface &                     ! scratch, correction factor, 2tke
,sfacg &                     ! scratch, correction factor, rrw
,sfacq &                     ! scratch, correction factor, spec. hum.
,sfacw &                     ! scratch, correction factor, condensate
,sumne &                     ! scratch, sum of negative changes, 2tke
,sumng &                     ! scratch, sum of negative changes, rrw
,sumnq &                     ! scratch, sum of negative changes, spec. hum.
,sumnw &                     ! scratch, sum of negative changes, cond.
,sumpe &                     ! scratch, sum of positive changes, 2tke
,sumpg &                     ! scratch, sum of positive changes, rrw
,sumpq &                     ! scratch, sum of positive changes, spec. hum.
,sumpw                       ! scratch, sum of positive changes, cond.

integer(kind=kint),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 la                          ! vertical index increment, next to departure pt.

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 afr &                       ! antifiltering factor,horizontal
,de &                        ! scratch, 2tke change
,dg &                        ! scratch, rrw change
,dq &                        ! scratch, specific humidity change
,dw                          ! scratch, condensate change
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      addt=dt*float(idtad)
!-----------------------------------------------------------------------
      do j=jts,jte
        do i=its,ite
          e2(i,j,1)=(q2(i,j,1)+epsq2)*0.5
        enddo
      enddo
!jaa      do l=2,lm
!jaa        do j=jts,jte
!jaa          do i=its,ite
!jaa            e2(i,j,l)=(q2(i,j,l-1)+q2(i,j,l))*0.5
!jaa          enddo
!jaa        enddo
!jaa      enddo
!.......................................................................
!$omp parallel 
!$omp do private(l,j,i)
!.......................................................................
      do l=1,lm
       if (l .gt.1) then
        do j=jts,jte
         do i=its,ite
           e2(i,j,l)=(q2(i,j,l-1)+q2(i,j,l))*0.5
         enddo
        enddo
       endif
        do j=jts,jte
          do i=its,ite
            q  (i,j,l)=max(q  (i,j,l),epsq)
            cw (i,j,l)=max(cw (i,j,l),epsq)
            rrw(i,j,l)=max(rrw(i,j,l),epsq)
            e2 (i,j,l)=max(e2 (i,j,l),epsq2)
            q1 (i,j,l)=q  (i,j,l)
            w1 (i,j,l)=cw (i,j,l)
            g1 (i,j,l)=rrw(i,j,l)
            e1 (i,j,l)=e2 (i,j,l)
          enddo
        enddo
      enddo
!.......................................................................
!$omp end do
!.......................................................................
!
!-----------------------------------------------------------------------
!-----------------vertical advection------------------------------------
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp do private (l,j,i,psgdtp,rr,lap,llap)
!.......................................................................
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            if(traditional) then
              if(l.eq.1) then
                psgdtp=psgdt(i,j,1)*0.5
              elseif(l.eq.lm) then
                psgdtp=psgdt(i,j,lm-1)*0.5
              else
                psgdtp=(psgdt(i,j,l-1)+psgdt(i,j,l))*0.5
              endif
            else
              if(l.eq.1) then
                psgdtp=(psgdt(i,j-1,1   )+psgdt(i-1,j,1   ) &
                       +psgdt(i+1,j,1   )+psgdt(i,j+1,1   ) &
                       +psgdt(i ,j ,1   )*4.)*0.0625
              elseif(l.eq.lm) then
                psgdtp=(psgdt(i,j-1,lm-1)+psgdt(i-1,j,lm-1) &
                       +psgdt(i+1,j,lm-1)+psgdt(i,j+1,lm-1) &
                       +psgdt(i ,j ,lm-1)*4.)*0.0625
              else
                psgdtp=(psgdt(i,j-1,l-1 )+psgdt(i-1,j,l-1 ) &
                       +psgdt(i+1,j,l-1 )+psgdt(i,j+1,l-1 ) &
                       +psgdt(i ,j ,l-1 )*4. &
                       +psgdt(i,j-1,l   )+psgdt(i-1,j,l   ) &
                       +psgdt(i+1,j,l   )+psgdt(i,j+1,l   ) &
                       +psgdt(i ,j ,l   )*4.)*0.0625
              endif
            endif
!
            rr=psgdtp*(-addt)
            if(rr.lt.0.) then
              lap=-1
            else
              lap=1
            endif
!
            la(i,j,l)=lap
            llap=l+lap
!
            if(llap.gt.0.and.llap.lt.lm+1) then ! internal and outflow points
              rr=abs(rr &
                /((sgml2(llap)-sgml2(l))*pd(i,j) &
                 +(psgml1(llap)-psgml1(l))))
              if(rr.gt.0.999) rr=0.999
              afr(i,j,l)=(((ff4*rr+ff3)*rr+ff2)*rr+ff1)*rr
              dq(i,j,l)=(q  (i,j,llap)-q  (i,j,l))*rr
              dw(i,j,l)=(cw (i,j,llap)-cw (i,j,l))*rr
              dg(i,j,l)=(rrw(i,j,llap)-rrw(i,j,l))*rr
              de(i,j,l)=(e2 (i,j,llap)-e2 (i,j,l))*rr
            elseif(llap.eq.lm+1) then
              bot(i,j)=.true.
              rr=abs(rr &
                /((1.-sgml2(l))*pd(i,j)))
              if(rr.gt.0.999) rr=0.999
              afr(i,j,l)=0.
!              dq(i,j,l)=-q  (i,j,l)*rr
!              dw(i,j,l)=-cw (i,j,l)*rr
!              dg(i,j,l)=-rrw(i,j,l)*rr
!              de(i,j,l)=-e2 (i,j,l)*rr
              dq(i,j,l)=0.
              dw(i,j,l)=0.
              dg(i,j,l)=0.
              de(i,j,l)=0.
            else
              rr=abs(rr &
                /(pdsg1(l)*0.5))
              if(rr.gt.0.999) rr=0.999
              afr(i,j,l)=0.
!              dq(i,j,l)=-q  (i,j,l)*rr
!              dw(i,j,l)=-cw (i,j,l)*rr
!              dg(i,j,l)=-rrw(i,j,l)*rr
!              de(i,j,l)=-e2 (i,j,l)*rr
              dq(i,j,l)=0.
              dw(i,j,l)=0.
              dg(i,j,l)=0.
              de(i,j,l)=0.
            endif
          enddo
        enddo
 !jaa     enddo
!
!jaa      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            q1(i,j,l)=q  (i,j,l)+dq(i,j,l)
            w1(i,j,l)=cw (i,j,l)+dw(i,j,l)
            g1(i,j,l)=rrw(i,j,l)+dg(i,j,l)
            e1(i,j,l)=e2 (i,j,l)+de(i,j,l)
          enddo
        enddo
      enddo
!.......................................................................
!$omp end do
!$omp end parallel
!.......................................................................
!
!----------------------------------------------------------------------
!--------------anti-filtering diffusion and limiters-------------------
!----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do &
!$omp private (l,j,i,q1p,w1p,g1p,e1p,lap,dpdn,dpup,afrp,rdpdn,rdpup, &
!$omp          d2pqq,d2pqw,d2pqg,d2pqe,qp,wp,gp,ep,q00,qp0,w00,wp0,  &
!$omp          g00,gp0,e00,ep0)
!.......................................................................
      do l=2,lm-1
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            dq(i,j,l)=0.
            dw(i,j,l)=0.
            dg(i,j,l)=0.
            de(i,j,l)=0.
!
            q1p=q1(i,j,l)
            w1p=w1(i,j,l)
            g1p=g1(i,j,l)
            e1p=e1(i,j,l)
!
            lap=la(i,j,l)
!
              dpdn=(sgml2(l+lap)-sgml2(l))*pd(i,j) &
                  +(psgml1(l+lap)-psgml1(l))
              dpup=(sgml2(l)-sgml2(l-lap))*pd(i,j) &
                  +(psgml1(l)-psgml1(l-lap))
!
              rdpdn=1./dpdn
              rdpup=1./dpup
!
              afrp=afr(i,j,l)*(dsg2(l)*pd(i,j)+pdsg1(l))
!
              d2pqq=((q1(i,j,l+lap)-q1p)*rdpdn &
                    -(q1p-q1(i,j,l-lap))*rdpup)*afrp
              d2pqw=((w1(i,j,l+lap)-w1p)*rdpdn &
                    -(w1p-w1(i,j,l-lap))*rdpup)*afrp
              d2pqg=((g1(i,j,l+lap)-g1p)*rdpdn &
                    -(g1p-g1(i,j,l-lap))*rdpup)*afrp
              d2pqe=((e1(i,j,l+lap)-e1p)*rdpdn &
                    -(e1p-e1(i,j,l-lap))*rdpup)*afrp
!
              qp=q1p-d2pqq
              wp=w1p-d2pqw
              gp=g1p-d2pqg
              ep=e1p-d2pqe
!
              q00=q(i,j,l)
              qp0=q(i,j,l+lap)
!
              w00=cw(i,j,l)
              wp0=cw(i,j,l+lap)
!
              g00=rrw(i,j,l)
              gp0=rrw(i,j,l+lap)
!
              e00=e2(i,j,l)
              ep0=e2(i,j,l+lap)
!
              qp=max(qp,min(q00,qp0))
              qp=min(qp,max(q00,qp0))
              wp=max(wp,min(w00,wp0))
              wp=min(wp,max(w00,wp0))
              gp=max(gp,min(g00,gp0))
              gp=min(gp,max(g00,gp0))
              ep=max(ep,min(e00,ep0))
              ep=min(ep,max(e00,ep0))
!
              dq(i,j,l)=qp-q1p
              dw(i,j,l)=wp-w1p
              dg(i,j,l)=gp-g1p
              de(i,j,l)=ep-e1p
!
          enddo
        enddo
      enddo
!.......................................................................
!$omp end parallel do
!.......................................................................
!
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
          dq(i,j,1 )=0.
          dw(i,j,1 )=0.
          dg(i,j,1 )=0.
          de(i,j,1 )=0.
!
          dq(i,j,lm)=0.
          dw(i,j,lm)=0.
          dg(i,j,lm)=0.
          de(i,j,lm)=0.
        enddo
      enddo
!--------------compensate + & - antifiltering changes-------------------
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
          sumpq(i,j)=0.
          sumnq(i,j)=0.
          sumpw(i,j)=0.
          sumnw(i,j)=0.
          sumpg(i,j)=0.
          sumng(i,j)=0.
          sumpe(i,j)=0.
          sumne(i,j)=0.
        enddo
      enddo
      do l=2,lm-1
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            pdsg=(dsg2(l)*pd(i,j)+pdsg1(l))
            if(dq(i,j,l).gt.0.) then
              sumpq(i,j)=dq(i,j,l)*pdsg+sumpq(i,j)
            else
              sumnq(i,j)=dq(i,j,l)*pdsg+sumnq(i,j)
            endif
            if(dw(i,j,l).gt.0.) then
              sumpw(i,j)=dw(i,j,l)*pdsg+sumpw(i,j)
            else
              sumnw(i,j)=dw(i,j,l)*pdsg+sumnw(i,j)
            endif
            if(dg(i,j,l).gt.0.) then
              sumpg(i,j)=dg(i,j,l)*pdsg+sumpg(i,j)
            else
              sumng(i,j)=dg(i,j,l)*pdsg+sumng(i,j)
            endif
            if(de(i,j,l).gt.0.) then
              sumpe(i,j)=de(i,j,l)*pdsg+sumpe(i,j)
            else
              sumne(i,j)=de(i,j,l)*pdsg+sumne(i,j)
            endif
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!--------------first moment conserving factor---------------------------
!-----------------------------------------------------------------------
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
!          if(sumpq(i,j).gt.1.e-9)    then
          if(sumpq(i,j)*(-sumnq(i,j)).gt.0.)    then
            sfacq(i,j)=-sumnq(i,j)/sumpq(i,j)
          else
            sfacq(i,j)=0.
          endif
!          if(sumpw(i,j).gt.1.e-9)    then
          if(sumpw(i,j)*(-sumnw(i,j)).gt.0.)    then
            sfacw(i,j)=-sumnw(i,j)/sumpw(i,j)
          else
            sfacw(i,j)=0.
          endif
!          if(sumpg(i,j).gt.1.e-9)    then
          if(sumpg(i,j)*(-sumng(i,j)).gt.0.)    then
            sfacg(i,j)=-sumng(i,j)/sumpg(i,j)
          else
            sfacg(i,j)=0.
          endif
!          if(sumpe(i,j).gt.1.e-9)    then
          if(sumpe(i,j)*(-sumne(i,j)).gt.0.)    then
            sface(i,j)=-sumne(i,j)/sumpe(i,j)
          else
            sface(i,j)=0.
          endif
!
        enddo
      enddo
!
!.......................................................................
!$omp parallel do &
!$omp private (l,j,i,dqp,dwp,dgp,dep)
!.......................................................................
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
!
            dqp=dq(i,j,l)
            if(sfacq(i,j).gt.0.) then
              if(sfacq(i,j).ge.1.) then
                if(dqp.lt.0.) dqp=dqp/sfacq(i,j)
              else
                if(dqp.gt.0.) dqp=dqp*sfacq(i,j)
              endif
            else
              dqp=0.
            endif
!
            q  (i,j,l)=q1(i,j,l)+dqp
!
            dwp=dw(i,j,l)
            if(sfacw(i,j).gt.0.) then
              if(sfacw(i,j).ge.1.) then
                if(dwp.lt.0.) dwp=dwp/sfacw(i,j)
              else
                if(dwp.gt.0.) dwp=dwp*sfacw(i,j)
              endif
            else
              dwp=0.
            endif
!
            cw (i,j,l)=w1(i,j,l)+dwp
!
            dgp=dg(i,j,l)
            if(sfacg(i,j).gt.0.) then
              if(sfacg(i,j).ge.1.) then
                if(dgp.lt.0.) dgp=dgp/sfacg(i,j)
              else
                if(dgp.gt.0.) dgp=dgp*sfacg(i,j)
              endif
            else
              dgp=0.
            endif
!
            rrw(i,j,l)=g1(i,j,l)+dgp
!
            dep=de(i,j,l)
            if(sface(i,j).gt.0.) then
              if(sface(i,j).ge.1.) then
                if(dep.lt.0.) dep=dep/sface(i,j)
              else
                if(dep.gt.0.) dep=dep*sface(i,j)
              endif
            else
              dep=0.
            endif
!
            e2 (i,j,l)=e1(i,j,l)+dep
!
          enddo
        enddo
!jaa      enddo
!-----------------------------------------------------------------------
!jaa      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            q (i,j,l)=max(q (i,j,l),epsq)
            cw(i,j,l)=max(cw(i,j,l),epsq)
            rrw(i,j,l)=max(rrw(i,j,l),epsq)
            e2(i,j,l)=max(e2(i,j,l),epsq2)
          enddo
        enddo
      enddo
!.......................................................................
!$omp end parallel do
!.......................................................................
!-----------------------------------------------------------------------
!
                        endsubroutine vadv2
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine hadv2 &
(global,ntsd,inpes,jnpes &
,lm,idtad &
,dt,rdyh &
,dsg2,pdsg1,psgml1,sgml2 &
,dare,rdxh &
,pd &
,u,v &
,cw,q,q2,rrw &
!temporary argument passing
,e2)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),parameter:: &
 nsmud=0                     ! number of smoothing iterations

real(kind=kfpt),parameter:: &
 conserve_max=1.1 &          ! max limit on conservation ratios
,conserve_min=0.9 &          ! min limit on conservation ratios
,epsq=1.e-20 &               ! floor value for specific humidity
,epsq2=0.02 &                ! floor value for 2tke
!,ff1=0.52500 &               ! antifiltering weighting factor
,ff1=0.50000 &               ! antifiltering weighting factor
,ff2=-0.64813 &              ! antifiltering weighting factor
,ff3=0.24520 &               ! antifiltering weighting factor
,ff4=-0.12189                ! antifiltering weighting factor
!-----------------------------------------------------------------------
logical(kind=klog),intent(in):: &
 global                      ! global or regional

integer(kind=kint),intent(in):: &
 idtad &                     ! timestep factor
,inpes &                     ! tasks in the x direction
,jnpes &                     ! tasks in y direction
,lm &                        ! total # of levels
,ntsd                        ! timestep

real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,rdyh                        ! 1/deltay, h points

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigma
,pdsg1 &                     ! delta pressure
,psgml1 &                    ! pressure at midlevels
,sgml2                       ! sigma at midlevels

real(kind=kfpt),dimension(jds:jde),intent(in):: &
 dare &                      ! grid box area, h points
,rdxh                        ! 1/dx, h points

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd                          ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 u &                         ! u wind component
,v                           ! v wind component

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 cw &                        ! condensate
,q &                         ! specific humidity
,q2 &                        ! 2tke
,rrw                         ! rt/p
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
!
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 e2                          ! 2TKE in the layers
!
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
!
!*** This next group of four arrays are dimensioned with haloes
!*** since they are used as primary arrays
!
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm) :: &
 e1 &                        ! scratch, 2tke
,g1 &                        ! scratch, rrw
,q1 &                        ! scratch, specific humidity
,w1                          ! scratch, condensate

integer(kind=kint):: &
 i &                         ! index in x direction
,iap &                       !
,irecv &                     !
,j &                         ! index in y direction
,jap &                       !
,l &                         ! index in p direction
,lngth                       !

real(kind=kfpt):: &
 addt &                      ! dt*idtad, time step
,app &                       !
,aqq &                       !
,darep &                     ! grid box area at the point
,dvolp &                     ! grid box volume at the point
,emhp &                      !
,enhp &                      !
,qfc &                       !
,q1p &                       !
,w1p &                       !
,g1p &                       !
,e1p &                       !
,dqp &                       !
,dwp &                       !
,dgp &                       !
,dep &                       !
,d2pqq &                     !
,d2pqw &                     !
,d2pqg &                     !
,d2pqe &                     !
,ep &                        !
,e00 &                       !
,e0q &                       !
,ep0 &                       !
,gp &                        !
,g00 &                       !
,g0q &                       !
,gp0 &                       !
,qp &                        !
,q00 &                       !
,q0q &                       !
,qp0 &                       !
,wp &                        !
,w00 &                       !
,w0q &                       !
,wp0 &                       !
,pp &                        !
,qq &                        !
,rdx &                       !
,rdy &                       !
,sumnel &                    ! sum of negative changes, 2tke
,sumpel &                    ! sum of positive changes, 2tke
,sumngl &                    ! sum of negative changes, rrw
,sumpgl &                    ! sum of positive changes, rrw
,sumnql &                    ! sum of negative changes, spec. hum.
,sumpql &                    ! sum of positive changes, spec. hum.
,sumnwl &                    ! sum of negative changes, condensate
,sumpwl &                    ! sum of positive changes, condensate
,up4 &                       !
,vp4                         !

real(kind=kfpt),dimension(1:lm):: &
 sfacep &                    ! correction factor, 2tke
,sfacgp &                    ! correction factor, rrw
,sfacqp &                    ! correction factor, spec. hum.
,sfacwp                      ! correction factor, condensate

real(kind=kdbl):: &
 xsump,gsump,vgsums

real(kind=kdbl),dimension(1:lm):: &
 xsumr &                     ! sum of neg/pos changes all local fields
,gsumr                       ! sum of neg/pos changes all global fields

real(kind=kdbl),dimension(8,1:lm):: &
 xsums &                     ! sum of neg/pos changes all local fields
,gsums                       ! sum of neg/pos changes all global fields

integer(kind=kint),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 ia &                        ! scratch, i index next to departure point
,ja                          ! scratch, j index next to departure point

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 afp &                       ! scratch, antifiltering weight, x direction
,afq                         ! scratch, antifiltering weight, y direction

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 de &                        ! scratch, 2tke change
,dg &                        ! scratch, rrw change
,dq &                        ! scratch, specific humidity change
,dw &                        ! scratch, condensate change
,dvol                        !

integer(kind=kint) :: ierr,istat
!!!integer(kind=kint),save :: iunit
logical(kind=klog) :: opened
logical(kind=klog),save :: sum_file_is_open=.false.
character(10) :: fstatus
real(kind=kfpt),dimension(8,1:lm) :: gsums_single
!-----------------
!-----------------
integer(kind=kint) :: &
 jstart &
,jstop &
,nth &
,omp_get_num_threads &
,omp_get_thread_num &
,tid
!-----------------
!-----------------
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
      addt=dt*real(idtad)
!
!-----------------------------------------------------------------------
!
      rdy=rdyh
      enhp=-addt*rdyh*0.25
!
!.......................................................................
!$omp parallel do private (i,j,l,rdx,emhp,up4,vp4,pp,qq,iap,app,jap, &
!$omp                      aqq,qfc,dqp,dwp,dgp,dep,darep,dvolp)
!.......................................................................
      do l=1,lm
        do j=jts_h1,jte_h1
          do i=its_h1,ite_h1
            q  (i,j,l)=max(q  (i,j,l),epsq)
            cw (i,j,l)=max(cw (i,j,l),epsq)
            rrw(i,j,l)=max(rrw(i,j,l),epsq)
            e2 (i,j,l)=max(e2 (i,j,l),epsq2)
            q1 (i,j,l)=q  (i,j,l)
            w1 (i,j,l)=cw (i,j,l)
            g1 (i,j,l)=rrw(i,j,l)
            e1 (i,j,l)=e2 (i,j,l)
          enddo
        enddo
!
!-----------------------------------------------------------------------
!--------------horizontal advection-------------------------------------
!-----------------------------------------------------------------------
!
!
        enhp=-addt*rdyh*0.25
        do j=jts_b1,jte_b1
          rdx=rdxh(j)
          emhp=-addt*rdxh(j)*0.25
          do i=its_b1,ite_b1
            up4=u(i-1,j-1,l)+u(i,j-1,l)+u(i-1,j,l)+u(i,j,l)
            vp4=v(i-1,j-1,l)+v(i,j-1,l)+v(i-1,j,l)+v(i,j,l)
!
            pp=up4*emhp
            qq=vp4*enhp
!
            if(pp.lt.0.) then
              iap=-1
            else
              iap= 1
            endif
!
            app=abs(pp)
!
            if(qq.lt.0.) then
              jap=-1
            else
              jap=1
            endif
!
            ia(i,j,l)=iap
            ja(i,j,l)=jap
!
            aqq=abs(qq)
!
            if(app.le.1.) then
              afp(i,j,l)=(((ff4*app+ff3)*app+ff2)*app+ff1)*app
              qfc=pp*qq*0.25
            else
              afp(i,j,l)=0.
!              qfc=0.
              qfc=pp*qq*0.25
            endif
            afq(i,j,l)=(((ff4*aqq+ff3)*aqq+ff2)*aqq+ff1)*aqq
!
            dq(i,j,l)=(q  (i+iap,j,l)-q  (i,j,l))*app &
                     +(q  (i,j+jap,l)-q  (i,j,l))*aqq &
                     +(q  (i+1,j+1,l)-q  (i-1,j+1,l) &
                      -q  (i+1,j-1,l)+q  (i-1,j-1,l))*qfc
            dw(i,j,l)=(cw (i+iap,j,l)-cw (i,j,l))*app &
                     +(cw (i,j+jap,l)-cw (i,j,l))*aqq &
                     +(cw (i+1,j+1,l)-cw (i-1,j+1,l) &
                      -cw (i+1,j-1,l)+cw (i-1,j-1,l))*qfc
            dg(i,j,l)=(rrw(i+iap,j,l)-rrw(i,j,l))*app &
                     +(rrw(i,j+jap,l)-rrw(i,j,l))*aqq &
                     +(rrw(i+1,j+1,l)-rrw(i-1,j+1,l) &
                      -rrw(i+1,j-1,l)+rrw(i-1,j-1,l))*qfc
            de(i,j,l)=(e2 (i+iap,j,l)-e2 (i,j,l))*app &
                     +(e2 (i,j+jap,l)-e2 (i,j,l))*aqq &
                     +(e2 (i+1,j+1,l)-e2 (i-1,j+1,l) &
                      -e2 (i+1,j-1,l)+e2 (i-1,j-1,l))*qfc
!
          enddo
        enddo
!-----------------------------------------------------------------------
        if(global) then
!-----------------------------------------------------------------------
          xsums(1,l)=0.
          xsums(2,l)=0.
          xsums(3,l)=0.
          xsums(4,l)=0.
          xsums(5,l)=0.
          xsums(6,l)=0.
          xsums(7,l)=0.
          xsums(8,l)=0.
!-----------------------------------------------------------------------
          do j=jts_b1,jte_b1
            darep=dare(j)
            do i=its_b1,ite_b1
              dvolp=(dsg2(l)*pd(i,j)+pdsg1(l))*darep
!
              dqp=dq(i,j,l)*dvolp
              dwp=dw(i,j,l)*dvolp
              dgp=dg(i,j,l)*dvolp
              dep=de(i,j,l)*dvolp
!
              if(dqp.gt.0.) then
                xsums(1,l)=xsums(1,l)+dqp
              else
                xsums(2,l)=xsums(2,l)+dqp
              endif
!
              if(dwp.gt.0.) then
                xsums(3,l)=xsums(3,l)+dwp
              else
                xsums(4,l)=xsums(4,l)+dwp
              endif
!
              if(dgp.gt.0.) then
                xsums(5,l)=xsums(5,l)+dgp
              else
                xsums(6,l)=xsums(6,l)+dgp
              endif
!
              if(dep.gt.0.) then
                xsums(7,l)=xsums(7,l)+dep
              else
                xsums(8,l)=xsums(8,l)+dep
              endif
          enddo
        enddo
!-----------------------------------------------------------------------
        endif
!-----------------------------------------------------------------------
!
      enddo
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
      if(global) then
!-----------------------------------------------------------------------
!***  Global reductions
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Skip computing the global reduction if they are to be read in
!***  from another run to check bit reproducibility.
!-----------------------------------------------------------------------
!
        if(.not.read_global_sums)then
          call mpi_allreduce(xsums,gsums,8*lm,mpi_double_precision &
                            ,mpi_sum,mpi_comm_comp,irecv)
!
          do l=1,lm
            gsums_single(1,l)=gsums(1,l)
            gsums_single(2,l)=gsums(2,l)
            gsums_single(3,l)=gsums(3,l)
            gsums_single(4,l)=gsums(4,l)
            gsums_single(5,l)=gsums(5,l)
            gsums_single(6,l)=gsums(6,l)
            gsums_single(7,l)=gsums(7,l)
            gsums_single(8,l)=gsums(8,l)
          enddo
!
        endif
!-----------------------------------------------------------------------
!***  For bit reproducibility, read/write global sums.
!-----------------------------------------------------------------------
!
        bitsad: if(read_global_sums.or.write_global_sums)then   !<--- NEVER SET BOTH READ AND WRITE TO .TRUE.
!!!       if(ntsd==0.and..not.sum_file_is_open)then
          if(.not.sum_file_is_open.and.mype==0)then
            open_unit_ad: do l=51,59
              inquire(l,opened=opened)
              if(.not.opened)then
                iunit_advec_sums=l
                if(read_global_sums)fstatus='OLD'
                if(write_global_sums)fstatus='REPLACE'
                open(unit=iunit_advec_sums,file='global_sums',status=fstatus &
                    ,form='UNFORMATTED',iostat=istat)
                sum_file_is_open=.true.
                exit open_unit_ad
              endif
            enddo open_unit_ad
            write(0,*)' hadv2 opened iunit_advec_sums=',iunit_advec_sums
          endif
!
          if(write_global_sums.and.mype==0)then
            write(0,*)' hadv2 writing to iunit_advec_sums=',iunit_advec_sums
            do l=1,lm
              write(iunit_advec_sums)gsums_single(1,l),gsums_single(2,l) &
                                    ,gsums_single(3,l),gsums_single(4,l) &
                                    ,gsums_single(5,l),gsums_single(6,l) &
                                    ,gsums_single(7,l),gsums_single(8,l)
            enddo
          endif
!
          if(read_global_sums)then
            if(mype==0)then
              do l=1,lm
                read(iunit_advec_sums)gsums_single(1,l),gsums_single(2,l) &
                                     ,gsums_single(3,l),gsums_single(4,l) &
                                     ,gsums_single(5,l),gsums_single(6,l) &
                                     ,gsums_single(7,l),gsums_single(8,l)
              enddo
            endif
!
            call mpi_bcast(gsums_single,8*lm,mpi_real,0,mpi_comm_comp,ierr)
          endif
!
        endif bitsad
!
!-----------------------------------------------------------------------
!--------------first moment conserving factor---------------------------
!-----------------------------------------------------------------------
        do l=1,lm
          sumpql=gsums_single(1,l)
          sumnql=gsums_single(2,l)
          sumpwl=gsums_single(3,l)
          sumnwl=gsums_single(4,l)
          sumpgl=gsums_single(5,l)
          sumngl=gsums_single(6,l)
          sumpel=gsums_single(7,l)
          sumnel=gsums_single(8,l)
!
          if(sumpql*(-sumnql).gt.1.)    then
            sfacqp(l)=-sumnql/sumpql
          else
            sfacqp(l)=0.
          endif
!
          if(sumpwl*(-sumnwl).gt.1.)    then
            sfacwp(l)=-sumnwl/sumpwl
          else
            sfacwp(l)=0.
          endif
!
          if(sumpgl*(-sumngl).gt.1.)    then
            sfacgp(l)=-sumngl/sumpgl
          else
            sfacgp(l)=0.
          endif
!
          if(sumpel*(-sumnel).gt.1.)    then
            sfacep(l)=-sumnel/sumpel
          else
            sfacep(l)=0.
          endif
!
        enddo
!
!-----------------------------------------------------------------------
!
      endif ! global
!
!-----------------------------------------------------------------------
!--------------impose conservation on global advection------------------
!-----------------------------------------------------------------------
      do l=1,lm
!-----------------------------------------------------------------------
        if(global) then
!-----------------------------------------------------------------------
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              dqp=dq(i,j,l)
              if(sfacqp(l).eq.0.) then
                dqp=0.
              elseif(sfacqp(l).ge.1.) then
                if(dqp.lt.0.) dqp=dqp/sfacqp(l)
              else
                if(dqp.gt.0.) dqp=dqp*sfacqp(l)
              endif
              q1(i,j,l)=q(i,j,l)+dqp
!
              dwp=dw(i,j,l)
              if(sfacwp(l).eq.0.) then
                dwp=0.
              elseif(sfacwp(l).ge.1.) then
                if(dwp.lt.0.) dwp=dwp/sfacwp(l)
              else
                if(dwp.gt.0.) dwp=dwp*sfacwp(l)
              endif
              w1(i,j,l)=cw(i,j,l)+dwp
!
              dgp=dg(i,j,l)
              if(sfacgp(l).eq.0.) then
                dgp=0.
              elseif(sfacgp(l).ge.1.) then
                if(dgp.lt.0.) dgp=dgp/sfacgp(l)
              else
                if(dgp.gt.0.) dgp=dgp*sfacgp(l)
              endif
              g1(i,j,l)=rrw(i,j,l)+dgp
!
              dep=de(i,j,l)
              if(sfacep(l).eq.0.) then
                dep=0.
              elseif(sfacep(l).ge.1.) then
                if(dep.lt.0.) dep=dep/sfacep(l)
              else
                if(dep.gt.0.) dep=dep*sfacep(l)
              endif
                e1(i,j,l)=e2(i,j,l)+dep
!
            enddo
          enddo
!
          btim=timef()
          call swaphn(q1(ims,jms,l),ims,ime,jms,jme,1,inpes)
          call swaphn(w1(ims,jms,l),ims,ime,jms,jme,1,inpes)
          call swaphn(g1(ims,jms,l),ims,ime,jms,jme,1,inpes)
          call swaphn(e1(ims,jms,l),ims,ime,jms,jme,1,inpes)
          swaphn_tim=swaphn_tim+timef()-btim
!
          btim=timef()
          call polehn(q1(ims,jms,l),ims,ime,jms,jme,1,inpes,jnpes)
          call polehn(w1(ims,jms,l),ims,ime,jms,jme,1,inpes,jnpes)
          call polehn(g1(ims,jms,l),ims,ime,jms,jme,1,inpes,jnpes)
          call polehn(e1(ims,jms,l),ims,ime,jms,jme,1,inpes,jnpes)
          polehn_tim=polehn_tim+timef()-btim
!-----------------------------------------------------------------------
        else
!-----------------------------------------------------------------------
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              q1(i,j,l)=q  (i,j,l)+dq(i,j,l)
              w1(i,j,l)=cw (i,j,l)+dw(i,j,l)
              g1(i,j,l)=rrw(i,j,l)+dg(i,j,l)
              e1(i,j,l)=e2 (i,j,l)+de(i,j,l)
            enddo
          enddo
!-----------------------------------------------------------------------
        endif !global
!-----------------------------------------------------------------------
      enddo
!-----------------------------------------------------------------------
!
      btim=timef()
      call halo_exch(q1,lm,w1,lm,g1,lm,e1,lm,1,1)
      exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!--------------anti-filtering limiters----------------------------------
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private (q1p,w1p,g1p,e1p,iap,jap,d2pqq,d2pqw,d2pqg,d2pqe, &
!$omp                      qp,wp,gp,ep,q00,qp0,q0q,w00,wp0,w0q,g00,gp0,g0q, &
!$omp                      e00,ep0,e0q,darep,dvolp,dqp,dwp,dgp,dep)
!.......................................................................
       do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            q1p=q1(i,j,l)
            w1p=w1(i,j,l)
            g1p=g1(i,j,l)
            e1p=e1(i,j,l)
!
            iap=ia(i,j,l)
            jap=ja(i,j,l)
!
            d2pqq=(q1(i+1,j,l)+q1(i-1,j,l)-q1p-q1p) &
                 *afp(i,j,l) &
                 +(q1(i,j+1,l)+q1(i,j-1,l)-q1p-q1p) &
                 *afq(i,j,l)
            d2pqw=(w1(i+1,j,l)+w1(i-1,j,l)-w1p-w1p) &
                 *afp(i,j,l) &
                 +(w1(i,j+1,l)+w1(i,j-1,l)-w1p-w1p) &
                 *afq(i,j,l)
            d2pqg=(g1(i+1,j,l)+g1(i-1,j,l)-g1p-g1p) &
                 *afp(i,j,l) &
                 +(g1(i,j+1,l)+g1(i,j-1,l)-g1p-g1p) &
                 *afq(i,j,l)
            d2pqe=(e1(i+1,j,l)+e1(i-1,j,l)-e1p-e1p) &
                 *afp(i,j,l) &
                 +(e1(i,j+1,l)+e1(i,j-1,l)-e1p-e1p) &
                 *afq(i,j,l)
!
            qp=q1p-d2pqq
            wp=w1p-d2pqw
            gp=g1p-d2pqg
            ep=e1p-d2pqe
!
            q00=q(i,j,l)
            qp0=q(i+iap,j,l)
            q0q=q(i,j+jap,l)
!
            w00=cw(i,j,l)
            wp0=cw(i+iap,j,l)
            w0q=cw(i,j+jap,l)
!
            g00=rrw(i,j,l)
            gp0=rrw(i+iap,j,l)
            g0q=rrw(i,j+jap,l)
!
            e00=e2(i,j,l)
            ep0=e2(i+iap,j,l)
            e0q=e2(i,j+jap,l)
!
            qp=max(qp,min(q00,qp0,q0q))
            qp=min(qp,max(q00,qp0,q0q))
            wp=max(wp,min(w00,wp0,w0q))
            wp=min(wp,max(w00,wp0,w0q))
            gp=max(gp,min(g00,gp0,g0q))
            gp=min(gp,max(g00,gp0,g0q))
            ep=max(ep,min(e00,ep0,e0q))
            ep=min(ep,max(e00,ep0,e0q))
!
!            dq(i,j,l)=qp-q00
!            dw(i,j,l)=wp-w00
!            dg(i,j,l)=gp-g00
!            de(i,j,l)=ep-e00
            dq(i,j,l)=qp-q1p
            dw(i,j,l)=wp-w1p
            dg(i,j,l)=gp-g1p
            de(i,j,l)=ep-e1p
          enddo
        enddo
!-----------------------------------------------------------------------
        xsums(1,l)=0.
        xsums(2,l)=0.
        xsums(3,l)=0.
        xsums(4,l)=0.
        xsums(5,l)=0.
        xsums(6,l)=0.
        xsums(7,l)=0.
        xsums(8,l)=0.
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1
          darep=dare(j)
          do i=its_b1,ite_b1
            dvolp=(dsg2(l)*pd(i,j)+pdsg1(l))*darep
            dvol(i,j,l)=dvolp
!
            dqp=dq(i,j,l)*dvolp
            dwp=dw(i,j,l)*dvolp
            dgp=dg(i,j,l)*dvolp
            dep=de(i,j,l)*dvolp
!
            if(dqp.gt.0.) then
              xsums(1,l)=xsums(1,l)+dqp
            else
              xsums(2,l)=xsums(2,l)+dqp
            endif
!
            if(dwp.gt.0.) then
              xsums(3,l)=xsums(3,l)+dwp
            else
              xsums(4,l)=xsums(4,l)+dwp
            endif
!
            if(dgp.gt.0.) then
              xsums(5,l)=xsums(5,l)+dgp
            else
              xsums(6,l)=xsums(6,l)+dgp
            endif
!
            if(dep.gt.0.) then
              xsums(7,l)=xsums(7,l)+dep
            else
              xsums(8,l)=xsums(8,l)+dep
            endif
          enddo
        enddo
      enddo
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  Global reductions
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Skip computing the global reduction if they are to be read in
!***  from another run to check bit reproducibility.
!-----------------------------------------------------------------------
!
      if(.not.read_global_sums)then
        call mpi_allreduce(xsums,gsums,8*lm,mpi_double_precision &
                          ,mpi_sum,mpi_comm_comp,irecv)
!
        do l=1,lm
          gsums_single(1,l)=gsums(1,l)
          gsums_single(2,l)=gsums(2,l)
          gsums_single(3,l)=gsums(3,l)
          gsums_single(4,l)=gsums(4,l)
          gsums_single(5,l)=gsums(5,l)
          gsums_single(6,l)=gsums(6,l)
          gsums_single(7,l)=gsums(7,l)
          gsums_single(8,l)=gsums(8,l)
        enddo
!
      endif
!-----------------------------------------------------------------------
!***  For bit reproducibility, read/write global sums.
!-----------------------------------------------------------------------
!
      bitsaf: if(read_global_sums.or.write_global_sums)then   !<--- NEVER SET BOTH READ AND WRITE TO .TRUE.
!
        if(write_global_sums.and.mype==0)then
          write(0,*)' hadv2 2 writing to iunit_advec_sums=',iunit_advec_sums
          do l=1,lm
            write(iunit_advec_sums)gsums_single(1,l),gsums_single(2,l) &
                                  ,gsums_single(3,l),gsums_single(4,l) &
                                  ,gsums_single(5,l),gsums_single(6,l) &
                                  ,gsums_single(7,l),gsums_single(8,l)
          enddo
        endif
!
        if(read_global_sums)then
          if(mype==0)then
            do l=1,lm
              read(iunit_advec_sums)gsums_single(1,l),gsums_single(2,l) &
                                   ,gsums_single(3,l),gsums_single(4,l) &
                                   ,gsums_single(5,l),gsums_single(6,l) &
                                   ,gsums_single(7,l),gsums_single(8,l)
            enddo
          endif
!
          call mpi_bcast(gsums_single,8*lm,mpi_real,0,mpi_comm_comp,ierr)
!
        endif
!
      endif bitsaf
!
!-----------------------------------------------------------------------
!--------------first moment conserving factor---------------------------
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private(sumpql,sumnql,sumpwl,sumnwl,sumpgl,sumngl,sumpel, &
!$omp                     sumnel,dqp,dwp,dgp,dep)
!.......................................................................
      do l=1,lm
        sumpql=gsums_single(1,l)
        sumnql=gsums_single(2,l)
        sumpwl=gsums_single(3,l)
        sumnwl=gsums_single(4,l)
        sumpgl=gsums_single(5,l)
        sumngl=gsums_single(6,l)
        sumpel=gsums_single(7,l)
        sumnel=gsums_single(8,l)
!
        if(sumpql*(-sumnql).gt.1.)    then
          sfacqp(l)=-sumnql/sumpql
        else
          sfacqp(l)=0.
        endif
!
        if(sumpwl*(-sumnwl).gt.1.)    then
          sfacwp(l)=-sumnwl/sumpwl
        else
          sfacwp(l)=0.
        endif
!
        if(sumpgl*(-sumngl).gt.1.)    then
          sfacgp(l)=-sumngl/sumpgl
        else
          sfacgp(l)=0.
        endif
!
        if(sumpel*(-sumnel).gt.1.)    then
          sfacep(l)=-sumnel/sumpel
        else
          sfacep(l)=0.
        endif
!
!        if(sfacqp(l).lt.conserve_min.or.sfacqp(l).gt.conserve_max) &
!           sfacqp(l)=1.
!        if(sfacwp(l).lt.conserve_min.or.sfacwp(l).gt.conserve_max) &
!           sfacwp(l)=1.
!        if(sfacgp(l).lt.conserve_min.or.sfacgp(l).gt.conserve_max) &
!           sfacgp(l)=1.
!        if(sfacep(l).lt.conserve_min.or.sfacep(l).gt.conserve_max) &
!           sfacep(l)=1.
!
!-----------------------------------------------------------------------
!--------------impose conservation on anti-filtering--------------------
!-----------------------------------------------------------------------
!
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
!
            dqp=dq(i,j,l)
            if(sfacqp(l).eq.0.) then
              dqp=0.
            elseif(sfacqp(l).ge.1.) then
              if(dqp.lt.0.) dqp=dqp/sfacqp(l)
            else
              if(dqp.gt.0.) dqp=dqp*sfacqp(l)
            endif
!            q(i,j,l)=q(i,j,l)+dqp
            q(i,j,l)=q1(i,j,l)+dqp
!
            dwp=dw(i,j,l)
            if(sfacwp(l).eq.0.) then
              dwp=0.
            elseif(sfacwp(l).ge.1.) then
              if(dwp.lt.0.) dwp=dwp/sfacwp(l)
            else
              if(dwp.gt.0.) dwp=dwp*sfacwp(l)
            endif
!            cw(i,j,l)=cw(i,j,l)+dwp
            cw(i,j,l)=w1(i,j,l)+dwp
!
            dgp=dg(i,j,l)
            if(sfacgp(l).eq.0.) then
              dgp=0.
            elseif(sfacgp(l).ge.1.) then
              if(dgp.lt.0.) dgp=dgp/sfacgp(l)
            else
              if(dgp.gt.0.) dgp=dgp*sfacgp(l)
            endif
!            rrw(i,j,l)=rrw(i,j,l)+dgp
            rrw(i,j,l)=g1(i,j,l)+dgp
!
            dep=de(i,j,l)
            if(sfacep(l).eq.0.) then
              dep=0.
            elseif(sfacep(l).ge.1.) then
              if(dep.lt.0.) dep=dep/sfacep(l)
            else
              if(dep.gt.0.) dep=dep*sfacep(l)
            endif
!              e2(i,j,l)=e2(i,j,l)+dep
              e2(i,j,l)=e1(i,j,l)+dep
!
          enddo
        enddo
!
!-----------------------------------------------------------------------
!
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            q (i,j,l)=max(q (i,j,l),epsq)
            cw(i,j,l)=max(cw(i,j,l),epsq)
            rrw(i,j,l)=max(rrw(i,j,l),epsq)
            e2(i,j,l)=max(e2(i,j,l),epsq2)
          enddo
        enddo
!-----------------------------------------------------------------------
!
      enddo
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!-----------------
!-----------------
!.......................................................................
!$omp parallel &
!$omp private(nth,tid,i,j,jstart,jstop)
!.......................................................................
        nth = omp_get_num_threads()
        tid = omp_get_thread_num()
        call looplimits(tid, nth,jts_b1,jte_b1,jstart,jstop) 
!-----------------
!-----------------
      do j=jstart,jstop
        do i=its_b1,ite_b1
          e2(i,j,1)=e2(i,j,1)-(epsq2+q2(i,j,1))*0.5
        enddo
      enddo
!
      do l=2,lm
        do j=jstart,jstop
          do i=its_b1,ite_b1
            e2(i,j,l)=e2(i,j,l)-(q2(i,j,l-1)+q2(i,j,l))*0.5
          enddo
        enddo
      enddo
!
      do l=1,lm-1
        do j=jstart,jstop
          do i=its_b1,ite_b1
            q2(i,j,l)=(e2(i,j,l)*(dsg2(l)*pd(i,j)+pdsg1(l)) &
                      +e2(i,j,l+1)*(dsg2(l+1)*pd(i,j)+pdsg1(l+1)))*0.5 &
                    /((sgml2(l+1)-sgml2(l))*pd(i,j) &
                     +(psgml1(l+1)-psgml1(l))) &
                    +q2(i,j,l)
            q2(i,j,l)=max(q2(i,j,l),epsq2)
          enddo
        enddo
      enddo
!-----------------
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
!-----------------
!zjwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
!-----------------------------------------------------------------------
      do l=1,lm
        xsumr(l)=0.
        gsumr(l)=0.
!
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            xsumr(l)=rrw(i,j,l)*dvol(i,j,l)+xsumr(l)
          enddo
        enddo
      enddo
!
      xsump=0.
      gsump=0.
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
          xsump=pd(i,j)*dare(j)+xsump
        enddo
      enddo
!-----------------------------------------------------------------------
!***  GLOBAL REDUCTION
!-----------------------------------------------------------------------
      lngth=2
      call mpi_allreduce(xsump,gsump,lngth &
                        ,mpi_double_precision &
                        ,mpi_sum,mpi_comm_comp,irecv)
!-----------------------------------------------------------------------
!***  END OF GLOBAL REDUCTION
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  GLOBAL REDUCTION
!-----------------------------------------------------------------------
      lngth=2*lm
      call mpi_allreduce(xsumr,gsumr,lngth &
                        ,mpi_double_precision &
                        ,mpi_sum,mpi_comm_comp,irecv)
!-----------------------------------------------------------------------
!***  END OF GLOBAL REDUCTION
!-----------------------------------------------------------------------
        vgsums=0.
!
        do l=1,lm
          vgsums=gsumr(l)+vgsums
        enddo
!
      if(mype.eq.0) then
        write(0,1000) vgsums,gsump
      endif
 1000 format('global vol sums ',10d13.5)
!zjmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

!-----------------------------------------------------------------------
!
                        endsubroutine hadv2
!
!----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine vadv2_scal &
(lm,idtad &
,dt &
,dsg2,pdsg1,psgml1,sgml2 &
,pd &
,psgdt &
,scal &
,num_scal,indx_start)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
logical(kind=klog),parameter:: &
 traditional=.true.          !

integer(kind=kint),parameter:: &
 nsmud=0                     ! number of smoothing iterations

real(kind=kfpt),parameter:: &
 conserve_max=1.5 &          ! max limit on conservation ratios
,conserve_min=0.5 &          ! min limit on conservation ratios
,epsq=1.e-20 &               ! floor value for scalar
,ff1=0.52500 &               ! antifiltering weighting factor
,ff2=-0.64813 &              ! antifiltering weighting factor
,ff3=0.24520 &               ! antifiltering weighting factor
,ff4=-0.12189                ! antifiltering weighting factor
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,idtad &                     ! timestep factor
,indx_start &                ! starting index in 4-D scalar array
,num_scal                    ! vertical extent of scalar array

real(kind=kfpt),intent(in):: &
 dt                          ! dynamics time step

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigma
,pdsg1 &                     ! delta pressure
,psgml1 &                    ! pressure at midlevels
,sgml2                       ! sigma at midlevels

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd                          ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(in):: &
 psgdt                       ! vertical mass flux

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,1:num_scal),intent(inout):: &
 scal                        ! scalar
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
!
!*** This working array is dimensioned with haloes
!*** since they are used as primary arrays.
!
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm) :: &
 s1                          ! scratch, scalar
!
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l &                         ! index in p direction
,lap &                       ! l increment next to departure point
,llap &                      ! vertical index next to departure point
,n                           ! general index for scalar array

real(kind=kfpt):: &
 addt &                      ! dt*idtad, time step
,afrp &                      !
,s1p &                       !
,dsp &                       !
,dpdn &                      !
,dpup &                      !
,rdpdn &                     !
,rdpup &                     !
,d2pqs &                     !
,sp &                        !
,s00 &                       !
,sp0 &                       !
,pdsg &                      !
,psgdtp &                    ! vertical mass flux
,rfc &                       !
,rr                          !

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1):: &
 sfacs &                     ! scratch, correction factor, scalar
,sumns &                     ! scratch, sum of negative changes, scalar
,sumps                       ! scratch, sum of positive changes, scalar

integer(kind=kint),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 la                          ! vertical index increment, next to departure pt.

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 afr &                       ! antifiltering factor,horizontal
,ds                          ! scratch, scalar change
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      addt=dt*float(idtad)
!-----------------------------------------------------------------------
!
      scalar_loop: do n=indx_start,num_scal
!
!-----------------------------------------------------------------------
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            scal(i,j,l,n)=max(scal(i,j,l,n),epsq)
            s1  (i,j,l)=scal(i,j,l,n)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!-----------------vertical advection------------------------------------
!-----------------------------------------------------------------------
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            if(traditional) then
              if(l.eq.1) then
                psgdtp=psgdt(i,j,1)*0.5
              elseif(l.eq.lm) then
                psgdtp=psgdt(i,j,lm-1)*0.5
              else
                psgdtp=(psgdt(i,j,l-1)+psgdt(i,j,l))*0.5
              endif
            else
              if(l.eq.1) then
                psgdtp=(psgdt(i,j-1,1   )+psgdt(i-1,j,1   ) &
                       +psgdt(i+1,j,1   )+psgdt(i,j+1,1   ) &
                       +psgdt(i ,j ,1   )*4.)*0.0625
              elseif(l.eq.lm) then
                psgdtp=(psgdt(i,j-1,lm-1)+psgdt(i-1,j,lm-1) &
                       +psgdt(i+1,j,lm-1)+psgdt(i,j+1,lm-1) &
                       +psgdt(i ,j ,lm-1)*4.)*0.0625
              else
                psgdtp=(psgdt(i,j-1,l-1 )+psgdt(i-1,j,l-1 ) &
                       +psgdt(i+1,j,l-1 )+psgdt(i,j+1,l-1 ) &
                       +psgdt(i ,j ,l-1 )*4. &
                       +psgdt(i,j-1,l   )+psgdt(i-1,j,l   ) &
                       +psgdt(i+1,j,l   )+psgdt(i,j+1,l   ) &
                       +psgdt(i ,j ,l   )*4.)*0.0625
              endif
            endif
!
            rr=psgdtp*(-addt)
            if(rr.lt.0.) then
              lap=-1
            else
              lap=1
            endif
!
            la(i,j,l)=lap
            llap=l+lap
!
            if(llap.gt.0.and.llap.lt.lm+1) then ! internal and outflow points
              rr=abs(rr &
                /((sgml2(llap)-sgml2(l))*pd(i,j) &
                 +(psgml1(llap)-psgml1(l))))
              if(rr.gt.0.999) rr=0.999
              afr(i,j,l)=(((ff4*rr+ff3)*rr+ff2)*rr+ff1)*rr
              ds(i,j,l)=(scal(i,j,llap,n)-scal(i,j,l,n))*rr
            elseif(llap.eq.lm+1) then
              rr=abs(rr &
                /((1.-sgml2(l))*pd(i,j)))
              if(rr.gt.0.999) rr=0.999
              afr(i,j,l)=0.
!              ds(i,j,l)=-scal(i,j,l,n)*rr
              ds(i,j,l)=0.
!
            else
              rr=abs(rr &
                /(pdsg1(l)*0.5))
              if(rr.gt.0.999) rr=0.999
              afr(i,j,l)=0.
              ds(i,j,l)=-scal(i,j,l,n)*rr
            endif
          enddo
        enddo
      enddo
!
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            s1(i,j,l)=scal(i,j,l,n)+ds(i,j,l)
          enddo
        enddo
      enddo
!----------------------------------------------------------------------
!--------------anti-filtering diffusion and limiters-------------------
!----------------------------------------------------------------------
      do l=2,lm-1
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            ds(i,j,l)=0.
!
            s1p=s1(i,j,l)
!
            lap=la(i,j,l)
!
              dpdn=(sgml2(l+lap)-sgml2(l))*pd(i,j) &
                  +(psgml1(l+lap)-psgml1(l))
              dpup=(sgml2(l)-sgml2(l-lap))*pd(i,j) &
                  +(psgml1(l)-psgml1(l-lap))
!
              rdpdn=1./dpdn
              rdpup=1./dpup
!
              afrp=afr(i,j,l)*(dsg2(l)*pd(i,j)+pdsg1(l))
!
              d2pqs=((s1(i,j,l+lap)-s1p)*rdpdn &
                    -(s1p-s1(i,j,l-lap))*rdpup)*afrp
!
              sp=s1p-d2pqs
!
              s00=scal(i,j,l,n)
              sp0=scal(i,j,l+lap,n)
!
              sp=max(sp,min(s00,sp0))
              sp=min(sp,max(s00,sp0))
!
              ds(i,j,l)=sp-s1p
!
          enddo
        enddo
      enddo
!
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
          ds(i,j,1 )=0.
!
          ds(i,j,lm)=0.
        enddo
      enddo
!--------------compensate + & - antifiltering changes-------------------
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
          sumps(i,j)=0.
          sumns(i,j)=0.
        enddo
      enddo
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            pdsg=(dsg2(l)*pd(i,j)+pdsg1(l))
            if(ds(i,j,l).gt.0.) then
              sumps(i,j)=ds(i,j,l)*pdsg+sumps(i,j)
            else
              sumns(i,j)=ds(i,j,l)*pdsg+sumns(i,j)
            endif
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!--------------first moment conserving factor---------------------------
!-----------------------------------------------------------------------
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
!          if(sumps(i,j).gt.1.e-9)    then
          if(sumps(i,j)*(-sumns(i,j)).gt.0.)    then
            sfacs(i,j)=-sumns(i,j)/sumps(i,j)
          else
            sfacs(i,j)=0.
          endif
!
        enddo
      enddo
!
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
!
            dsp=ds(i,j,l)
            if(sfacs(i,j).gt.0.) then
              if(sfacs(i,j).ge.1.) then
                if(dsp.lt.0.) dsp=dsp/sfacs(i,j)
              else
                if(dsp.gt.0.) dsp=dsp*sfacs(i,j)
              endif
            else
              dsp=0.
            endif
            scal(i,j,l,n)=s1(i,j,l)+dsp
!
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            scal(i,j,l,n)=max(scal(i,j,l,n),epsq)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!
      enddo scalar_loop
!
!-----------------------------------------------------------------------
!
                        endsubroutine vadv2_scal
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine hadv2_scal &
(global,ntsd,inpes,jnpes &
,lm,idtad &
,dt,rdyh &
,dsg2,pdsg1,psgml1,sgml2 &
,dare,rdxh &
,pd &
,u,v &
,scal &
,num_scal,indx_start)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),parameter:: &
 nsmud=0                     ! number of smoothing iterations

real(kind=kfpt),parameter:: &
 conserve_max=1.1 &          ! max limit on conservation ratios
,conserve_min=0.9 &          ! min limit on conservation ratios
,epsq=1.e-20 &               ! floor value for specific humidity
,epsq2=0.02 &                ! floor value for 2tke
!,ff1=0.52500 &               ! antifiltering weighting factor
,ff1=0.50000 &               ! antifiltering weighting factor
,ff2=-0.64813 &              ! antifiltering weighting factor
,ff3=0.24520 &               ! antifiltering weighting factor
,ff4=-0.12189                ! antifiltering weighting factor
!-----------------------------------------------------------------------
logical(kind=klog),intent(in):: &
 global                      ! global or regional

integer(kind=kint),intent(in):: &
 idtad &                     ! timestep factor
,indx_start &                ! starting index for general scalar
,inpes &                     ! tasks in the x direction
,jnpes &                     ! tasks in y direction
,lm &                        ! total # of levels
,ntsd &                      ! timestep
,num_scal                    ! vertical extent of scalar array

real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,rdyh                        ! 1/deltay, h points

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigma
,pdsg1 &                     ! delta pressure
,psgml1 &                    ! pressure at midlevels
,sgml2                       ! sigma at midlevels

real(kind=kfpt),dimension(jds:jde),intent(in):: &
 dare &                      ! grid box area, h points
,rdxh                        ! 1/dx, h points

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd                          ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 u &                         ! u wind component
,v                           ! v wind component

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,num_scal),intent(inout):: &
 scal                        ! general scalar

!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
!
!*** This working array is dimensioned with haloes
!*** since it is used as a primary array
!
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm) :: &
 s1                          ! scratch, scalar

integer(kind=kint):: &
 i &                         ! index in x direction
,iap &                       !
,irecv &                     !
,j &                         ! index in y direction
,jap &                       !
,l &                         ! index in p direction
,n

real(kind=kfpt):: &
 addt &                      ! dt*idtad, time step
,app &                       !
,aqq &                       !
,darep &                     ! grid box area at the point
,dvolp &                     ! grid box volume at the point
,emhp &                      !
,enhp &                      !
,qfc &                       !
,s1p &                       !
,dsp &                       !
,d2pqs &                     !
,sp &                        !
,s00 &                       !
,s0q &                       !
,sp0 &                       !
,pp &                        !
,qq &                        !
,rdx &                       !
,rdy &                       !
,sumnql &                    ! sum of negative changes, spec. hum.
,sumpql &                    ! sum of positive changes, spec. hum.
,sumnsl &                    ! sum of negative changes, scalar
,sumpsl &                    ! sum of positive changes, scalar
,up4 &                       !
,vp4                         !

real(kind=kfpt),dimension(1:lm):: &
 sfacsp                      ! correction factor, scalar

real(kind=kdbl),dimension(2,1:lm):: &
 xsums &                     ! sum of neg/pos changes all local fields
,gsums                       ! sum of neg/pos changes all global fields

integer(kind=kint),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 ia &                        ! scratch, i index next to departure point
,ja                          ! scratch, j index next to departure point

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 afp &                       ! scratch, antifiltering weight, x direction
,afq                         ! scratch, antifiltering weight, y direction

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 ds                          ! scratch, scalar

integer(kind=kint) :: ierr,istat
!!!integer(kind=kint),save :: iunit
logical(kind=klog) :: opened
logical(kind=klog),save :: sum_file_is_open=.false.
character(10) :: fstatus
real(kind=kfpt),dimension(8,1:lm) :: gsums_single
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      mype=mype_share
      addt=dt*real(idtad)
!-----------------------------------------------------------------------
!
      scalar_loop: do n=indx_start,num_scal
!
!-----------------------------------------------------------------------
      do l=1,lm
        do j=jts_h1,jte_h1
          do i=its_h1,ite_h1
            scal(i,j,l,n)=max(scal(i,j,l,n),epsq)
            s1  (i,j,l)=scal(i,j,l,n)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!--------------horizontal advection-------------------------------------
!-----------------------------------------------------------------------
      rdy=rdyh
!
      do l=1,lm
        enhp=-addt*rdyh*0.25
        do j=jts_b1,jte_b1
          rdx=rdxh(j)
          emhp=-addt*rdxh(j)*0.25
          do i=its_b1,ite_b1
            up4=u(i-1,j-1,l)+u(i,j-1,l)+u(i-1,j,l)+u(i,j,l)
            vp4=v(i-1,j-1,l)+v(i,j-1,l)+v(i-1,j,l)+v(i,j,l)
!
            pp=up4*emhp
            qq=vp4*enhp
!
            if(pp.lt.0.) then
              iap=-1
            else
              iap= 1
            endif
!
            app=abs(pp)
!
            if(qq.lt.0.) then
              jap=-1
            else
              jap=1
            endif
!
            ia(i,j,l)=iap
            ja(i,j,l)=jap
!
            aqq=abs(qq)
!
            if(app.le.1.) then
              afp(i,j,l)=(((ff4*app+ff3)*app+ff2)*app+ff1)*app
              qfc=pp*qq*0.25
            else
              afp(i,j,l)=0.
!zjtest              qfc=0.
              qfc=pp*qq*0.25
            endif
            afq(i,j,l)=(((ff4*aqq+ff3)*aqq+ff2)*aqq+ff1)*aqq
!
            ds(i,j,l)=(scal(i+iap,j,l,n)-scal(i,j,l,n))*app &
                     +(scal(i,j+jap,l,n)-scal(i,j,l,n))*aqq &
                     +(scal(i+1,j+1,l,n)-scal(i-1,j+1,l,n) &
                      -scal(i+1,j-1,l,n)+scal(i-1,j-1,l,n))*qfc
!
          enddo
        enddo
!-----------------------------------------------------------------------
        if(global) then
!-----------------------------------------------------------------------
          xsums(1,l)=0.
          xsums(2,l)=0.
!-----------------------------------------------------------------------
          do j=jts_b1,jte_b1
            darep=dare(j)
            do i=its_b1,ite_b1
              dvolp=(dsg2(l)*pd(i,j)+pdsg1(l))*darep
!
              dsp=ds(i,j,l)*dvolp
!
              if(dsp.gt.0.) then
                xsums(1,l)=xsums(1,l)+dsp
              else
                xsums(2,l)=xsums(2,l)+dsp
              endif
!
          enddo
        enddo
!-----------------------------------------------------------------------
        endif
!-----------------------------------------------------------------------
      enddo
!-----------------------------------------------------------------------
      if(global) then
!-----------------------------------------------------------------------
!***  Global reductions
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Skip computing the global reduction if they are to be read in
!***  from another run to check bit reproducibility.
!-----------------------------------------------------------------------
!
        if(.not.read_global_sums)then
          call mpi_allreduce(xsums,gsums,8*lm,mpi_double_precision &
                            ,mpi_sum,mpi_comm_comp,irecv)
!
          do l=1,lm
            gsums_single(1,l)=gsums(1,l)
            gsums_single(2,l)=gsums(2,l)
          enddo
!
        endif
!-----------------------------------------------------------------------
!***  For bit reproducibility, read/write global sums.
!-----------------------------------------------------------------------
!
        bitsad: if(read_global_sums.or.write_global_sums)then   !<--- NEVER SET BOTH READ AND WRITE TO .TRUE.
!!!       if(ntsd==0.and..not.sum_file_is_open)then
          if(.not.sum_file_is_open.and.mype==0)then
            open_unit_ad: do l=51,59
              inquire(l,opened=opened)
              if(.not.opened)then
                iunit_advec_sums=l
                if(read_global_sums)fstatus='OLD'
                if(write_global_sums)fstatus='REPLACE'
                open(unit=iunit_advec_sums,file='global_sums',status=fstatus &
                    ,form='UNFORMATTED',iostat=istat)
                sum_file_is_open=.true.
                exit open_unit_ad
              endif
            enddo open_unit_ad
            write(0,*)' hadv2_scal opened iunit_advec_sums=',iunit_advec_sums
          endif
!
          if(write_global_sums.and.mype==0)then
            write(0,*)' hadv2_scal writing to iunit_advec_sums=',iunit_advec_sums
            do l=1,lm
              write(iunit_advec_sums)gsums_single(1,l),gsums_single(2,l)
            enddo
          endif
!
          if(read_global_sums)then
            if(mype==0)then
              do l=1,lm
                read(iunit_advec_sums)gsums_single(1,l),gsums_single(2,l)
              enddo
            endif
!
            call mpi_bcast(gsums_single,2*lm,mpi_real,0,mpi_comm_comp,ierr)
          endif
!
        endif bitsad
!
!-----------------------------------------------------------------------
!--------------first moment conserving factor---------------------------
!-----------------------------------------------------------------------
        do l=1,lm
          sumpsl=gsums_single(1,l)
          sumnsl=gsums_single(2,l)
!
          if(sumpsl*(-sumnsl).gt.1.)    then
            sfacsp(l)=-sumnsl/sumpsl
          else
            sfacsp(l)=0.
          endif
!
        enddo
!
!-----------------------------------------------------------------------
!
      endif ! global
!
!-----------------------------------------------------------------------
!--------------impose conservation on global advection------------------
!-----------------------------------------------------------------------
      do l=1,lm
!-----------------------------------------------------------------------
        if(global) then
!-----------------------------------------------------------------------
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              dsp=ds(i,j,l)
              if(sfacsp(l).eq.0.) then
                dsp=0.
              elseif(sfacsp(l).ge.1.) then
                if(dsp.lt.0.) dsp=dsp/sfacsp(l)
              else
                if(dsp.gt.0.) dsp=dsp*sfacsp(l)
              endif
              s1(i,j,l)=scal(i,j,l,n)+dsp
!
            enddo
          enddo
!
          btim=timef()
          call swaphn(s1(ims,jms,l),ims,ime,jms,jme,1,inpes)
          swaphn_tim=swaphn_tim+timef()-btim
!
          btim=timef()
          call polehn(s1(ims,jms,l),ims,ime,jms,jme,1,inpes,jnpes)
          polehn_tim=polehn_tim+timef()-btim
!-----------------------------------------------------------------------
        else
!-----------------------------------------------------------------------
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              s1(i,j,l)=scal(i,j,l,n)+ds(i,j,l)
            enddo
          enddo
!-----------------------------------------------------------------------
        endif !global
!-----------------------------------------------------------------------
      enddo
!-----------------------------------------------------------------------
!
      btim=timef()
      call halo_exch(s1,lm,1,1)
      exch_dyn_tim=exch_dyn_tim+timef()-btim
!
!-----------------------------------------------------------------------
!--------------anti-filtering limiters----------------------------------
!-----------------------------------------------------------------------
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            s1p=s1(i,j,l)
!
            iap=ia(i,j,l)
            jap=ja(i,j,l)
!
            d2pqs=(s1(i+1,j,l)+s1(i-1,j,l)-s1p-s1p) &
                 *afp(i,j,l) &
                 +(s1(i,j+1,l)+s1(i,j-1,l)-s1p-s1p) &
                 *afq(i,j,l)
!
            sp=s1p-d2pqs
!
            s00=scal(i,j,l,n)
            sp0=scal(i+iap,j,l,n)
            s0q=scal(i,j+jap,l,n)
!
            sp=max(sp,min(s00,sp0,s0q))
            sp=min(sp,max(s00,sp0,s0q))
!
!            ds(i,j,l)=sp-s00
            ds(i,j,l)=sp-s1p
          enddo
        enddo
!-----------------------------------------------------------------------
        xsums(1,l)=0.
        xsums(2,l)=0.
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1
          darep=dare(j)
          do i=its_b1,ite_b1
            dvolp=(dsg2(l)*pd(i,j)+pdsg1(l))*darep
!
            dsp=ds(i,j,l)*dvolp
!
            if(dsp.gt.0.) then
              xsums(1,l)=xsums(1,l)+dsp
            else
              xsums(2,l)=xsums(2,l)+dsp
            endif
!
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!***  Global reductions
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Skip computing the global reduction if they are to be read in
!***  from another run to check bit reproducibility.
!-----------------------------------------------------------------------
!
      if(.not.read_global_sums)then
        call mpi_allreduce(xsums,gsums,8*lm,mpi_double_precision &
                          ,mpi_sum,mpi_comm_comp,irecv)
!
        do l=1,lm
          gsums_single(1,l)=gsums(1,l)
          gsums_single(2,l)=gsums(2,l)
        enddo
!
      endif
!-----------------------------------------------------------------------
!***  For bit reproducibility, read/write global sums.
!-----------------------------------------------------------------------
!
      bitsaf: if(read_global_sums.or.write_global_sums)then   !<--- NEVER SET BOTH READ AND WRITE TO .TRUE.
!
        if(write_global_sums.and.mype==0)then
          write(0,*)' hadv2_scal 2 writing to iunit_advec_sums=',iunit_advec_sums
          do l=1,lm
            write(iunit_advec_sums)gsums_single(1,l),gsums_single(2,l)
          enddo
        endif
!
        if(read_global_sums)then
          if(mype==0)then
            do l=1,lm
              read(iunit_advec_sums)gsums_single(1,l),gsums_single(2,l)
            enddo
          endif
!
          call mpi_bcast(gsums_single,2*lm,mpi_real,0,mpi_comm_comp,ierr)
!
        endif
!
      endif bitsaf
!
!-----------------------------------------------------------------------
!--------------first moment conserving factor---------------------------
!-----------------------------------------------------------------------
      do l=1,lm
        sumpql=gsums_single(1,l)
        sumnql=gsums_single(2,l)
!
        if(sumpsl*(-sumnsl).gt.1.)    then
          sfacsp(l)=-sumnsl/sumpsl
        else
          sfacsp(l)=0.
        endif
!
!        if(sfacsp(l).lt.conserve_min.or.sfacsp(l).gt.conserve_max) &
!           sfacsp(l)=1.
!
      enddo
!-----------------------------------------------------------------------
!--------------impose conservation on anti-filtering--------------------
!-----------------------------------------------------------------------
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
!
            dsp=ds(i,j,l)
            if(sfacsp(l).eq.0.) then
              dsp=0.
            elseif(sfacsp(l).ge.1.) then
              if(dsp.lt.0.) dsp=dsp/sfacsp(l)
            else
              if(dsp.gt.0.) dsp=dsp*sfacsp(l)
            endif
!           scal(i,j,l,n)=scal(i,j,l)+dsp
            scal(i,j,l,n)=s1(i,j,l)+dsp
!
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            scal(i,j,l,n)=max(scal(i,j,l,n),epsq)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!
      enddo scalar_loop
!
!-----------------------------------------------------------------------
!
                        endsubroutine hadv2_scal
!
!----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
!
                        end module module_dynamics_routines
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
