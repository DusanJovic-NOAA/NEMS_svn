!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
    subroutine post_alctvars(imi,jmi,lmi,mypei,nwtlpes,mpicomp,mygridtype,  &
               mymaptype,post_gribversion,mynsoil,lead_write,jts,jte,jtsgrp,jtegrp)
!
!-----------------------------------------------------------------------
!*** allocate post variables
!-----------------------------------------------------------------------
!
      use vrbls4d
      use vrbls3d
      use vrbls2d
      use soil
      use masks
      use ctlblk_mod
      use params_mod
      use gridspec_mod
      use lookup_mod
!
!-----------------------------------------------------------------------
!
      implicit none
!
      include 'mpif.h'
!
!-----------------------------------------------------------------------
!
      integer,intent(in)            :: imi,jmi,lmi,mypei,nwtlpes,mpicomp
      character(1),intent(in)       :: mygridtype
      character(5),intent(in)       :: post_gribversion
      integer,intent(in)            :: mymaptype,mynsoil
      integer,intent(in)            :: lead_write
      integer,intent(in)            :: jts,jte
      integer,intent(in)            :: jtsgrp(nwtlpes),jtegrp(nwtlpes)
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      integer ii,jj,i,j,l,NPOS_END,NPOS_START,indx_2d,indx_num,nfield
      integer ierr,last_write_task
      integer indx,ll,mype
!
      REAL,allocatable :: SLDPTH2(:),dxh(:),dummy(:,:)
      REAL FACT,tsph,tstart
      REAL RINC(5)
!
!-----------------------------------------------------------------------
!*** get dims from int_state
!-----------------------------------------------------------------------
!
      print *,'in post_alctvars,im=',imi,'jm=',jmi,'lm=',lmi,'grib=',post_gribversion
      im=imi
      jm=jmi
      lm=lmi
      im_jm=im*jm
      lp1=lm + 1
      grib=trim(post_gribversion)
! set ndegr
      if(grib=='grib1') then
        gdsdegr=1000.
      else if (grib=='grib2') then
        gdsdegr=1000000.
      endif
      IOFORM='grib'
      mype=mypei
      me=mype-lead_write
      last_write_task=lead_write+nwtlpes-1
      MPI_COMM_COMP=mpicomp
      num_servers=nwtlpes
      NUM_PROCS=nwtlpes
      NUM_SERVERS=0
      GRIDTYPE=mygridtype
      MAPTYPE=mymaptype
      NSOIL=mynsoil
      print *,'grib=',grib,'ioform=',ioform,'mype=',mype,'me=',me, &
         'lead_write=',lead_write,'last_write_task=',last_write_task, &
         'num_servers=',num_servers,'NUM_PROCS=',NUM_PROCS,'GRIDTYPE=', &
          GRIDTYPE,'maptype=',maptype,'nsoil=',nsoil,'gdsdegr=',gdsdegr
!
!
      allocate(dxh(jm))
      allocate(SLDPTH2(nsoil))
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE ARRAYS OF THE POST.
!-----------------------------------------------------------------------
!
      jsta=jts
      jend=jte
      jsta_m  = jsta
      jsta_m2 = jsta
      jend_m  = jend
      jend_m2 = jend
      if ( mype .eq. lead_write ) then
         jsta_m  = 2
         jsta_m2 = 3
      end if
      if ( mype .eq. last_write_task ) then
         jend_m  = jm - 1
         jend_m2 = jm - 2
      end if
!** neighbors
!jw      iup = me + 1
!jw      idn = me - 1
      iup = mype + 1 - lead_write
      idn = mype - 1 - lead_write
      if ( mype .eq. lead_write ) then
         idn = MPI_PROC_NULL
      end if
      if ( mype .eq. last_write_task ) then
         iup = MPI_PROC_NULL
      end if
      print *,'lead_write_task=',lead_write,'last taks=',last_write_task, &
        'idn=',idn,'iup=',iup,'MPI_PROC_NULL=',MPI_PROC_NULL,'jsta=',jsta,'jend=',jend
!
!     counts, disps for gatherv and scatterv
!
      do i = 1, NUM_PROCS
       icnt(i-1) = (jtegrp(i)-jtsgrp(i)+1)*im
       idsp(i-1) = (jtsgrp(i)-1)*im
       if ( mype .eq. lead_write ) then
           print *, ' i, icnt(i),idsp(i) = ',i-1,icnt(i-1),idsp(i-1)
       end if
      enddo
!
!     extraction limits -- set to two rows
!
      jsta_2l = max(jsta - 2,  1 )
      jend_2u = min(jend + 2, jm )
! special for c-grid v
      jvend_2u = min(jend + 2, jm+1 )
      print *,'im=',im,'jsta_2l=',jsta_2l,'jend_2u=',jend_2u,'lm=',lm
!
!-----------------------------------------------------------------------
!***  3D vars
!-----------------------------------------------------------------------
!
!jw check memory
      allocate(u(im+1,jsta_2l:jend_2u,lm))
      allocate(v(im,jsta_2l:jvend_2u,lm))
      allocate(t(im,jsta_2l:jend_2u,lm))
      allocate(q(im,jsta_2l:jend_2u,lm))
      allocate(uh(im,jsta_2l:jend_2u,lm))
      allocate(vh(im,jsta_2l:jend_2u,lm))
      allocate(wh(im,jsta_2l:jend_2u,lm))
      allocate(pmid(im,jsta_2l:jend_2u,lm))
      allocate(pmidv(im,jsta_2l:jend_2u,lm))
      allocate(pint(im,jsta_2l:jend_2u,lp1))
      allocate(alpint(im,jsta_2l:jend_2u,lp1))
      allocate(zmid(im,jsta_2l:jend_2u,lm))
      allocate(zint(im,jsta_2l:jend_2u,lp1))
      allocate(q2(im,jsta_2l:jend_2u,lm))
      allocate(omga(im,jsta_2l:jend_2u,lm))
      allocate(T_ADJ(im,jsta_2l:jend_2u,lm))
      allocate(ttnd(im,jsta_2l:jend_2u,lm))
      allocate(rswtt(im,jsta_2l:jend_2u,lm))
      allocate(rlwtt(im,jsta_2l:jend_2u,lm))
      allocate(exch_h(im,jsta_2l:jend_2u,lm))
      allocate(train(im,jsta_2l:jend_2u,lm))
      allocate(tcucn(im,jsta_2l:jend_2u,lm))
      allocate(el_myj(im,jsta_2l:jend_2u,lm))
      write(0,*)'in set_varposts, after 3D fields 1'
!     MP FIELD
      allocate(cwm(im,jsta_2l:jend_2u,lm))
      allocate(F_ice(im,jsta_2l:jend_2u,lm))
      allocate(F_rain(im,jsta_2l:jend_2u,lm))
      allocate(F_RimeF(im,jsta_2l:jend_2u,lm))
      allocate(QQW(im,jsta_2l:jend_2u,lm))
      allocate(QQI(im,jsta_2l:jend_2u,lm))
      allocate(QQR(im,jsta_2l:jend_2u,lm))
      allocate(QQS(im,jsta_2l:jend_2u,lm))
      allocate(QQG(im,jsta_2l:jend_2u,lm))
      allocate(EXTCOF55(im,jsta_2l:jend_2u,lm))
      allocate(CFR(im,jsta_2l:jend_2u,lm))
      allocate(DBZ(im,jsta_2l:jend_2u,lm))
      allocate(DBZR(im,jsta_2l:jend_2u,lm))
      allocate(DBZI(im,jsta_2l:jend_2u,lm))
      allocate(DBZC(im,jsta_2l:jend_2u,lm))
      allocate(mcvg(im,jsta_2l:jend_2u,lm))
      allocate(NLICE(im,jsta_2l:jend_2u,lm))

!GFS FIELD
      allocate(o3(im,jsta_2l:jend_2u,lm))
! Add GFS d3d fields
      allocate(vdifftt(im,jsta_2l:jend_2u,lm))
      allocate(tcucns(im,jsta_2l:jend_2u,lm))
      allocate(vdiffmois(im,jsta_2l:jend_2u,lm))
      allocate(dconvmois(im,jsta_2l:jend_2u,lm))
      allocate(sconvmois(im,jsta_2l:jend_2u,lm))
      allocate(nradtt(im,jsta_2l:jend_2u,lm))
      allocate(o3vdiff(im,jsta_2l:jend_2u,lm))
      allocate(o3prod(im,jsta_2l:jend_2u,lm))
      allocate(o3tndy(im,jsta_2l:jend_2u,lm))
      allocate(mwpv(im,jsta_2l:jend_2u,lm))
      allocate(unknown(im,jsta_2l:jend_2u,lm))
      allocate(vdiffzacce(im,jsta_2l:jend_2u,lm))
      allocate(zgdrag(im,jsta_2l:jend_2u,lm))
      allocate(cnvctummixing(im,jsta_2l:jend_2u,lm))
      allocate(vdiffmacce(im,jsta_2l:jend_2u,lm))
      allocate(mgdrag(im,jsta_2l:jend_2u,lm))
      allocate(cnvctvmmixing(im,jsta_2l:jend_2u,lm))
      allocate(ncnvctcfrac(im,jsta_2l:jend_2u,lm))
      allocate(cnvctumflx(im,jsta_2l:jend_2u,lm))
      allocate(cnvctdmflx(im,jsta_2l:jend_2u,lm))
      allocate(cnvctdetmflx(im,jsta_2l:jend_2u,lm))
      allocate(cnvctzgdrag(im,jsta_2l:jend_2u,lm))
      allocate(cnvctmgdrag(im,jsta_2l:jend_2u,lm))
      write(0,*)'in set_varposts, after 3D fields 3'
!
      allocate(htm(im,jsta_2l:jend_2u,lm))
      allocate(vtm(im,jsta_2l:jend_2u,lm))
! add GFIP ICING
      allocate(icing_gfip(im,jsta_2l:jend_2u,lm))
!
!-----------------------------------------------------------------------
!***  soil vars
!-----------------------------------------------------------------------
!
      allocate(stc(im,jsta_2l:jend_2u,nsoil))
      allocate(smc(im,jsta_2l:jend_2u,nsoil))
      allocate(sh2o(im,jsta_2l:jend_2u,nsoil))
      allocate(sldpth(nsoil))
      allocate(sllevel(nsoil))
      allocate(rtdpth(nsoil))
      write(0,*)'in set_varposts, after soil fields'
!
!-----------------------------------------------------------------------
!***  2D vars
!-----------------------------------------------------------------------
! SRD
      allocate(wspd10max(im,jsta_2l:jend_2u))
      allocate(w_up_max(im,jsta_2l:jend_2u))
      allocate(w_dn_max(im,jsta_2l:jend_2u))
      allocate(w_mean(im,jsta_2l:jend_2u))
      allocate(refd_max(im,jsta_2l:jend_2u))
      allocate(up_heli_max(im,jsta_2l:jend_2u))
      allocate(grpl_max(im,jsta_2l:jend_2u))
! SRD
      allocate(u10(im,jsta_2l:jend_2u))
      allocate(v10(im,jsta_2l:jend_2u))
      allocate(tshltr(im,jsta_2l:jend_2u))
      allocate(qshltr(im,jsta_2l:jend_2u))
      allocate(smstav(im,jsta_2l:jend_2u))
      allocate(ssroff(im,jsta_2l:jend_2u))
      allocate(bgroff(im,jsta_2l:jend_2u))
      allocate(vegfrc(im,jsta_2l:jend_2u))
      allocate(acsnow(im,jsta_2l:jend_2u))
      allocate(acsnom(im,jsta_2l:jend_2u))
      allocate(cmc(im,jsta_2l:jend_2u))
      allocate(sst(im,jsta_2l:jend_2u))
      allocate(qz0(im,jsta_2l:jend_2u))
      allocate(thz0(im,jsta_2l:jend_2u))
      allocate(uz0(im,jsta_2l:jend_2u))
      allocate(vz0(im,jsta_2l:jend_2u))
      allocate(qs(im,jsta_2l:jend_2u))
      allocate(ths(im,jsta_2l:jend_2u))
      allocate(sno(im,jsta_2l:jend_2u))
!NAMstart
      allocate(snoavg(im,jsta_2l:jend_2u))
      allocate(psfcavg(im,jsta_2l:jend_2u))
      allocate(t10m(im,jsta_2l:jend_2u))
      allocate(t10avg(im,jsta_2l:jend_2u))
      allocate(akmsavg(im,jsta_2l:jend_2u))
      allocate(akhsavg(im,jsta_2l:jend_2u))
      allocate(u10max(im,jsta_2l:jend_2u))
      allocate(v10max(im,jsta_2l:jend_2u))
!NAMend
      allocate(akms(im,jsta_2l:jend_2u))
      allocate(akhs(im,jsta_2l:jend_2u))
      allocate(cuprec(im,jsta_2l:jend_2u))
      allocate(acprec(im,jsta_2l:jend_2u))
      allocate(ancprc(im,jsta_2l:jend_2u))
      allocate(cuppt(im,jsta_2l:jend_2u))
! GSDstart
      allocate(rainc_bucket(im,jsta_2l:jend_2u))
      allocate(rainnc_bucket(im,jsta_2l:jend_2u))
      allocate(pcp_bucket(im,jsta_2l:jend_2u))
      allocate(snow_bucket(im,jsta_2l:jend_2u))
      allocate(qrmax(im,jsta_2l:jend_2u))
      allocate(tmax(im,jsta_2l:jend_2u))
      allocate(snownc(im,jsta_2l:jend_2u))
      allocate(graupelnc(im,jsta_2l:jend_2u))
! GSDend
      allocate(rswin(im,jsta_2l:jend_2u))
      allocate(rlwin(im,jsta_2l:jend_2u))
      allocate(rlwtoa(im,jsta_2l:jend_2u))
      allocate(tg(im,jsta_2l:jend_2u))
      allocate(sfcshx(im,jsta_2l:jend_2u))
      allocate(sfclhx(im,jsta_2l:jend_2u))
      allocate(fis(im,jsta_2l:jend_2u))
      allocate(t500(im,jsta_2l:jend_2u))
      allocate(cfracl(im,jsta_2l:jend_2u))
      allocate(cfracm(im,jsta_2l:jend_2u))
      allocate(cfrach(im,jsta_2l:jend_2u))
      allocate(acfrst(im,jsta_2l:jend_2u))
      allocate(acfrcv(im,jsta_2l:jend_2u))
      allocate(hbot(im,jsta_2l:jend_2u))
      allocate(htop(im,jsta_2l:jend_2u))
      allocate(aswin(im,jsta_2l:jend_2u))
      allocate(alwin(im,jsta_2l:jend_2u))
      allocate(aswout(im,jsta_2l:jend_2u))
      allocate(alwout(im,jsta_2l:jend_2u))
      allocate(aswtoa(im,jsta_2l:jend_2u))
      allocate(alwtoa(im,jsta_2l:jend_2u))
      allocate(czen(im,jsta_2l:jend_2u))
      allocate(czmean(im,jsta_2l:jend_2u))
      allocate(sigt4(im,jsta_2l:jend_2u))
      allocate(rswout(im,jsta_2l:jend_2u))
      allocate(radot(im,jsta_2l:jend_2u))
      allocate(ncfrst(im,jsta_2l:jend_2u))  ! real
      allocate(ncfrcv(im,jsta_2l:jend_2u))  ! real
      allocate(smstot(im,jsta_2l:jend_2u))
      allocate(pctsno(im,jsta_2l:jend_2u))
      allocate(pshltr(im,jsta_2l:jend_2u))
      allocate(th10(im,jsta_2l:jend_2u))
      allocate(q10(im,jsta_2l:jend_2u))
      allocate(sr(im,jsta_2l:jend_2u))
      allocate(prec(im,jsta_2l:jend_2u))
      allocate(subshx(im,jsta_2l:jend_2u))
      allocate(snopcx(im,jsta_2l:jend_2u))
      allocate(sfcuvx(im,jsta_2l:jend_2u))
      allocate(sfcevp(im,jsta_2l:jend_2u))
      allocate(potevp(im,jsta_2l:jend_2u))
      allocate(z0(im,jsta_2l:jend_2u))
      allocate(ustar(im,jsta_2l:jend_2u))
      allocate(pblh(im,jsta_2l:jend_2u))
      allocate(mixht(im,jsta_2l:jend_2u))
      allocate(twbs(im,jsta_2l:jend_2u))
      allocate(qwbs(im,jsta_2l:jend_2u))
      allocate(sfcexc(im,jsta_2l:jend_2u))
      allocate(grnflx(im,jsta_2l:jend_2u))
      allocate(soiltb(im,jsta_2l:jend_2u))
      allocate(z1000(im,jsta_2l:jend_2u))
      allocate(slp(im,jsta_2l:jend_2u))
      allocate(pslp(im,jsta_2l:jend_2u))
      allocate(f(im,jsta_2l:jend_2u))
      allocate(albedo(im,jsta_2l:jend_2u))
      allocate(albase(im,jsta_2l:jend_2u))
      allocate(cldfra(im,jsta_2l:jend_2u))
      allocate(cprate(im,jsta_2l:jend_2u))
      allocate(cnvcfr(im,jsta_2l:jend_2u))
      allocate(ivgtyp(im,jsta_2l:jend_2u))
      allocate(isltyp(im,jsta_2l:jend_2u))
      allocate(hbotd(im,jsta_2l:jend_2u))
      allocate(htopd(im,jsta_2l:jend_2u))
      allocate(hbots(im,jsta_2l:jend_2u))
      allocate(htops(im,jsta_2l:jend_2u))
      allocate(cldefi(im,jsta_2l:jend_2u))
      allocate(islope(im,jsta_2l:jend_2u))
      allocate(si(im,jsta_2l:jend_2u))
      allocate(lspa(im,jsta_2l:jend_2u))
      allocate(rswinc(im,jsta_2l:jend_2u))
      allocate(vis(im,jsta_2l:jend_2u))
      allocate(pd(im,jsta_2l:jend_2u))
      allocate(mxsnal(im,jsta_2l:jend_2u))
      write(0,*)'in set_varposts, after 2D fields 1'
! add GFS fields
      allocate(sfcux(im,jsta_2l:jend_2u))
      allocate(sfcvx(im,jsta_2l:jend_2u))
      allocate(avgalbedo(im,jsta_2l:jend_2u))
      allocate(avgcprate(im,jsta_2l:jend_2u))
      allocate(avgprec(im,jsta_2l:jend_2u))
      allocate(ptop(im,jsta_2l:jend_2u))
      allocate(pbot(im,jsta_2l:jend_2u))
      allocate(avgcfrach(im,jsta_2l:jend_2u))
      allocate(avgcfracm(im,jsta_2l:jend_2u))
      allocate(avgcfracl(im,jsta_2l:jend_2u))
      allocate(avgtcdc(im,jsta_2l:jend_2u))
      allocate(auvbin(im,jsta_2l:jend_2u))
      allocate(auvbinc(im,jsta_2l:jend_2u))
      allocate(ptopl(im,jsta_2l:jend_2u))
      allocate(pbotl(im,jsta_2l:jend_2u))
      allocate(ttopl(im,jsta_2l:jend_2u))
      allocate(ptopm(im,jsta_2l:jend_2u))
      allocate(pbotm(im,jsta_2l:jend_2u))
      allocate(ttopm(im,jsta_2l:jend_2u))
      allocate(ptoph(im,jsta_2l:jend_2u))
      allocate(pboth(im,jsta_2l:jend_2u))
      allocate(ttoph(im,jsta_2l:jend_2u))
      allocate(sfcugs(im,jsta_2l:jend_2u))
      allocate(sfcvgs(im,jsta_2l:jend_2u))
      allocate(pblcfr(im,jsta_2l:jend_2u))
      allocate(cldwork(im,jsta_2l:jend_2u))
      allocate(gtaux(im,jsta_2l:jend_2u))
      allocate(gtauy(im,jsta_2l:jend_2u))
      allocate(runoff(im,jsta_2l:jend_2u))
      allocate(maxtshltr(im,jsta_2l:jend_2u))
      allocate(mintshltr(im,jsta_2l:jend_2u))
      allocate(maxrhshltr(im,jsta_2l:jend_2u))
      allocate(minrhshltr(im,jsta_2l:jend_2u))
      allocate(dzice(im,jsta_2l:jend_2u))
      allocate(alwinc(im,jsta_2l:jend_2u))
      allocate(alwoutc(im,jsta_2l:jend_2u))
      allocate(alwtoac(im,jsta_2l:jend_2u))
      allocate(aswinc(im,jsta_2l:jend_2u))
      allocate(aswoutc(im,jsta_2l:jend_2u))
      allocate(aswtoac(im,jsta_2l:jend_2u))
      allocate(aswintoa(im,jsta_2l:jend_2u))
      allocate(smcwlt(im,jsta_2l:jend_2u))
      allocate(suntime(im,jsta_2l:jend_2u))
      allocate(fieldcapa(im,jsta_2l:jend_2u))
      allocate(avisbeamswin(im,jsta_2l:jend_2u))
      allocate(avisdiffswin(im,jsta_2l:jend_2u))
      allocate(airbeamswin(im,jsta_2l:jend_2u))
      allocate(airdiffswin(im,jsta_2l:jend_2u))
      allocate(snowfall(im,jsta_2l:jend_2u))

      write(0,*)'in set_varposts, after 2D fields 2'

!-----------------------------------------------------------------------
!*** FROM MASKS
!-----------------------------------------------------------------------
!
      allocate(hbm2(im,jsta_2l:jend_2u))
      allocate(sm(im,jsta_2l:jend_2u))
      allocate(sice(im,jsta_2l:jend_2u))
      allocate(lmh(im,jsta_2l:jend_2u))  ! real
      allocate(lmv(im,jsta_2l:jend_2u))  ! real
      allocate(gdlat(im,jsta_2l:jend_2u))
      allocate(gdlon(im,jsta_2l:jend_2u))
      allocate(dx(im,jsta_2l:jend_2u))
      allocate(dy(im,jsta_2l:jend_2u))
! vrbls4d
      allocate(dust(im,jsta_2l:jend_2u,lm,5))
!output from fixed
!jw      allocate(lsmask(im,jsta_2l:jend_2u))
!jw      allocate(sstsm(im,jsta_2l:jend_2u))
      write(0,*)'in set_varposts, after masks fields '
!
!***
! LMH always = LM for sigma-type vert coord
! LMV always = LM for sigma-type vert coord

       do j = jsta_2l, jend_2u
        do i = 1, im
            LMV ( i, j ) = lm
            LMH ( i, j ) = lm
        end do
       end do
!
! HTM VTM all 1 for sigma-type vert coord

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            HTM ( i, j, l ) = 1.0
            VTM ( i, j, l ) = 1.0
        end do
       end do
      end do
    end subroutine post_alctvars
!
!---------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!---------------------------------------------------------------------
!
  subroutine read_postnmlt(kpo,kth,kpv,po,th,pv)
!
      use ctlblk_mod, only : komax,fileNameD3D,lsm,lsmp1,SPL,SPLDEF,  &
                              lsmdef,ALSL,me
!
      implicit none
!---
      integer :: kpo,kth,kpv
      real,dimension(komax) :: po,th,pv
      namelist/nampgb/kpo,po,kth,th,kpv,pv
      integer l,k,iret
!---------------------------------------------------------------------
!
      print *,'in read_postnmlt'
!
! set default for kpo, kth, th, kpv, pv
      kpo=0
      po=0
      kth=1
      th=(/320.,(0.,k=kth+1,komax)/) ! isentropic level to output
      kpv=8
      pv=(/0.5,-0.5,1.0,-1.0,1.5,-1.5,2.0,-2.0,(0.,k=kpv+1,komax)/)
      read(5,nampgb,iostat=iret,end=118)
 118  continue
      if(me==0)print*,'komax,iret for nampgb= ',komax,iret
      if(me==0)print*,'komax,kpo,kth,th,kpv,pv= ',komax,kpo            &
     &  ,kth,th(1:kth),kpv,pv(1:kpv)
       fileNameD3D='/dev/null'
 119  continue
! set up pressure level from POSTGPVARS or DEFAULT
      if(kpo == 0)then
! use default pressure levels
        print*,'using default pressure levels,spldef=',(spldef(l),l=1,lsmdef)
        lsm=lsmdef
        do l=1,lsm
         spl(l)=spldef(l)
        end do
      else
! use POSTGPVARS
        print*,'using pressure levels from POSTGPVARS'
        lsm=kpo
        if(po(lsm)<po(1))then ! post logic assumes asscending
         do l=1,lsm
          spl(l)=po(lsm-l+1)*100.
         end do
        else
         do l=1,lsm
          spl(l)=po(l)*100.
         end do
        end if
      end if
      print*,'LSM, SPL = ',lsm,spl(1:lsm)
      lsmp1=lsm+1
!
!     COMPUTE DERIVED MAP OUTPUT CONSTANTS.
      DO L = 1,LSM
!jw real4         ALSL(L) = ALOG(SPL(L))
         ALSL(L) = LOG(SPL(L))
      END DO
      write(0,*)' after ALSL'
!
1000  continue

      end subroutine read_postnmlt
