craa********************************************************************
c Starting May 27, 2008: modified by RAA to impose horizontal physical
C     dissipation on top of numerical diffusion in all layers.
c Jan 4, 2007: modified by Rashid Akmaev based on DELDIFS_hyb to do
c horizontal viscosity, thermal conduction, and diffusion of major 
c species (O and O2) with global mean coefficients.
craa********************************************************************
      SUBROUTINE idea_deldifs(RTE,WE,QME,XE,YE,TEME,
     X                   RTO,WO,QMO,XO,YO,TEMO,DELTIM,SL,
     X                   LS_NODE,COEF00,K_LEVEL,hybrid,gen_coord_hybrid,
     $     visc,cond,diff)
!
! Apr 06 2012   Henry Juang, copy from deldif and modify for idea
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def				! hmhj
      use gfs_dyn_deldifs_def
      use gfs_dyn_physcons, rerth => con_rerth
     &                    ,  rd => con_rd, cp => con_cp
      
      IMPLICIT NONE
!
      logical hybrid, gen_coord_hybrid
      REAL(KIND=KIND_EVOD) RTE(LEN_TRIE_LS,2,LEVS,NTRAC)
     &,                     WE(LEN_TRIE_LS,2)
     &,                    QME(LEN_TRIE_LS,2)
     &,                     XE(LEN_TRIE_LS,2)
     &,                     YE(LEN_TRIE_LS,2)
     &,                   TEME(LEN_TRIE_LS,2)
     &,                     PE(LEN_TRIE_LS,2)
!
      REAL(KIND=KIND_EVOD) RTO(LEN_TRIO_LS,2,LEVS,NTRAC)
     &,                     WO(LEN_TRIO_LS,2)
     &,                    QMO(LEN_TRIO_LS,2)
     &,                     XO(LEN_TRIO_LS,2)
     &,                     YO(LEN_TRIO_LS,2)
     &,                   TEMO(LEN_TRIO_LS,2)
     &,                     PO(LEN_TRIO_LS,2)
!
      REAL(KIND=KIND_EVOD) DELTIM, SL(LEVS)
!
      INTEGER              LS_NODE(LS_DIM,3)
!
!CMR  LS_NODE(1,1) ... LS_NODE(LS_MAX_NODE,1) : VALUES OF L
!CMR  LS_NODE(1,2) ... LS_NODE(LS_MAX_NODE,2) : VALUES OF JBASEV
!CMR  LS_NODE(1,3) ... LS_NODE(LS_MAX_NODE,3) : VALUES OF JBASOD
!
      REAL(KIND=KIND_EVOD) COEF00(LEVS,NTRAC)
!
      INTEGER              K_LEVEL
!
      INTEGER              I,IS,IT,JDEL,JDELH,K,KD,KU
craa********************************************************************
c idea change1
c IDEA-related changes
      INTEGER              L,LOCL,N,N0,ND,NP,NPD

      INTEGER              N00

!
      INTEGER              INDEV
      INTEGER              INDOD
      integer              indev1,indev2
      integer              indod1,indod2
      real(kind=kind_evod), parameter :: rkappa = cp / rd
      REAL(KIND=KIND_EVOD) DN1,REALVAL,RTNP,DF_DK,FACT
     &,                    SLRD0,FTRD1,RFACT,RFACTRD,RTRD1,FSHK, tem
c INPUT
c - global mean coefficients of viscosity, conductivity, and diffusion
c
      REAL(KIND=KIND_EVOD),intent(in):: visc(levs),cond(levs),diff(levs)
c
c LOCAL VARIABLES
c Physical coefficients DN*, FACT*, RFACT*, etc., are all different for
C     temperature, dynamics, and tracer transport: dn2, dn3, dn4, fact2,
C     fact3, fact4, etc.
c
      REAL(KIND=KIND_EVOD) dn2,dn3,dn4,fact2,fact3,fact4,
     $     rfact2,rfact3,rfact4
c
c - these arrays are analogous to the original DNE and DNO, kept in 
c     module deldifs_hyb_def, but are different for conduction,
C     viscosity, and diffusion. Although recalculated at every time
c     step, it appears they need to be SAVEd as allocatable arrays
! hmhj all saved are moving to module to have thread safe !!!!
c
!     REAL(KIND=KIND_EVOD),allocatable,save:: dnce(:),dnco(:),dnve(:),
!    $     dnvo(:),dnde(:),dndo(:)
CC
c idea change1 end
craa********************************************************************
!
      REAL(KIND=KIND_EVOD), parameter :: CONS0=0.0, CONS1=1.0, CONS2=2.0
!
      INTEGER              INDLSEV,JBASEV
      INTEGER              INDLSOD,JBASOD
!
      INCLUDE 'function2'
!
!     print *,' enter idea_deldifs ' 				! hmhj
!......................................................................
!
craa********************************************************************
craa********************************************************************
      IF(K_LEVEL.EQ.0) THEN
c Begin initialization
craa********************************************************************
craa********************************************************************
!
        CALL COUNTPERF(0,15,0.)
!!
        allocate(RTRD(LEVS),RTHK(LEVS),sf(levs))
        ALLOCATE ( DNE(LEN_TRIE_LS) )
        ALLOCATE ( DNO(LEN_TRIO_LS) )
        ALLOCATE ( BKLY(levs) )        					! hmhj
        ALLOCATE ( CKLY(levs) )        					! hmhj
craa********************************************************************
c idea change2
c IDEA-related changes
         ALLOCATE (dnce(LEN_TRIE_LS) )
         ALLOCATE (dnco(LEN_TRIO_LS) )
         ALLOCATE (dnve(LEN_TRIE_LS) )
         ALLOCATE (dnvo(LEN_TRIO_LS) )
         ALLOCATE (dnde(LEN_TRIE_LS) )
         ALLOCATE (dndo(LEN_TRIO_LS) )
c idea change2 end
craa********************************************************************
        BKLY(:) = 1.0
        CKLY(:) = 0.0
        if (gen_coord_hybrid) then					! hmhj
          DO  k=1,LEVS							! hmhj
! hmhj ak5, bk5, ck5 in gen_coord_hybrid is the same order as model index
            BKLY(k)=0.5*(bk5(k)+bk5(k+1))				! hmhj
            CKLY(k)=0.5*(ck5(k)+ck5(k+1))*rkappa/thref(k)	        ! hmhj
            if( me.eq.0 )						! hmhj
     &         print*,'sl bkly ckly  in deldif=',k,sl(k),bkly(k),ckly(k)! hmhj
          enddo								! hmhj
        else if (hybrid) then						! hmhj
          DO  k=1,LEVS
! hmhj   sl(k) go bottom to top but bk(k) go top to bottom
            BKLY(k)=0.5*(bk5(levs-k+1)+bk5(levs-k+2))/SL(k)
            if( me.eq.0 ) print*,'sl bkly in deldif=',k,sl(k),bkly(k)
          enddo
        endif
!
        IF(JCAP.GT.170) THEN
!         RECIPROCAL OF TIME SCALE OF DIFFUSION AT REFERENCE WAVENUMBER NP
          RTNP=(JCAP/170.)**4*1.1/3600
          NP=JCAP
          N0=0             ! MAXIMUM WAVENUMBER FOR ZERO DIFFUSION
          JDEL=8           ! ORDER OF DIFFUSION (EVEN POWER TO RAISE DEL)
          FSHK=2.2         ! EXTRA HEIGHT-DEPENDENT DIFFUSION FACTOR PER SCALE HEIGHT
        ELSEIF(JCAP.EQ.170) THEN
!         RECIPROCAL OF TIME SCALE OF DIFFUSION AT REFERENCE WAVENUMBER NP
          RTNP=4*3.E15/(RERTH**4)*FLOAT(80*81)**2
          NP=JCAP
          N0=0.55*JCAP     ! MAXIMUM WAVENUMBER FOR ZERO DIFFUSION
          JDEL=2           ! ORDER OF DIFFUSION (EVEN POWER TO RAISE DEL)
          FSHK=1.0         ! EXTRA HEIGHT-DEPENDENT DIFFUSION FACTOR PER SCALE HEIGHT
        ELSEIF(JCAP.EQ.126) THEN					! hmhj
!         BELOW HAS BEEN TESTED IN SIHMA-THETA FOR 2 YEAR CFS RUN	! hmhj
          RTNP=4*3.E15/(RERTH**4)*FLOAT(80*81)**2			! hmhj
          NP=JCAP							! hmhj
          N0=0.0           						! hmhj
          JDEL=4           						! hmhj
          FSHK=1.0         						! hmhj
        ELSE
!         RECIPROCAL OF TIME SCALE OF DIFFUSION AT REFERENCE WAVENUMBER NP
          RTNP=1*3.E15/(RERTH**4)*FLOAT(80*81)**2
          NP=JCAP
c idea change3
          N0=0.55*JCAP     ! MAXIMUM WAVENUMBER FOR ZERO DIFFUSION
          JDEL=2           ! ORDER OF DIFFUSION (EVEN POWER TO RAISE DEL)
          FSHK=1.0         ! EXTRA HEIGHT-DEPENDENT DIFFUSION FACTOR PER SCALE HEIGHT
c         N0=0.     ! MAXIMUM WAVENUMBER FOR ZERO DIFFUSION
c         JDEL=4           ! ORDER OF DIFFUSION (EVEN POWER TO RAISE DEL)
c         FSHK=1.5         ! EXTRA HEIGHT-DEPENDENT DIFFUSION FACTOR PER SCALE HEIGHT
        ENDIF
!
        N00=N0
!
        SLRD0=0.002     ! SIGMA LEVEL AT WHICH TO BEGIN RAYLEIGH DAMPING
!       RTRD1=1./(10.*86400.) ! RECIPROCAL OF TIME SCALE PER SCALE HEIGHT
        RTRD1=0.
!       RTRD1=1./(5.*86400.) ! RECIPROCAL OF TIME SCALE PER SCALE HEIGHT
!       RTRD1=1./(2.*86400) ! RECIPROCAL OF TIME SCALE PER SCALE HEIGHT
                    !  ABOVE BEGINNING SIGMA LEVEL FOR RAYLEIGH DAMPING

        IF (ME.EQ.0) THEN
          PRINT 6,RTNP,NP,N0,JDEL
    6     FORMAT(' HORIZONTAL DIFFUSION PARAMETERS'/
     &  '   EFFECTIVE ',6PF10.3,' MICROHERTZ AT WAVENUMBER ',I4/
     &  '   MAXIMUM WAVENUMBER FOR ZERO DIFFUSION ',I4/
     &  '   ORDER OF DIFFUSION ',I2)

         print *, '***IDEA*** Using physical diffusion in all layers'
         print *,JCAP,N0,FSHK,rtrd1

        ENDIF
c idea change3 end
craa********************************************************************
!
        DO K=1,LEVS
          IF(SL(K).LT.SLRD0) THEN
            if (k .gt. levr) then
! idea
              RTRD(K)=RTRD1*LOG(SLRD0/SL(K)) ** 2
!             RTRD(K)=RTRD1*LOG(SLRD0/SL(K)) ** 3
            else
              RTRD(K)=RTRD1*LOG(SLRD0/SL(K))
            endif
c idea
c             RTRD(K)=min(1.e-5,RTRD(K))
              RTRD(K)=min(RTRD1,RTRD(K))
          ELSE
            RTRD(K)=0
          ENDIF
          RTHK(K)=(SL(K))**LOG(1/FSHK)
        ENDDO
!
        JDELH=JDEL/2
        NPD=MAX(NP-N0,0)
        REALVAL=NPD*(NPD+1)
        DN1=CONS2*RTNP/REALVAL**JDELH
!
!......................................................................
!
        DO LOCL=1,LS_MAX_NODE
               L=LS_NODE(LOCL,1)
          JBASEV=LS_NODE(LOCL,2)
          INDEV=INDLSEV(L,L)
          DO N=L,JCAP,2
            ND=MAX(N-N0,0)
            REALVAL=ND*(ND+1)
            DNE(INDEV)=DN1*REALVAL**JDELH
            INDEV=INDEV+1
          ENDDO
        ENDDO
!
!......................................................................
!
        DO LOCL=1,LS_MAX_NODE
               L=LS_NODE(LOCL,1)
          JBASEV=LS_NODE(LOCL,2)
          if (mod(L,2).eq.mod(jcap+1,2)) then
            DNE(INDLSEV(JCAP+1,L))=CONS0 ! SET THE EVEN (N-L) TERMS OF THE TOP ROW TO ZERO
          ENDIF
        ENDDO
!
!......................................................................
!
        DO LOCL=1,LS_MAX_NODE
               L=LS_NODE(LOCL,1)
          JBASOD=LS_NODE(LOCL,3)
          INDOD=INDLSOD(L+1,L)
          DO N=L+1,JCAP,2
            ND=MAX(N-N0,0)
            REALVAL=ND*(ND+1)
            DNO(INDOD)=DN1*REALVAL**JDELH
            INDOD=INDOD+1
          ENDDO
        ENDDO
!
!......................................................................
!
        DO LOCL=1,LS_MAX_NODE
               L=LS_NODE(LOCL,1)
          JBASOD=LS_NODE(LOCL,3)
          if (mod(L,2).ne.mod(jcap+1,2)) then
            DNO(INDLSOD(JCAP+1,L))=CONS0 ! SET THE ODD (N-L) TERMS OF THE TOP ROW TO ZERO
          ENDIF
        ENDDO
!
!......................................................................
!
        DO K=1,LEVS
          KD=MAX(K-1,1)
          KU=MIN(K+1,LEVS)
          SF(K)=SL(K)/(SL(KU)-SL(KD))/SQRT(CONS2)     !CONSTANT
        ENDDO
!
        CALL COUNTPERF(1,15,0.)
!!
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     if(me.eq.0) then 
c        realval=real((JCAP-N00)*(JCAP-N00+1))
c        fact=real(JCAP*(JCAP+1))
c        print '(a,i5,5es11.3)','www5',k,rtrd(k),dn1*rthk(k)*realval,
c    $        dn2*fact,dn3*fact,dn4*fact
c     endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        RETURN
craa********************************************************************
craa********************************************************************
c End initialization
      ENDIF ! K_LEVEL=00
craa********************************************************************
craa********************************************************************
!
!......................................................................
!
      CALL COUNTPERF(0,13,0.)
!!
      K=K_LEVEL
!
!     TEM = COEF00(K,1) * BKLY(K)
!!
craa********************************************************************
c idea change4
c
c Global mean horizontal transport of heat, momentum and tracers is
c incorporated in all layers based on following assumptions
c     1) physical diffusion acts on all scales (N0=0) and
c     2) it is quadratic (JDEL=2)
c     3) diffusion coefficients for all tracers are the same as for
c       pairs O-N2 and O-O2, which are equal (D12=D13), which should be
c       fine for all tracers as of May 2008, because those different
c       from these three are not present where horizontal diffusion is 
c       important
c
c Precalculate arrays of coefficients dnce, dnve, dnde, etc., analogous
c to the coefficients dne and dno for numerical diffusion. Has to be 
c done every timestep. The only reason these coefficients are 
c precalculated is to keep the original structure of this subroutine.
c Otherwise all the relevant calculations could have been done 
c compactly and probably more efficiently in one loop.
c
c First calculate scalar coefficients (for a given vertical index k)
c dn* for conductivity, viscosity, and diffusion 
c
      dn2=cons2*cond(k)/(rerth**2)
      dn3=cons2*visc(k)/(rerth**2)
      dn4=cons2*diff(k)/(rerth**2)
c
c (I think cons2 is here because the time step is 2.*deltim)
c
c -Evens
c
      DO LOCL=1,LS_MAX_NODE
              L=LS_NODE(LOCL,1)
         JBASEV=LS_NODE(LOCL,2)
         INDEV=INDLSEV(L,L)
         DO N=L,JCAP,2
            REALVAL=real(N*(N+1))
            DNcE(INDEV)=DN2*REALVAL
            DNvE(INDEV)=DN3*REALVAL
            DNdE(INDEV)=DN4*REALVAL
            INDEV=INDEV+1
         ENDDO
      ENDDO
c     print*,'ok1'
c     stop
CC
      DO LOCL=1,LS_MAX_NODE
              L=LS_NODE(LOCL,1)
         JBASEV=LS_NODE(LOCL,2)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            DNcE(INDLSEV(JCAP+1,L))=CONS0     !CONSTANT
            DNvE(INDLSEV(JCAP+1,L))=CONS0     !CONSTANT
            DNdE(INDLSEV(JCAP+1,L))=CONS0     !CONSTANT
         ENDIF
      ENDDO
c     print*,'ok2'
c
c -Odds
c
      DO LOCL=1,LS_MAX_NODE
              L=LS_NODE(LOCL,1)
         JBASOD=LS_NODE(LOCL,3)
         INDOD=INDLSOD(L+1,L)
         DO N=L+1,JCAP,2
            REALVAL=real(N*(N+1))
            DNcO(INDOD)=DN2*REALVAL
            DNvO(INDOD)=DN3*REALVAL
            DNdO(INDOD)=DN4*REALVAL
            INDOD=INDOD+1
         ENDDO
      ENDDO
c     print*,'ok3'
CC
      DO LOCL=1,LS_MAX_NODE
              L=LS_NODE(LOCL,1)
         JBASOD=LS_NODE(LOCL,3)
         if (mod(L,2).ne.mod(jcap+1,2)) then
            DNcO(INDLSOD(JCAP+1,L))=CONS0     !CONSTANT
            DNvO(INDLSOD(JCAP+1,L))=CONS0     !CONSTANT
            DNdO(INDLSOD(JCAP+1,L))=CONS0     !CONSTANT
         ENDIF
      ENDDO
c     print*,'ok4'
CC
CC......................................................................
CC
c The following do-loops originally did only numerical diffusion and
c Rayleigh friction. Now physical diffusion is inserted as well.
c
c -Evens
c
      DO LOCL=1,LS_MAX_NODE
              L=LS_NODE(LOCL,1)
         JBASEV=LS_NODE(LOCL,2)
         IF (L.EQ.0) THEN
            N0=2
         ELSE
            N0=L
         ENDIF
         indev1 = indlsev(N0,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
           indev2 = indlsev(jcap+1,L)
         else
           indev2 = indlsev(jcap  ,L)
         endif
!!       DO N = N0, JCAP+1, 2
         DO INDEV = indev1 , indev2
c
c fact* and rfact* are for conductivity, viscosity, and diffusion
c
            FACT             = DELTIM*DNE(INDEV)*RTHK(K)
            fact2=deltim*dnce(indev)            
            fact3=deltim*dnve(indev)
            fact4=deltim*dnde(indev)
c            RFACT            = CONS1/(CONS1+FACT)
c            RFACTRD          = CONS1/(CONS1+FACT+DELTIM*RTRD(K))
            rfact2  =CONS1/(CONS1+FACT+fact2) 
            rfact3  =CONS1/(CONS1+FACT+DELTIM*RTRD(K)+fact3) 
            rfact4  =CONS1/(CONS1+FACT+fact4) 
            
c            WE(INDEV,1)      = WE(INDEV,1)*RFACTRD
c            WE(INDEV,2)      = WE(INDEV,2)*RFACTRD
c            XE(INDEV,1)      = XE(INDEV,1)*RFACTRD
c            XE(INDEV,2)      = XE(INDEV,2)*RFACTRD
            WE(INDEV,1)      = WE(INDEV,1)*rfact3
            WE(INDEV,2)      = WE(INDEV,2)*rfact3
            XE(INDEV,1)      = XE(INDEV,1)*rfact3
            XE(INDEV,2)      = XE(INDEV,2)*rfact3
            
c            RTE(INDEV,1,1,1) = RTE(INDEV,1,1,1)*RFACT
c            RTE(INDEV,2,1,1) = RTE(INDEV,2,1,1)*RFACT
            RTE(INDEV,1,1,1) = RTE(INDEV,1,1,1)*rfact4
            RTE(INDEV,2,1,1) = RTE(INDEV,2,1,1)*rfact4

            PE(INDEV,1)=BKLY(K)*QME(INDEV,1)+CKLY(K)*TEME(INDEV,1) ! hmhj
            PE(INDEV,2)=BKLY(K)*QME(INDEV,2)+CKLY(K)*TEME(INDEV,2) ! hmhj
c            YE(INDEV,1)      =  ( YE(INDEV,1)+
c     X           FACT*COEF00(K,1)* PE(INDEV,1) )*RFACT ! hmhj
c            YE(INDEV,2)      = ( YE(INDEV,2) +
c     X           FACT*COEF00(K,1)* PE(INDEV,2) )*RFACT ! hmhj
            YE(INDEV,1)      =  ( YE(INDEV,1)+
     X           (FACT+fact2)*COEF00(K,1)* PE(INDEV,1) )*rfact2
            YE(INDEV,2)      = ( YE(INDEV,2) +
     X           (FACT+fact2)*COEF00(K,1)* PE(INDEV,2) )*rfact2
 
         ENDDO
      ENDDO
c     print*,'ok5'
c
c -Odds
c
      DO LOCL=1,LS_MAX_NODE
              L=LS_NODE(LOCL,1)
         JBASOD=LS_NODE(LOCL,3)
         indod1 = indlsod(L+1,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indod2 = indlsod(jcap  ,L)
         else
            indod2 = indlsod(jcap+1,L)
         endif
!        DO N = L+1, JCAP+1, 2
         DO INDOD = indod1 , indod2
 
            FACT             = DELTIM*DNO(INDOD)*RTHK(K)
            fact2=deltim*dnco(indod)
            fact3=deltim*dnvo(indod)
            fact4=deltim*dndo(indod)
c            RFACT            = CONS1/(CONS1+FACT)
c            RFACTRD          = CONS1/(CONS1+FACT+DELTIM*RTRD(K))
            rfact2  =CONS1/(CONS1+FACT+fact2)
            rfact3  =CONS1/(CONS1+FACT+DELTIM*RTRD(K)+fact3)
            rfact4  =CONS1/(CONS1+FACT+fact4)
            
c            WO(INDOD,1)      = WO(INDOD,1)*RFACTRD
c            WO(INDOD,2)      = WO(INDOD,2)*RFACTRD
c            XO(INDOD,1)      = XO(INDOD,1)*RFACTRD
c            XO(INDOD,2)      = XO(INDOD,2)*RFACTRD
            WO(INDOD,1)      = WO(INDOD,1)*rfact3
            WO(INDOD,2)      = WO(INDOD,2)*rfact3
            XO(INDOD,1)      = XO(INDOD,1)*rfact3
            XO(INDOD,2)      = XO(INDOD,2)*rfact3
            
c            RTO(INDOD,1,1,1) = RTO(INDOD,1,1,1)*RFACT
c            RTO(INDOD,2,1,1) = RTO(INDOD,2,1,1)*RFACT
            RTO(INDOD,1,1,1) = RTO(INDOD,1,1,1)*rfact4
            RTO(INDOD,2,1,1) = RTO(INDOD,2,1,1)*rfact4
            
            PO(INDOD,1)=BKLY(K)*QMO(INDOD,1)+CKLY(K)*TEMO(INDOD,1) ! hmhj
            PO(INDOD,2)=BKLY(K)*QMO(INDOD,2)+CKLY(K)*TEMO(INDOD,2) ! hmhj
c            YO(INDOD,1)      = ( YO(INDOD,1)+
c     X           FACT*COEF00(K,1)* PO(INDOD,1) )*RFACT ! hmhj
c            YO(INDOD,2)      = ( YO(INDOD,2)+
c     X           FACT*COEF00(K,1)* PO(INDOD,2) )*RFACT ! hmhj
            YO(INDOD,1)      = ( YO(INDOD,1)+
     X           (FACT+fact2)*COEF00(K,1)* PO(INDOD,1) )*rfact2
            YO(INDOD,2)      = ( YO(INDOD,2)+
     X           (FACT+fact2)*COEF00(K,1)* PO(INDOD,2) )*rfact2
            
         ENDDO
      ENDDO
c     print*,'ok6'
!
!......................................................................
!
c Now do tracers 2-ntrac. Apply the same diffusion coefficient to all.
c This version also does correction for these tracers similar to the one
c done for enthalpy above. The correction may be gone in the future.
c
      if (ntrac .gt. 1) then
         do it=2,ntrac
c
c -Evens
c
            DO LOCL=1,LS_MAX_NODE
                  L=LS_NODE(LOCL,1)
               JBASEV=LS_NODE(LOCL,2)
               IF (L.EQ.0) THEN
                  N0=2
               ELSE
                  N0=L
               ENDIF
               indev1 = indlsev(N0,L)
               if (mod(L,2).eq.mod(jcap+1,2)) then
                  indev2 = indlsev(jcap+1,L)
               else
                  indev2 = indlsev(jcap  ,L)
               endif
               DO INDEV = indev1 , indev2
 
                  FACT              = DELTIM*DNE(INDEV)*RTHK(K)
                  fact4=deltim*dnde(indev)
c                  RFACT             = CONS1/(CONS1+FACT)
                  rfact4  =CONS1/(CONS1+FACT+fact4)
 
                  PE(INDEV,1)=BKLY(K)*QME(INDEV,1)+CKLY(K)*TEME(INDEV,1) ! hmhj
                  PE(INDEV,2)=BKLY(K)*QME(INDEV,2)+CKLY(K)*TEME(INDEV,2) ! hmhj

c                  RTE(INDEV,1,1,IT) = ( RTE(INDEV,1,1,IT)+
c     X                 FACT*COEF00(K,it)* PE(INDEV,1) )*RFACT ! hmhj
c                  RTE(INDEV,2,1,IT) = ( RTE(INDEV,2,1,IT)+
c     X                 FACT*COEF00(K,it)* PE(INDEV,2) )*RFACT ! hmhj
                  RTE(INDEV,1,1,IT) = ( RTE(INDEV,1,1,IT)+
     X                 (FACT+fact4)*COEF00(K,it)*PE(INDEV,1))*rfact4
                  RTE(INDEV,2,1,IT) = ( RTE(INDEV,2,1,IT)+
     X                 (FACT+fact4)*COEF00(K,it)* PE(INDEV,2))*rfact4
 
               ENDDO
            ENDDO
c     print*,'ok7'
!
!......................................................................
!
c
c -Odds
c
            DO LOCL=1,LS_MAX_NODE
                  L=LS_NODE(LOCL,1)
               JBASOD=LS_NODE(LOCL,3)
               indod1 = indlsod(L+1,L)
               if (mod(L,2).eq.mod(jcap+1,2)) then
                  indod2 = indlsod(jcap  ,L)
               else
                  indod2 = indlsod(jcap+1,L)
               endif

c                 if(me.eq.0)print*,'www7',indod1,indod2,K,LEN_TRIO_LS
               DO INDOD = indod1 , indod2
 
c                 if(me.eq.0)print*,'www7',INDOD,K
                  FACT              = DELTIM*DNO(INDOD)*RTHK(K)
c                 fact4=deltim*dndo(indev)
                  fact4=deltim*dndo(indod)
c                  RFACT             = CONS1/(CONS1+FACT)
                  rfact4  =CONS1/(CONS1+FACT+fact4)
 
                  PO(INDOD,1)=BKLY(K)*QMO(INDOD,1)+CKLY(K)*TEMO(INDOD,1) ! hmhj
                  PO(INDOD,2)=BKLY(K)*QMO(INDOD,2)+CKLY(K)*TEMO(INDOD,2) ! hmhj

c                  RTO(INDOD,1,1,IT) = ( RTO(INDOD,1,1,IT)+
c     X                 FACT*COEF00(K,it)* PO(INDOD,1) )*RFACT ! hmhj
c                  RTO(INDOD,2,1,IT) = ( RTO(INDOD,2,1,IT)+
c     X                 FACT*COEF00(K,it)* PO(INDOD,2) )*RFACT ! hmhj
                  RTO(INDOD,1,1,IT) = ( RTO(INDOD,1,1,IT)+
     X                 (FACT+fact4)*COEF00(K,it)* PO(INDOD,1) )*rfact4
                  RTO(INDOD,2,1,IT) = ( RTO(INDOD,2,1,IT)+
     X                 (FACT+fact4)*COEF00(K,it)* PO(INDOD,2) )*rfact4
 
               ENDDO
            ENDDO
c     print*,'ok8'
         enddo                  ! Tracer do loop end
      endif                     ! If ntrac > 1
c idea change4 end
craa********************************************************************
!
      CALL COUNTPERF(1,13,0.)

!     print *,' leave idea_deldifs ' 	! hmhj
!!
      RETURN
      END
