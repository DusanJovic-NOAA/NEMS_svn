      SUBROUTINE SASCNV(IM,IX,KM,JCAP,DELT,DEL,PRSL,PS,PHIL,QL,
!     SUBROUTINE SASCNV(IM,IX,KM,JCAP,DELT,DEL,PRSL,PHIL,QL,
     &       Q1,T1,U1,V1,RCS,CLDWRK,RN,KBOT,KTOP,KUO,SLIMSK,
     &       DOT,XKT2,ncloud)
!  for cloud water version
!     parameter(ncloud=0)
!     SUBROUTINE SASCNV(KM,JCAP,DELT,DEL,SL,SLK,PS,QL,
!    &       Q1,T1,U1,V1,RCS,CLDWRK,RN,KBOT,KTOP,KUO,SLIMSK,
!    &       DOT,xkt2,ncloud)
!
      USE MACHINE , ONLY : kind_phys
      USE FUNCPHYS , ONLY : fpvs
      USE PHYSCONS, grav => con_g, CP => con_CP, HVAP => con_HVAP
     &,             RV => con_RV, FV => con_fvirt, T0C => con_T0C
     &,             CVAP => con_CVAP, CLIQ => con_CLIQ
     &,             EPS => con_eps, EPSM1 => con_epsm1
      implicit none
!
!
      integer            IM, IX,  KM, JCAP, ncloud,
     &                   KBOT(IM), KTOP(IM), KUO(IM)
      real(kind=kind_phys) DELT
      real(kind=kind_phys) PS(IM),     DEL(IX,KM),  PRSL(IX,KM),
!     real(kind=kind_phys)             DEL(IX,KM),  PRSL(IX,KM),
     &                     QL(IX,KM,2),Q1(IX,KM),   T1(IX,KM),
     &                     U1(IX,KM),  V1(IX,KM),   RCS(IM),
     &                     CLDWRK(IM), RN(IM),      SLIMSK(IM),
     &                     DOT(IX,KM), XKT2(IM),    PHIL(IX,KM)
!
      integer              I, INDX, jmn, k, knumb, latd, lond, km1
!
      real(kind=kind_phys) adw,     alpha,   alphal,  alphas,
     &                     aup,     beta,    betal,   betas,
     &                     c0,      cpoel,   dellat,  delta,
     &                     desdt,   deta,    detad,   dg,
     &                     dh,      dhh,     dlnsig,  dp,
     &                     dq,      dqsdp,   dqsdt,   dt,
     &                     dt2,     dtmax,   dtmin,   dv1,
     &                     dv1q,    dv2,     dv2q,    dv1u,
     &                     dv1v,    dv2u,    dv2v,    dv3u,
     &                     dv3v,    dv3,     dv3q,    dvq1,
     &                     dz,      dz1,     e1,      edtmax,
     &                     edtmaxl, edtmaxs, el2orc,  elocp,
     &                     es,      etah,
     &                     evef,    evfact,  evfactl, fact1,
     &                     fact2,   factor,  fjcap,   fkm,
     &                     fuv,     g,       gamma,   onemf,
     &                     onemfu,  pdetrn,  pdpdwn,  pprime,
     &                     qc,      qlk,     qrch,    qs,
     &                     rain,    rfact,   shear,   tem1,
     &                     tem2,    terr,    val,     val1,
     &                     val2,    w1,      w1l,     w1s,
     &                     w2,      w2l,     w2s,     w3,
     &                     w3l,     w3s,     w4,      w4l,
     &                     w4s,     xdby,    xpw,     xpwd,
     &                     xqc,     xqrch,   xlambu,  mbdt,
     &                     tem
!
!
      integer              JMIN(IM), KB(IM), KBCON(IM), KBDTR(IM),
     &                     KT2(IM),  KTCON(IM), LMIN(IM),
     &                     kbm(IM),  kbmax(IM), kmax(IM)
!
      real(kind=kind_phys) AA1(IM),     ACRT(IM),   ACRTFCT(IM),
     &                     DELHBAR(IM), DELQ(IM),   DELQ2(IM),
     &                     DELQBAR(IM), DELQEV(IM), DELTBAR(IM),
     &                     DELTV(IM),   DTCONV(IM), EDT(IM),
     &                     EDTO(IM),    EDTX(IM),   FLD(IM),
     &                     HCDO(IM),    HKBO(IM),   HMAX(IM),
     &                     HMIN(IM),    HSBAR(IM),  UCDO(IM),
     &                     UKBO(IM),    VCDO(IM),   VKBO(IM),
     &                     PBCDIF(IM),  PDOT(IM),   PO(IM,KM),
     &                                  PWAVO(IM),  PWEVO(IM),
!    &                     PSFC(IM),    PWAVO(IM),  PWEVO(IM),
     &                     QCDO(IM),    QCOND(IM),  QEVAP(IM),
     &                     QKBO(IM),    RNTOT(IM),  VSHEAR(IM),
     &                     XAA0(IM),    XHCD(IM),   XHKB(IM),
     &                     XK(IM),      XLAMB(IM),  XLAMD(IM),
     &                     XMB(IM),     XMBMAX(IM), XPWAV(IM),
     &                     XPWEV(IM),   XQCD(IM),   XQKB(IM)
cc
C  PHYSICAL PARAMETERS
      PARAMETER(G=grav)
      PARAMETER(CPOEL=CP/HVAP,ELOCP=HVAP/CP,
     &          EL2ORC=HVAP*HVAP/(RV*CP))
      PARAMETER(TERR=0.,C0=.002,DELTA=fv)
      PARAMETER(FACT1=(CVAP-CLIQ)/RV,FACT2=HVAP/RV-FACT1*T0C)
C  LOCAL VARIABLES AND ARRAYS
      real(kind=kind_phys) PFLD(IM,KM),    TO(IM,KM),     QO(IM,KM),
     &                     UO(IM,KM),      VO(IM,KM),     QESO(IM,KM)
c  cloud water
      real(kind=kind_phys) QLKO_KTCON(IM), DELLAL(IM),    TVO(IM,KM),
     &                     DBYO(IM,KM),    ZO(IM,KM),     SUMZ(IM,KM),
     &                     SUMH(IM,KM),    HEO(IM,KM),    HESO(IM,KM),
     &                     QRCD(IM,KM),    DELLAH(IM,KM), DELLAQ(IM,KM),
     &                     DELLAU(IM,KM),  DELLAV(IM,KM), HCKO(IM,KM),
     &                     UCKO(IM,KM),    VCKO(IM,KM),   QCKO(IM,KM),
     &                     ETA(IM,KM),     ETAU(IM,KM),   ETAD(IM,KM),
     &                     QRCDO(IM,KM),   PWO(IM,KM),    PWDO(IM,KM),
     &                     RHBAR(IM),      TX1(IM)
!
      LOGICAL TOTFLG, CNVFLG(IM), DWNFLG(IM), DWNFLG2(IM), FLG(IM)
!
      real(kind=kind_phys) PCRIT(15), ACRITT(15), ACRIT(15)
cmy      SAVE PCRIT, ACRITT
      DATA PCRIT/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,
     &           350.,300.,250.,200.,150./
      DATA ACRITT/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,
     &           .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
C  GDAS DERIVED ACRIT
C     DATA ACRITT/.203,.515,.521,.566,.625,.665,.659,.688,
C    &            .743,.813,.886,.947,1.138,1.377,1.896/
cc
      real(kind=kind_phys) tf, tcr, tcrf
      parameter (TF=233.16, TCR=263.16, TCRF=1.0/(TCR-TF)) ! From Lord(1978)
!
!     parameter (tf=258.16, tcr=273.16, tcrf=1.0/(tcr-tf))
!
      real(kind=kind_phys), parameter :: cons_0=0.0
c
c--------------------------------------------------------------------
!
      km1 = km - 1
C  INITIALIZE ARRAYS
C
      DO I=1,IM
        RN(I)=0.
        KBOT(I)=KM+1
        KTOP(I)=0
!       KUO(I)=0
        CNVFLG(I) = .TRUE.
        DTCONV(I) = 3600.
        CLDWRK(I) = 0.
        PDOT(I) = 0.
        KT2(I) = 0
        QLKO_KTCON(I) = 0.
        DELLAL(I) = 0.
      ENDDO
!!
      DO K = 1, 15
        ACRIT(K) = ACRITT(K) * (975. - PCRIT(K))
      ENDDO
      DT2 = DELT
      val   =         1200.
      dtmin = max(dt2, val )
      val   =         3600.
      dtmax = max(dt2, val )
C  MODEL TUNABLE PARAMETERS ARE ALL HERE
      MBDT    = 10.
      EDTMAXl = .3
      EDTMAXs = .3
      ALPHAl  = .5
      ALPHAs  = .5
      BETAl   = .15
      betas   = .15
      BETAl   = .05
      betas   = .05
c     EVEF    = 0.07
      evfact  = 0.3
      evfactl = 0.3
      PDPDWN  = 0.
      PDETRN  = 200.
      xlambu  = 1.e-4
      fjcap   = (float(jcap) / 126.) ** 2
      val     =           1.
      fjcap   = max(fjcap,val)
      fkm     = (float(km) / 28.) ** 2
      fkm     = max(fkm,val)
      W1l     = -8.E-3
      W2l     = -4.E-2
      W3l     = -5.E-3
      W4l     = -5.E-4
      W1s     = -2.E-4
      W2s     = -2.E-3
      W3s     = -1.E-3
      W4s     = -2.E-5
CCCCC IF(IM.EQ.384) THEN
        LATD  = 92
        lond  = 189
CCCCC ELSEIF(IM.EQ.768) THEN
CCCCC   LATD = 80
CCCCC ELSE
CCCCC   LATD = 0
CCCCC ENDIF
C
C  DEFINE TOP LAYER FOR SEARCH OF THE DOWNDRAFT ORIGINATING LAYER
C  AND THE MAXIMUM THETAE FOR UPDRAFT
C
      DO I=1,IM
        KBMAX(I) = KM
        KBM(I)   = KM
        KMAX(I)  = KM
        TX1(I)   = 1.0 / PS(I)
      ENDDO
!
      DO K = 1, KM
        DO I=1,IM
          IF (prSL(I,K)*tx1(I) .GT. 0.45) KBMAX(I) = K + 1
          IF (prSL(I,K)*tx1(I) .GT. 0.70) KBM(I)   = K + 1
          IF (prSL(I,K)*tx1(I) .GT. 0.04) KMAX(I)  = K + 1
        ENDDO
      ENDDO
      DO I=1,IM
        KBMAX(I) = MIN(KBMAX(I),KMAX(I))
        KBM(I)   = MIN(KBM(I),KMAX(I))
      ENDDO
C
C   CONVERT SURFACE PRESSURE TO MB FROM CB
C
!!
      DO K = 1, KM
        DO I=1,IM
          if (K .le. kmax(i)) then
            PFLD(I,k) = PRSL(I,K) * 10.0
            PWO(I,k)  = 0.
            PWDO(I,k) = 0.
            TO(I,k)   = T1(I,k)
            QO(I,k)   = Q1(I,k)
            UO(I,k)   = U1(I,k)
            VO(I,k)   = V1(I,k)
            DBYO(I,k) = 0.
            SUMZ(I,k) = 0.
            SUMH(I,k) = 0.
          endif
        ENDDO
      ENDDO
C
C  COLUMN VARIABLES
C  P IS PRESSURE OF THE LAYER (MB)
C  T IS TEMPERATURE AT T-DT (K)..TN
C  Q IS MIXING RATIO AT T-DT (KG/KG)..QN
C  TO IS TEMPERATURE AT T+DT (K)... THIS IS AFTER ADVECTION AND TURBULAN
C  QO IS MIXING RATIO AT T+DT (KG/KG)..Q1
C
      DO K = 1, KM
        DO I=1,IM
          if (k .le. kmax(i)) then
!jfe        QESO(I,k) = 10. * FPVS(T1(I,k))
!
            QESO(I,k) = 0.01 * fpvs(T1(I,K))      ! fpvs is in Pa
!
            QESO(I,k) = EPS * QESO(I,k) / (PFLD(I,k) + EPSM1*QESO(I,k))
            val1      =             1.E-8
            QESO(I,k) = MAX(QESO(I,k), val1)
            val2      =           1.e-10
            QO(I,k)   = max(QO(I,k), val2 )
c           QO(I,k)   = MIN(QO(I,k),QESO(I,k))
            TVO(I,k)  = TO(I,k) + DELTA * TO(I,k) * QO(I,k)
          endif
        ENDDO
      ENDDO
C
C  HYDROSTATIC HEIGHT ASSUME ZERO TERR
C
      DO K = 1, KM
        DO I=1,IM
          ZO(I,k) = PHIL(I,k) / G
        ENDDO
      ENDDO
C  COMPUTE MOIST STATIC ENERGY
      DO K = 1, KM
        DO I=1,IM
          if (K .le. kmax(i)) then
!           tem       = G * ZO(I,k) + CP * TO(I,k)
            tem       = PHIL(I,k) + CP * TO(I,k)
            HEO(I,k)  = tem  + HVAP * QO(I,k)
            HESO(I,k) = tem  + HVAP * QESO(I,k)
C           HEO(I,k)  = MIN(HEO(I,k),HESO(I,k))
          endif
        ENDDO
      ENDDO
C
C  DETERMINE LEVEL WITH LARGEST MOIST STATIC ENERGY
C  THIS IS THE LEVEL WHERE UPDRAFT STARTS
C
      DO I=1,IM
        HMAX(I) = HEO(I,1)
        KB(I) = 1
      ENDDO
!!
      DO K = 2, KM
        DO I=1,IM
          if (k .le. kbm(i)) then
            IF(HEO(I,k).GT.HMAX(I).AND.CNVFLG(I)) THEN
              KB(I)   = K
              HMAX(I) = HEO(I,k)
            ENDIF
          endif
        ENDDO
      ENDDO
C     DO K = 1, KMAX - 1
C         TOL(k) = .5 * (TO(I,k) + TO(I,k+1))
C         QOL(k) = .5 * (QO(I,k) + QO(I,k+1))
C         QESOL(I,k) = .5 * (QESO(I,k) + QESO(I,k+1))
C         HEOL(I,k) = .5 * (HEO(I,k) + HEO(I,k+1))
C         HESOL(I,k) = .5 * (HESO(I,k) + HESO(I,k+1))
C     ENDDO
      DO K = 1, KM1
        DO I=1,IM
          if (k .le. kmax(i)-1) then
            DZ      = .5 * (ZO(I,k+1) - ZO(I,k))
            DP      = .5 * (PFLD(I,k+1) - PFLD(I,k))
!jfe        ES      = 10. * FPVS(TO(I,k+1))
!
            ES      = 0.01 * fpvs(TO(I,K+1))      ! fpvs is in Pa
!
            PPRIME  = PFLD(I,k+1) + EPSM1 * ES
            QS      = EPS * ES / PPRIME
            DQSDP   = - QS / PPRIME
            DESDT   = ES * (FACT1 / TO(I,k+1) + FACT2 / (TO(I,k+1)**2))
            DQSDT   = QS * PFLD(I,k+1) * DESDT / (ES * PPRIME)
            GAMMA   = EL2ORC * QESO(I,k+1) / (TO(I,k+1)**2)
            DT      = (G * DZ + HVAP * DQSDP * DP) / (CP * (1. + GAMMA))
            DQ      = DQSDT * DT + DQSDP * DP
            TO(I,k) = TO(I,k+1) + DT
            QO(I,k) = QO(I,k+1) + DQ
            PO(I,k) = .5 * (PFLD(I,k) + PFLD(I,k+1))
          endif
        ENDDO
      ENDDO
!
      DO K = 1, KM1
        DO I=1,IM
          if (k .le. kmax(I)-1) then
!jfe        QESO(I,k) = 10. * FPVS(TO(I,k))
!
            QESO(I,k) = 0.01 * fpvs(TO(I,K))      ! fpvs is in Pa
!
            QESO(I,k) = EPS * QESO(I,k) / (PO(I,k) + EPSM1*QESO(I,k))
            val1      =             1.E-8
            QESO(I,k) = MAX(QESO(I,k), val1)
            val2      =           1.e-10
            QO(I,k)   = max(QO(I,k), val2 )
c           QO(I,k)   = MIN(QO(I,k),QESO(I,k))
            HEO(I,k)  = .5 * G * (ZO(I,k) + ZO(I,k+1)) +
     &                  CP * TO(I,k) + HVAP * QO(I,k)
            HESO(I,k) = .5 * G * (ZO(I,k) + ZO(I,k+1)) +
     &                  CP * TO(I,k) + HVAP * QESO(I,k)
            UO(I,k)   = .5 * (UO(I,k) + UO(I,k+1))
            VO(I,k)   = .5 * (VO(I,k) + VO(I,k+1))
          endif
        ENDDO
      ENDDO
c     k = kmax
c       HEO(I,k) = HEO(I,k)
c       hesol(k) = HESO(I,k)
c      IF(LAT.EQ.LATD.AND.lon.eq.lond.and.CNVFLG(I)) THEN
c        PRINT *, '   HEO ='
c        PRINT 6001, (HEO(I,K),K=1,KMAX)
c        PRINT *, '   HESO ='
c        PRINT 6001, (HESO(I,K),K=1,KMAX)
c        PRINT *, '   TO ='
c        PRINT 6002, (TO(I,K)-273.16,K=1,KMAX)
c        PRINT *, '   QO ='
c        PRINT 6003, (QO(I,K),K=1,KMAX)
c        PRINT *, '   QSO ='
c        PRINT 6003, (QESO(I,K),K=1,KMAX)
c      ENDIF
C
C  LOOK FOR CONVECTIVE CLOUD BASE AS THE LEVEL OF FREE CONVECTION
C
      DO I=1,IM
        IF(CNVFLG(I)) THEN
          INDX    = KB(I)
          HKBO(I) = HEO(I,INDX)
          QKBO(I) = QO(I,INDX)
          UKBO(I) = UO(I,INDX)
          VKBO(I) = VO(I,INDX)
        ENDIF
        FLG(I)    = CNVFLG(I)
        KBCON(I)  = KMAX(I)
      ENDDO
!!
      DO K = 1, KM
        DO I=1,IM
          if (k .le. kbmax(i)) then
            IF(FLG(I).AND.K.GT.KB(I)) THEN
              HSBAR(I)   = HESO(I,k)
              IF(HKBO(I).GT.HSBAR(I)) THEN
                FLG(I)   = .FALSE.
                KBCON(I) = K
              ENDIF
            ENDIF
          endif
        ENDDO
      ENDDO
      DO I=1,IM
        IF(CNVFLG(I)) THEN
          PBCDIF(I) = -PFLD(I,KBCON(I)) + PFLD(I,KB(I))
          PDOT(I)   = 10.* DOT(I,KBCON(I))
          IF(PBCDIF(I).GT.150.)    CNVFLG(I) = .FALSE.
          IF(KBCON(I).EQ.KMAX(I))  CNVFLG(I) = .FALSE.
        ENDIF
      ENDDO
!!
      TOTFLG = .TRUE.
      DO I=1,IM
        TOTFLG = TOTFLG .AND. (.NOT. CNVFLG(I))
      ENDDO
      IF(TOTFLG) RETURN
C  FOUND LFC, CAN DEFINE REST OF VARIABLES
 6001 FORMAT(2X,-2P10F12.2)
 6002 FORMAT(2X,10F12.2)
 6003 FORMAT(2X,3P10F12.2)
C
C  DETERMINE ENTRAINMENT RATE BETWEEN KB AND KBCON
C
      DO I = 1, IM
        alpha = alphas
        if(SLIMSK(I).eq.1.) alpha = alphal
        IF(CNVFLG(I)) THEN
          IF(KB(I).EQ.1) THEN
            DZ = .5 * (ZO(I,KBCON(I)) + ZO(I,KBCON(I)-1)) - ZO(I,1)
          ELSE
            DZ = .5 * (ZO(I,KBCON(I)) + ZO(I,KBCON(I)-1))
     &         - .5 * (ZO(I,KB(I)) + ZO(I,KB(I)-1))
          ENDIF
          IF(KBCON(I).NE.KB(I)) THEN
            XLAMB(I) = - LOG(ALPHA) / DZ
          ELSE
            XLAMB(I) = 0.
          ENDIF
        ENDIF
      ENDDO
C  DETERMINE UPDRAFT MASS FLUX
      DO K = 1, KM
        DO I = 1, IM
          if (k .le. kmax(i) .and. CNVFLG(I)) then
            ETA(I,k)  = 1.
            ETAU(I,k) = 1.
          ENDIF
        ENDDO
      ENDDO
      DO K = KM1, 2, -1
        DO I = 1, IM
          if (k .le. kbmax(i)) then
            IF(CNVFLG(I).AND.K.LT.KBCON(I).AND.K.GE.KB(I)) THEN
              DZ        = .5 * (ZO(I,k+1) - ZO(I,k-1))
              ETA(I,k)  = ETA(I,k+1) * EXP(-XLAMB(I) * DZ)
              ETAU(I,k) = ETA(I,k)
            ENDIF
          endif
        ENDDO
      ENDDO
      DO I = 1, IM
        IF(CNVFLG(I).AND.KB(I).EQ.1.AND.KBCON(I).GT.1) THEN
          DZ = .5 * (ZO(I,2) - ZO(I,1))
          ETA(I,1) = ETA(I,2) * EXP(-XLAMB(I) * DZ)
          ETAU(I,1) = ETA(I,1)
        ENDIF
      ENDDO
C
C  WORK UP UPDRAFT CLOUD PROPERTIES
C
      DO I = 1, IM
        IF(CNVFLG(I)) THEN
          INDX         = KB(I)
          HCKO(I,INDX) = HKBO(I)
          QCKO(I,INDX) = QKBO(I)
          UCKO(I,INDX) = UKBO(I)
          VCKO(I,INDX) = VKBO(I)
          PWAVO(I)     = 0.
        ENDIF
      ENDDO
C
C  CLOUD PROPERTY BELOW CLOUD BASE IS MODIFIED BY THE ENTRAINMENT PROCES
C
      DO K = 2, KM1
        DO I = 1, IM
          if (k .le. kmax(i)-1) then
            IF(CNVFLG(I).AND.K.GT.KB(I).AND.K.LE.KBCON(I)) THEN
              FACTOR = ETA(I,k-1) / ETA(I,k)
              ONEMF = 1. - FACTOR
              HCKO(I,k) = FACTOR * HCKO(I,k-1) + ONEMF *
     &                    .5 * (HEO(I,k) + HEO(I,k+1))
              UCKO(I,k) = FACTOR * UCKO(I,k-1) + ONEMF *
     &                    .5 * (UO(I,k) + UO(I,k+1))
              VCKO(I,k) = FACTOR * VCKO(I,k-1) + ONEMF *
     &                    .5 * (VO(I,k) + VO(I,k+1))
              DBYO(I,k) = HCKO(I,k) - HESO(I,k)
            ENDIF
            IF(CNVFLG(I).AND.K.GT.KBCON(I)) THEN
              HCKO(I,k) = HCKO(I,k-1)
              UCKO(I,k) = UCKO(I,k-1)
              VCKO(I,k) = VCKO(I,k-1)
              DBYO(I,k) = HCKO(I,k) - HESO(I,k)
            ENDIF
          endif
        ENDDO
      ENDDO
C  DETERMINE CLOUD TOP
      DO I = 1, IM
        FLG(I) = CNVFLG(I)
        KTCON(I) = 1
      ENDDO
C     DO K = 2, KMAX
C       KK = KMAX - K + 1
C         IF(DBYO(I,kK).GE.0..AND.FLG(I).AND.KK.GT.KBCON(I)) THEN
C           KTCON(I) = KK + 1
C           FLG(I) = .FALSE.
C         ENDIF
C     ENDDO
      DO K = 2, KM
        DO I = 1, IM
          if (k .le. kmax(i)) then
            IF(DBYO(I,k).LT.0..AND.FLG(I).AND.K.GT.KBCON(I)) THEN
              KTCON(I) = K
              FLG(I) = .FALSE.
            ENDIF
          endif
        ENDDO
      ENDDO
      DO I = 1, IM
        IF(CNVFLG(I).AND.(PFLD(I,KBCON(I)) - PFLD(I,KTCON(I))).LT.150.)
     &  CNVFLG(I) = .FALSE.
      ENDDO
      TOTFLG = .TRUE.
      DO I = 1, IM
        TOTFLG = TOTFLG .AND. (.NOT. CNVFLG(I))
      ENDDO
      IF(TOTFLG) RETURN
C
C  SEARCH FOR DOWNDRAFT ORIGINATING LEVEL ABOVE THETA-E MINIMUM
C
      DO I = 1, IM
        HMIN(I) = HEO(I,KBCON(I))
        LMIN(I) = KBMAX(I)
        JMIN(I) = KBMAX(I)
      ENDDO
      DO I = 1, IM
        DO K = KBCON(I), KBMAX(I)
          IF(HEO(I,k).LT.HMIN(I).AND.CNVFLG(I)) THEN
            LMIN(I) = K + 1
            HMIN(I) = HEO(I,k)
          ENDIF
        ENDDO
      ENDDO
C
C  Make sure that JMIN(I) is within the cloud
C
      DO I = 1, IM
        IF(CNVFLG(I)) THEN
          JMIN(I) = MIN(LMIN(I),KTCON(I)-1)
          XMBMAX(I) = .1
          JMIN(I) = MAX(JMIN(I),KBCON(I)+1)
        ENDIF
      ENDDO
C
C  ENTRAINING CLOUD
C
      do k = 2, km1
        DO I = 1, IM
          if (k .le. kmax(i)-1) then
            if(CNVFLG(I).and.k.gt.JMIN(I).and.k.le.KTCON(I)) THEN
              SUMZ(I,k) = SUMZ(I,k-1) + .5 * (ZO(I,k+1) - ZO(I,k-1))
              SUMH(I,k) = SUMH(I,k-1) + .5 * (ZO(I,k+1) - ZO(I,k-1))
     &                  * HEO(I,k)
            ENDIF
          endif
        enddo
      enddo
!!
      DO I = 1, IM
        IF(CNVFLG(I)) THEN
c         call random_number(XKT2)
c         call srand(fhour)
c         XKT2(I) = rand()
          KT2(I) = nint(XKT2(I)*float(KTCON(I)-JMIN(I))-.5)+JMIN(I)+1
!         KT2(I) = nint(sqrt(XKT2(I))*float(KTCON(I)-JMIN(I))-.5) + JMIN(I) + 1
c         KT2(I) = nint(ranf() *float(KTCON(I)-JMIN(I))-.5) + JMIN(I) + 1
          tem1 = (HCKO(I,JMIN(I)) - HESO(I,KT2(I)))
          tem2 = (SUMZ(I,KT2(I)) * HESO(I,KT2(I)) - SUMH(I,KT2(I)))
          if (abs(tem2) .gt. 0.000001) THEN
            XLAMB(I) = tem1 / tem2
          else
            CNVFLG(I) = .false.
          ENDIF
!         XLAMB(I) = (HCKO(I,JMIN(I)) - HESO(I,KT2(I)))
!    &          / (SUMZ(I,KT2(I)) * HESO(I,KT2(I)) - SUMH(I,KT2(I)))
          XLAMB(I) = max(XLAMB(I),cons_0)
          XLAMB(I) = min(XLAMB(I),2.3/SUMZ(I,KT2(I)))
        ENDIF
      ENDDO
!!
      DO I = 1, IM
       DWNFLG(I)  = CNVFLG(I)
       DWNFLG2(I) = CNVFLG(I)
       IF(CNVFLG(I)) THEN
        if(KT2(I).ge.KTCON(I)) DWNFLG(I) = .false.
      if(XLAMB(I).le.1.e-30.or.HCKO(I,JMIN(I))-HESO(I,KT2(I)).le.1.e-30)
     &  DWNFLG(I) = .false.
        do k = JMIN(I), KT2(I)
          if(DWNFLG(I).and.HEO(I,k).gt.HESO(I,KT2(I))) DWNFLG(I)=.false.
        enddo
c       IF(CNVFLG(I).AND.(PFLD(KBCON(I))-PFLD(KTCON(I))).GT.PDETRN)
c    &     DWNFLG(I)=.FALSE.
        IF(CNVFLG(I).AND.(PFLD(I,KBCON(I))-PFLD(I,KTCON(I))).LT.PDPDWN)
     &     DWNFLG2(I)=.FALSE.
       ENDIF
      ENDDO
!!
      DO K = 2, KM1
        DO I = 1, IM
          if (k .le. kmax(i)-1) then
            IF(DWNFLG(I).AND.K.GT.JMIN(I).AND.K.LE.KT2(I)) THEN
              DZ        = .5 * (ZO(I,k+1) - ZO(I,k-1))
c             ETA(I,k)  = ETA(I,k-1) * EXP( XLAMB(I) * DZ)
c  to simplify matter, we will take the linear approach here
c
              ETA(I,k)  = ETA(I,k-1) * (1. + XLAMB(I) * dz)
              ETAU(I,k) = ETAU(I,k-1) * (1. + (XLAMB(I)+xlambu) * dz)
            ENDIF
          endif
        ENDDO
      ENDDO
!!
      DO K = 2, KM1
        DO I = 1, IM
          if (k .le. kmax(i)-1) then
c           IF(.NOT.DWNFLG(I).AND.K.GT.JMIN(I).AND.K.LE.KT2(I)) THEN
            IF(.NOT.DWNFLG(I).AND.K.GT.JMIN(I).AND.K.LE.KTCON(I)) THEN
              DZ        = .5 * (ZO(I,k+1) - ZO(I,k-1))
              ETAU(I,k) = ETAU(I,k-1) * (1. + xlambu * dz)
            ENDIF
          endif
        ENDDO
      ENDDO
c      IF(LAT.EQ.LATD.AND.lon.eq.lond.and.CNVFLG(I)) THEN
c        PRINT *, ' LMIN(I), KT2(I)=', LMIN(I), KT2(I)
c        PRINT *, ' KBOT, KTOP, JMIN(I) =', KBCON(I), KTCON(I), JMIN(I)
c      ENDIF
c     IF(LAT.EQ.LATD.AND.lon.eq.lond) THEN
c       print *, ' xlamb =', xlamb
c       print *, ' eta =', (eta(k),k=1,KT2(I))
c       print *, ' ETAU =', (ETAU(I,k),k=1,KT2(I))
c       print *, ' HCKO =', (HCKO(I,k),k=1,KT2(I))
c       print *, ' SUMZ =', (SUMZ(I,k),k=1,KT2(I))
c       print *, ' SUMH =', (SUMH(I,k),k=1,KT2(I))
c     ENDIF
      DO I = 1, IM
        if(DWNFLG(I)) THEN
          KTCON(I) = KT2(I)
        ENDIF
      ENDDO
C
C  CLOUD PROPERTY ABOVE CLOUD Base IS MODIFIED BY THE DETRAINMENT PROCESS
C
      DO K = 2, KM1
        DO I = 1, IM
          if (k .le. kmax(i)-1) then
cjfe
            IF(CNVFLG(I).AND.K.GT.KBCON(I).AND.K.LE.KTCON(I)) THEN
cjfe      IF(K.GT.KBCON(I).AND.K.LE.KTCON(I)) THEN
              FACTOR    = ETA(I,k-1) / ETA(I,k)
              ONEMF     = 1. - FACTOR
              fuv       = ETAU(I,k-1) / ETAU(I,k)
              onemfu    = 1. - fuv
              HCKO(I,k) = FACTOR * HCKO(I,k-1) + ONEMF *
     &                    .5 * (HEO(I,k) + HEO(I,k+1))
              UCKO(I,k) = fuv * UCKO(I,k-1) + ONEMFu *
     &                    .5 * (UO(I,k) + UO(I,k+1))
              VCKO(I,k) = fuv * VCKO(I,k-1) + ONEMFu *
     &                    .5 * (VO(I,k) + VO(I,k+1))
              DBYO(I,k) = HCKO(I,k) - HESO(I,k)
            ENDIF
          endif
        ENDDO
      ENDDO
c      IF(LAT.EQ.LATD.AND.lon.eq.lond.and.CNVFLG(I)) THEN
c        PRINT *, ' UCKO=', (UCKO(I,k),k=KBCON(I)+1,KTCON(I))
c        PRINT *, ' uenv=', (.5*(UO(I,k)+UO(I,k-1)),k=KBCON(I)+1,KTCON(I))
c      ENDIF
      DO I = 1, IM
        if(CNVFLG(I).and.DWNFLG2(I).and.JMIN(I).le.KBCON(I))
     &     THEN
          CNVFLG(I) = .false.
          DWNFLG(I) = .false.
          DWNFLG2(I) = .false.
        ENDIF
      ENDDO
!!
      TOTFLG = .TRUE.
      DO I = 1, IM
        TOTFLG = TOTFLG .AND. (.NOT. CNVFLG(I))
      ENDDO
      IF(TOTFLG) RETURN
!!
C
C  COMPUTE CLOUD MOISTURE PROPERTY AND PRECIPITATION
C
      DO I = 1, IM
          AA1(I) = 0.
          RHBAR(I) = 0.
      ENDDO
      DO K = 1, KM
        DO I = 1, IM
          if (k .le. kmax(i)) then
            IF(CNVFLG(I).AND.K.GT.KB(I).AND.K.LT.KTCON(I)) THEN
              DZ = .5 * (ZO(I,k+1) - ZO(I,k-1))
              DZ1 = (ZO(I,k) - ZO(I,k-1))
              GAMMA = EL2ORC * QESO(I,k) / (TO(I,k)**2)
              QRCH = QESO(I,k)
     &             + GAMMA * DBYO(I,k) / (HVAP * (1. + GAMMA))
              FACTOR = ETA(I,k-1) / ETA(I,k)
              ONEMF = 1. - FACTOR
              QCKO(I,k) = FACTOR * QCKO(I,k-1) + ONEMF *
     &                    .5 * (QO(I,k) + QO(I,k+1))
              DQ = ETA(I,k) * QCKO(I,k) - ETA(I,k) * QRCH
              RHBAR(I) = RHBAR(I) + QO(I,k) / QESO(I,k)
C
C  BELOW LFC CHECK IF THERE IS EXCESS MOISTURE TO RELEASE LATENT HEAT
C
              IF(DQ.GT.0.) THEN
                ETAH = .5 * (ETA(I,k) + ETA(I,k-1))
                QLK = DQ / (ETA(I,k) + ETAH * C0 * DZ)
                AA1(I) = AA1(I) - DZ1 * G * QLK
                QC = QLK + QRCH
                PWO(I,k) = ETAH * C0 * DZ * QLK
                QCKO(I,k) = QC
                PWAVO(I) = PWAVO(I) + PWO(I,k)
              ENDIF
            ENDIF
          endif
        ENDDO
      ENDDO
      DO I = 1, IM
        RHBAR(I) = RHBAR(I) / float(KTCON(I) - KB(I) - 1)
      ENDDO
c
c  this section is ready for cloud water
c
      if(ncloud.gt.0) THEN
c
c  compute liquid and vapor separation at cloud top
c
      DO I = 1, IM
        k = KTCON(I)
        IF(CNVFLG(I)) THEN
          GAMMA = EL2ORC * QESO(I,K) / (TO(I,K)**2)
          QRCH = QESO(I,K)
     &         + GAMMA * DBYO(I,K) / (HVAP * (1. + GAMMA))
          DQ = QCKO(I,K-1) - QRCH
C
C  CHECK IF THERE IS EXCESS MOISTURE TO RELEASE LATENT HEAT
C
          IF(DQ.GT.0.) THEN
            QLKO_KTCON(I) = dq
            QCKO(I,K-1) = QRCH
          ENDIF
        ENDIF
      ENDDO
      ENDIF
C
C  CALCULATE CLOUD WORK FUNCTION AT T+DT
C
      DO K = 1, KM
        DO I = 1, IM
          if (k .le. kmax(i)) then
            IF(CNVFLG(I).AND.K.GT.KBCON(I).AND.K.LE.KTCON(I)) THEN
              DZ1 = ZO(I,k) - ZO(I,k-1)
              GAMMA = EL2ORC * QESO(I,k-1) / (TO(I,k-1)**2)
              RFACT =  1. + DELTA * CP * GAMMA
     &                 * TO(I,k-1) / HVAP
              AA1(I) = AA1(I) +
     &                 DZ1 * (G / (CP * TO(I,k-1)))
     &                 * DBYO(I,k-1) / (1. + GAMMA)
     &                 * RFACT
              val = 0.
              AA1(I)=AA1(I)+
     &                 DZ1 * G * DELTA *
     &                 MAX(val,(QESO(I,k-1) - QO(I,k-1)))
            ENDIF
          endif
        ENDDO
      ENDDO
      DO I = 1, IM
        IF(CNVFLG(I).AND.AA1(I).LE.0.) DWNFLG(I)  = .FALSE.
        IF(CNVFLG(I).AND.AA1(I).LE.0.) DWNFLG2(I) = .FALSE.
        IF(CNVFLG(I).AND.AA1(I).LE.0.) CNVFLG(I)  = .FALSE.
      ENDDO
!!
      TOTFLG = .TRUE.
      DO I = 1, IM
        TOTFLG = TOTFLG .AND. (.NOT. CNVFLG(I))
      ENDDO
      IF(TOTFLG) RETURN
!!
ccccc IF(LAT.EQ.LATD.AND.lon.eq.lond.and.CNVFLG(I)) THEN
ccccc   PRINT *, ' AA1(I) BEFORE DWNDRFT =', AA1(I)
ccccc ENDIF
C
C------- DOWNDRAFT CALCULATIONS
C
C
C--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
C
      DO I = 1, IM
        IF(CNVFLG(I)) THEN
          VSHEAR(I) = 0.
        ENDIF
      ENDDO
      DO K = 1, KM
        DO I = 1, IM
          if (k .le. kmax(i)) then
            IF(K.GE.KB(I).AND.K.LE.KTCON(I).AND.CNVFLG(I)) THEN
              shear=rcs(I) * sqrt((UO(I,k+1)-UO(I,k)) ** 2
     &                          + (VO(I,k+1)-VO(I,k)) ** 2)
              VSHEAR(I) = VSHEAR(I) + SHEAR
            ENDIF
          endif
        ENDDO
      ENDDO
      DO I = 1, IM
        EDT(I) = 0.
        IF(CNVFLG(I)) THEN
          KNUMB = KTCON(I) - KB(I) + 1
          KNUMB = MAX(KNUMB,1)
          VSHEAR(I) = 1.E3 * VSHEAR(I) / (ZO(I,KTCON(I))-ZO(I,KB(I)))
          E1=1.591-.639*VSHEAR(I)
     &       +.0953*(VSHEAR(I)**2)-.00496*(VSHEAR(I)**3)
          EDT(I)=1.-E1
          val =         .9
          EDT(I) = MIN(EDT(I),val)
          val =         .0
          EDT(I) = MAX(EDT(I),val)
          EDTO(I)=EDT(I)
          EDTX(I)=EDT(I)
        ENDIF
      ENDDO
C  DETERMINE DETRAINMENT RATE BETWEEN 1 AND KBDTR
      DO I = 1, IM
        KBDTR(I) = KBCON(I)
        beta = betas
        if(SLIMSK(I).eq.1.) beta = betal
        IF(CNVFLG(I)) THEN
          KBDTR(I) = KBCON(I)
          KBDTR(I) = MAX(KBDTR(I),1)
          XLAMD(I) = 0.
          IF(KBDTR(I).GT.1) THEN
            DZ = .5 * ZO(I,KBDTR(I)) + .5 * ZO(I,KBDTR(I)-1)
     &         - ZO(I,1)
            XLAMD(I) =  LOG(BETA) / DZ
          ENDIF
        ENDIF
      ENDDO
C  DETERMINE DOWNDRAFT MASS FLUX
      DO K = 1, KM
        DO I = 1, IM
          IF(k .le. kmax(i)) then
            IF(CNVFLG(I)) THEN
              ETAD(I,k) = 1.
            ENDIF
            QRCDO(I,k) = 0.
          endif
        ENDDO
      ENDDO
      DO K = KM1, 2, -1
        DO I = 1, IM
          if (k .le. kbmax(i)) then
            IF(CNVFLG(I).AND.K.LT.KBDTR(I)) THEN
              DZ        = .5 * (ZO(I,k+1) - ZO(I,k-1))
              ETAD(I,k) = ETAD(I,k+1) * EXP(XLAMD(I) * DZ)
            ENDIF
          endif
        ENDDO
      ENDDO
      K = 1
      DO I = 1, IM
        IF(CNVFLG(I).AND.KBDTR(I).GT.1) THEN
          DZ = .5 * (ZO(I,2) - ZO(I,1))
          ETAD(I,k) = ETAD(I,k+1) * EXP(XLAMD(I) * DZ)
        ENDIF
      ENDDO
C
C--- DOWNDRAFT MOISTURE PROPERTIES
C
      DO I = 1, IM
        PWEVO(I) = 0.
        FLG(I) = CNVFLG(I)
      ENDDO
      DO I = 1, IM
        IF(CNVFLG(I)) THEN
          JMN = JMIN(I)
          HCDO(I) = HEO(I,JMN)
          QCDO(I) = QO(I,JMN)
          QRCDO(I,JMN) = QESO(I,JMN)
          UCDO(I) = UO(I,JMN)
          VCDO(I) = VO(I,JMN)
        ENDIF
      ENDDO
      DO K = KM1, 1, -1
        DO I = 1, IM
          if (k .le. kmax(i)-1) then
            IF(CNVFLG(I).AND.K.LT.JMIN(I)) THEN
              DQ = QESO(I,k)
              DT = TO(I,k)
              GAMMA      = EL2ORC * DQ / DT**2
              DH         = HCDO(I) - HESO(I,k)
              QRCDO(I,k) = DQ+(1./HVAP)*(GAMMA/(1.+GAMMA))*DH
              DETAD      = ETAD(I,k+1) - ETAD(I,k)
              PWDO(I,k)  = ETAD(I,k+1) * QCDO(I) -
     &                     ETAD(I,k) * QRCDO(I,k)
              PWDO(I,k)  = PWDO(I,k) - DETAD *
     &                    .5 * (QRCDO(I,k) + QRCDO(I,k+1))
              QCDO(I)    = QRCDO(I,k)
              PWEVO(I)   = PWEVO(I) + PWDO(I,k)
            ENDIF
          endif
        ENDDO
      ENDDO
C     IF(LAT.EQ.LATD.AND.lon.eq.lond.and.DWNFLG(I)) THEN
C       PRINT *, ' PWAVO(I), PWEVO(I) =', PWAVO(I), PWEVO(I)
C     ENDIF
C
C--- FINAL DOWNDRAFT STRENGTH DEPENDENT ON PRECIP
C--- EFFICIENCY (EDT), NORMALIZED CONDENSATE (PWAV), AND
C--- EVAPORATE (PWEV)
C
      DO I = 1, IM
        edtmax = edtmaxl
        if(SLIMSK(I).eq.0.) edtmax = edtmaxs
        IF(DWNFLG2(I)) THEN
          IF(PWEVO(I).LT.0.) THEN
            EDTO(I) = -EDTO(I) * PWAVO(I) / PWEVO(I)
            EDTO(I) = MIN(EDTO(I),EDTMAX)
          ELSE
            EDTO(I) = 0.
          ENDIF
        ELSE
          EDTO(I) = 0.
        ENDIF
      ENDDO
C
C
C--- DOWNDRAFT CLOUDWORK FUNCTIONS
C
C
      DO K = KM1, 1, -1
        DO I = 1, IM
          if (k .le. kmax(i)-1) then
            IF(DWNFLG2(I).AND.K.LT.JMIN(I)) THEN
              GAMMA = EL2ORC * QESO(I,k+1) / TO(I,k+1)**2
              DHH=HCDO(I)
              DT=TO(I,k+1)
              DG=GAMMA
              DH=HESO(I,k+1)
              DZ=-1.*(ZO(I,k+1)-ZO(I,k))
              AA1(I)=AA1(I)+EDTO(I)*DZ*(G/(CP*DT))*((DHH-DH)/(1.+DG))
     &               *(1.+DELTA*CP*DG*DT/HVAP)
              val=0.
              AA1(I)=AA1(I)+EDTO(I)*
     &        DZ*G*DELTA*MAX(val,(QESO(I,k+1)-QO(I,k+1)))
            ENDIF
          endif
        ENDDO
      ENDDO
ccccc IF(LAT.EQ.LATD.AND.lon.eq.lond.and.DWNFLG2(I)) THEN
ccccc   PRINT *, '  AA1(I) AFTER DWNDRFT =', AA1(I)
ccccc ENDIF
      DO I = 1, IM
        IF(AA1(I).LE.0.) CNVFLG(I)  = .FALSE.
        IF(AA1(I).LE.0.) DWNFLG(I)  = .FALSE.
        IF(AA1(I).LE.0.) DWNFLG2(I) = .FALSE.
      ENDDO
!!
      TOTFLG = .TRUE.
      DO I = 1, IM
        TOTFLG = TOTFLG .AND. (.NOT. CNVFLG(I))
      ENDDO
      IF(TOTFLG) RETURN
!!
C
C
C--- WHAT WOULD THE CHANGE BE, THAT A CLOUD WITH UNIT MASS
C--- WILL DO TO THE ENVIRONMENT?
C
      DO K = 1, KM
        DO I = 1, IM
          IF(k .le. kmax(i) .and. CNVFLG(I)) THEN
            DELLAH(I,k) = 0.
            DELLAQ(I,k) = 0.
            DELLAU(I,k) = 0.
            DELLAV(I,k) = 0.
          ENDIF
        ENDDO
      ENDDO
      DO I = 1, IM
        IF(CNVFLG(I)) THEN
          DP = 1000. * DEL(I,1)
          DELLAH(I,1) = EDTO(I) * ETAD(I,1) * (HCDO(I)
     &                - HEO(I,1)) * G / DP
          DELLAQ(I,1) = EDTO(I) * ETAD(I,1) * (QCDO(I)
     &                - QO(I,1)) * G / DP
          DELLAU(I,1) = EDTO(I) * ETAD(I,1) * (UCDO(I)
     &                - UO(I,1)) * G / DP
          DELLAV(I,1) = EDTO(I) * ETAD(I,1) * (VCDO(I)
     &                - VO(I,1)) * G / DP
        ENDIF
      ENDDO
C
C--- CHANGED DUE TO SUBSIDENCE AND ENTRAINMENT
C
      DO K = 2, KM1
        DO I = 1, IM
          if (k .le. kmax(i)-1) then
            IF(CNVFLG(I).AND.K.LT.KTCON(I)) THEN
              AUP = 1.
              IF(K.LE.KB(I)) AUP = 0.
              ADW = 1.
              IF(K.GT.JMIN(I)) ADW = 0.
              DV1= HEO(I,k)
              DV2 = .5 * (HEO(I,k) + HEO(I,k+1))
              DV3= HEO(I,k-1)
              DV1Q= QO(I,k)
              DV2Q = .5 * (QO(I,k) + QO(I,k+1))
              DV3Q= QO(I,k-1)
              DV1U= UO(I,k)
              DV2U = .5 * (UO(I,k) + UO(I,k+1))
              DV3U= UO(I,k-1)
              DV1V= VO(I,k)
              DV2V = .5 * (VO(I,k) + VO(I,k+1))
              DV3V= VO(I,k-1)
              DP = 1000. * DEL(I,K)
              DZ = .5 * (ZO(I,k+1) - ZO(I,k-1))
              DETA = ETA(I,k) - ETA(I,k-1)
              DETAD = ETAD(I,k) - ETAD(I,k-1)
              DELLAH(I,k) = DELLAH(I,k) +
     &            ((AUP * ETA(I,k) - ADW * EDTO(I) * ETAD(I,k)) * DV1
     &        - (AUP * ETA(I,k-1) - ADW * EDTO(I) * ETAD(I,k-1))* DV3
     &                    - AUP * DETA * DV2
     &                    + ADW * EDTO(I) * DETAD * HCDO(I)) * G / DP
              DELLAQ(I,k) = DELLAQ(I,k) +
     &            ((AUP * ETA(I,k) - ADW * EDTO(I) * ETAD(I,k)) * DV1Q
     &        - (AUP * ETA(I,k-1) - ADW * EDTO(I) * ETAD(I,k-1))* DV3Q
     &                    - AUP * DETA * DV2Q
     &       +ADW*EDTO(I)*DETAD*.5*(QRCDO(I,k)+QRCDO(I,k-1))) * G / DP
              DELLAU(I,k) = DELLAU(I,k) +
     &            ((AUP * ETA(I,k) - ADW * EDTO(I) * ETAD(I,k)) * DV1U
     &        - (AUP * ETA(I,k-1) - ADW * EDTO(I) * ETAD(I,k-1))* DV3U
     &                     - AUP * DETA * DV2U
     &                    + ADW * EDTO(I) * DETAD * UCDO(I)
     &                    ) * G / DP
              DELLAV(I,k) = DELLAV(I,k) +
     &            ((AUP * ETA(I,k) - ADW * EDTO(I) * ETAD(I,k)) * DV1V
     &        - (AUP * ETA(I,k-1) - ADW * EDTO(I) * ETAD(I,k-1))* DV3V
     &                     - AUP * DETA * DV2V
     &                    + ADW * EDTO(I) * DETAD * VCDO(I)
     &                    ) * G / DP
            ENDIF
          endif
        ENDDO
      ENDDO
C
C------- CLOUD TOP
C
      DO I = 1, IM
        IF(CNVFLG(I)) THEN
          INDX = KTCON(I)
          DP = 1000. * DEL(I,INDX)
          DV1 = HEO(I,INDX-1)
          DELLAH(I,INDX) = ETA(I,INDX-1) *
     &                     (HCKO(I,INDX-1) - DV1) * G / DP
          DVQ1 = QO(I,INDX-1)
          DELLAQ(I,INDX) = ETA(I,INDX-1) *
     &                     (QCKO(I,INDX-1) - DVQ1) * G / DP
          DV1U = UO(I,INDX-1)
          DELLAU(I,INDX) = ETA(I,INDX-1) *
     &                     (UCKO(I,INDX-1) - DV1U) * G / DP
          DV1V = VO(I,INDX-1)
          DELLAV(I,INDX) = ETA(I,INDX-1) *
     &                     (VCKO(I,INDX-1) - DV1V) * G / DP
c
c  cloud water
c
          DELLAL(I) = ETA(I,INDX-1) * QLKO_KTCON(I) * g / dp
        ENDIF
      ENDDO
C
C------- FINAL CHANGED VARIABLE PER UNIT MASS FLUX
C
      DO K = 1, KM
        DO I = 1, IM
          if (k .le. kmax(i)) then
            IF(CNVFLG(I).and.k.gt.KTCON(I)) THEN
              QO(I,k) = Q1(I,k)
              TO(I,k) = T1(I,k)
              UO(I,k) = U1(I,k)
              VO(I,k) = V1(I,k)
            ENDIF
            IF(CNVFLG(I).AND.K.LE.KTCON(I)) THEN
              QO(I,k) = DELLAQ(I,k) * MBDT + Q1(I,k)
              DELLAT = (DELLAH(I,k) - HVAP * DELLAQ(I,k)) / CP
              TO(I,k) = DELLAT * MBDT + T1(I,k)
              val   =           1.e-10
              QO(I,k) = max(QO(I,k), val  )
            ENDIF
          endif
        ENDDO
      ENDDO
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C--- THE ABOVE CHANGED ENVIRONMENT IS NOW USED TO CALULATE THE
C--- EFFECT THE ARBITRARY CLOUD (WITH UNIT MASS FLUX)
C--- WOULD HAVE ON THE STABILITY,
C--- WHICH THEN IS USED TO CALCULATE THE REAL MASS FLUX,
C--- NECESSARY TO KEEP THIS CHANGE IN BALANCE WITH THE LARGE-SCALE
C--- DESTABILIZATION.
C
C--- ENVIRONMENTAL CONDITIONS AGAIN, FIRST HEIGHTS
C
      DO K = 1, KM
        DO I = 1, IM
          IF(k .le. kmax(i) .and. CNVFLG(I)) THEN
!jfe        QESO(I,k) = 10. * FPVS(TO(I,k))
!
            QESO(I,k) = 0.01 * fpvs(TO(I,K))      ! fpvs is in Pa
!
            QESO(I,k) = EPS * QESO(I,k) / (PFLD(I,k)+EPSM1*QESO(I,k))
            val       =             1.E-8
            QESO(I,k) = MAX(QESO(I,k), val )
            TVO(I,k)  = TO(I,k) + DELTA * TO(I,k) * QO(I,k)
          ENDIF
        ENDDO
      ENDDO
      DO I = 1, IM
        IF(CNVFLG(I)) THEN
          XAA0(I) = 0.
          XPWAV(I) = 0.
        ENDIF
      ENDDO
C
C  HYDROSTATIC HEIGHT ASSUME ZERO TERR
C
!     DO I = 1, IM
!       IF(CNVFLG(I)) THEN
!         DLNSIG =  LOG(PRSL(I,1)/PS(I))
!         ZO(I,1) = TERR - DLNSIG * RD / G * TVO(I,1)
!       ENDIF
!     ENDDO
!     DO K = 2, KM
!       DO I = 1, IM
!         IF(k .le. kmax(i) .and. CNVFLG(I)) THEN
!           DLNSIG =  LOG(PRSL(I,K) / PRSL(I,K-1))
!           ZO(I,k) = ZO(I,k-1) - DLNSIG * RD / G
!    &             * .5 * (TVO(I,k) + TVO(I,k-1))
!         ENDIF
!       ENDDO
!     ENDDO
C
C--- MOIST STATIC ENERGY
C
      DO K = 1, KM1
        DO I = 1, IM
          IF(k .le. kmax(i)-1 .and. CNVFLG(I)) THEN
            DZ = .5 * (ZO(I,k+1) - ZO(I,k))
            DP = .5 * (PFLD(I,k+1) - PFLD(I,k))
cjfe        ES = 10. * FPVS(TO(I,k+1))
!
            ES = 0.01 * fpvs(TO(I,K+1))      ! fpvs is in Pa
!
            PPRIME = PFLD(I,k+1) + EPSM1 * ES
            QS = EPS * ES / PPRIME
            DQSDP = - QS / PPRIME
            DESDT = ES * (FACT1 / TO(I,k+1) + FACT2 / (TO(I,k+1)**2))
            DQSDT = QS * PFLD(I,k+1) * DESDT / (ES * PPRIME)
            GAMMA = EL2ORC * QESO(I,k+1) / (TO(I,k+1)**2)
            DT = (G * DZ + HVAP * DQSDP * DP) / (CP * (1. + GAMMA))
            DQ = DQSDT * DT + DQSDP * DP
            TO(I,k) = TO(I,k+1) + DT
            QO(I,k) = QO(I,k+1) + DQ
            PO(I,k) = .5 * (PFLD(I,k) + PFLD(I,k+1))
          ENDIF
        ENDDO
      ENDDO
      DO K = 1, KM1
        DO I = 1, IM
          IF(k .le. kmax(i)-1 .and. CNVFLG(I)) THEN
cjfe        QESO(I,k) = 10. * FPVS(TO(I,k))
!
            QESO(I,k) = 0.01 * fpvs(TO(I,K))      ! fpvs is in Pa
!
            QESO(I,k) = EPS * QESO(I,k) / (PO(I,k) + EPSM1 * QESO(I,k))
            val1      =             1.E-8
            QESO(I,k) = MAX(QESO(I,k), val1)
            val2      =           1.e-10
            QO(I,k)   = max(QO(I,k), val2 )
c           QO(I,k)   = MIN(QO(I,k),QESO(I,k))
            HEO(I,k)   = .5 * G * (ZO(I,k) + ZO(I,k+1)) +
     &                    CP * TO(I,k) + HVAP * QO(I,k)
            HESO(I,k) = .5 * G * (ZO(I,k) + ZO(I,k+1)) +
     &                  CP * TO(I,k) + HVAP * QESO(I,k)
          ENDIF
        ENDDO
      ENDDO
      DO I = 1, IM
        k = kmax(i)
        IF(CNVFLG(I)) THEN
          HEO(I,k) = G * ZO(I,k) + CP * TO(I,k) + HVAP * QO(I,k)
          HESO(I,k) = G * ZO(I,k) + CP * TO(I,k) + HVAP * QESO(I,k)
c         HEO(I,k) = MIN(HEO(I,k),HESO(I,k))
        ENDIF
      ENDDO
      DO I = 1, IM
        IF(CNVFLG(I)) THEN
          INDX = KB(I)
          XHKB(I) = HEO(I,INDX)
          XQKB(I) = QO(I,INDX)
          HCKO(I,INDX) = XHKB(I)
          QCKO(I,INDX) = XQKB(I)
        ENDIF
      ENDDO
C
C
C**************************** STATIC CONTROL
C
C
C------- MOISTURE AND CLOUD WORK FUNCTIONS
C
      DO K = 2, KM1
        DO I = 1, IM
          if (k .le. kmax(i)-1) then
C           IF(CNVFLG(I).AND.K.GT.KB(I).AND.K.LE.KBCON(I)) THEN
            IF(CNVFLG(I).AND.K.GT.KB(I).AND.K.LE.KTCON(I)) THEN
              FACTOR = ETA(I,k-1) / ETA(I,k)
              ONEMF = 1. - FACTOR
              HCKO(I,k) = FACTOR * HCKO(I,k-1) + ONEMF *
     &                    .5 * (HEO(I,k) + HEO(I,k+1))
            ENDIF
C           IF(CNVFLG(I).AND.K.GT.KBCON(I)) THEN
C             HEO(I,k) = HEO(I,k-1)
C           ENDIF
          endif
        ENDDO
      ENDDO
      DO K = 2, KM1
        DO I = 1, IM
          if (k .le. kmax(i)-1) then
            IF(CNVFLG(I).AND.K.GT.KB(I).AND.K.LT.KTCON(I)) THEN
              DZ = .5 * (ZO(I,k+1) - ZO(I,k-1))
              GAMMA = EL2ORC * QESO(I,k) / (TO(I,k)**2)
              XDBY = HCKO(I,k) - HESO(I,k)
              val  =          0.
              XDBY = MAX(XDBY,val)
              XQRCH = QESO(I,k)
     &              + GAMMA * XDBY / (HVAP * (1. + GAMMA))
              FACTOR = ETA(I,k-1) / ETA(I,k)
              ONEMF = 1. - FACTOR
              QCKO(I,k) = FACTOR * QCKO(I,k-1) + ONEMF *
     &                    .5 * (QO(I,k) + QO(I,k+1))
              DQ = ETA(I,k) * QCKO(I,k) - ETA(I,k) * XQRCH
              IF(DQ.GT.0.) THEN
                ETAH = .5 * (ETA(I,k) + ETA(I,k-1))
                QLK = DQ / (ETA(I,k) + ETAH * C0 * DZ)
                XAA0(I) = XAA0(I) - (ZO(I,k) - ZO(I,k-1)) * G * QLK
                XQC = QLK + XQRCH
                XPW = ETAH * C0 * DZ * QLK
                QCKO(I,k) = XQC
                XPWAV(I) = XPWAV(I) + XPW
              ENDIF
            ENDIF
c           IF(CNVFLG(I).AND.K.GT.KBCON(I).AND.K.LT.KTCON(I)) THEN
            IF(CNVFLG(I).AND.K.GT.KBCON(I).AND.K.LE.KTCON(I)) THEN
              DZ1 = ZO(I,k) - ZO(I,k-1)
              GAMMA = EL2ORC * QESO(I,k-1) / (TO(I,k-1)**2)
              RFACT =  1. + DELTA * CP * GAMMA
     &                 * TO(I,k-1) / HVAP
              XDBY = HCKO(I,k-1) - HESO(I,k-1)
              XAA0(I) = XAA0(I)
     &                + DZ1 * (G / (CP * TO(I,k-1)))
     &                * XDBY / (1. + GAMMA)
     &                * RFACT
              val=0.
              XAA0(I)=XAA0(I)+
     &                 DZ1 * G * DELTA *
     &                 MAX(val,(QESO(I,k-1) - QO(I,k-1)))
            ENDIF
          endif
        ENDDO
      ENDDO
ccccc IF(LAT.EQ.LATD.AND.lon.eq.lond.and.CNVFLG(I)) THEN
ccccc   PRINT *, ' XAA BEFORE DWNDRFT =', XAA0(I)
ccccc ENDIF
C
C------- DOWNDRAFT CALCULATIONS
C
C
C--- DOWNDRAFT MOISTURE PROPERTIES
C
      DO I = 1, IM
        XPWEV(I) = 0.
      ENDDO
      DO I = 1, IM
        IF(DWNFLG2(I)) THEN
          JMN = JMIN(I)
          XHCD(I) = HEO(I,JMN)
          XQCD(I) = QO(I,JMN)
          QRCD(I,JMN) = QESO(I,JMN)
        ENDIF
      ENDDO
      DO K = KM1, 1, -1
        DO I = 1, IM
          if (k .le. kmax(i)-1) then
            IF(DWNFLG2(I).AND.K.LT.JMIN(I)) THEN
              DQ = QESO(I,k)
              DT = TO(I,k)
              GAMMA    = EL2ORC * DQ / DT**2
              DH       = XHCD(I) - HESO(I,k)
              QRCD(I,k)=DQ+(1./HVAP)*(GAMMA/(1.+GAMMA))*DH
              DETAD    = ETAD(I,k+1) - ETAD(I,k)
              XPWD     = ETAD(I,k+1) * QRCD(I,k+1) -
     &                   ETAD(I,k) * QRCD(I,k)
              XPWD     = XPWD - DETAD *
     &                 .5 * (QRCD(I,k) + QRCD(I,k+1))
              XPWEV(I) = XPWEV(I) + XPWD
            ENDIF
          endif
        ENDDO
      ENDDO
C
      DO I = 1, IM
        edtmax = edtmaxl
        if(SLIMSK(I).eq.0.) edtmax = edtmaxs
        IF(DWNFLG2(I)) THEN
          IF(XPWEV(I).GE.0.) THEN
            EDTX(I) = 0.
          ELSE
            EDTX(I) = -EDTX(I) * XPWAV(I) / XPWEV(I)
            EDTX(I) = MIN(EDTX(I),EDTMAX)
          ENDIF
        ELSE
          EDTX(I) = 0.
        ENDIF
      ENDDO
C
C
C
C--- DOWNDRAFT CLOUDWORK FUNCTIONS
C
C
      DO K = KM1, 1, -1
        DO I = 1, IM
          if (k .le. kmax(i)-1) then
            IF(DWNFLG2(I).AND.K.LT.JMIN(I)) THEN
              GAMMA = EL2ORC * QESO(I,k+1) / TO(I,k+1)**2
              DHH=XHCD(I)
              DT= TO(I,k+1)
              DG= GAMMA
              DH= HESO(I,k+1)
              DZ=-1.*(ZO(I,k+1)-ZO(I,k))
              XAA0(I)=XAA0(I)+EDTX(I)*DZ*(G/(CP*DT))*((DHH-DH)/(1.+DG))
     &                *(1.+DELTA*CP*DG*DT/HVAP)
              val=0.
              XAA0(I)=XAA0(I)+EDTX(I)*
     &        DZ*G*DELTA*MAX(val,(QESO(I,k+1)-QO(I,k+1)))
            ENDIF
          endif
        ENDDO
      ENDDO
ccccc IF(LAT.EQ.LATD.AND.lon.eq.lond.and.DWNFLG2(I)) THEN
ccccc   PRINT *, '  XAA AFTER DWNDRFT =', XAA0(I)
ccccc ENDIF
C
C  CALCULATE CRITICAL CLOUD WORK FUNCTION
C
      DO I = 1, IM
        ACRT(I) = 0.
        IF(CNVFLG(I)) THEN
C       IF(CNVFLG(I).AND.SLIMSK(I).NE.1.) THEN
          IF(PFLD(I,KTCON(I)).LT.PCRIT(15))THEN
            ACRT(I)=ACRIT(15)*(975.-PFLD(I,KTCON(I)))
     &              /(975.-PCRIT(15))
          ELSE IF(PFLD(I,KTCON(I)).GT.PCRIT(1))THEN
            ACRT(I)=ACRIT(1)
          ELSE
            K =  int((850. - PFLD(I,KTCON(I)))/50.) + 2
            K = MIN(K,15)
            K = MAX(K,2)
            ACRT(I)=ACRIT(K)+(ACRIT(K-1)-ACRIT(K))*
     *           (PFLD(I,KTCON(I))-PCRIT(K))/(PCRIT(K-1)-PCRIT(K))
           ENDIF
C        ELSE
C          ACRT(I) = .5 * (PFLD(I,KBCON(I)) - PFLD(I,KTCON(I)))
         ENDIF
      ENDDO
      DO I = 1, IM
        ACRTFCT(I) = 1.
        IF(CNVFLG(I)) THEN
          if(SLIMSK(I).eq.1.) THEN
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          ENDIF
C       IF(CNVFLG(I).AND.SLIMSK(I).EQ.1.) THEN
C         ACRTFCT(I) = PDOT(I) / W3
c
c  modify critical cloud workfunction by cloud base vertical velocity
c
          IF(PDOT(I).LE.W4) THEN
            ACRTFCT(I) = (PDOT(I) - W4) / (W3 - W4)
          ELSEIF(PDOT(I).GE.-W4) THEN
            ACRTFCT(I) = - (PDOT(I) + W4) / (W4 - W3)
          ELSE
            ACRTFCT(I) = 0.
          ENDIF
          val1    =             -1.
          ACRTFCT(I) = MAX(ACRTFCT(I),val1)
          val2    =             1.
          ACRTFCT(I) = MIN(ACRTFCT(I),val2)
          ACRTFCT(I) = 1. - ACRTFCT(I)
c
c  modify ACRTFCT(I) by colume mean rh if RHBAR(I) is greater than 80 percent
c
c         if(RHBAR(I).ge..8) THEN
c           ACRTFCT(I) = ACRTFCT(I) * (.9 - min(RHBAR(I),.9)) * 10.
c         ENDIF
c
c  modify adjustment time scale by cloud base vertical velocity
c
          DTCONV(I) = DT2 + max((1800. - DT2),cons_0) *
     &                (PDOT(I) - W2) / (W1 - W2)
c         DTCONV(I) = MAX(DTCONV(I), DT2)
c         DTCONV(I) = 1800. * (PDOT(I) - w2) / (w1 - w2)
          DTCONV(I) = max(DTCONV(I),dtmin)
          DTCONV(I) = min(DTCONV(I),dtmax)
 
        ENDIF
      ENDDO
C
C--- LARGE SCALE FORCING
C
      DO I= 1, IM
        FLG(I) = CNVFLG(I)
        IF(CNVFLG(I)) THEN
C         F = AA1(I) / DTCONV(I)
          FLD(I) = (AA1(I) - ACRT(I) * ACRTFCT(I)) / DTCONV(I)
          IF(FLD(I).LE.0.) FLG(I) = .FALSE.
        ENDIF
        CNVFLG(I) = FLG(I)
        IF(CNVFLG(I)) THEN
C         XAA0(I) = MAX(XAA0(I),0.)
          XK(I) = (XAA0(I) - AA1(I)) / MBDT
          IF(XK(I).GE.0.) FLG(I) = .FALSE.
        ENDIF
C
C--- KERNEL, CLOUD BASE MASS FLUX
C
        CNVFLG(I) = FLG(I)
        IF(CNVFLG(I)) THEN
          XMB(I) = -FLD(I) / XK(I)
          XMB(I) = MIN(XMB(I),XMBMAX(I))
        ENDIF
      ENDDO
c      IF(LAT.EQ.LATD.AND.lon.eq.lond.and.CNVFLG(I)) THEN
c        print *, ' RHBAR(I), ACRTFCT(I) =', RHBAR(I), ACRTFCT(I)
c        PRINT *, '  A1, XA =', AA1(I), XAA0(I)
c        PRINT *, ' XMB(I), ACRT =', XMB(I), ACRT
c      ENDIF
      TOTFLG = .TRUE.
      DO I = 1, IM
        TOTFLG = TOTFLG .AND. (.NOT. CNVFLG(I))
      ENDDO
      IF(TOTFLG) RETURN
c
c  restore t0 and QO to t1 and q1 in case convection stops
c
      do k = 1, km
        DO I = 1, IM
          if (k .le. kmax(i)) then
            TO(I,k) = T1(I,k)
            QO(I,k) = Q1(I,k)
!jfe        QESO(I,k) = 10. * FPVS(T1(I,k))
!
            QESO(I,k) = 0.01 * fpvs(T1(I,K))      ! fpvs is in Pa
!
            QESO(I,k) = EPS * QESO(I,k) / (PFLD(I,k) + EPSM1*QESO(I,k))
            val     =             1.E-8
            QESO(I,k) = MAX(QESO(I,k), val )
          endif
        enddo
      enddo
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C--- FEEDBACK: SIMPLY THE CHANGES FROM THE CLOUD WITH UNIT MASS FLUX
C---           MULTIPLIED BY  THE MASS FLUX NECESSARY TO KEEP THE
C---           EQUILIBRIUM WITH THE LARGER-SCALE.
C
      DO I = 1, IM
        DELHBAR(I) = 0.
        DELQBAR(I) = 0.
        DELTBAR(I) = 0.
        QCOND(I) = 0.
      ENDDO
      DO K = 1, KM
        DO I = 1, IM
          if (k .le. kmax(i)) then
            IF(CNVFLG(I).AND.K.LE.KTCON(I)) THEN
              AUP = 1.
              IF(K.Le.KB(I)) AUP = 0.
              ADW = 1.
              IF(K.GT.JMIN(I)) ADW = 0.
              DELLAT = (DELLAH(I,k) - HVAP * DELLAQ(I,k)) / CP
              T1(I,k) = T1(I,k) + DELLAT * XMB(I) * DT2
              Q1(I,k) = Q1(I,k) + DELLAQ(I,k) * XMB(I) * DT2
              U1(I,k) = U1(I,k) + DELLAU(I,k) * XMB(I) * DT2
              V1(I,k) = V1(I,k) + DELLAV(I,k) * XMB(I) * DT2
              DP = 1000. * DEL(I,K)
              DELHBAR(I) = DELHBAR(I) + DELLAH(I,k)*XMB(I)*DP/G
              DELQBAR(I) = DELQBAR(I) + DELLAQ(I,k)*XMB(I)*DP/G
              DELTBAR(I) = DELTBAR(I) + DELLAT*XMB(I)*DP/G
            ENDIF
          endif
        ENDDO
      ENDDO
      DO K = 1, KM
        DO I = 1, IM
          if (k .le. kmax(i)) then
            IF(CNVFLG(I).AND.K.LE.KTCON(I)) THEN
!jfe          QESO(I,k) = 10. * FPVS(T1(I,k))
!
              QESO(I,k) = 0.01 * fpvs(T1(I,K))      ! fpvs is in Pa
!
              QESO(I,k) = EPS * QESO(I,k)/(PFLD(I,k) + EPSM1*QESO(I,k))
              val     =             1.E-8
              QESO(I,k) = MAX(QESO(I,k), val )
c
c  cloud water
c
              if(ncloud.gt.0.and.cnvflg(i).and.k.eq.ktcon(i)) then
                tem  = dellal(i) * xmb(i) * dt2
                tem1 = max(0.0, min(1.0, (tcr-t1(i,k))*tcrf))
                if (ql(i,k,2) .gt. -999.0) then
                  ql(i,k,1) = ql(i,k,1) + tem * tem1            ! ice
                  ql(i,k,2) = ql(i,k,2) + tem *(1.0-tem1)       ! water
                else
                  ql(i,k,1) = ql(i,k,1) + tem
                endif
                dp = 1000. * del(i,k)
                dellal(i) = dellal(i) * xmb(i) * dp / g
              endif
!
!             if(ncloud.gt.0.and.CNVFLG(I).and.k.eq.KTCON(I)) THEN
!               QL(I,k) = QL(I,k) + DELLAL(I) * XMB(I) * dt2
!               dp = 1000. * del(i,k)
!               DELLAL(I) = DELLAL(I) * XMB(I) * dp / g
!             ENDIF
!
            ENDIF
          endif
        ENDDO
      ENDDO
c     IF(LAT.EQ.LATD.AND.lon.eq.lond.and.CNVFLG(I) ) THEN
c       PRINT *, ' DELHBAR, DELQBAR, DELTBAR ='
c       PRINT *, DELHBAR, HVAP*DELQBAR, CP*DELTBAR
c       PRINT *, '   DELLBAR ='
c       PRINT 6003,  HVAP*DELLbar
c       PRINT *, '   DELLAQ ='
c       PRINT 6003, (HVAP*DELLAQ(I,k)*XMB(I),K=1,KMAX)
c       PRINT *, '   DELLAT ='
c       PRINT 6003, (DELLAH(i,k)*XMB(I)-HVAP*DELLAQ(I,k)*XMB(I),
c    &               K=1,KMAX)
c     ENDIF
      DO I = 1, IM
        RNTOT(I) = 0.
        DELQEV(I) = 0.
        DELQ2(I) = 0.
        FLG(I) = CNVFLG(I)
      ENDDO
      DO K = KM, 1, -1
        DO I = 1, IM
          if (k .le. kmax(i)) then
            IF(CNVFLG(I).AND.K.LE.KTCON(I)) THEN
              AUP = 1.
              IF(K.Le.KB(I)) AUP = 0.
              ADW = 1.
              IF(K.GT.JMIN(I)) ADW = 0.
              rain =  AUP * PWO(I,k) + ADW * EDTO(I) * PWDO(I,k)
              RNTOT(I) = RNTOT(I) + rain * XMB(I) * .001 * dt2
            ENDIF
          endif
        ENDDO
      ENDDO
      DO K = KM, 1, -1
        DO I = 1, IM
          if (k .le. kmax(i)) then
            DELTV(I) = 0.
            DELQ(I) = 0.
            QEVAP(I) = 0.
            IF(CNVFLG(I).AND.K.LE.KTCON(I)) THEN
              AUP = 1.
              IF(K.Le.KB(I)) AUP = 0.
              ADW = 1.
              IF(K.GT.JMIN(I)) ADW = 0.
              rain =  AUP * PWO(I,k) + ADW * EDTO(I) * PWDO(I,k)
              RN(I) = RN(I) + rain * XMB(I) * .001 * dt2
            ENDIF
            IF(FLG(I).AND.K.LE.KTCON(I)) THEN
              evef = EDT(I) * evfact
              if(SLIMSK(I).eq.1.) evef=EDT(I) * evfactl
!             if(SLIMSK(I).eq.1.) evef=.07
c             if(SLIMSK(I).ne.1.) evef = 0.
              QCOND(I) = EVEF * (Q1(I,k) - QESO(I,k))
     &                 / (1. + EL2ORC * QESO(I,k) / T1(I,k)**2)
              DP = 1000. * DEL(I,K)
              IF(RN(I).GT.0..AND.QCOND(I).LT.0.) THEN
                QEVAP(I) = -QCOND(I) * (1.-EXP(-.32*SQRT(DT2*RN(I))))
                QEVAP(I) = MIN(QEVAP(I), RN(I)*1000.*G/DP)
                DELQ2(I) = DELQEV(I) + .001 * QEVAP(I) * dp / g
              ENDIF
              if(RN(I).gt.0..and.QCOND(I).LT.0..and.
     &           DELQ2(I).gt.RNTOT(I)) THEN
                QEVAP(I) = 1000.* g * (RNTOT(I) - DELQEV(I)) / dp
                FLG(I) = .false.
              ENDIF
              IF(RN(I).GT.0..AND.QEVAP(I).gt.0.) THEN
                Q1(I,k) = Q1(I,k) + QEVAP(I)
                T1(I,k) = T1(I,k) - ELOCP * QEVAP(I)
                RN(I) = RN(I) - .001 * QEVAP(I) * DP / G
                DELTV(I) = - ELOCP*QEVAP(I)/DT2
                DELQ(I) =  + QEVAP(I)/DT2
                DELQEV(I) = DELQEV(I) + .001*dp*QEVAP(I)/g
              ENDIF
              DELLAQ(I,k) = DELLAQ(I,k) + DELQ(I) / XMB(I)
              DELQBAR(I) = DELQBAR(I) + DELQ(I)*DP/G
              DELTBAR(I) = DELTBAR(I) + DELTV(I)*DP/G
            ENDIF
          endif
        ENDDO
      ENDDO
c      IF(LAT.EQ.LATD.AND.lon.eq.lond.and.CNVFLG(I) ) THEN
c        PRINT *, '   DELLAH ='
c        PRINT 6003, (DELLAH(k)*XMB(I),K=1,KMAX)
c        PRINT *, '   DELLAQ ='
c        PRINT 6003, (HVAP*DELLAQ(I,k)*XMB(I),K=1,KMAX)
c        PRINT *, ' DELHBAR, DELQBAR, DELTBAR ='
c        PRINT *, DELHBAR, HVAP*DELQBAR, CP*DELTBAR
c        PRINT *, ' PRECIP =', HVAP*RN(I)*1000./DT2
CCCCC   PRINT *, '   DELLBAR ='
CCCCC   PRINT *,  HVAP*DELLbar
c      ENDIF
C
C  PRECIPITATION RATE CONVERTED TO ACTUAL PRECIP
C  IN UNIT OF M INSTEAD OF KG
C
      DO I = 1, IM
        IF(CNVFLG(I)) THEN
C
C  IN THE EVENT OF UPPER LEVEL RAIN EVAPORATION AND LOWER LEVEL DOWNDRAF
C    MOISTENING, RN CAN BECOME NEGATIVE, IN THIS CASE, WE BACK OUT OF TH
C    HEATING AND THE MOISTENING
C
          if(RN(I).lt.0..and..not.FLG(I)) RN(I) = 0.
          IF(RN(I).LE.0.) THEN
            RN(I) = 0.
          ELSE
            KTOP(I) = KTCON(I)
            KBOT(I) = KBCON(I)
            KUO(I) = 1
            CLDWRK(I) = AA1(I)
          ENDIF
        ENDIF
      ENDDO
      DO K = 1, KM
        DO I = 1, IM
          if (k .le. kmax(i)) then
            IF(CNVFLG(I).AND.RN(I).LE.0.) THEN
              T1(I,k) = TO(I,k)
              Q1(I,k) = QO(I,k)
            ENDIF
          endif
        ENDDO
      ENDDO
!!
      RETURN
      END
