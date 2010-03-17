!----------------------------------------------------------------------
! !MODULE: GEFS_CplComp_ESMFMod
!        --- ESMF coupler gridded component of the GFS Ensemble 
!            Forecast Operational system. 
!
! !DESCRIPTION: GFS coupler gridded component main module.
!
! !REVISION HISTORY:
!
!  April    2006     Weiyu Yang Initial code.
!                           
!
!  March 2007         Dingchen Hou added Stochatic Perturbation Combination Coefficient array.
!  January to November 2007      Dingchen and Weiyu Yang
!                     Added Broadcasting procedure and Global variables/arrays
!  November 2007      Dingchen, added minimum documentation, mainly for the arrays added during 2007
!  March    2009      Weiyu yang, modified for the NEMS model.

! !INTERFACE:
!
 MODULE GEFS_CplComp_ESMFMod
 
!!USES:
!------

! Define the ESMF internal state and all routines related to run 
! the GFS ensemble coupler grid component.
!---------------------------------------------------------------
 USE GEFS_Cpl_ESMFMod

 IMPLICIT none

#include "ESMF_LogMacros.inc"

 PRIVATE   ! By default data is private to this module
!
! !PUBLIC TYPES:
!---------------

 PUBLIC GEFS_CplCompSetServices

!EOP
!-------------------------------------------------------------------------


 CONTAINS


!----------------------------------------------------------------------
!BOP
!
! !ROUTINE: GEFS_CplCompSetServices --- Set services for GFS Ensemble 
!                                       Coupler Gridded Component.
! 
! !INTERFACE:
!
 SUBROUTINE GEFS_CplCompSetServices(CplGEFS, rc)

! !ARGUMENTS:
!------------

 TYPE(ESMF_CplComp),  INTENT(inout) :: CplGEFS ! gridded component
 INTEGER,             INTENT(out)   :: rc      ! return code
     
! !DESCRIPTION: Set services (register) for the GFS Ensemble Coupler
!               Grid Component.
!         
!EOP         
!----------------------------------------------------------------------
  
 INTEGER                            :: rc1     = ESMF_SUCCESS
 rc = ESMF_SUCCESS

! REGISTER SERVICES FOR THIS COMPONENT
! ------------------------------------

     CALL ESMF_LogWrite("Set Entry Point for Cpl Initialize",             &
                        ESMF_LOG_INFO, rc = rc1)

 CALL ESMF_CplCompSetEntryPoint (CplGEFS, ESMF_SETINIT,  Cpl_Initialize, &
                                 ESMF_SINGLEPHASE, rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Set Entry Point for Cpl Initialize")) THEN
         rc = ESMF_FAILURE
         PRINT*, 'Error Happened When Setting the Entry Point for Cpl Initialize, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

     CALL ESMF_LogWrite("Set Entry Point for Cpl Run",                    &
                        ESMF_LOG_INFO, rc = rc1)

 CALL ESMF_CplCompSetEntryPoint (CplGEFS, ESMF_SETRUN,   Cpl_Run,        &
                                 ESMF_SINGLEPHASE, rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Set Entry Point for Cpl Run")) THEN
         rc = ESMF_FAILURE
         PRINT*, 'Error Happened When Setting the Entry Point for Cpl Run, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

     CALL ESMF_LogWrite("Set Entry Point for Cpl Finalize",               &
                        ESMF_LOG_INFO, rc = rc1)

 CALL ESMF_CplCompSetEntryPoint (CplGEFS, ESMF_SETFINAL, Cpl_Finalize,   &
                                 ESMF_SINGLEPHASE, rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Set Entry Point for Cpl Finalize")) THEN
         rc = ESMF_FAILURE
         PRINT*, 'Error Happened When Setting the Entry Point for Cpl Finalize, rc = ', rc1
     END IF

 IF(rc == ESMF_SUCCESS) THEN
     PRINT*, "PASS: GEFS_CplCompSetServices."
 ELSE
     PRINT*, "FAIL: GEFS_CplCompSetServices."
 END IF

 END SUBROUTINE GEFS_CplCompSetServices





!----------------------------------------------------------------------
!BOP
! !ROUTINE:  Cpl_Initialize --- initialize routine to initialize 
!                               and set up the GFS ensemble coupler.
!
! !DESCRIPTION: This subroutine initializes the GFS ensemble coupler
!               before the main running routine.
!
!
! !REVISION HISTORY:
!
!  April    2006     Weiyu Yang Initial code.
!
! !INTERFACE:
!

 SUBROUTINE Cpl_Initialize(CplGEFS, impGEFS, expGEFS, clock, rcfinal)

 USE GEFS_GetParameterFromStateMod
 
!
! !INPUT/OUTPUT VARIABLES AND PARAMETERS:
!----------------------------------------

 TYPE(ESMF_CplComp), INTENT(inout)     :: CplGEFS 
 TYPE(ESMF_State),   INTENT(inout)     :: impGEFS
 TYPE(ESMF_State),   INTENT(inout)     :: expGEFS
 TYPE(ESMF_Clock),   INTENT(inout)     :: clock

!
! !OUTPUT VARIABLES AND PARAMETERS:
!----------------------------------

 INTEGER,             INTENT(out)       :: rcfinal

! !EOP
!------------------------------------------------------------------------- 
 
! !WORKING ARRAYS AND LOCAL PARAMETERS.
!--------------------------------------
 TYPE(ESMF_VM)                          :: vm
 TYPE(GEFS_Cpl_wrap)                    :: wrap
 TYPE(GEFS_Cpl_InternalState), POINTER  :: Cpl_Int_State
 TYPE(ESMF_Config)                      :: Cf
 TYPE(ESMF_State)                       :: imp_atm
 TYPE(ESMF_State)                       :: exp_atm
 TYPE(ESMF_State)                       :: imp_driver
 TYPE(ESMF_State)                       :: exp_driver
 TYPE(ESMF_State)                       :: imp_gfs_dyn
 TYPE(ESMF_State)                       :: exp_gfs_dyn
 CHARACTER(ESMF_MAXSTR)                 :: Cf_fname
 CHARACTER(12)                          :: PELAB

 INTEGER                                :: i, j, k, l
 INTEGER                                :: rc1 

 rc1     = ESMF_SUCCESS
 rcfinal = ESMF_SUCCESS

!These are the standard ESMF internal state lines.
!-------------------------------------------------
     CALL ESMF_LogWrite("Allocate the Cpl Internal State",                &
                    ESMF_LOG_INFO, rc = rc1)

 ALLOCATE(Cpl_Int_State, stat = rc1)

     IF(ESMF_LogMsgFoundAllocError(rc1, " - Allocate the Cpl Internal State")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Allocating the Cpl Internal State, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 wrap%Cpl_Int_State => Cpl_Int_State

! Attach internal state to the GFS Ensemble Coupler grid component.
!------------------------------------------------------------------
     CALL ESMF_LogWrite("Set Up the GFS Ensemble coupler Internal State", &
                        ESMF_LOG_INFO, rc = rc1)

 CALL ESMF_CplCompSetInternalState(CplGEFS, wrap, rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Set Up the  GFS Ensemble coupler Internal State")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Setting Up the  GFS Ensemble coupler Internal State, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

! Get the coupler grid component vm.
!-----------------------------------
     CALL ESMF_LogWrite("Get the GLobal VM", ESMF_LOG_INFO, rc = rc1)

 CALL ESMF_VMGetGlobal(vm, rc = rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Get the VM")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting the VM, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

! Set up parameters of MPI communications.
! Use ESMF utility to get PE identification and total number of PEs.
!-------------------------------------------------------------------
     CALL ESMF_LogWrite("Get me and NODES from VM", ESMF_LOG_INFO, rc = rc1)

 CALL ESMF_VMGet(vm, localPet = Cpl_Int_State%me,    &
                     petCount = Cpl_Int_State%nodes, &
                     rc       = rc1)
 Cpl_Int_State%mm1 = Cpl_Int_State%me + 1

     IF(ESMF_LogMsgFoundError(rc1, "Get me and NODES from VM")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting me and NODES from VM, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

! Set up parameters from the configure file.
!-------------------------------------------
     CALL ESMF_LogWrite("Get parameters from the configure file", &
                         ESMF_LOG_INFO, rc = rc1)

 Cf       = ESMF_ConfigCreate(rc = rc1)
 Cf_fname = 'atm_namelist.rc'

 CALL ESMF_ConfigLoadFile(Cf, Cf_fname, rc = rc1)

 CALL ESMF_ConfigGetAttribute(Cf,                         &
                              Cpl_Int_State%Total_member, &
                              label = 'total_member:',    &
                              rc    = rc1)
 CALL ESMF_ConfigGetAttribute(Cf,                         &
                              Cpl_Int_State%hh_increase,  &
                              label = 'HH_INCREASE:',     &
                              rc    = rc1)
 CALL ESMF_ConfigGetAttribute(Cf,                         &
                              Cpl_Int_State%hh_start,  &
                              label = 'HH_START:',     &
                              rc    = rc1)
 CALL ESMF_ConfigGetAttribute(Cf,                         &
                              Cpl_Int_State%hh_final,  &
                              label = 'HH_FINAL:',     &
                              rc    = rc1)

 CALL ESMF_ConfigDestroy(Cf, rc = rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Get parameters from the configure file")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting parameters from the configure file, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

! Set up the pointer links for the GEFS export state.
!----------------------------------------------------
 CALL ESMF_StateGet(impGEFS,    "ATM Driver Export State", exp_driver,  rc = rc1)
 CALL ESMF_StateGet(exp_driver, "atm export state",        exp_atm,     rc = rc1)
 CALL ESMF_StateGet(exp_atm,    "dynamics export",         exp_gfs_dyn, rc = rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Get the nested state from the import ESMF state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting the nested state from the import ESMF state, rc = ', rc1
     END IF

 CALL ESMF_StateGet(expGEFS,    "ATM Driver Import State", imp_driver,  rc = rc1)
 CALL ESMF_StateGet(imp_driver, "atm import state",        imp_atm,     rc = rc1)
 CALL ESMF_StateGet(imp_atm,    "dynamics import",         imp_gfs_dyn, rc = rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Get the nested state from the export ESMF state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting the nested state from the export ESMF state, rc = ', rc1
     END IF

! Get required parameters from the ESMF import state to do the spectral 
! transform processes for the second step of the stochastic perturbation 
! scheme.
!-----------------------------------------------------------------------
 NULLIFY (Cpl_Int_State%t   )
 NULLIFY (Cpl_Int_State%u   )
 NULLIFY (Cpl_Int_State%v   )
 NULLIFY (Cpl_Int_State%ps  )

 NULLIFY (Cpl_Int_State%tw  )
 NULLIFY (Cpl_Int_State%uw  )
 NULLIFY (Cpl_Int_State%vw  )
 NULLIFY (Cpl_Int_State%psw )

 NULLIFY (Cpl_Int_State%tm  )
 NULLIFY (Cpl_Int_State%um  )
 NULLIFY (Cpl_Int_State%vm  )
 NULLIFY (Cpl_Int_State%psm )

 NULLIFY (Cpl_Int_State%twm )
 NULLIFY (Cpl_Int_State%uwm )
 NULLIFY (Cpl_Int_State%vwm )
 NULLIFY (Cpl_Int_State%pswm)

 NULLIFY (Cpl_Int_State%t6  )
 NULLIFY (Cpl_Int_State%u6  )
 NULLIFY (Cpl_Int_State%v6  )
 NULLIFY (Cpl_Int_State%ps6 )

 NULLIFY (Cpl_Int_State%t6m )
 NULLIFY (Cpl_Int_State%u6m )
 NULLIFY (Cpl_Int_State%v6m )
 NULLIFY (Cpl_Int_State%ps6m)

 NULLIFY (Cpl_Int_State%tracer  )
 NULLIFY (Cpl_Int_State%tracerw )
 NULLIFY (Cpl_Int_State%tracerm )
 NULLIFY (Cpl_Int_State%tracerwm)
 NULLIFY (Cpl_Int_State%tracer6 )
 NULLIFY (Cpl_Int_State%tracer6m)

 CALL GEFS_GetParameterFromState(exp_gfs_dyn, Cpl_Int_State, rc = rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Get the parameters from the ESMF state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting the parameters from the ESMF state, rc = ', rc1
     END IF

! Currently ntrac = 3.
!---------------------
 ALLOCATE(Cpl_Int_State%vname(Cpl_Int_State%ntrac, 5))

 Cpl_Int_State%vname(1, 1) = 'spfh'
 Cpl_Int_State%vname(1, 2) = 'spfh_q'
 Cpl_Int_State%vname(1, 3) = 'spfh_m'
 Cpl_Int_State%vname(1, 4) = 'spfh_q6'
 Cpl_Int_State%vname(1, 5) = 'spfh_m6'
 Cpl_Int_State%vname(2, 1) = 'o3mr'
 Cpl_Int_State%vname(2, 2) = 'o3mr_q'
 Cpl_Int_State%vname(2, 3) = 'o3mr_m'
 Cpl_Int_State%vname(2, 4) = 'o3mr_q6'
 Cpl_Int_State%vname(2, 5) = 'o3mr_m6'
 Cpl_Int_State%vname(3, 1) = 'clwmr'
 Cpl_Int_State%vname(3, 2) = 'clwmr_q'
 Cpl_Int_State%vname(3, 3) = 'clwmr_m'
 Cpl_Int_State%vname(3, 4) = 'clwmr_q6'
 Cpl_Int_State%vname(3, 5) = 'clwmr_m6'

     CALL ESMF_LogWrite("Cpl Internal State to ESMF Export State", &
                        ESMF_LOG_INFO, rc = rc1)

 CALL GEFS_Cpl_InternalState2ESMFExportState(exp_gfs_dyn, imp_gfs_dyn, Cpl_Int_State, rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Cpl Internal State to ESMF Export State")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Cpl Setting up the ESMF Export State, rc = ', rc1
     END IF

!DHOU, 11/13/2007 assignement of latmax (Number of Gaussian latitudes following
! setlonsgg subroutine in GFS_Initialize_ESMFMod.f)
!------------------------------------------------------------------------------
 IF(Cpl_Int_State%jcap == 62)THEN
     Cpl_Int_State%latmax = 94
 END IF
 IF(Cpl_Int_State%jcap == 126)THEN
     Cpl_Int_State%latmax = 190
 END IF
 IF(Cpl_Int_State%jcap == 170)THEN
     Cpl_Int_State%latmax = 256
 END IF
 IF(Cpl_Int_State%jcap == 190)THEN
     Cpl_Int_State%latmax = 288
 END IF
 IF(Cpl_Int_State%jcap == 254)THEN
     Cpl_Int_State%latmax = 384
 END IF
 IF(Cpl_Int_State%jcap == 382)THEN
     Cpl_Int_State%latmax = 576
 END IF
 IF(Cpl_Int_State%jcap == 510)THEN
     Cpl_Int_State%latmax = 766
 END IF

!DHOU, 09/25/2007 set up the Gaussian latitude (sin of it) slat_work and wlat_work  array
 ALLOCATE(Cpl_Int_State%slat_work(Cpl_Int_State%latmax))
 ALLOCATE(Cpl_Int_State%wlat_work(Cpl_Int_State%latmax))
 CALL SPLAT(4, Cpl_Int_State%latmax, Cpl_Int_State%slat_work, Cpl_Int_State%wlat_work)

!Set up the Cpl_Int_State%factor1_work array and initialize it with values of 1.0.
 Cpl_Int_State%nreg=3
 ALLOCATE(Cpl_Int_State%factor1_work(Cpl_Int_State%nreg,Cpl_Int_State%Total_member))
 DO i=1,Cpl_Int_State%nreg
   DO j=1,Cpl_Int_State%Total_member
        Cpl_Int_State%factor1_work(i,j)=1.0
   ENDDO
 ENDDO

! Set up the member_id global arraycw and populate it by local assignments.
! This is a GLOBAL array,specifying the ens_member for processed by each cpu
!----------------------
 ALLOCATE(Cpl_Int_State%member_id(Cpl_Int_State%nodes))
 i = Cpl_Int_State%nodes / Cpl_Int_State%Total_member
 DO k = 1, Cpl_Int_State%Total_member
     DO j = 1, i
         Cpl_Int_State%member_id((k - 1) * i + j) = k
     END DO
 END DO

! Initialize the number of times of calling the coupler run routine and the combination coefficients.
!-------------------------------------------------------------------
 Cpl_Int_State%Cpl_Run_Calling_Number = Cpl_Int_State%hh_start / Cpl_Int_State%hh_increase + 1
 Cpl_Int_State%Cpl_Run_Calling_Start  = Cpl_Int_State%hh_start / Cpl_Int_State%hh_increase + 1
 Cpl_Int_State%Cpl_Run_Calling_Final  = Cpl_Int_State%hh_final / Cpl_Int_State%hh_increase
 ALLOCATE(Cpl_Int_State%Sto_Coef(Cpl_Int_State%Total_member-1, Cpl_Int_State%Total_member-1))
 Cpl_Int_State%Sto_Coef = 0.0

!DHOU, 09/25/2007
! Set up the global array for the lats_node local variable (lats_node_global, demension 90 for 90cpu)
! this GLOBAL array contains the number of latitudes in each cpu
! The values are 31 or 32 for (T126, 15 members, 90 cpus 190/(90/15)=31.666)
!----------------------------------------------------------------------------------------------------
 ALLOCATE(Cpl_Int_State%lats_node_global(Cpl_Int_State%nodes))

 DO i = 1, Cpl_Int_State%nodes
     IF(Cpl_Int_State%mm1 == i) Cpl_Int_State%lats_node_global(i) = Cpl_Int_State%lats_node_a
     CALL GEFS_bcst_global_i4(Cpl_Int_State%lats_node_global(i), i - 1, rc1)
 END DO

! Set up the global lats_global array. (lats_global (32,90) for 90cpu, 15 members)
! this GLOBAL array contains the index of the latitudes in each latitude of each cpu
! (in global order, value between 1 and 190 for T126)
! Or the other hand, the lats_node_global(:) array takes value between 1 and 31/32.
!------------------------------------------------------------------------------------
 k = 0
 DO i = 1, Cpl_Int_State%nodes
     k = MAX(k, Cpl_Int_State%lats_node_global(i))
 END DO
 ALLOCATE(Cpl_Int_State%lats_global(k, Cpl_Int_State%nodes))

!populate the Cpl_Int_State%lats_global array by assignment at one cpu and BROADCASTING to all cpus.
!---------------------------------------------------------------------------------------------------
 DO i = 1, Cpl_Int_State%nodes
      DO j = 1, Cpl_Int_State%lats_node_global(i)
          IF(Cpl_Int_State%mm1 == i) Cpl_Int_State%lats_global(j, i) = &
              Cpl_Int_State%global_lats_a(Cpl_Int_State%ipt_lats_node_a - 1 + j)
          CALL GEFS_bcst_global_i4(Cpl_Int_State%lats_global(j, i), i - 1, rc1)
      END DO
 END DO

 IF(rcfinal == ESMF_SUCCESS) THEN
     PRINT*, "PASS: Cpl_Initialize."
 ELSE
     PRINT*, "FAIL: Cpl_Initialize."
 END IF

 END SUBROUTINE Cpl_Initialize





!----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Cpl_Run --- Main grid component routine to run the GFS 
!                       ensemble coupler.
!
! !DESCRIPTION: This subroutine will run the most part computations 
!               of the GFS ensemble coupler.
!
! !REVISION HISTORY:
!
!  April    2006     Weiyu Yang Initial code.
!
! !INTERFACE:
!

 SUBROUTINE Cpl_Run(CplGEFS, impGEFS, expGEFS, clock, rcfinal)

!
! !INPUT VARIABLES AND PARAMETERS:
!---------------------------------
 TYPE(ESMF_CplComp), INTENT(inout)     :: CplGEFS   
 TYPE(ESMF_State),   INTENT(in)        :: impGEFS 
 
! !OUTPUT VARIABLES AND PARAMETERS:
!----------------------------------
 TYPE(ESMF_Clock),   INTENT(inout)     :: clock
 TYPE(ESMF_State),   INTENT(inout)     :: expGEFS
 INTEGER,            INTENT(out)       :: rcfinal 
!
!EOP
!-------------------------------------------------------------------------

!
! !WORKING ARRAYS AND LOCAL PARAMETERS.
!--------------------------------------
 TYPE(ESMF_State)                       :: imp_atm
 TYPE(ESMF_State)                       :: exp_atm
 TYPE(ESMF_State)                       :: imp_driver
 TYPE(ESMF_State)                       :: exp_driver
 TYPE(ESMF_State)                       :: imp_gfs_dyn
 TYPE(ESMF_State)                       :: exp_gfs_dyn
 TYPE(ESMF_VM)                         :: vm
 TYPE(GEFS_Cpl_wrap)                   :: wrap
 TYPE(GEFS_Cpl_InternalState), POINTER :: Cpl_Int_State
 INTEGER                               :: rc1

 rc1     = ESMF_SUCCESS
 rcfinal = ESMF_SUCCESS

! Retrieve the ESMF internal state.
!---------------------------------- 
     CALL ESMF_LogWrite("Get the Internal State in the Cpl Run Routine", &
                        ESMF_LOG_INFO, rc = rc1)
 CALL ESMF_CplCompGetInternalState(CplGEFS, wrap, rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Get the Internal State in the Cpl Run Routine")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting the Internal State in Cpl Run Routine, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

 Cpl_Int_State => wrap%Cpl_Int_State

 CALL ESMF_VMGetGlobal(vm, rc = rc1)
 CALL ESMF_VMBarrier  (vm, rc = rc1)

! Transfer the GFS export fields to the working arrays in the Cpl internal state.
!--------------------------------------------------------------------------------
     CALL ESMF_LogWrite("ESMF import State to the Cpl Internal State", &
                        ESMF_LOG_INFO, rc = rc1)

 CALL ESMF_StateGet(impGEFS,    "ATM Driver Export State", exp_driver,  rc = rc1)
 CALL ESMF_StateGet(exp_driver, "atm export state",        exp_atm,     rc = rc1)
 CALL ESMF_StateGet(exp_atm,    "dynamics export",         exp_gfs_dyn, rc = rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Get the nested state from the import ESMF state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting the nested state from the import ESMF state, rc = ', rc1
     END IF

 CALL ESMF_StateGet(expGEFS,    "ATM Driver Import State", imp_driver,  rc = rc1)
 CALL ESMF_StateGet(imp_driver, "atm import state",        imp_atm,     rc = rc1)
 CALL ESMF_StateGet(imp_atm,    "dynamics import",         imp_gfs_dyn, rc = rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Get the nested state from the export ESMF state.")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Getting the nested state from the export ESMF state, rc = ', rc1
     END IF

 CALL GEFS_Cpl_ESMFImportState2InternalState(exp_gfs_dyn, Cpl_Int_State, rc = rc1)

     IF(ESMF_LogMsgFoundError(rc1, "ESMF import State to the Cpl Internal State")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Moving the ESMF import State to the Cpl Internal State, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

!Run the GEFS_Cpl_Run, to create new initial conditions for the next ensemble run.
!---------------------------------------------------------------------------------
     CALL ESMF_LogWrite("Run the GEFS_Cpl_Run", ESMF_LOG_INFO, rc = rc1)

 CALL ESMF_VMBarrier(vm, rc = rc1)
 CALL GEFS_Cpl_Run(exp_gfs_dyn, clock, Cpl_Int_State,  rc = rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Run the GEFS_Cpl_Run")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Running the GEFS_Cpl_Run, rc = ', rc1
         rc1 = ESMF_SUCCESS
     END IF

! Update the number of times of calling the coupler run routine.
!---------------------------------------------------------------
 Cpl_Int_State%Cpl_Run_Calling_Number = Cpl_Int_State%Cpl_Run_Calling_Number + 1

! Transfer the new initial conditions to the ESMF Cpl exprot state.
!------------------------------------------------------------------
     CALL ESMF_LogWrite("Cpl Internal State to ESMF Export State", &
                        ESMF_LOG_INFO, rc = rc1)

 CALL ESMF_VMBarrier(vm, rc = rc1)
 CALL GEFS_Cpl_InternalState2ESMFExportState(exp_gfs_dyn, imp_gfs_dyn, Cpl_Int_State, rc = rc1)

     IF(ESMF_LogMsgFoundError(rc1, "Cpl Internal State to ESMF Export State")) THEN
         rcfinal = ESMF_FAILURE
         PRINT*, 'Error Happened When Cpl Setting up the ESMF Export State, rc = ', rc1
     END IF

!End of the grid component run.
!------------------------------
 IF(rcfinal == ESMF_SUCCESS) THEN
     PRINT*, "PASS: Cpl_Run."
 ELSE
     PRINT*, "FAIL: Cpl_Run."
 END IF

 END SUBROUTINE Cpl_Run





!----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Cpl_Finalize --- finalizing routine to finish the 
!                            GFS ensemble coupler.
!
! !DESCRIPTION: This subroutine will finish the GFS ensemble coupler
! !             and will release the memory space.
!
! !REVISION HISTORY:
!
!  April    2006     Weiyu Yang Initial code.
!
! !INTERFACE:

 SUBROUTINE Cpl_Finalize(CplGEFS, impGEFS, expGEFS, clock, rcfinal)

!
! !INPUT VARIABLES AND PARAMETERS:
!---------------------------------
 TYPE(ESMF_CplComp), INTENT(inout)  :: CplGEFS
 TYPE(ESMF_State),   INTENT(inout)  :: impGEFS
 TYPE(ESMF_State),   INTENT(inout)  :: expGEFS
 TYPE(ESMF_Clock),   INTENT(inout)  :: clock

! !OUTPUT VARIABLES AND PARAMETERS:
!----------------------------------
 INTEGER,            INTENT(out)    :: rcfinal

! !WORKING ARRAYS AND LOCAL PARAMETERS.
!--------------------------------------
 TYPE(GEFS_Cpl_wrap)                   :: wrap
 TYPE(GEFS_Cpl_InternalState), POINTER :: Cpl_Int_State
 INTEGER                               :: rc1

!EOP
!-------------------------------------------------------------------------

 rc1     = ESMF_SUCCESS
 rcfinal = ESMF_SUCCESS

! Retrieve the ESMF internal state.
!----------------------------------
!    CALL ESMF_LogWrite("Get the Internal State in the Cpl Finalize Routine", &
!                       ESMF_LOG_INFO, rc = rc1)

!CALL ESMF_CplCompGetInternalState(CplGEFS, wrap, rc1)

!    IF(ESMF_LogMsgFoundError(rc1, "Get the Internal State in the Cpl Finalize Routine")) THEN
!        rcfinal = ESMF_FAILURE
!        PRINT*, 'Error Happened When Getting the Internal State in Cpl FInalize Routine, rc = ', rc1
!    END IF

!Cpl_Int_State => wrap%Cpl_Int_State

 IF(rcfinal == ESMF_SUCCESS) THEN
     PRINT*, "PASS: Cpl_Finalize."
 ELSE
     PRINT*, "FAIL: Cpl_Finalize."
 END IF

 END SUBROUTINE Cpl_Finalize

 END MODULE GEFS_CplComp_ESMFMod