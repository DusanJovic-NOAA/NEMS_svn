!-----------------------------------------------------------------------
!
      MODULE MODULE_WRITE_ROUTINES
!
!-----------------------------------------------------------------------
!***  THIS MODULE CONTAINS ROUTINES NEEDED BY THE RUN STEP OF THE
!***  WRITE GRIDDED COMPONENT IN WHICH HISTORY OUTPUT DATA FROM
!***  THE FORECAST TASKS ARE ASSEMBLED AND WRITTEN TO HISTORY FILES
!***  BY THE WRITE TASKS.
!-----------------------------------------------------------------------
!***
!***  HISTORY   
!***
!       06 Sep 2007:  T. Black - Created module.
!                                Moved the FIRST block from the old
!                                write component here for clarity.
!       15 Aug 2008:  J. Wang  - Revised for NEMS-IO
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
      USE MODULE_WRITE_INTERNAL_STATE
!
      USE MODULE_ATM_INTERNAL_STATE,ONLY: ATM_INTERNAL_STATE
!
      USE MODULE_DM_PARALLEL,ONLY : PARA_RANGE,MPI_COMM_INTER_ARRAY     &
                                   ,LOCAL_ISTART,LOCAL_IEND             &
                                   ,LOCAL_JSTART,LOCAL_JEND
!
      USE MODULE_CONTROL,ONLY : TIMEF
!
      USE MODULE_INCLUDE
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
      USE NEMSIO_MODULE
!
      USE MODULE_CONSTANTS,ONLY : A,PI,G
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------

      PRIVATE
!
      PUBLIC :: FIRST_PASS                                              &
               ,WRITE_ASYNC                                             &
               ,WRITE_INIT                                              &
               ,WRITE_NEMSIO_RUNHISTORY                                 &
               ,WRITE_NEMSIOCTL                                         &
               ,WRITE_RUNHISTORY
!
!-----------------------------------------------------------------------
!
      INTEGER :: MYPE       !<-- My MPI task ID
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE FIRST_PASS(IMP_STATE_WRITE                             &
                           ,HISTORY_BUNDLE                              &
                           ,WRT_INT_STATE                               &
                           ,NTASKS                                      &
                           ,MYPE                                        &
                           ,NCURRENT_GROUP                              &
                            )
! 
!-----------------------------------------------------------------------
!***  EACH TIME A NEW GROUP OF WRITE TASKS IS INVOKED FOR THE FIRST
!***  TIME THIS ROUTINE WILL PERFORM CERTAIN COMPUTATIONS AND TASKS
!***  THAT ONLY NEED TO BE DONE ONCE FOR EACH WRITE GROUP.
!***  THE ROUTINE WILL UNLOAD SOME VALUES FROM THE WRITE COMPONENT'S
!***  IMPORT STATE.  BECAUSE THESE QUANTITIES DO NOT CHANGE WITH
!***  FORECAST TIME, THE ROUTINE IS NOT NEEDED FOR SUBSEQUENT USE
!***  OF EACH WRITE GROUP.  THE DATA BEING UNLOADED CONSISTS OF
!***  EVERYTHING EXCEPT THE 2D/3D GRIDDED FORECAST ARRAYS.
!***  ALSO BASIC INFORMATION IS PROVIDED TO FORECAST AND WRITE TASKS
!***  THAT WILL BE NECESSARY IN THE QUILTING OF LOCAL 2-D GRIDDED
!***  HISTORY DATA INTO FULL DOMAIN 2-D FIELDS.
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: NTASKS
      INTEGER,INTENT(IN) :: MYPE
      INTEGER,INTENT(IN) :: NCURRENT_GROUP
!
      TYPE(ESMF_State)          ,INTENT(INOUT) :: IMP_STATE_WRITE
      TYPE(ESMF_FieldBundle)    ,INTENT(INOUT) :: HISTORY_BUNDLE
      TYPE(WRITE_INTERNAL_STATE),INTENT(INOUT) :: WRT_INT_STATE
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER,SAVE                 :: ITS,ITE,JTS,JTE
!
      INTEGER                      :: I,IERR,IM,J,JM,L                  &
                                     ,N,NN,NUM_ATTRIB,NWTPG             &
                                     ,RC,RC_WRT
!
      INTEGER,DIMENSION(:),POINTER :: INPES,JNPES
      INTEGER,DIMENSION(:),POINTER :: IHALO,JHALO
!
      INTEGER                      :: JROW_FIRST,JROW_LAST,JROWS
!
      INTEGER                      :: LAST_FCST_TASK                    &
                                     ,LEAD_WRITE_TASK                   &
                                     ,LAST_WRITE_TASK
!
      INTEGER,SAVE                 :: NCHAR_I1D                         &
                                     ,NCHAR_R1D                         &
                                     ,NCHAR_LOG
!
      INTEGER                      :: JEND_WRITE,JSTA_WRITE
!
      INTEGER                      :: KOUNT_I1D                         &
                                     ,KOUNT_I2D                         &
                                     ,KOUNT_R1D                         &
                                     ,KOUNT_R2D                         &
                                     ,KOUNT_LOG
!
      INTEGER                      :: LENGTH                            &
                                     ,LENGTH_SUM_I1D                    &
                                     ,LENGTH_SUM_R1D                    &
                                     ,LENGTH_SUM_LOG
!
      INTEGER                      :: NPOSN_START,NPOSN_END
!
      INTEGER                      :: NUM_FIELD_NAMES                   &
                                     ,NUM_PES_FCST
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      INTEGER,DIMENSION(:),POINTER :: LOCAL_ISTART                      &
                                     ,LOCAL_IEND                        &
                                     ,LOCAL_JSTART                      &
                                     ,LOCAL_JEND
!
      INTEGER,DIMENSION(:),POINTER :: NCHAR_I2D                         &
                                     ,NCHAR_R2D
!
      INTEGER,DIMENSION(:),POINTER :: WORK_ARRAY_I1D
!
      REAL(4),DIMENSION(:),POINTER :: WORK_ARRAY_R1D
!
      CHARACTER(ESMF_MAXSTR)       :: ATTRIB_NAME
!
      TYPE(ESMF_TypeKind)          :: DATATYPE
!
      TYPE(ESMF_Field)             :: FIELD_WORK1
!
      TYPE(ESMF_Logical)              :: WORK_LOGICAL
      TYPE(ESMF_Logical),DIMENSION(1) :: NO_FIELDS
!
      TYPE(ESMF_VM)                :: VM
!
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT)              :: btim,btim0
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  FIRST WE NEED THE NUMBER OF WRITE TASKS IN EACH GROUP.
!-----------------------------------------------------------------------
!
      NWTPG=wrt_int_state%WRITE_TASKS_PER_GROUP
!
      LAST_FCST_TASK =NTASKS-NWTPG-1
      LEAD_WRITE_TASK=LAST_FCST_TASK+1
      LAST_WRITE_TASK=NTASKS-1
!
!-----------------------------------------------------------------------
!***  THE INTEGER QUANTITIES 'INPES' AND 'JNPES' MUST BE POINTERS
!***  OF LENGTH 1 SINCE THEY WILL BE PASSED THROUGH ESMF Send/Recv.
!***  LIKEWISE WITH THE TOTAL LENGTH OF ALL NAMES OF 2D DATA 
!***  AND THE HALO DEPTHS.
!-----------------------------------------------------------------------
!
      ALLOCATE(INPES(1))
      ALLOCATE(JNPES(1))
      ALLOCATE(IHALO(1))
      ALLOCATE(JHALO(1))
      ALLOCATE(NCHAR_I2D(1))
      ALLOCATE(NCHAR_R2D(1))
!
!-----------------------------------------------------------------------
!***  EXTRACT THE FULL DOMAIN LIMITS FROM THE COMPONENT'S IMPORT 
!***  STATE.  THESE ARE NEEDED FOR ALLOCATING THE WORKING ARRAYS 
!***  THAT WILL MOVE THE 2-D AND 3-D FIELDS FROM THE IMPORT TO THE
!***  EXPORT STATE.
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      domain_limits: IF(MYPE<=LAST_FCST_TASK)THEN                          !<-- This selects only forecast tasks to do extractions
                                                                           !    since only they know what is in the import state
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE ARRAYS THAT HOLD ALL START AND END POINTS IN I,J
!***  FOR ALL FORECAST TASKS.
!-----------------------------------------------------------------------
!
        ALLOCATE(LOCAL_ISTART(0:LAST_FCST_TASK),stat=RC)
        ALLOCATE(LOCAL_IEND  (0:LAST_FCST_TASK),stat=RC)
        ALLOCATE(LOCAL_JSTART(0:LAST_FCST_TASK),stat=RC)
        ALLOCATE(LOCAL_JEND  (0:LAST_FCST_TASK),stat=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Global Parameters from History Bundle" 
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(bundle    =HISTORY_BUNDLE          &  !<-- The Bundle of history data
                              ,name      ='IM'                    &  !<-- Name of the Attribute to extract
                              ,count     =1                       &  !<-- Length of Attribute
                              ,valueList =wrt_int_state%IM        &  !<-- Extract this Attribute from History Bundle
                              ,rc        =RC)
!
        CALL ESMF_AttributeGet(bundle    =HISTORY_BUNDLE          &  !<-- The Bundle of history data
                              ,name      ='JM'                    &  !<-- Name of the Attribute to extract
                              ,count     =1                       &  !<-- Length of Attribute
                              ,valueList =wrt_int_state%JM        &  !<-- Extract this Attribute from History Bundle
                              ,rc        =RC)
!
        CALL ESMF_AttributeGet(bundle    =HISTORY_BUNDLE          &  !<-- The Bundle of history data
                              ,name      ='LM'                    &  !<-- Name of the Attribute to extract
                              ,count     =1                       &  !<-- Length of Attribute
                              ,valueList =wrt_int_state%LM        &  !<-- Extract this Attribute from History Bundle
                              ,rc        =RC)
!
        CALL ESMF_AttributeGet(bundle =HISTORY_BUNDLE             &  !<-- The Bundle of history data
                              ,name   ='GLOBAL'                   &  !<-- Name of the Attribute to extract
                              ,value  =wrt_int_state%GLOBAL       &  !<-- Extract this Attribute from History Bundle
                              ,rc     =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(wrt_int_state%GLOBAL==ESMF_TRUE)THEN                            !<-- Increase lateral dimensions by 2 for global runs
          wrt_int_state%IM(1)=wrt_int_state%IM(1)+2
          wrt_int_state%JM(1)=wrt_int_state%JM(1)+2
        ENDIF
!
!-----------------------------------------------------------------------
!***  NOW EXTRACT LOCAL SUBDOMAIN LIMITS.
!***  THESE WILL BE USED TO ALLOCATE THE WORKING ARRAY TO HOLD FIELDS
!***  ON EACH SUBDOMAIN PRIOR TO QUILTING THEM TOGETHER.
!***  WE FIRST NEED THE NUMBER OF FORECAST TASKS SINCE THAT
!***  DETERMINES THE SIZE OF THE ARRAYS HOLDING THE LOCAL
!***  SUBDOMAIN LIMITS.
!
!***  ALSO EXTRACT THE HALO DEPTHS SINCE THEY ARE NEEDED FOR
!***  EXCLUDING HALO POINTS FROM THE FINAL HISTORY DATA.
!
!***  THESE VALUES ARE NOT TO BE WRITTEN TO THE HISTORY FILES SO
!***  THEY WERE NOT INSERTED INTO THE HISTORY DATA Bundle INSIDE
!***  THE WRITE COMPONENT'S IMPORT STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Local Quilting Info from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state =IMP_STATE_WRITE             &  !<-- The Write component's import state
                              ,name  ='INPES'                     &  !<-- Name of the Attribute to extract
                              ,value =INPES(1)                    &  !<-- Extract this Attribute from import state
                              ,rc    =RC)
!
        CALL ESMF_AttributeGet(state =IMP_STATE_WRITE             &  !<-- The Write component's import state
                              ,name  ='JNPES'                     &  !<-- Name of the Attribute to extract
                              ,value =JNPES(1)                    &  !<-- Extract this Attribute from import state
                              ,rc    =RC)
!
        wrt_int_state%INPES=INPES(1)                                 !<-- Place in internal state for later use
        wrt_int_state%JNPES=JNPES(1)                                 !<-- Place in internal state for later use
!
        NUM_PES_FCST=INPES(1)*JNPES(1)                               !<-- Number of fcst tasks
!
        CALL ESMF_AttributeGet(state =IMP_STATE_WRITE             &  !<-- The Write component's import state
                              ,name  ='IHALO'                     &  !<-- Name of the Attribute to extract
                              ,value =IHALO(1)                    &  !<-- Extract this Attribute from import state
                              ,rc    =RC)
!
        CALL ESMF_AttributeGet(state =IMP_STATE_WRITE             &  !<-- The Write component's import state
                              ,name  ='JHALO'                     &  !<-- Name of the Attribute to extract
                              ,value =JHALO(1)                    &  !<-- Extract this Attribute from import state
                              ,rc    =RC)
!
        wrt_int_state%IHALO=IHALO(1)                                 !<-- Place in internal state for later use
        wrt_int_state%JHALO=JHALO(1)                                 !<-- Place in internal state for later use
!
        CALL ESMF_AttributeGet(state     =IMP_STATE_WRITE         &  !<-- The Write component's import state
                              ,name      ='LOCAL_ISTART'          &  !<-- Name of the Attribute to extract
                              ,count     =NUM_PES_FCST            &  !<-- Length of Attribute
                              ,valueList =LOCAL_ISTART            &  !<-- Extract local subdomain starting I's
                              ,rc=RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE          &  !<-- The Write component's import state
                              ,name     ='LOCAL_IEND'             &  !<-- Name of the Attribute to extract
                              ,count    =NUM_PES_FCST             &  !<-- Length of Attribute
                              ,valueList= LOCAL_IEND              &  !<-- Extract local subdomain ending I's
                              ,rc=RC)
!
        CALL ESMF_AttributeGet(state     =IMP_STATE_WRITE         &  !<-- The Write component's import state
                              ,name      ='LOCAL_JSTART'          &  !<-- Name of the Attribute to extract
                              ,count     =NUM_PES_FCST            &  !<-- Length of Attribute
                              ,valueList =LOCAL_JSTART            &  !<-- Extract local subdomain starting J's
                              ,rc=RC)
!
        CALL ESMF_AttributeGet(state     =IMP_STATE_WRITE         &  !<-- The Write component's import state
                              ,name      ='LOCAL_JEND'            &  !<-- Name of the Attribute to extract
                              ,count     =NUM_PES_FCST            &  !<-- Length of Attribute
                              ,valueList =LOCAL_JEND              &  !<-- Extract local subdomain ending J's
                              ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        DO N=0,LAST_FCST_TASK
          wrt_int_state%LOCAL_ISTART(N)=LOCAL_ISTART(N)
          wrt_int_state%LOCAL_IEND  (N)=LOCAL_IEND(N)
          wrt_int_state%LOCAL_JSTART(N)=LOCAL_JSTART(N)
          wrt_int_state%LOCAL_JEND  (N)=LOCAL_JEND(N)
        ENDDO
!
        ITS=LOCAL_ISTART(MYPE)
        ITE=LOCAL_IEND(MYPE)
        JTS=LOCAL_JSTART(MYPE)
        JTE=LOCAL_JEND(MYPE)
!
        DEALLOCATE(LOCAL_ISTART)
        DEALLOCATE(LOCAL_IEND  )
        DEALLOCATE(LOCAL_JSTART)
        DEALLOCATE(LOCAL_JEND  )
!
!-----------------------------------------------------------------------
!
      ENDIF domain_limits
!
!-----------------------------------------------------------------------
!***  FORECAST TASK 0 SENDS THE DOMAIN SIZE INFORMATION
!***  TO THE FIRST WRITE TASK IN EACH WRITE GROUP BECAUSE 
!***  THE WRITE TASKS NEED TO KNOW THIS TO ASSEMBLE THE
!***  FINAL GRIDDED DATA.
!***  FIRST WE NEED THE VM.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get the Current VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGetCurrent(vm=VM                                      &  !<-- The ESMF virtual machine for this group of tasks
                            ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!                            -- IM --
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN                                                      !<-- Forecast task 0 sends
        DO N=0,NWTPG-1
          CALL MPI_SEND(wrt_int_state%IM                                &  !<-- Send this data
                       ,1                                               &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,N                                               &  !<-- Send to each of the write tasks (local IDs)
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)            &  !<-- MPI communicator
                       ,IERR)
!
          IF(IERR/=0)WRITE(0,*)' Failed to send IM from fcst task0 to write tasks'
!
!       CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
!                       ,sendData=wrt_int_state%IM                      &  !<-- Send this data
!                       ,count   =1                                     &  !<-- Words sent
!                       ,dst     =LEAD_WRITE_TASK                       &  !<-- Receiving task (1st write task in group)
!                       ,rc      =RC)
        ENDDO
      ENDIF 

      IF(MYPE>=LEAD_WRITE_TASK)THEN                                        !<-- Write tasks in this group receive
        CALL MPI_RECV(wrt_int_state%IM                                  &  !<-- Recv this data
                     ,1                                                 &  !<-- Words received
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,0                                                 &  !<-- Recv from fcst 0
                     ,0                                                 &  !<-- An MPI tag
                     ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)              &  !<-- MPI communicator
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR)
!
        IF(IERR/=0)WRITE(0,*)' Write tasks failed to receive IM from fcst task0'
!
!       CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
!                       ,recvData=wrt_int_state%IM                      &  !<-- Recv this data
!                       ,count   =1                                     &  !<-- Words received
!                       ,src     =0                                     &  !<-- Sending task (forecast task 0)
!                       ,rc      =RC)
      ENDIF
!
!-----------------------------------------------------------------------
!                          -- JM --
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN                                                      !<-- Forecast task 0 sends
        DO N=0,NWTPG-1                                            
          CALL MPI_SEND(wrt_int_state%JM                                &  !<-- Send this data
                       ,1                                               &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,N                                               &  !<-- Send to each of the write tasks (local IDs)
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)            &  !<-- MPI communicator
                       ,IERR)
!
          IF(IERR/=0)WRITE(0,*)' Failed to send JM from fcst task0 to write tasks'
!
        ENDDO
      ENDIF 
!
      IF(MYPE>=LEAD_WRITE_TASK)THEN                                        !<-- Write tasks in this group receive
        CALL MPI_RECV(wrt_int_state%JM                                  &  !<-- Recv this data
                     ,1                                                 &  !<-- Words received
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,0                                                 &  !<-- Recv from fcst 0
                     ,0                                                 &  !<-- An MPI tag
                     ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)              &  !<-- MPI communicator
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR)
!
        IF(IERR/=0)WRITE(0,*)' Write tasks failed to receive JM from fcst task0'
!
      ENDIF
!
!-----------------------------------------------------------------------
!                          -- LM --
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN                                                      !<-- Forecast task 0 sends
        DO N=0,NWTPG-1                                            
          CALL MPI_SEND(wrt_int_state%LM                                &  !<-- Send this data
                       ,1                                               &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,N                                               &  !<-- Send to each of the write tasks (local IDs)
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)            &  !<-- MPI communicator
                       ,IERR)
!
          IF(IERR/=0)WRITE(0,*)' Failed to send LM from fcst task0 to write tasks'
!
        ENDDO
      ENDIF 
!
      IF(MYPE>=LEAD_WRITE_TASK)THEN                                        !<-- Write tasks in this group receive
        CALL MPI_RECV(wrt_int_state%LM                                  &  !<-- Recv this data
                     ,1                                                 &  !<-- Words received
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,0                                                 &  !<-- Recv from fcst 0
                     ,0                                                 &  !<-- An MPI tag
                     ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)              &  !<-- MPI communicator
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR)
!
        IF(IERR/=0)WRITE(0,*)' Write tasks failed to receive LM from fcst task0'
!
      ENDIF
!
      IM=wrt_int_state%IM(1)
      JM=wrt_int_state%JM(1)
!
!-----------------------------------------------------------------------
!***  THE NUMBER OF Attributes (FOR SCALARS AND 1D ARRAYS) AND
!***  Fields (FOR GRIDDED 2D ARRAYS) IN THE WRITE COMPONENT'S
!***  IMPORT STATE ARE NOT KNOWN A PRIORI.  IN ORDER TO TRANSFER
!***  THEM TO THE WRITE TASKS, EXTRACT THE NUMBER OF EACH OF
!***  THEM ALONG WITH THEIR NAMES.  THE SCALARS CAN BE LUMPED IN
!***  WITH THE 1D ARRAYS AT THIS POINT.
!
!***  EVEN THOUGH THESE COUNTS ARE JUST SCALAR INTEGERS THEIR
!***  POINTERS WERE ALLOCATED IN WRT_INIT TO LENGTH 1 SINCE THEY
!***  WILL BE USED IN ESMF_Send/Recv WHICH REQUIRE THEM TO BE 
!***  CONTIGUOUS DATA ARRAYS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ALL INTEGER QUANTITIES (AS 1D ARRAYS) AND 1D AND 2D REAL
!***  QUANTITIES WILL BE STRUNG TOGETHER IN SINGLE ARRAYS OF 
!***  EACH PARTICULAR TYPE.  ARRAYS THAT WILL HOLD THE LENGTH OF  
!***  EACH OF THE QUANTITIES IN THESE 'STRINGS' WERE ALLOCATED
!***  IN WRT_INIT.
!-----------------------------------------------------------------------
!
      KOUNT_I1D=0
      KOUNT_R1D=0
      KOUNT_LOG=0
!
      LENGTH_SUM_I1D=0
      LENGTH_SUM_R1D=0
      LENGTH_SUM_LOG=0
!
      NCHAR_I1D=0
      NCHAR_R1D=0
      NCHAR_LOG=0
!
!-----------------------------------------------------------------------
!
      fcst_tasks: IF(MYPE<=LAST_FCST_TASK)THEN                           !<-- Only forecast tasks will extract output information
                                                                         !    from the import state because only they participated
                                                                         !    in filling the import state in the Dynamics/Physics
                                                                         !    components.
!
!-----------------------------------------------------------------------
!***  FIRST FIND THE NUMBER OF Attributes IN THE HISTORY DATA Bundle
!***  IN THE IMPORT STATE AND THEN FIND THEIR NAMES, LENGTHS, AND
!***  DATATYPES.
!***  EXTRACT THE INTEGER AND REAL DATA AND PACK IT INTO INTEGER
!***  AND REAL BUFFERS.  LATER THE BUFFERS WILL BE SENT FROM THE
!***  FORECAST TASKS (THE ONLY ONES WHO CAN SEE THE ORIGINAL
!***  DATA) TO THE WRITE TASKS.
!
!***  THE FACT THAT THE Attribute HISTORY DATA IS BEING COLLECTED
!***  HERE IN A BLOCK THAT EXECUTES ONLY ONCE PER WRITE GROUP
!***  IMPLIES THE ASSUMPTION THAT ONLY THE 2D/3D DATA
!***  ASSOCIATED WITH THE FORECAST GRID CAN CHANGE WITH TIME.
!***  IF ANY SCALAR/1D Attribute DATA CHANGE WITH TIME THEN
!***  THIS MUST BE MOVED OUT OF THIS 'FIRST' BLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Attribute Count from History Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(bundle =HISTORY_BUNDLE        &  !<-- The write component's history data Bundle
                              ,count  =NUM_ATTRIB            &  !<-- # of Attributes in the history data Bundle
                              ,rc     =RC)       
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
        attribute_loop: DO N=1,NUM_ATTRIB                                  !<-- Loop through all the Attributes
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Get Attribute Names, Datatypes, Lengths"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(bundle         =HISTORY_BUNDLE  &  !<-- The write component's history data Bundle
                                ,attributeIndex =N               &  !<-- Index of each Attribute
                                ,name           =ATTRIB_NAME     &  !<-- Each Attribute's name
                                ,typekind       =DATATYPE        &  !<-- Each Attribute's ESMF Datatype
                                ,count          =LENGTH          &  !<-- Each Attribute's length
                                ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!                 -- SCALAR AND 1D INTEGER HISTORY DATA --
!-----------------------------------------------------------------------
!
          IF(DATATYPE==ESMF_TYPEKIND_I4)THEN                               !<-- Extract integer data with rank <2
!
            ALLOCATE(WORK_ARRAY_I1D(LENGTH),stat=RC)                       !<-- This length is from the preceding call 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Scalar/1-D Integer Data from History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(bundle    =HISTORY_BUNDLE      &  !<-- The write component's history data Bundle
                                  ,name      =ATTRIB_NAME         &  !<-- Name of the Attribute to extract
                                  ,count     =LENGTH              &  !<-- Length of Attribute
                                  ,valueList =WORK_ARRAY_I1D      &  !<-- Place the Attribute here
                                  ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_I1D=KOUNT_I1D+1                                          !<-- Count # of integer Attributes
!
            NPOSN_END=KOUNT_I1D*ESMF_MAXSTR
            NPOSN_START=NPOSN_END-ESMF_MAXSTR+1     
            wrt_int_state%NAMES_I1D_STRING(NPOSN_START:NPOSN_END)=ATTRIB_NAME  !<-- Save the 1D integer names
            NCHAR_I1D=NCHAR_I1D+ESMF_MAXSTR                                !<-- Save #of characters in all scalar/1D integer names
                                                                           !    Note that each name is being given
                                                                           !    EMSF_MAXSTR total spaces
!
            DO L=1,LENGTH
              wrt_int_state%ALL_DATA_I1D(LENGTH_SUM_I1D+L)=WORK_ARRAY_I1D(L)  !<-- String together the integer data
            ENDDO
!
            LENGTH_SUM_I1D=LENGTH_SUM_I1D+LENGTH                           !<-- Total word sum of scalar/1D integer data
            wrt_int_state%LENGTH_DATA_I1D(KOUNT_I1D)=LENGTH                !<-- Store length of each individual integer variable
!
            DEALLOCATE(WORK_ARRAY_I1D)
!
!-----------------------------------------------------------------------
!                  -- SCALAR AND 1D REAL HISTORY DATA --
!-----------------------------------------------------------------------
!
          ELSEIF(DATATYPE==ESMF_TYPEKIND_R4)THEN                           ! <-- Extract real data with rank <2
!
            ALLOCATE(WORK_ARRAY_R1D(LENGTH),stat=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Scalar/1-D Real Data from History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(bundle    =HISTORY_BUNDLE            &  !<-- The write component's history data Bundle
                                  ,name      =ATTRIB_NAME               &  !<-- Name of the Attribute to extract
                                  ,count     =LENGTH                    &  !<-- Length of Attribute
                                  ,valueList =WORK_ARRAY_R1D            &  !<-- Place the Attribute here
                                  ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_R1D=KOUNT_R1D+1                                          !<-- Count # of real Attributes
!
            NPOSN_END=KOUNT_R1D*ESMF_MAXSTR
            NPOSN_START=NPOSN_END-ESMF_MAXSTR+1     
            wrt_int_state%NAMES_R1D_STRING(NPOSN_START:NPOSN_END)=ATTRIB_NAME  !<-- Save the scalar/1D real names
            NCHAR_R1D=NCHAR_R1D+ESMF_MAXSTR                                !<-- Save #of characters in all scalar/1D real names
                                                                           !    Note that each name is being given
                                                                           !    EMSF_MAXSTR total spaces
!
            DO L=1,LENGTH
              wrt_int_state%ALL_DATA_R1D(LENGTH_SUM_R1D+L)=WORK_ARRAY_R1D(L)  !<-- String together the real data
            ENDDO
!
            LENGTH_SUM_R1D=LENGTH_SUM_R1D+LENGTH                           !<-- Total word sum of scalar/1D real data
            wrt_int_state%LENGTH_DATA_R1D(KOUNT_R1D)=LENGTH                !<-- Store length of each individual real variable
!
            DEALLOCATE(WORK_ARRAY_R1D)
!
!-----------------------------------------------------------------------
!                          -- LOGICAL DATA --                       
!-----------------------------------------------------------------------
!
!ratko    ELSEIF(DATATYPE==ESMF_DATA_LOGICAL)THEN                          ! <-- Extract logical data
! --- nothing else than I4 and R4; should FIX later
          ELSE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Logical Data from History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(bundle =HISTORY_BUNDLE               &  !<-- The write component's history data Bundle
                                  ,name   =ATTRIB_NAME                  &  !<-- Name of the Attribute to extract
                                  ,value  =WORK_LOGICAL                 &  !<-- Place the Attribute here
                                  ,rc     =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_LOG=KOUNT_LOG+1                                          !<-- Count # of logical Attributes
!
            NPOSN_END=KOUNT_LOG*ESMF_MAXSTR
            NPOSN_START=NPOSN_END-ESMF_MAXSTR+1     
            wrt_int_state%NAMES_LOG_STRING(NPOSN_START:NPOSN_END)=ATTRIB_NAME  !<-- Save the logical names
            NCHAR_LOG=NCHAR_LOG+ESMF_MAXSTR                                !<-- Save #of characters in all logical names
                                                                           !    Note that each name is being given
                                                                           !    EMSF_MAXSTR total spaces
!
            wrt_int_state%ALL_DATA_LOG(KOUNT_LOG)=WORK_LOGICAL             !<-- String together the logical data
!
            LENGTH_SUM_LOG=LENGTH_SUM_LOG+1                                !<-- Total length of all logical data variables
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO attribute_loop
!
!-----------------------------------------------------------------------
!***  INSERT NUMBER AND LENGTHS OF SCALAR/1D INTEGER AND REAL QUANTITIES
!***  AND LOGICALS INTO THE WRITE COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
        wrt_int_state%KOUNT_I1D(1)=KOUNT_I1D
        wrt_int_state%KOUNT_R1D(1)=KOUNT_R1D
        wrt_int_state%KOUNT_LOG(1)=KOUNT_LOG
!
        wrt_int_state%LENGTH_SUM_I1D(1)=LENGTH_SUM_I1D 
        wrt_int_state%LENGTH_SUM_R1D(1)=LENGTH_SUM_R1D
        wrt_int_state%LENGTH_SUM_LOG(1)=LENGTH_SUM_LOG
!
!-----------------------------------------------------------------------
!***  NOW EXTRACT THE NUMBER OF ESMF Fields IN THE HISTORY DATA Bundle
!***  WRITE COMPONENT'S IMPORT STATE ALONG WITH THEIR NAMES.
!***  SAVE THE Field INFORMATION OF THE 2D HISTORY DATA SINCE 
!***  IT WILL BE NEEDED FOR DATA EXTRACTION FROM THE IMPORT STATE
!***  AND THE TRANSFER TO THE WRITE TASKS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  FIND OUT THE NUMBER OF ESMF Fields IN THE HISTORY DATA Bundle.
!***  IT WAS Fields INTO WHICH THE 2D GRIDDED HISTORY DATA WAS PLACED.
!
!***  THIS INFORMATION WILL BE SAVED IN THE INTERNAL STATE
!***  FOR ALL OUTPUT TIMES AND NOT BE RETRIEVED OVER AND OVER AGAIN.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Field Count from History Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldBundleGet(bundle    =HISTORY_BUNDLE                   &  !<-- The write component's history data Bundle
                                ,fieldCount=wrt_int_state%NCOUNT_FIELDS(1)   &  !<-- Get total # of Fields in the history data Bundle
                                ,rc        =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  NOW EXTRACT THE NAMES OF ALL THE Fields IN THE BUNDLE.  
!***  ALSO, THE NUMBER OF Field NAMES RETURNED SHOULD EQUAL 
!***  THE FIELD COUNT IN THE PRECEDING CALL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Field Names from History Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldBundleGet(bundle    =HISTORY_BUNDLE           &  !<-- The write component's history data Bundle
                                ,nameList  =wrt_int_state%FIELD_NAME &  !<-- Array of ESMF Field names in the Bundle
                                ,nameCount =NUM_FIELD_NAMES          &  !<-- Number of Field names in the Bundle
                                ,rc        =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(NUM_FIELD_NAMES/=wrt_int_state%NCOUNT_FIELDS(1))THEN
          WRITE(0,*)' WARNING: Number of Fields in Bundle of history'    &
                   ,' output does not equal the number of Field names'
          WRITE(0,*)' They are ',NUM_FIELD_NAMES,' and '                 &
                   ,wrt_int_state%NCOUNT_FIELDS(1),', respectively'
        ENDIF
!
!-----------------------------------------------------------------------
!***  DO A PRELIMINARY EXTRACTION OF THE Fields THEMSELVES IN ORDER TO
!***  COUNT THE NUMBER OF REAL AND INTEGER 2D ARRAYS.  
!-----------------------------------------------------------------------
!
        KOUNT_R2D=0
        KOUNT_I2D=0
        NCHAR_I2D(1)=0
        NCHAR_R2D(1)=0
!
        DO N=1,wrt_int_state%NCOUNT_FIELDS(1)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Fields from History Bundle for Counting"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleGet(bundle=HISTORY_BUNDLE                &  !<-- The write component's history data Bundle
                                  ,name  =wrt_int_state%FIELD_NAME(N)   &  !<-- The ESMF Field's name
                                  ,field =FIELD_WORK1                   &  !<-- The ESMF Field taken from the Bundle
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Get Datatype of Fields for Counting Real/Integer"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field=FIELD_WORK1                          &  !<-- The ESMF 2D Field
                            ,typekind=DATATYPE                          &  !<-- The Field's ESMF Datatype
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(DATATYPE==ESMF_TYPEKIND_I4)THEN
            KOUNT_I2D=KOUNT_I2D+1                                          !<-- Add up the total number of integer 2D Fields
!
            IF(KOUNT_I2D>MAX_DATA_I2D)THEN
              WRITE(0,*)' FATAL: YOU HAVE EXCEEDED MAX NUMBER OF INTEGER 2D FIELDS FOR OUTPUT'
              WRITE(0,*)' YOU MUST INCREASE VALUE OF MAX_DATA_I2D WHICH NOW EQUALS ',MAX_DATA_I2D
            ENDIF
!
            NPOSN_END  =KOUNT_I2D*ESMF_MAXSTR
            NPOSN_START=NPOSN_END-ESMF_MAXSTR+1
            wrt_int_state%NAMES_I2D_STRING(NPOSN_START:NPOSN_END)=wrt_int_state%FIELD_NAME(N) !<-- Save the 2D integer Field names 
                                                                                              !<-- in one long string
            NCHAR_I2D(1)=NCHAR_I2D(1)+ESMF_MAXSTR
!
          ELSEIF(DATATYPE==ESMF_TYPEKIND_R4)THEN
            KOUNT_R2D=KOUNT_R2D+1                                          !<-- Add up the total number of real 2D Fields
!
            IF(KOUNT_R2D>MAX_DATA_R2D)THEN
              WRITE(0,*)' FATAL: YOU HAVE EXCEEDED MAX NUMBER OF REAL 2D FIELDS FOR OUTPUT'
              WRITE(0,*)' KOUNT_R2D=',KOUNT_R2D
              WRITE(0,*)' YOU MUST INCREASE VALUE OF MAX_DATA_R2D WHICH NOW EQUALS ',MAX_DATA_R2D
            ENDIF
!
            NPOSN_END  =KOUNT_R2D*ESMF_MAXSTR
            NPOSN_START=NPOSN_END-ESMF_MAXSTR+1
            wrt_int_state%NAMES_R2D_STRING(NPOSN_START:NPOSN_END)=wrt_int_state%FIELD_NAME(N) !<-- Save the 2D real Field names 
                                                                                              !<-- in one long string
            NCHAR_R2D(1)=NCHAR_R2D(1)+ESMF_MAXSTR
!
          ENDIF
!
        ENDDO
!
        wrt_int_state%KOUNT_R2D(1)=KOUNT_R2D
        wrt_int_state%KOUNT_I2D(1)=KOUNT_I2D
!
!-----------------------------------------------------------------------
!
      ENDIF fcst_tasks
!
!-----------------------------------------------------------------------
!***  IF THERE ARE NO QUANTITIES SPECIFIED FOR HISTORY OUTPUT,
!***  FORECAST TASK 0 WILL INFORM THE WRITE TASKS AND THEN
!***  EVERYONE WILL RETURN.
!-----------------------------------------------------------------------
!
      NO_FIELDS(1)=ESMF_FALSE
!
      IF(MYPE==0)THEN
        IF(wrt_int_state%NCOUNT_FIELDS(1)==0)NO_FIELDS(1)=ESMF_TRUE          !<-- Reset flag saying there are no history quantities
        LAST_WRITE_TASK=NTASKS-1                                             !<-- The last write task in this group
!
        DO N=LEAD_WRITE_TASK,LAST_WRITE_TASK                                 !<-- Loop through all the write tasks in the write group
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Fcst Task0 Informs All That There Are No Fields"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                          ,sendData=NO_FIELDS                             &  !<-- Send this data
                          ,count   =1                                     &  !<-- Words sent
                          ,dst     =N                                     &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDDO
      ENDIF
!
      IF(MYPE>=LEAD_WRITE_TASK)THEN                                        !<-- All write tasks in this group receive
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Tasks Told By Fcst Task0 There Are No Fields"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=NO_FIELDS                             &  !<-- Recv this data
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      IF(NO_FIELDS(1)==ESMF_TRUE)THEN
        IF(MYPE==0)THEN
          WRITE(6,*)'WARNING: No Import ESMF quantities for the Write Component'
          WRITE(0,*)'WARNING: No Import ESMF quantities for the Write Component'
        ENDIF
!
        RETURN                                                             !<-- All tasks return if there is no history output
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  FORECAST TASK 0 SENDS ALL THE WRITE TASKS THE NUMBER OF
!***  REAL AND INTEGER 2D GRIDDED QUANTITIES PLUS ALL OF THE
!***  LOCAL HORIZONTAL DOMAIN LIMITS IN PREPARATION FOR THE
!***  WRITE TASKS' RECEIVING AND ASSEMBLING THE LOCAL HISTORY
!***  DATA THEY RECEIVE FROM THE FORECAST TASKS.
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN                                                        !<-- Forecast task 0 sends
!
        LAST_WRITE_TASK=NTASKS-1                                             !<-- The last write task in this group
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Fcst Task0 Sends Write Tasks Info for Quilting"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        DO N=LEAD_WRITE_TASK,LAST_WRITE_TASK                                 !<-- Loop through all the write tasks in the write group
!
          CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                          ,sendData=INPES                                 &  !<-- Send this data
                          ,count   =1                                     &  !<-- Words sent
                          ,dst     =N                                     &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                          ,sendData=JNPES                                 &  !<-- Send this data
                          ,count   =1                                     &  !<-- Words sent
                          ,dst     =N                                     &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                          ,sendData=IHALO                                 &  !<-- Send this data
                          ,count   =1                                     &  !<-- Words sent
                          ,dst     =N                                     &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                          ,sendData=JHALO                                 &  !<-- Send this data
                          ,count   =1                                     &  !<-- Words sent
                          ,dst     =N                                     &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                          ,sendData=wrt_int_state%NCOUNT_FIELDS           &  !<-- Send this data
                          ,count   =1                                     &  !<-- Words sent
                          ,dst     =N                                     &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                          ,sendData=wrt_int_state%KOUNT_R2D               &  !<-- Send this data
                          ,count   =1                                     &  !<-- Words sent
                          ,dst     =N                                     &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                          ,sendData=wrt_int_state%KOUNT_I2D               &  !<-- Send this data
                          ,count   =1                                     &  !<-- Words sent
                          ,dst     =N                                     &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                          ,sendData=wrt_int_state%LOCAL_ISTART            &  !<-- Send this data
                          ,count   =NUM_PES_FCST                          &  !<-- Words sent
                          ,dst     =N                                     &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                          ,sendData=wrt_int_state%LOCAL_IEND              &  !<-- Send this data
                          ,count   =NUM_PES_FCST                          &  !<-- Words sent
                          ,dst     =N                                     &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                          ,sendData=wrt_int_state%LOCAL_JSTART            &  !<-- Send this data
                          ,count   =NUM_PES_FCST                          &  !<-- Words sent
                          ,dst     =N                                     &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                          ,sendData=wrt_int_state%LOCAL_JEND              &  !<-- Send this data
                          ,count   =NUM_PES_FCST                          &  !<-- Words sent
                          ,dst     =N                                     &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
        ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      IF(MYPE>=LEAD_WRITE_TASK)THEN                                        !<-- All write tasks in this group receive
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Tasks Recv Quilting Info From Fcst Task0"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=INPES                                 &  !<-- Recv this data
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=JNPES                                 &  !<-- Recv this data
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        NUM_PES_FCST=INPES(1)*JNPES(1)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=IHALO                                 &  !<-- Recv this data
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=JHALO                                 &  !<-- Recv this data
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        wrt_int_state%IHALO=IHALO(1)
        wrt_int_state%JHALO=JHALO(1)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%NCOUNT_FIELDS           &  !<-- Recv this data
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%KOUNT_R2D               &  !<-- Recv this data
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%KOUNT_I2D               &  !<-- Recv this data
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%LOCAL_ISTART            &  !<-- Recv this data
                        ,count   =NUM_PES_FCST                          &  !<-- Words received (#of fcst tasks)
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%LOCAL_IEND              &  !<-- Recv this data
                        ,count   =NUM_PES_FCST                          &  !<-- Words received (#of fcst tasks)
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%LOCAL_JSTART            &  !<-- Recv this data
                        ,count   =NUM_PES_FCST                          &  !<-- Words received (#of fcst tasks)
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%LOCAL_JEND              &  !<-- Recv this data
                        ,count   =NUM_PES_FCST                          &  !<-- Words received (#of fcst tasks)
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  FORECAST TASK 0 SENDS THE 2D DATA NAMES TO THE LEAD WRITE TASK.
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN                                                      !<-- Fcst task0 alone can send write task preliminary info
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Fcst Task0 Sends Write Tasks 2D Data Names"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=NCHAR_I2D                             &  !<-- Send total length of the names of 2D integer data
                        ,count   =1                                     &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                       &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%NAMES_I2D_STRING        &  !<-- Send names of 2D integer history variables
                        ,count   =NCHAR_I2D(1)                          &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                       &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=NCHAR_R2D                             &  !<-- Send total length of the names of 2D real data
                        ,count   =1                                     &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                       &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%NAMES_R2D_STRING        &  !<-- Send names of 2D real history variables
                        ,count   =NCHAR_R2D(1)                          &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                       &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
      ELSEIF(MYPE==NTASKS-NWTPG)THEN                                       !<-- 1st write task receives 2D preliminary info
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=NCHAR_I2D                             &  !<-- Recv total length of the names of 2D integer data
                        ,count   =1                                     &  !<-- Words sent
                        ,src     =0                                     &  !<-- Sending task (fcst task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%NAMES_I2D_STRING        &  !<-- Recv names of 2D integer history variables
                        ,count   =NCHAR_I2D(1)                          &  !<-- Words sent
                        ,src     =0                                     &  !<-- Sending task (fcst task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=NCHAR_R2D                             &  !<-- Recv total length of the names of 2D gridded data
                        ,count   =1                                     &  !<-- Words sent
                        ,src     =0                                     &  !<-- Sending task (fcst task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%NAMES_R2D_STRING        &  !<-- Recv names of 2D real history variables
                        ,count   =NCHAR_R2D(1)                          &  !<-- Words sent
                        ,src     =0                                     &  !<-- Sending task (fcst task 0)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!
!-----------------------------------------------------------------------
!***  EACH WRITE TASK MUST KNOW THE IDs OF THE FORECAST TASKS
!***  FROM WHICH IT WILL RECEIVE 2D GRIDDED HISTORY DATA.
!-----------------------------------------------------------------------
!
      IF(MYPE>=LEAD_WRITE_TASK)THEN                                        !<-- The write tasks
!
        ALLOCATE(wrt_int_state%ID_FTASK_RECV_STA(LEAD_WRITE_TASK:LAST_WRITE_TASK))
        ALLOCATE(wrt_int_state%ID_FTASK_RECV_END(LEAD_WRITE_TASK:LAST_WRITE_TASK))
!
        NN=0
        DO N=LEAD_WRITE_TASK,LAST_WRITE_TASK
          NN=NN+1
          CALL PARA_RANGE(JNPES(1),NWTPG,NN                             &  !<-- Find each write task's first and last rows of
                         ,JROW_FIRST,JROW_LAST)                            !<--   fcst tasks from which it will recv
!
          wrt_int_state%ID_FTASK_RECV_STA(N)=(JROW_FIRST-1)*INPES(1)       !<-- First fcst task that sends to this write task
          wrt_int_state%ID_FTASK_RECV_END(N)=JROW_LAST*INPES(1)-1          !<-- Last fcst task that sends to this write task
        ENDDO
!
!-----------------------------------------------------------------------
!***  EACH WRITE TASK ALSO MUST KNOW THE NORTH-SOUTH EXTENT OF THE
!***  FULL 2D DOMAIN THAT IT WILL HANDLE.  THIS IS DETERMINED BY
!***  THE COVERAGE OF THE FCST TASKS THAT SEND TO IT.
!-----------------------------------------------------------------------
!
        JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- JTS of 1st fcst task that sends to this write task
        JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- JTE of last fcst task that sends to this write task
!
!-----------------------------------------------------------------------
!***  NOW EACH WRITE TASK ALLOCATES ITS OWN SECTION OF THE 2D DOMAIN
!***  FOR ALL THE 2D VARIABLES IT WILL RECEIVE AND ITS 1D EQUIVALENT
!***  USED TO TRANSFER THE DATA TO THE LEAD WRITE TASK.
!-----------------------------------------------------------------------
!
        ALLOCATE(wrt_int_state%WRITE_SUBSET_I(1:IM,JSTA_WRITE:JEND_WRITE  &
                                             ,wrt_int_state%KOUNT_I2D(1)))
        LENGTH=IM*(JEND_WRITE-JSTA_WRITE+1)*wrt_int_state%KOUNT_I2D(1)
        ALLOCATE(wrt_int_state%BUFF_INT(LENGTH))
!
        ALLOCATE(wrt_int_state%WRITE_SUBSET_R(1:IM,JSTA_WRITE:JEND_WRITE  &
                                             ,wrt_int_state%KOUNT_R2D(1)))
        LENGTH=IM*(JEND_WRITE-JSTA_WRITE+1)*wrt_int_state%KOUNT_R2D(1)
        ALLOCATE(wrt_int_state%BUFF_REAL(LENGTH))
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  THE LEAD WRITE TASK ALLOCATES ITS WORKING ARRAYS INTO WHICH
!***  IT WILL ASSEMBLE EACH INDIVIDUAL 2D FIELD THAT WILL BE
!***  WRITTEN TO THE HISTORY FILES.
!-----------------------------------------------------------------------
!
      IF(MYPE==LEAD_WRITE_TASK)THEN
        ALLOCATE(wrt_int_state%OUTPUT_ARRAY_I2D(1:IM,1:JM,wrt_int_state%KOUNT_I2D(1)))
        ALLOCATE(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM,wrt_int_state%KOUNT_R2D(1)))
      ENDIF
!
!-----------------------------------------------------------------------
!***  SINCE ALL SCALAR/1D DATA IS IDENTICAL ON ALL FORECAST TASKS,
!***  TASK 0 ALONE CAN SEND THE INFORMATION TO THE LEAD WRITE TASK
!***  THAT WILL LATER WRITE IT TO THE HISTORY FILE.
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      task_0_sends: IF(MYPE==0)THEN                                      !<-- Forecast task 0 sends
!--------------------------------------------------------------------
!
!------------------------------------------------
!***  SEND SCALAR/1D INTEGER HISTORY INFORMATION.
!------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Fcst Task0 Sends Scalar/1D Integer History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%KOUNT_I1D             &  !<-- Send # of scalar/1D integer history variables
                        ,count   =1                                   &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%LENGTH_SUM_I1D        &  !<-- Send length of string of all such integer history variables
                        ,count   =1                                   &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%LENGTH_DATA_I1D       &  !<-- Send lengths of each scalar/1D integer history variable
                        ,count   =wrt_int_state%KOUNT_I1D(1)          &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%NAMES_I1D_STRING      &  !<-- Send names of each scalar/1D integer history variable
                        ,count   =NCHAR_I1D                           &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%ALL_DATA_I1D          &  !<-- Send the full string of all scalar/1D integer history data
                        ,count   =wrt_int_state%LENGTH_SUM_I1D(1)     &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------------------------------
!***  SEND SCALAR/1D REAL HISTORY INFORMATION.
!---------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Fcst Task0 Sends Scalar/1D Real History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%KOUNT_R1D             &  !<-- Send # of scalar/1D real history variables
                        ,count   =1                                   &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%LENGTH_SUM_R1D        &  !<-- Send length of string of all such real history variables
                        ,count   =1                                   &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%LENGTH_DATA_R1D       &  !<-- Send lengths of each scalar/1D real history variable
                        ,count   =wrt_int_state%KOUNT_R1D(1)          &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%NAMES_R1D_STRING      &  !<-- Send names of each scalar/1D real history variable
                        ,count   =NCHAR_R1D                           &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%ALL_DATA_R1D          &  !<-- Send the full string of all scalar/1D real history data
                        ,count   =wrt_int_state%LENGTH_SUM_R1D(1)     &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------
!***  SEND LOGICAL HISTORY INFORMATION.
!--------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Fcst Task0 Sends Logical History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%KOUNT_LOG             &  !<-- Send # of logical history variables
                        ,count   =1                                   &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%LENGTH_SUM_LOG        &  !<-- Send length of string of all logical variables
                        ,count   =1                                   &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%NAMES_LOG_STRING      &  !<-- Send names of each logical history variable
                        ,count   =NCHAR_LOG                           &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%ALL_DATA_LOG          &  !<-- Send the full string of all logical history data
                        ,count   =wrt_int_state%LENGTH_SUM_LOG(1)     &  !<-- Words sent
                        ,dst     =LEAD_WRITE_TASK                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
      ENDIF task_0_sends
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      write_task_recvs: IF(MYPE==LEAD_WRITE_TASK)THEN                      !<-- 1st write task in this group receives
                                                                           !    all of the data just sent to it by
                                                                           !    fcst task 0
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RECEIVE SCALAR/1D INTEGER HISTORY INFORMATION
!***  FROM FORECAST TASK 0.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Tasks Recv Scalar/1D Integer History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%KOUNT_I1D               &  !<-- Recv # of integer history variables
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%LENGTH_SUM_I1D          &  !<-- Recv length of string of all integer history variables
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%LENGTH_DATA_I1D         &  !<-- Recv lengths of each integer history variable
                        ,count   =wrt_int_state%KOUNT_I1D(1)            &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        NCHAR_I1D=wrt_int_state%KOUNT_I1D(1)*ESMF_MAXSTR
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%NAMES_I1D_STRING        &  !<-- Recv names of integer history variables
                        ,count   =NCHAR_I1D                             &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%ALL_DATA_I1D            &  !<-- Recv the string of integer history data
                        ,count   =wrt_int_state%LENGTH_SUM_I1D(1)       &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RECEIVE SCALAR/1D REAL HISTORY INFORMATION.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Tasks Recv Scalar/1D Real History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%KOUNT_R1D               &  !<-- Recv # of scalar/1D real history variables
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%LENGTH_SUM_R1D          &  !<-- Recv length of string of all such real history variables
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%LENGTH_DATA_R1D         &  !<-- Recv lengths of each scalar/1D real history variable
                        ,count   =wrt_int_state%KOUNT_R1D(1)            &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        NCHAR_R1D=wrt_int_state%KOUNT_R1D(1)*ESMF_MAXSTR
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%NAMES_R1D_STRING        &  !<-- Recv names of scalar/1D real history variables
                        ,count   =NCHAR_R1D                             &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%ALL_DATA_R1D            &  !<-- Recv the string of all scalar/1D real history data
                        ,count   =wrt_int_state%LENGTH_SUM_R1D(1)       &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RECEIVE LOGICAL HISTORY INFORMATION.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Tasks Recv Logical Real History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%KOUNT_LOG               &  !<-- Recv # of logical history variables
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%LENGTH_SUM_LOG          &  !<-- Recv length of string of all logical history variables
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        NCHAR_LOG=wrt_int_state%KOUNT_LOG(1)*ESMF_MAXSTR
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%NAMES_LOG_STRING        &  !<-- Recv names of logical history variables
                        ,count   =NCHAR_LOG                             &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%ALL_DATA_LOG            &  !<-- Recv the string of all logical history data
                        ,count   =wrt_int_state%LENGTH_SUM_LOG(1)       &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDIF write_task_recvs
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(INPES,JNPES)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE FIRST_PASS
!
!-----------------------------------------------------------------------
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-----------------------------------------------------------------------
!
      SUBROUTINE OPEN_DISK_FILE(WRT_INT_STATE)
!
!-----------------------------------------------------------------------
!***  OPEN A DISK FILE.
!-----------------------------------------------------------------------
!
      TYPE(WRITE_INTERNAL_STATE),INTENT(INOUT) :: WRT_INT_STATE          !<-- The I/O component's internal state
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: IO_UNIT,N,RC
!
      LOGICAL :: OPENED
!
      CHARACTER(ESMF_MAXSTR) :: FILENAME
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  SPECIFYING 'DEFERRED' AS THE FILENAME IN THE CONFIGURE FILE
!***  MEANS THAT WE WANT TO CONSTRUCT THE OUTPUT FILENAME FROM
!***  FILENAME_BASE (FROM THE CONFIG FILE) APPENDED WITH THE
!***  FORECAST HOUR.
!-----------------------------------------------------------------------
!
      write(0,*)' OPEN_DISK_FILE wrt_int_state%IO_FILE=',trim(wrt_int_state%IO_FILE)
!
      IF(wrt_int_state%IO_FILE=='DEFERRED')THEN
        WRITE(FILENAME,100)wrt_int_state%FILENAME_BASE,wrt_int_state%NFHOUR
  100   FORMAT(A14,I3.3)
!!!     WRITE(0,*)' Created filename=',filename,' FILENAME_BASE=',wrt_int_state%FILENAME_BASE
      ELSE
        FILENAME=wrt_int_state%IO_FILE
      ENDIF
!
!-----------------------------------------------------------------------
!***  FIND AN UNOPENED UNIT NUMBER IF ONE WAS NOT DESIGNATED IN
!***  THE CONFIGURE FILE.
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%IO_UNIT==-999)THEN
        DO N=51,99
          INQUIRE(N,opened=OPENED)
          IF(.NOT.OPENED)THEN
            IO_UNIT=N
            EXIT
          ENDIF
        ENDDO
        wrt_int_state%IO_UNIT=IO_UNIT
!
      ELSE
        IO_UNIT=wrt_int_state%IO_UNIT
      ENDIF
!
!-----------------------------------------------------------------------
!***  OPEN THE FILE NOW.
!-----------------------------------------------------------------------
!
      OPEN(unit  =IO_UNIT                                               &
          ,file  =FILENAME                                              &
          ,status=wrt_int_state%IO_STATUS                               &
          ,access=wrt_int_state%IO_ACCESS                               &
          ,form  =wrt_int_state%IO_FORM                                 &
          ,iostat=RC)
!
      IF(RC==0)THEN
        WRITE(0,*)' Opened IO_UNIT=',IO_UNIT,' for history'
        write(0,*)' iostat=',rc,' file=',trim(filename)
        write(0,*)' status=',trim(wrt_int_state%IO_STATUS), &
                  ' access=',trim(wrt_int_state%IO_ACCESS), &
                  ' form='  ,trim(wrt_int_state%IO_FORM)
      ELSE
        WRITE(0,*)' Failed to OPEN IO_UNIT=',IO_UNIT,' for history'
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE OPEN_DISK_FILE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_INIT(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM)
! 
!-----------------------------------------------------------------------
!***  EXECUTE THE INITIALIZE STEP OF THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ATM_INTERNAL_STATE),INTENT(INOUT) :: ATM_INT_STATE             !<-- The ATM Internal State
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_ATM                 !<-- The ATM Component's ESMF Clock
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config)      :: CF                                        !<-- The config object
      TYPE(ESMF_VM)          :: VM                                        !<-- The ESMF virtual machine.
!
      INTEGER :: I,INPES,J,JNPES,NUM_PES_FCST,RC,RC_INIT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE CONFIG OBJECT CF FROM THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Config Object from ATM Component in Write Init"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  HOW MANY FORECAST TASKS DO WE HAVE?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Write Init: Get INPES/JNPES from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ESMF configure object
                                  ,value =INPES                         &  !<-- # of fcst tasks in I direction
                                  ,label ='inpes:'                      &  !<-- Give the value of this label to INPES
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ESMF configure object
                                  ,value =JNPES                         &  !<-- # of fcst tasks in J direction
                                  ,label ='jnpes:'                      &  !<-- Give the value of this label to JNPES
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NUM_PES_FCST=INPES*JNPES                                             !<-- Total number of forecast tasks
!
!-----------------------------------------------------------------------
!***  EXECUTE THE INITIALIZE STEP FOR THE WRITE COMPONENTS.
!***  THESE ARE THE INITIALIZE SUBROUTINES SPECIFIED IN THE
!***  REGISTER ROUTINES CALLED IN ESMF_GridCompSetServices.
!-----------------------------------------------------------------------
!
      DO J=1,atm_int_state%WRITE_GROUPS
!
!!!!    CALL ESMF_VMBarrier(VM,rc=RC)    ! Insert barrier since fcst tasks are involved in each iteration of write groups
!
        DO I=1,NUM_PES_FCST+atm_int_state%WRITE_TASKS_PER_GROUP
          IF(MYPE==atm_int_state%PETLIST_WRITE(I,J))THEN                   !<--  Forecast tasks plus the Write tasks in each write group
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Execute Initialize Step of Write Component"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_GridCompInitialize(atm_int_state%WRT_COMPS(J)                 &  !<-- The Write gridded components
                                        ,importstate=atm_int_state%IMP_STATE_WRITE  &  !<-- The Write import state
                                        ,exportstate=atm_int_state%EXP_STATE_WRITE  &  !<-- The Write export state
                                        ,clock      =CLOCK_ATM                      &  !<-- The ESMF clock of the ATM component
                                        ,phase      =ESMF_SINGLEPHASE               &
                                        ,rc         =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDIF
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  SET THE FIRST WRITE GROUP AS THE FIRST ONE TO ACT.
!-----------------------------------------------------------------------
!
      atm_int_state%WRITE_GROUP_READY_TO_GO=1
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_INIT
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_ASYNC(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM,MYPE)
!
!-----------------------------------------------------------------------
!***  WRITE OUT A HISTORY FILE USING THE ASYNCHRONOUS QUILTING.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp)     ,INTENT(INOUT) :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ATM_INTERNAL_STATE),INTENT(INOUT) :: ATM_INT_STATE             !<-- The ATM Internal State
      TYPE(ESMF_Clock)        ,INTENT(INOUT) :: CLOCK_ATM                 !<-- The ATM Component's ESMF Clock
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config) :: CF                                             !<-- The configure object (~namelist)
      TYPE(ESMF_Time)   :: CURRTIME                                       !<-- The current forecast time (ESMF)
      INTEGER,INTENT(IN) :: MYPE
!
      INTEGER :: YY,MM,DD,H,M,S                                           !<-- Year, Month, Day, Hour, Minute, Second (integer)
!
      INTEGER :: I,INPES,JNPES,N_GROUP,NUM_PES_FCST                     &
                ,WRITE_GROUPS,WRITE_TASKS_PER_GROUP                     &
                ,RC,RC_ASYNC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  EXTRACT THE CONFIGURE OBJECT IN ORDER TO KNOW THE NUMBER
!***  OF FORECAST TASKS, WRITE GROUPS, AND WRITE TASKS PER GROUP.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Extract Config Object from ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Get General Info from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =INPES                         &  !<-- # of fcst tasks in I direction
                                  ,label ='inpes:'                      &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =JNPES                         &  !<-- # of fcst tasks in J direction
                                  ,label ='jnpes:'                      &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =WRITE_GROUPS                  &  !<-- Number of write groups
                                  ,label ='write_groups:'               &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =WRITE_TASKS_PER_GROUP         &  !<-- Number of write tasks per group
                                  ,label ='write_tasks_per_group:'      &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NUM_PES_FCST=INPES*JNPES                                             !<-- Number of forecast tasks
!
!-----------------------------------------------------------------------
!***  WHAT IS THE CURRENT FORECAST TIME?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Get Current Time from ATM Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock   =CLOCK_ATM                             &  !<-- The ATM component's ESMF Clock
                        ,currTime=CURRTIME                              &  !<-- The current forecast time (ESMF)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Convert ESMF Time to Real Time"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeGet (time=CURRTIME                                  &  !<-- The current forecast time (ESMF)
                        ,yy  =YY                                        &  !<-- The current year (integer)
                        ,mm  =MM                                        &  !<-- The current month (integer)
                        ,dd  =DD                                        &  !<-- The current day (integer)
                        ,h   =H                                         &  !<-- The current hour (integer)
                        ,m   =M                                         &  !<-- The current minute (integer)
                        ,s   =S                                         &  !<-- The current second (integer)
                        ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE EXPORT STATE OF THE DYNAMICS COMPONENT LIES WITHIN THE
!***  INTERNAL STATE OF THE ATM GRIDDED COMPONENT AND HOLDS THE
!***  IMPORT STATE OF THE WRITE COMPONENT.
!***  EXTRACT THAT WRITE COMPONENT'S IMPORT STATE SINCE WE ARE 
!***  ABOUT TO EXECUTE THE RUN STEP OF THE WRITE COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Extract Write Import State from Dyn Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state          =atm_int_state%EXP_STATE_DYN    &  !<-- The Dyn component's export state
                        ,itemName       ="Write Import State"           &  !<-- Name of state to be extracted
                        ,nestedState    =atm_int_state%IMP_STATE_WRITE  &  !<-- The extracted state
                        ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  ALL FORECAST TASKS PLUS THOSE WRITE TASKS IN THE APPROPRIATE
!***  WRITE GROUP EXECUTE THE RUN STEP OF A WRITE COMPONENT.
!-----------------------------------------------------------------------
!
      N_GROUP=atm_int_state%WRITE_GROUP_READY_TO_GO                          !<-- The active write group
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Execute Run Step of Write Components" 
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I=1,NUM_PES_FCST+WRITE_TASKS_PER_GROUP

        IF(MYPE==atm_int_state%PETLIST_WRITE(I,N_GROUP))THEN
!
          CALL ESMF_GridCompRun(atm_int_state%WRT_COMPS(N_GROUP)          &  !<-- The write gridded component
                               ,importState=atm_int_state%IMP_STATE_WRITE &  !<-- Its import state
                               ,exportState=atm_int_state%EXP_STATE_WRITE &  !<-- Its export state
                               ,clock      =CLOCK_ATM                     &  !<-- The ATM Clock
                               ,phase      =ESMF_SINGLEPHASE              &
                               ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(I==NUM_PES_FCST+1)THEN                                          !<-- The first write task tells us the history output time
            WRITE(0,101)YY,MM,DD,H,M,S
  101       FORMAT(' Wrote History File at ',I4.4,'_',I2.2,'_',I2.2,'_',I2.2,':',I2.2,':',I2.2)
          ENDIF
!
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  PREPARE TO USE THE NEXT WRITE GROUP AT THE NEXT OUTPUT TIME.
!***  RETURN TO THE 1ST GROUP IF WE HAVE CYCLED THROUGH ALL OF THEM.
!-----------------------------------------------------------------------
!
      IF(atm_int_state%WRITE_GROUP_READY_TO_GO==WRITE_GROUPS)THEN
        atm_int_state%WRITE_GROUP_READY_TO_GO=1
      ELSE
        atm_int_state%WRITE_GROUP_READY_TO_GO=atm_int_state%WRITE_GROUP_READY_TO_GO+1
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_ASYNC
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_RUNHISTORY(WRT_INT_STATE                         &
                                 ,IYEAR_FCST                            &
                                 ,IMONTH_FCST                           &
                                 ,IDAY_FCST                             &
                                 ,IHOUR_FCST                            &
                                 ,IMINUTE_FCST                          &
                                 ,SECOND_FCST                           &
                                 ,NF_HOURS                              &
                                 ,NF_MINUTES                            &
                                 ,NF_SECONDS                            &
                                 ,LEAD_WRITE_TASK)
!
!-----------------------------------------------------------------------
!***  WRITE OUT A BINARY RUN HISTORY FILE.
!-----------------------------------------------------------------------
!
      TYPE(WRITE_INTERNAL_STATE),INTENT(INOUT) :: WRT_INT_STATE            !<-- The Write component's internal state
!
      INTEGER,INTENT(IN) :: IYEAR_FCST                                  &
                           ,IMONTH_FCST                                 &
                           ,IDAY_FCST                                   &
                           ,IHOUR_FCST                                  &
                           ,IMINUTE_FCST                                &
                           ,NF_HOURS                                    &
                           ,NF_MINUTES                                  &
                           ,LEAD_WRITE_TASK
!
      REAL,INTENT(IN) :: NF_SECONDS,SECOND_FCST
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: N,N1,N2,NPOSN_1,NPOSN_2,LENGTH
      INTEGER :: NFIELD,RC
      CHARACTER(ESMF_MAXSTR)                :: NAME
      INTEGER,DIMENSION(:),POINTER          :: WORK_ARRAY_I1D
      REAL(4),DIMENSION(:),POINTER          :: WORK_ARRAY_R1D
      LOGICAL                               :: WRITE_LOGICAL
!
      TYPE(ESMF_Logical)                    :: WORK_LOGICAL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  OPEN THE HISTORY FILE AND WRITE THE CURRENT FORECAST TIME
!***  AND ELAPSED TIME.
!-----------------------------------------------------------------------
!
      CALL OPEN_DISK_FILE(WRT_INT_STATE)
      WRITE(0,*)' Opened unit=',wrt_int_state%IO_UNIT,' for history output'
!
      WRITE(wrt_int_state%IO_UNIT,iostat=RC)IYEAR_FCST
      WRITE(wrt_int_state%IO_UNIT,iostat=RC)IMONTH_FCST
      WRITE(wrt_int_state%IO_UNIT,iostat=RC)IDAY_FCST
      WRITE(wrt_int_state%IO_UNIT,iostat=RC)IHOUR_FCST
      WRITE(wrt_int_state%IO_UNIT,iostat=RC)IMINUTE_FCST
      WRITE(wrt_int_state%IO_UNIT,iostat=RC)SECOND_FCST
      WRITE(wrt_int_state%IO_UNIT,iostat=RC)NF_HOURS
      WRITE(wrt_int_state%IO_UNIT,iostat=RC)NF_MINUTES
      WRITE(wrt_int_state%IO_UNIT,iostat=RC)NF_SECONDS
!
      IF(wrt_int_state%NFHOUR==0)THEN
        WRITE(0,*)' Wrote IYEAR_FCST to history file unit ',wrt_int_state%IO_UNIT
        WRITE(0,*)' Wrote IMONTH_FCST to history file unit ',wrt_int_state%IO_UNIT
        WRITE(0,*)' Wrote IDAY_FCST to history file unit ',wrt_int_state%IO_UNIT
        WRITE(0,*)' Wrote IHOUR_FCST to history file unit ',wrt_int_state%IO_UNIT
        WRITE(0,*)' Wrote IMINUTE_FCST to history file unit ',wrt_int_state%IO_UNIT
        WRITE(0,*)' Wrote SECOND_FCST to history file unit ',wrt_int_state%IO_UNIT
        WRITE(0,*)' Wrote NF_HOURS to history file unit ',wrt_int_state%IO_UNIT
        WRITE(0,*)' Wrote NF_MINUTES to history file unit ',wrt_int_state%IO_UNIT
        WRITE(0,*)' Wrote NF_SECONDS to history file unit ',wrt_int_state%IO_UNIT
      ENDIF
!
!-----------------------------------------------------------------------
!***  INTEGER SCALAR/1D HISTORY VARIABLES
!-----------------------------------------------------------------------
!
      N2=0                                                                !<-- Word counter for full string of integer scalar/1D data
!
      DO N=1,wrt_int_state%KOUNT_I1D(1)                                   !<-- Loop through all scalar/1D integer data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_I1D_STRING(NPOSN_1:NPOSN_2)              !<-- The variable's name
        LENGTH=wrt_int_state%LENGTH_DATA_I1D(N)                           !<-- The variable's length in words
        ALLOCATE(WORK_ARRAY_I1D(LENGTH),stat=RC)
!
        DO N1=1,LENGTH
          N2=N2+1
          WORK_ARRAY_I1D(N1)=wrt_int_state%ALL_DATA_I1D(N2)               !<-- Extract the individual data from the data string
        ENDDO
!
        WRITE(wrt_int_state%IO_UNIT,iostat=RC)WORK_ARRAY_I1D              !<-- Write out the data
!
        IF(wrt_int_state%NFHOUR==0)THEN
          WRITE(0,*)'Wrote ',TRIM(NAME),' to history file unit ',wrt_int_state%IO_UNIT
        ENDIF
!
        DEALLOCATE(WORK_ARRAY_I1D)
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  REAL SCALAR/1D HISTORY VARIABLES
!-----------------------------------------------------------------------
!
      N2=0                                                                !<-- Word counter for full string of real scalar/1D data
!
      DO N=1,wrt_int_state%KOUNT_R1D(1)                                   !<-- Loop through all scalar/1D real data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_R1D_STRING(NPOSN_1:NPOSN_2)              !<-- The variable's name
        LENGTH=wrt_int_state%LENGTH_DATA_R1D(N)                           !<-- The variable's length
        ALLOCATE(WORK_ARRAY_R1D(LENGTH),stat=RC)
!
        DO N1=1,LENGTH
          N2=N2+1
          WORK_ARRAY_R1D(N1)=wrt_int_state%ALL_DATA_R1D(N2)               !<-- Extract the individual data from the data string
        ENDDO
!
        WRITE(wrt_int_state%IO_UNIT,iostat=RC)WORK_ARRAY_R1D              !<-- Write out the data
!
        IF(wrt_int_state%NFHOUR==0)THEN
          WRITE(0,*)'Wrote ',TRIM(NAME),' to history file unit ',wrt_int_state%IO_UNIT
        ENDIF
!
        DEALLOCATE(WORK_ARRAY_R1D)
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  LOGICAL HISTORY VARIABLES
!-----------------------------------------------------------------------
!
      N2=0                                                                !<-- Counter for full string of logical data
!
      DO N=1,wrt_int_state%KOUNT_LOG(1)                                   !<-- Loop through all logical data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_LOG_STRING(NPOSN_1:NPOSN_2)              !<-- The variable's name
!
        N2=N2+1
        WORK_LOGICAL=wrt_int_state%ALL_DATA_LOG(N2)                       !<-- Extract the individual data from the data string
        WRITE_LOGICAL=WORK_LOGICAL                                        !<-- Convert from ESMF_Logical to F90 logical
!
        WRITE(wrt_int_state%IO_UNIT,iostat=RC)WRITE_LOGICAL               !<-- Write out the data
!
        IF(wrt_int_state%NFHOUR==0)THEN
          WRITE(0,*)'Wrote ',TRIM(NAME),' to history file unit ',wrt_int_state%IO_UNIT
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  FIRST LOOP THROUGH ALL OF THE INTEGER Fields
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      field_loop_int: DO NFIELD=1,wrt_int_state%KOUNT_I2D(1)               !<-- Loop through all 2D integer gridded history data
!-----------------------------------------------------------------------
!
        NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
        NPOSN_2=NFIELD*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_I2D_STRING(NPOSN_1:NPOSN_2)               !<-- The variable's name

        WRITE(wrt_int_state%IO_UNIT,iostat=RC)                          &  !<-- Lead write task writes out the 2D real data
          wrt_int_state%OUTPUT_ARRAY_I2D(:,:,NFIELD) 
!
        IF(wrt_int_state%NFHOUR==0)THEN
          WRITE(0,*)'Wrote ',TRIM(NAME),' to history file unit ',wrt_int_state%IO_UNIT
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO field_loop_int
!
!-----------------------------------------------------------------------
!***  NOW LOOP THROUGH ALL THE REAL Fields
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      field_loop_real: DO NFIELD=1,wrt_int_state%KOUNT_R2D(1)              !<-- Loop through all 2D real gridded history data
!-----------------------------------------------------------------------
!
        NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
        NPOSN_2=NFIELD*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_R2D_STRING(NPOSN_1:NPOSN_2)               !<-- The name of this 2D real history quantity
!
        WRITE(wrt_int_state%IO_UNIT,iostat=RC)                          &  !<-- Lead write task writes out the 2D real data
           wrt_int_state%OUTPUT_ARRAY_R2D(:,:,NFIELD)
!
        IF(wrt_int_state%NFHOUR==0)THEN
          WRITE(0,*)'Wrote ',TRIM(NAME),' to history file unit ',wrt_int_state%IO_UNIT,NFIELD
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO field_loop_real
!
      WRITE(0,*)'Wrote all the real array'
!
!-----------------------------------------------------------------------
!***  CLOSE THE DISK FILE IF NEEDED.
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%WRITE_FLAG.and.wrt_int_state%MYPE==LEAD_WRITE_TASK)THEN
        CLOSE(wrt_int_state%IO_UNIT)
        write(0,*)' Closed history file with unit=',wrt_int_state%IO_UNIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_RUNHISTORY
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_NEMSIO_RUNHISTORY(WRT_INT_STATE                  &
                                        ,IYEAR_FCST                     &
                                        ,IMONTH_FCST                    &
                                        ,IDAY_FCST                      &
                                        ,IHOUR_FCST                     &
                                        ,IMINUTE_FCST                   &
                                        ,SECOND_FCST                    &
                                        ,NF_HOURS                       &
                                        ,NF_MINUTES                     &
                                        ,NF_SECONDS                     &
                                        ,LEAD_WRITE_TASK)
!
!-----------------------------------------------------------------------
!***  WRITE OUT A NEMSIO BINARY RUN HISTORY FILE.
!-----------------------------------------------------------------------
!
      TYPE(WRITE_INTERNAL_STATE),INTENT(INOUT) :: WRT_INT_STATE             !<-- The Write component's internal state
!
      INTEGER,INTENT(IN)  :: IYEAR_FCST                                 &
                            ,IMONTH_FCST                                &
                            ,IDAY_FCST                                  &
                            ,IHOUR_FCST                                 &
                            ,IMINUTE_FCST                               &
                            ,NF_HOURS                                   &
                            ,NF_MINUTES                                 &
                            ,LEAD_WRITE_TASK
!
      REAL,INTENT(IN)     :: NF_SECONDS                                 &
                            ,SECOND_FCST
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,J,N,N1,N2,NPOSN_1,NPOSN_2,LENGTH,MAXLENGTH
!
      INTEGER :: DIM1,DIM2,FIELDSIZE,IM,JM,LM,IDATE(7),INDX_2D,IRET     &
                ,N2ISCALAR,N2IARY,N2RSCALAR,N2RARY,N2LSCALAR            &
                ,NBDR,NDYH,NDXH,NFRAME,NPT,NPDTOP,NREC                  &
                ,NSG1,NSG2,NSGML1,NSGML2,NSOIL,TLMETA,VLEV
!
      INTEGER :: NFIELD,RC
!
      INTEGER,DIMENSION(:),POINTER :: ARYILEN                           &
                                     ,ARYRLEN                           &
                                     ,RECLEV                            &
                                     ,VARIVAL
!
      INTEGER,DIMENSION(:,:),POINTER :: ARYIVAL
!
      REAL(4) :: DEGRAD,DXCTL,DYCTL,TPH0D,TLM0D
!
      REAL(4),DIMENSION(:),POINTER :: DX,DY,DXH                         &
                                     ,GLAT1D,GLON1D                     &
                                     ,TMP
!
      REAL(4),DIMENSION(:,:,:),POINTER :: VCOORD
!
      REAL(KIND=KFPT),DIMENSION(:)  ,POINTER :: VARRVAL
      REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: ARYRVAL
!
      LOGICAL                      :: GLOBAL
      LOGICAL,DIMENSION(:),POINTER :: VARLVAL
!
      CHARACTER(6)                       :: MODEL_LEVEL
      CHARACTER(16)                       :: VLEVTYP
!
      CHARACTER(8),DIMENSION(:) ,POINTER :: ARYINAME                    &
                                           ,ARYRNAME                    &
                                           ,RECNAME                     &
                                           ,VARINAME                    &
                                           ,VARRNAME                    &
                                           ,VARLNAME
!
      CHARACTER(16),DIMENSION(:),POINTER :: RECLEVTYP
!
      CHARACTER(ESMF_MAXSTR) :: NAME,FILENAME
!
      TYPE(NEMSIO_GFILE) :: NEMSIOFILE
!
      TYPE(ESMF_Logical) :: WORK_LOGICAL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
     write(0,*)'in nemsio_runhistory'
      IDATE(1)=IYEAR_FCST
      IDATE(2)=IMONTH_FCST
      IDATE(3)=IDAY_FCST
      IDATE(4)=IHOUR_FCST
      IDATE(5)=IMINUTE_FCST
      IDATE(6)=nint(SECOND_FCST*100.)
      IDATE(7)=100
!
!-----------------------------------------------------------------------
!***  INTEGER SCALAR/1D HISTORY VARIABLES
!-----------------------------------------------------------------------
!
!-------------------------------------------------------------
!*** Find out the total number of int scalars and int arrays.
!-------------------------------------------------------------
!
      N2ISCALAR=0
      N2IARY=0
      MAXLENGTH=1
!
      DO N=1,wrt_int_state%KOUNT_I1D(1)                                    !<-- Loop through all scalar/1D integer data
        LENGTH=wrt_int_state%LENGTH_DATA_I1D(N)
!
        IF(LENGTH==1)THEN
          N2ISCALAR=N2ISCALAR+1
        ELSE
          N2IARY=N2IARY+1
          MAXLENGTH=MAX(LENGTH,MAXLENGTH)
        ENDIF
!
      ENDDO
!
      ALLOCATE(VARINAME(N2ISCALAR),VARIVAL(N2ISCALAR))
      ALLOCATE(ARYINAME(N2IARY),ARYILEN(N2IARY),ARYIVAL(MAXLENGTH,N2IARY))
!
!---------------------------------------
!***  SET VALUE TO AVRIVAL and ARYIVAL.
!---------------------------------------
!
      N2=0                                                                 !<-- Word counter for full string of integer scalar/1D data
      N2ISCALAR=0
      N2IARY=0
!
      DO N=1,wrt_int_state%KOUNT_I1D(1)                                    !<-- Loop through all scalar/1D integer data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_I1D_STRING(NPOSN_1:NPOSN_2)               !<-- The variable's name
        LENGTH=wrt_int_state%LENGTH_DATA_I1D(N)                            !<-- The variable's length in words
!
        IF(LENGTH==1)THEN
          N2=N2+1
          N2ISCALAR=N2ISCALAR+1
          VARINAME(N2ISCALAR)=TRIM(NAME)
          VARIVAL(N2ISCALAR)=wrt_int_state%ALL_DATA_I1D(N2)
        ELSE
          N2IARY=N2IARY+1
          ARYINAME(N2IARY)=TRIM(NAME)
          ARYILEN(N2IARY)=LENGTH
!            write(0,*)'in I1D array,aryiname=',aryiname(N2IARY),'len=',aryilen(N2IARY),  &
!              wrt_int_state%ALL_DATA_I1D(N2+1:N2+length)
!
          DO N1=1,LENGTH
            N2=N2+1
            ARYIVAL(N1,N2IARY)=wrt_int_state%ALL_DATA_I1D(N2)              !<-- Extract the individual data from the data string
          ENDDO
!
            write(0,*)'in I1D array,aryival=',aryival(:,N2IARY)
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  REAL SCALAR/1D HISTORY VARIABLES
!-----------------------------------------------------------------------
!
!------------------------------------------------------------
!***  Find the total number of real scalars and real arrays.
!------------------------------------------------------------
!
      N2RSCALAR=0
      N2RARY=0
      MAXLENGTH=1
!
      DO N=1,wrt_int_state%KOUNT_R1D(1)                                    !<-- Loop through all scalar/1D real data
        LENGTH=wrt_int_state%LENGTH_DATA_R1D(N)                            !<-- The variable's length
        IF(LENGTH==1)THEN
           N2RSCALAR=N2RSCALAR+1
        ELSE
          N2RARY=N2RARY+1
          MAXLENGTH=MAX(LENGTH,MAXLENGTH)
        ENDIF
      ENDDO
!
      ALLOCATE(VARRNAME(N2RSCALAR),VARRVAL(N2RSCALAR))
      ALLOCATE(ARYRNAME(N2RARY),ARYRLEN(N2RARY),ARYRVAL(MAXLENGTH,N2RARY))
!
!------------------------------------------------------
!***  Set values for the real scalars and real arrays.
!------------------------------------------------------
      N2=0                                                                 !<-- Word counter for full string of real scalar/1D data
      N2RSCALAR=0
      N2RARY=0
!
      DO N=1,wrt_int_state%KOUNT_R1D(1)                                    !<-- Loop through all scalar/1D real data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_R1D_STRING(NPOSN_1:NPOSN_2)               !<-- The variable's name
        LENGTH=wrt_int_state%LENGTH_DATA_R1D(N)                            !<-- The variable's length
!
        IF(LENGTH==1)THEN
          N2=N2+1
          N2RSCALAR=N2RSCALAR+1
          VARRNAME(N2RSCALAR)=TRIM(NAME)
          VARRVAL(N2RSCALAR)=wrt_int_state%ALL_DATA_R1D(N2)
!
          IF( TRIM(NAME)=='PT') THEN
            NPT=N2RSCALAR
          ELSEIF ( TRIM(NAME)=='PDTOP') THEN
            NPDTOP=N2RSCALAR
          ELSEIF ( TRIM(NAME)=='DYH') THEN
            NDYH=N2RSCALAR
          ELSEIF ( TRIM(NAME)=='TPH0D') THEN
            TPH0D=VARRVAL(N2RSCALAR)
          ELSEIF ( trim(NAME)=='TLM0D') THEN
            TLM0D=VARRVAL(N2RSCALAR)
          ENDIF
!
        ELSE
          N2RARY=N2RARY+1
          ARYRNAME(N2RARY)=TRIM(NAME)
          ARYRLEN(N2RARY)=LENGTH
!
          DO N1=1,LENGTH
            N2=N2+1
            ARYRVAL(N1,N2RARY)=wrt_int_state%ALL_DATA_R1D(N2)              !<-- Extract the individual data from the data string
          ENDDO
!
          IF( TRIM(NAME)=='SG1') THEN
            NSG1=N2RARY
          ELSEIF ( TRIM(NAME)=='SG2') THEN
            NSG2=N2RARY
          ELSEIF ( TRIM(NAME)=='SGML1' ) THEN
            NSGML1=N2RARY
          ELSEIF (TRIM(NAME)=='SGML2') THEN
            NSGML2=N2RARY
          ELSEIF (TRIM(NAME)=='DXH') THEN
            NDXH=N2RARY
          ENDIF

        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  LOGICAL HISTORY VARIABLES
!-----------------------------------------------------------------------
!
      N2LSCALAR=wrt_int_state%KOUNT_LOG(1)                                 !<-- Counter for full string of logical data
!
      ALLOCATE(VARLNAME(N2LSCALAR),VARLVAL(N2LSCALAR))
      N2LSCALAR=0                         
!
      DO N=1,wrt_int_state%KOUNT_LOG(1)                                    !<-- Loop through all logical data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_LOG_STRING(NPOSN_1:NPOSN_2)               !<-- The variable's name
!
        N2LSCALAR=N2LSCALAR+1
        WORK_LOGICAL=wrt_int_state%ALL_DATA_LOG(N2LSCALAR)                 !<-- Extract the individual data from the data string
        VARLNAME(N2LSCALAR)=NAME
        VARLVAL(N2LSCALAR)=WORK_LOGICAL
        IF(TRIM(NAME)=='GLOBAL') GLOBAL=WORK_LOGICAL
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  NOW OPEN NEMSIO FILE
!-----------------------------------------------------------------------
!
      write(0,*)' OPEN_NEMSIO_FILE wrt_int_state%IO_NEMSIOFILE=',        &
          trim(wrt_int_state%IO_FILE)
!
      IF(wrt_int_state%IO_FILE=='DEFERRED')THEN
        N=LEN_TRIM(wrt_int_state%FILENAME_BASE)-1
        WRITE(FILENAME,100)wrt_int_state%FILENAME_BASE(1:n)//'_nemsio.'  &
                          ,wrt_int_state%NFHOUR
        write(0,*)'FILENAME=',trim(FILENAME),'n=',n
  100   FORMAT(A21,I3.3)
      ELSE
        FILENAME=wrt_int_state%IO_FILE//'_nemsio'
      ENDIF
!
!----------------------------------------------------
!***  Prepare variables needed by the nemsip header:
!----------------------------------------------------
!
!dimension
      IF(GLOBAL) THEN
!for global im/jm for data field
        IM=wrt_int_state%im(1)+2
        JM=wrt_int_state%jm(1)+2
        DIM1=wrt_int_state%im(1)-3
        DIM2=wrt_int_state%jm(1)-2
      ELSE
!for regional
        IM=wrt_int_state%im(1)
        JM=wrt_int_state%jm(1)
        DIM1=wrt_int_state%im(1)-2
        DIM2=wrt_int_state%jm(1)-2
      ENDIF
!
      LM=wrt_int_state%LM(1)
!
!for nmmb trimmed domain
      NBDR=1
      NFRAME=0
      FIELDSIZE=DIM1*DIM2
      NREC=wrt_int_state%kount_I2D(1)+wrt_int_state%kount_R2D(1)
!
!vcoord
      ALLOCATE(VCOORD(LM+1,3,2))
      VCOORD=0.
      VCOORD(1:LM+1,1,1)=0.1*(ARYRVAL(1:LM+1,NSG1)*VARRVAL(NPDTOP)      &
         -ARYRVAL(1:LM+1,NSG2)*(VARRVAL(NPDTOP)+VARRVAL(NPT))           &
         +VARRVAL(NPT) )
      VCOORD(1:LM+1,2,1)=ARYRVAL(1:LM+1,NSG2)
      VCOORD(1:LM+1,3,1)=0
      VCOORD(1:LM,1,2)=0.1*(ARYRVAL(1:LM,NSGML1)*VARRVAL(NPDTOP)        &
         -ARYRVAL(1:LM,NSGML2)*(VARRVAL(NPDTOP)+VARRVAL(NPT))           &
         +VARRVAL(NPT) )
      VCOORD(1:LM,2,2)=ARYRVAL(1:LM,NSGML2)
      VCOORD(1:LM,3,2)=0
      write(0,*)'after vcoord,count_I2d=',wrt_int_state%kount_I2D(1),'nrec=',nrec
!
!-----------------------------------------------------------------------
!***  Cut the output I2D array.
!-----------------------------------------------------------------------
!
      ALLOCATE(RECNAME(NREC),RECLEVTYP(NREC),RECLEV(NREC))
      NREC=0
!
      DO NFIELD=1,wrt_int_state%KOUNT_I2D(1)
!
        DO J=1,DIM2
        DO I=1,DIM1
           wrt_int_state%OUTPUT_ARRAY_I2D(I,J,NFIELD)=                  &
             wrt_int_state%OUTPUT_ARRAY_I2D(I+NBDR,J+NBDR,NFIELD)
        ENDDO
        ENDDO
!
        NREC=NREC+1
        NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
        NPOSN_2=NFIELD*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_I2D_STRING(NPOSN_1:NPOSN_2)               !<-- The name of this 2D integer history quantity
        INDX_2D=index(recname(nrec),"_2D")
!
        IF (INDX_2D > 0) THEN
          MODEL_LEVEL=NAME(INDX_2D-2:INDX_2D-1)
          RECLEV(NREC)=(ICHAR(MODEL_LEVEL(1:1))-48)*10+ICHAR(MODEL_LEVEL(2:2))-48
          RECNAME(NREC)=NAME(1:INDX_2D-4)
          RECLEVTYP(NREC)='mid_layer'
          IF (RECLEV(NREC)==LM+1) RECLEVTYP(NREC-LM:NREC)='layer'
        ELSE
          RECNAME(NREC)=TRIM(NAME)
          RECLEV(NREC)=1
          RECLEVTYP(NREC)='sfc'
        ENDIF
!
        IF (RECNAME(NREC)=='ISLTYP') RECNAME(NREC)='sltyp'
        IF (RECNAME(NREC)=='IVGTYP') RECNAME(NREC)='vgtyp'
        IF (RECNAME(NREC)=='NCFRCV') RECNAME(NREC)='cfrcv'
        IF (RECNAME(NREC)=='NCFRST') RECNAME(NREC)='cfrst'
        CALL LOWERCASE(RECNAME(NREC))
      ENDDO
!
      write(0,*)'after I2D,nrec=',nrec
!
!-----------------------------------------------------------------------
!*** Cut the output R2D array.
!-----------------------------------------------------------------------
!
!      NREC=21*lm+103
      NSOIL=0
!
      DO NFIELD=1,wrt_int_state%KOUNT_R2D(1)
!
        DO J=1,DIM2
        DO I=1,DIM1
           wrt_int_state%OUTPUT_ARRAY_R2D(I,J,NFIELD)=                  &
             wrt_int_state%OUTPUT_ARRAY_R2D(I+NBDR,J+NBDR,NFIELD)
        ENDDO
        ENDDO
!
        NREC=NREC+1
        NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
        NPOSN_2=NFIELD*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_R2D_STRING(NPOSN_1:NPOSN_2)  !<-- The name of this 2D integer history quantity
        INDX_2D=INDEX(NAME,"_2D")
!
        IF (INDX_2D > 0) THEN
          MODEL_LEVEL=NAME(INDX_2D-2:INDX_2D-1)
          RECLEV(NREC)=(ICHAR(MODEL_LEVEL(1:1))-48)*10+ICHAR(MODEL_LEVEL(2:2))-48
          RECNAME(NREC)=NAME(1:INDX_2D-4)
          RECLEVTYP(NREC)='mid layer'
          IF (RECLEV(NREC)==LM+1) RECLEVTYP(NREC-LM:NREC)='layer'
          IF (RECNAME(NREC)=='SMC') NSOIL=NSOIL+1
          IF(RECNAME(NREC)=='SMC'.OR.RECNAME(NREC)=='SH2O'.or.RECNAME(NREC)=='STC') &
             RECLEVTYP(NREC)='soil layer'
          IF (RECNAME(NREC)=='W') RECNAME(NREC)='vvel'
          IF (RECNAME(NREC)=='PINT') RECNAME(NREC)='pres'
          IF (RECNAME(NREC)=='CW') RECNAME(NREC)='clwmr'
          IF (RECNAME(NREC)=='U') RECNAME(NREC)='ugrd'
          IF (RECNAME(NREC)=='V') RECNAME(NREC)='vgrd'
          IF (RECNAME(NREC)=='T') RECNAME(NREC)='tmp'
          IF (RECNAME(NREC)=='Q') RECNAME(NREC)='spfh'
!
        ELSE
          RECLEV(NREC)=1
          RECNAME(NREC)=TRIM(NAME)
          RECLEVTYP(NREC)='sfc'
!
          IF (INDEX(RECNAME(NREC),"10")>0) RECLEVTYP(NREC)='10 m above gnd'
          IF (RECNAME(NREC)=='PD') THEN 
            RECNAME(NREC)='dpres'
            RECLEVTYP(NREC)='hybrid sig lev'
          ENDIF
!
          IF (RECNAME(NREC)=='SST') RECNAME(NREC)='tsea'
          IF (RECNAME(NREC)=='FIS') RECNAME(NREC)='hgt'
          IF (RECNAME(NREC)=='USTAR') RECNAME(NREC)='uustar'
          IF (RECNAME(NREC)=='Z0') RECNAME(NREC)='zorl'
        ENDIF
!
!change unit for 'FIS'
        IF(RECNAME(NREC)=='hgt') wrt_int_state%OUTPUT_ARRAY_R2D(:,:,NFIELD)= &
          wrt_int_state%OUTPUT_ARRAY_R2D(:,:,NFIELD)/G
        write(0,*)'nfield=',nfield,'recname=',recname(nfield+4)
        CALL LOWERCASE(RECNAME(NREC))
      ENDDO
!
      write(0,'(12A8)')'after recname=',recname
      write(0,'(12A16)')'after reclevtyp=',reclevtyp
      write(0,'(12I4)')'after reclev=',reclev
!
!glat1d and glon1d
      ALLOCATE(GLAT1D(FIELDSIZE),GLON1D(FIELDSIZE))
      DEGRAD=90./ASIN(1.)
      wrt_int_state%OUTPUT_ARRAY_R2D(1:DIM1,1:DIM2,2)=wrt_int_state%OUTPUT_ARRAY_R2D(1:DIM1,1:DIM2,2)*DEGRAD
      wrt_int_state%oUTPUT_ARRAY_R2D(1:DIM1,1:DIM2,3)=wrt_int_state%OUTPUT_ARRAY_R2D(1:DIM1,1:DIM2,3)*DEGRAD
      wrt_int_state%OUTPUT_ARRAY_R2D(1:DIM1,1:DIM2,5)=wrt_int_state%OUTPUT_ARRAY_R2D(1:DIM1,1:DIM2,5)*DEGRAD
      wrt_int_state%OUTPUT_ARRAY_R2D(1:DIM1,1:DIM2,6)=wrt_int_state%OUTPUT_ARRAY_R2D(1:DIM1,1:DIM2,6)*DEGRAD
      GLAT1D(1:FIELDSIZE)=RESHAPE(wrt_int_state%OUTPUT_ARRAY_R2D(1:DIM1,1:DIM2,2),(/FIELDSIZE/))
      GLON1D(1:FIELDSIZE)=RESHAPE(wrt_int_state%OUTPUT_ARRAY_R2D(1:DIM1,1:DIM2,3),(/FIELDSIZE/))
      write(0,*)'after glat1d=',maxval(glat1d),minval(glat1d),           &
        'glon1d=',maxval(glon1d),minval(glon1d)
!
!dx and dy
      ALLOCATE(DX(FIELDSIZE),DY(FIELDSIZE))
!
      DO I=1,FIELDSIZE
        DY(I)=VARRVAL(NDYH)
      ENDDO
!
      DO J=1,DIM2
      DO I=1,DIM1
        DX(I+(J-1)*DIM1)=ARYRVAL(J+NBDR,NDXH)
      ENDDO
      ENDDO
      write(0,*)'after dx=',maxval(dx),minval(dx),'dy=',maxval(dy),minval(dy)
!
!-----------------------------------------------------------------------
!                      SET UP NEMSIO WRITE
!-----------------------------------------------------------------------
!
      CALL NEMSIO_INIT(IRET=IRET)
      print *,'nemsio_init, iret=',iret

!-----------------------------------------------------------------------
!***  OPEN NEMSIO FILE
!-----------------------------------------------------------------------
!
      CALL NEMSIO_OPEN(nemsiofile,trim(FILENAME),'write',iret,           &
        modelname="NMMB", gdatatype="bin4", idate=IDATE,nfhour=NF_HOURS, &
        nfminute=NF_MINUTES,nfsecondn=nint(NF_SECONDS*100),              &
        nfsecondd=100,dimx=DIM1,dimy=DIM2,dimz=LM,nframe=NFRAME,         &
        nsoil=NSOIL,ntrac=3, ncldt=1,rlon_min=minval(glon1d),            &
        rlon_max=maxval(glon1d), rlat_max=maxval(glat1d),                &
        rlat_min=minval(glat1d),vcoord=vcoord,lon=glon1d,lat=glat1d,     &
        dx=dx,dy=dy,extrameta=.true.,nmetavari=N2ISCALAR,                &
        nmetavarr=N2RSCALAR,nmetavarl=N2LSCALAR,nmetaaryi=N2IARY,        &
        nmetaaryr=N2RARY,variname=VARINAME,varival=VARIVAL,              &
        varrname=VARRNAME,varrval=VARRVAL,varlname=VARLNAME,             &
        varlval=VARLVAL,aryiname=ARYINAME,aryilen=ARYILEN,               &
        aryival=ARYIVAL,aryrname=ARYRNAME,aryrlen=ARYRLEN,               &
        aryrval=ARYRVAL,recname=recname,reclevtyp=reclevtyp,reclev=reclev)
!
      ALLOCATE(TMP(FIELDSIZE))
!
!-----------------------------------------------------------------------
!***  GET VARIABLES NEEDED BY THE .ctl FILE.
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETFILEHEAD(NEMSIOFILE,TLMETA=TLMETA)
      DXCTL=MAXVAL(DX)*180./(A*PI)
      DYCTL=MAXVAL(DY)*180./(A*PI)
!
!-----------------------------------------------------------------------
!***  FIRST LOOP THROUGH ALL OF THE INTEGER Fields.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      field_loop_int: DO NFIELD=1,wrt_int_state%KOUNT_I2D(1)              !<-- Loop through all 2D integer gridded history data
!-----------------------------------------------------------------------
!
        CALL NEMSIO_GETRECHEAD(NEMSIOFILE,NFIELD,NAME,VLEVTYP,VLEV,IRET=IRET)
            write(0,*)'nfield=',nfield,trim(name),vlevtyp,vlev
        TMP=RESHAPE(wrt_int_state%OUTPUT_ARRAY_I2D(1:DIM1,1:DIM2,NFIELD),(/FIELDSIZE/))
        CALL NEMSIO_WRITEREC(NEMSIOFILE,NFIELD,TMP,IRET=IRET)             !<-- Lead write task writes out the 2D int data
!
!-----------------------------------------------------------------------
!
      ENDDO field_loop_int
!
!-----------------------------------------------------------------------
!***  NOW LOOP THROUGH ALL THE REAL Fields
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      field_loop_real: DO NFIELD=1,wrt_int_state%KOUNT_R2D(1)             !<-- Loop through all 2D real gridded history data
!-----------------------------------------------------------------------
!
        N=NFIELD+wrt_int_state%KOUNT_I2D(1)
        CALL NEMSIO_GETRECHEAD(NEMSIOFILE,N,NAME,VLEVTYP,VLEV,IRET=IRET)
            write(0,*)'nfield=',nfield,trim(name),vlevtyp,vlev
        TMP=RESHAPE(wrt_int_state%OUTPUT_ARRAY_R2D(1:DIM1,1:DIM2,NFIELD),(/FIELDSIZE/))
        CALL NEMSIO_WRITEREC(NEMSIOFILE,N,TMP,IRET=IRET)                  !<-- Lead write task writes out the 2D real data
!
!-----------------------------------------------------------------------
!
      ENDDO field_loop_real
!
!-----------------------------------------------------------------------
!***  CLOSE THE DISK FILE IF NEEDED.
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%WRITE_NEMSIOFLAG.and.wrt_int_state%MYPE==LEAD_WRITE_TASK)THEN
        CALL NEMSIO_CLOSE(NEMSIOFILE)
        WRITE(0,*)' Closed nemsio_history file'
        CALL NEMSIO_FINALIZE()
      ENDIF
!
!-----------------------------------------------------------------------
!***  WRITE OUT NEMSIO CTL FILE
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%WRITE_NEMSIOCTL.AND.wrt_int_state%MYPE==LEAD_WRITE_TASK)THEN
        CALL WRITE_NEMSIOCTL(GLOBAL,IHOUR_FCST,IDAY_FCST,IMONTH_FCST,       &
          IYEAR_FCST,FILENAME,TLMETA,DIM1,DIM2,LM,NSOIL,TLM0D,TPH0D,DXCTL,  & 
          DYCTL,NF_HOURS,NREC,RECNAME,RECLEVTYP)
      ENDIF
!
!-----------------------------------------------------------------------
!***  CLEAN UP
!-----------------------------------------------------------------------
!
      DEALLOCATE(VCOORD,GLAT1D,GLON1D,DX,DY,TMP)
      DEALLOCATE(VARINAME,VARIVAL,ARYINAME,ARYILEN,ARYIVAL)
      DEALLOCATE(VARRNAME,VARRVAL,ARYRNAME,ARYRLEN,ARYRVAL)
      DEALLOCATE(VARLNAME,VARLVAL)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_NEMSIO_RUNHISTORY
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_NEMSIOCTL(GLOBAL,IHOUR_FCST,IDAY_FCST,IMONTH_FCST, &
        IYEAR_FCST,FILENAME,TLMETA,DIM1,DIM2,LM,NSOIL,TLM0D,TPH0D,DXCTL,  &
        DYCTL,NF_HOURS,NREC,RECNAME,RECLEVTYP)
!
!-----------------------------------------------------------------------
!***  WRITE OUT CTL FILE
!-----------------------------------------------------------------------
!
!--------------
!*** Arguments
!--------------
!
      INTEGER,INTENT(IN) :: DIM1,DIM2                                   &
                           ,IHOUR_FCST,IDAY_FCST,IMONTH_FCST,IYEAR_FCST &
                           ,LM,NF_HOURS,NREC,NSOIL,TLMETA
!
      REAL,INTENT(IN)    :: dxctl,dyctl,tlm0d,tph0d
!
      LOGICAL,INTENT(IN) :: GLOBAL
!
      CHARACTER(*) ,INTENT(IN) :: FILENAME
      CHARACTER(8) ,INTENT(IN) :: RECNAME(:)
      CHARACTER(16),INTENT(IN) :: RECLEVTYP(:)
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER N,IO_UNIT
!
      CHARACTER(3)  CMON
      CHARACTER(32) DATE
!
      LOGICAL OPENED
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!---------------------
!***  Get unit number
!---------------------
!
      DO N=51,99
        INQUIRE(N,opened=OPENED)
          IF(.NOT.OPENED)THEN
            IO_UNIT=N
            EXIT
        ENDIF
      ENDDO
!
      CALL CMONTH(IMONTH_FCST,CMON)
      WRITE(DATE,'(I2.2,A,I2.2,A3,I4.4)')IHOUR_FCST,'Z',IDAY_FCST       &
           ,CMON,IYEAR_FCST
      OPEN(IO_UNIT,file=TRIM(FILENAME)//'.ctl',form='formatted')
!
      WRITE(IO_UNIT,105)TRIM(FILENAME)
      WRITE(IO_UNIT,106)
      WRITE(IO_UNIT,107)
      WRITE(IO_UNIT,108)TLMETA
      WRITE(IO_UNIT,109)
!
      IF (GLOBAL) THEN
        WRITE(IO_UNIT,121)DIM1,DXCTL
        WRITE(IO_UNIT,122)DIM2,DYCTl
      ELSE
        WRITE(IO_UNIT,110)DIM1,DIM2,TLM0D,TPH0D,DXCTL,DYCTL
        WRITE(IO_UNIT,111)
        WRITE(IO_UNIT,112)
      ENDIF
!
      WRITE(IO_UNIT,113)LM
      WRITE(IO_UNIT,114)1,TRIM(DATE)
!
 105  FORMAT('dset ^',A24)
 106  FORMAT('undef -9.E+20')
 107  FORMAT('options big_endian sequential')
 108  FORMAT('fileheader',I12.0)
 109  FORMAT('title EXP1')

 110  FORMAT('pdef ',I6,I6,' eta.u ',f8.1,f8.1,f12.6,f12.6)
 111  FORMAT('xdef  320 linear  -137.000  0.25')
 112  FORMAT('ydef  140 linear    19.000  0.25')
 121  FORMAT('xdef  ',I6,' linear  -180.000 ',f12.6)
 122  FORMAT('ydef  ',I6,' linear   -90.000  ',f12.6)

 113  FORMAT('zdef ',I6,' linear 1 1 ')
 114  FORMAT('tdef ',I6,' linear ',A12,' 6hr')
!
      WRITE(IO_UNIT,'(A)')'VARS 114'
      N=1
!
      DO WHILE (N<=NREC)
        IF(RECLEVTYP(N)=='mid layer') THEN
          WRITE(IO_UNIT,'(A8,I3,A)')RECNAME(N),LM,' 99 mid layer'
          N=N+LM
        ELSEIF(RECLEVTYP(N)=='layer') THEN
          WRITE(IO_UNIT,'(A8,I3,A)')RECNAME(N),LM+1,' 99 layer'
          N=N+LM+1
        ELSEIF(RECLEVTYP(N)=='soil layer') THEN
          WRITE(IO_UNIT,'(A8,I3,A)')RECNAME(N),NSOIL,' 99 soil layer'
          N=N+NSOIL
        ELSE
          WRITE(IO_UNIT,'(A8,A)')RECNAME(N),'  0 99 sfc'
          N=N+1
        ENDIF
      ENDDO
!
      WRITE(IO_UNIT,'(A8)')'endvars'
      CLOSE(IO_UNIT)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_NEMSIOCTL
!
!-----------------------------------------------------------------------
!
      elemental subroutine lowercase(word)
!
!-----------------------------------------------------------------------
!***  convert a word to lower case
!-----------------------------------------------------------------------
!
      character (len=*) , intent(inout) :: word
      integer :: i,ic,nlen
      nlen = len(word)
!
      do i=1,nlen
        ic = ichar(word(i:i))
        if (ic >= 65 .and. ic < 91) word(i:i) = char(ic+32)
      end do
!
!-----------------------------------------------------------------------
!
      end subroutine lowercase 
!
!-----------------------------------------------------------------------
!
      SUBROUTINE CMONTH(IMON,CMON)
!
!-----------------------------------------------------------------------
!***  CONVERT MONTH
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IMON
      CHARACTER(LEN=3)   :: CMON
!
!-----------------------------------------------------------------------
!
      SELECT CASE (IMON)
        CASE(1) 
            CMON='Jan'
        CASE(2) 
            CMON='Feb'
        CASE(3) 
            CMON='Mar'
        CASE(4) 
            CMON='Apr'
        CASE(5) 
            CMON='May'
        CASE(6) 
            CMON='Jun'
        CASE(7) 
            CMON='Jul'
        CASE(8) 
            CMON='Aug'
        CASE(9) 
            CMON='Sep'
        CASE(10) 
            CMON='Oct'
        CASE(11) 
            CMON='Nov'
        CASE(12) 
            CMON='Dec'
      END SELECT
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CMONTH
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_WRITE_ROUTINES
!
!-----------------------------------------------------------------------
