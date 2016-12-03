#!/bin/bash
#set -eu

mkdir -p ${RUNDIR}

####################################################################################################
# Make configure and run files
####################################################################################################

source atparse.bash

wrtdopost="${WRITE_DOPOST}"
postgrbvs="${POST_GRIBVERSION}"
SRCD="${PATHTR}"
RUND="${RUNDIR}"

atparse > nmm_run < nmm_conf/nmm_${GBRG}_run.IN

atparse < nmm_conf/nmm_${GBRG}_conf.IN > configure_file_01

if [ ${nems_configure}"x" == "x" ]; then
  nems_configure=atm_nostep
  atm_model=nmm
fi
atparse < nems.configure.${nems_configure}.IN > nems.configure
atparse < atmos.configure_nmm > atmos.configure

if [ $SCHEDULER = 'moab' ]; then

atparse < nmm_conf/nmm_msub.IN > nmm_msub

elif [ $SCHEDULER = 'pbs' ]; then

atparse < nmm_conf/nmm_qsub.IN > nmm_qsub

elif [ $SCHEDULER = 'lsf' ]; then

atparse < nmm_conf/nmm_bsub.IN > nmm_bsub

fi

if [ ${GBRG} = nests ]; then
  cat nmm_conf/nmm_nests_conf_02.IN | atparse > configure_file_02

  cat nmm_conf/nmm_nests_conf_03.IN | atparse > configure_file_03

  cat nmm_conf/nmm_nests_conf_04.IN | atparse > configure_file_04

fi

if [ ${GBRG} = mnests ]; then
  rm -f configure_file_02 configure_file_03 configure_file_04
  cat nmm_conf/nmm_mnests_conf_02.IN | atparse > configure_file_02
  cat nmm_conf/nmm_mnests_conf_03.IN | atparse > configure_file_03
  cat nmm_conf/nmm_mnests_conf_04.IN | atparse > configure_file_04
fi

if [ ${MODE} = 2-way  ]; then
  rm -f configure_file_02 configure_file_03 configure_file_04
  cat nmm_conf/nmm_mnests_2way_conf_02.IN | atparse > configure_file_02
  cat nmm_conf/nmm_mnests_2way_conf_03.IN | atparse > configure_file_03
  cat nmm_conf/nmm_mnests_2way_conf_04.IN | atparse > configure_file_04
fi

if [ ${GBRG} = fltr ]; then
  rm -f configure_file_02 configure_file_03 configure_file_04
  cp nmm_conf/nmm_fltr_conf_02 configure_file_02
  cp nmm_conf/nmm_fltr_conf_03 configure_file_03
fi

if [ ${GBRG} = fltr_zombie ]; then
  rm -f configure_file_02 configure_file_03 configure_file_04
  cp nmm_conf/nmm_fltr_conf_02 configure_file_02
  cp nmm_conf/nmm_fltr_zombie_conf_03 configure_file_03
fi

####################################################################################################
# Submit test
####################################################################################################

sh ./nmm_run

# wait for the job to enter the queue
count=0
job_running=0
until [ $job_running -eq 1 ]
do
echo "TEST is waiting to enter the queue"
if [ $SCHEDULER = 'moab' ]; then
  job_running=`showq -u ${USER} -n | grep ${JBNME} | wc -l`;sleep 5
elif [ $SCHEDULER = 'pbs' ]; then
  job_running=`qstat -u ${USER} -n | grep ${JBNME} | wc -l`;sleep 5
elif [ $SCHEDULER = 'lsf' ]; then
  job_running=`bjobs -u ${USER} -J ${JBNME} 2>/dev/null | grep ${QUEUE} | wc -l`;sleep 5
fi
(( count=count+1 )) ; if [ $count -eq 13 ] ; then echo "No job in queue after one minute, exiting..." ; exit 2 ; fi
done

# wait for the job to finish and compare results
job_running=1
n=1
until [ $job_running -eq 0 ]
do

sleep 30
if [ $SCHEDULER = 'moab' ]; then
  job_running=`showq -u ${USER} -n | grep ${JBNME} | wc -l`
elif [ $SCHEDULER = 'pbs' ]; then
  job_running=`qstat -u ${USER} -n | grep ${JBNME} | wc -l`
elif [ $SCHEDULER = 'lsf' ]; then
  job_running=`bjobs -u ${USER} -J ${JBNME} 2>/dev/null | wc -l`
fi

if [ $SCHEDULER = 'moab' ]; then

  status=`showq -u ${USER} -n | grep ${JBNME} | awk '{print $3}'` ; status=${status:--}
  if [ -f ${RUNDIR}/err ] ; then FnshHrs=`grep Finished ${RUNDIR}/err | tail -1 | awk '{ print $10 }'` ; fi
  FnshHrs=${FnshHrs:-0}
  if   [ $status = 'Idle' ];       then echo "$n/2min TEST ${TEST_NR} is waiting in a queue, Status: " $status
  elif [ $status = 'Running' ];    then echo "$n/2min TEST ${TEST_NR} is running,            Status: " $status  ", Finished " $FnshHrs "hours"
  elif [ $status = 'Starting' ];   then echo "$n/2min TEST ${TEST_NR} is ready to run,       Status: " $status  ", Finished " $FnshHrs "hours"
  elif [ $status = 'Completed' ];  then echo "$n/2min TEST ${TEST_NR} is finished,           Status: " $status ; job_running=0
  else                                  echo "$n/2min TEST ${TEST_NR} is finished,           Status: " $status  ", Finished " $FnshHrs "hours"
  fi

elif [ $SCHEDULER = 'pbs' ]; then

  status=`qstat -u ${USER} -n | grep ${JBNME} | awk '{print $"10"}'` ; status=${status:--}
  if [ -f ${RUNDIR}/err ] ; then FnshHrs=`tail -100 ${RUNDIR}/err | grep Finished | tail -1 | awk '{ print $10 }'` ; fi
  FnshHrs=${FnshHrs:-0}
  if   [ $status = 'Q' ];  then echo "$n/2min TEST ${TEST_NR} is waiting in a queue, Status: " $status
  elif [ $status = 'H' ];  then echo "$n/2min TEST ${TEST_NR} is held in a queue,    Status: " $status
  elif [ $status = 'R' ];  then echo "$n/2min TEST ${TEST_NR} is running,            Status: " $status  ", Finished " $FnshHrs "hours"
  elif [ $status = 'E' -o $status = 'C' ];  then
    jobid=`qstat -u ${USER} | grep ${JBNME} | awk '{print $1}'`
    exit_status=`qstat ${jobid} -f | grep exit_status | awk '{print $3}'`
    if [ $exit_status != 0 ]; then
      echo "Test ${TEST_NR} FAIL " >> ${REGRESSIONTEST_LOG}
      (echo;echo;echo)             >> ${REGRESSIONTEST_LOG}
      echo "Test ${TEST_NR} FAIL "
      (echo;echo;echo)
      echo $TEST_NAME >> fail_test
      exit 0
    fi
    echo "$n/2min TEST ${TEST_NR} is finished,           Status: " $status
    job_running=0
  elif [ $status = 'C' ];  then echo "$n/2min TEST ${TEST_NR} is finished,           Status: " $status ; job_running=0
  else                          echo "$n/2min TEST ${TEST_NR} is finished,           Status: " $status  ", Finished " $FnshHrs "hours"
  fi

elif [ $SCHEDULER = 'lsf' ]; then

  status=`bjobs -u ${USER} -J ${JBNME} 2>/dev/null | grep ${QUEUE} | awk '{print $3}'` ; status=${status:--}
#  if [ $status != '-' -a $status != 'PEND' ] ; then FnshHrs=`bpeek -J ${JBNME} | grep Finished | tail -1 | awk '{ print $10 }'` ; fi
  if [ -f ${RUNDIR}/err ] ; then FnshHrs=`grep Finished ${RUNDIR}/err | tail -1 | awk '{ print $10 }'` ; fi
  FnshHrs=${FnshHrs:-0}
  if   [ $status = 'PEND' ];  then echo "$n/2min TEST ${TEST_NR} is waiting in a queue, Status: " $status
  elif [ $status = 'RUN'  ];  then echo "$n/2min TEST ${TEST_NR} is running,            Status: " $status  ", Finished " $FnshHrs "hours"
  elif [ $status = 'EXIT' ];  then
    echo "Test ${TEST_NR} FAIL " >> ${REGRESSIONTEST_LOG}
    (echo;echo;echo)             >> ${REGRESSIONTEST_LOG}
    echo "Test ${TEST_NR} FAIL "
    (echo;echo;echo)
    echo $TEST_NAME >> fail_test
    exit 0
  else                             echo "$n/2min TEST ${TEST_NR} is finished,           Status: " $status  ", Finished " $FnshHrs "hours"
    exit_status=`bjobs -u ${USER} -J ${JBNME} -a 2>/dev/null | grep $QUEUE | awk '{print $3}'`
    if [ $exit_status = 'EXIT' ];  then
      echo "Test ${TEST_NR} FAIL " >> ${REGRESSIONTEST_LOG}
      (echo;echo;echo)             >> ${REGRESSIONTEST_LOG}
      echo "Test ${TEST_NR} FAIL "
      (echo;echo;echo)
      echo $TEST_NAME >> fail_test
      exit 0
    fi
  fi


fi
(( n=n+1 ))
done

####################################################################################################
# Check results
####################################################################################################

test_status='PASS'

# Give 10 seconds for data to show up on file system
sleep 10

(echo;echo;echo "baseline dir = ${RTPWD}/${CNTL_DIR}")  >> ${REGRESSIONTEST_LOG}
           echo "working dir  = ${RUNDIR}"              >> ${REGRESSIONTEST_LOG}
           echo "Checking test ${TEST_NR} results ...." >> ${REGRESSIONTEST_LOG}
(echo;echo;echo "baseline dir = ${RTPWD}/${CNTL_DIR}")
           echo "working dir  = ${RUNDIR}"
           echo "Checking test ${TEST_NR} results ...."

#
     if [ ${CREATE_BASELINE} = false ]; then
#
# --- regression test comparison ----
#

for i in ${LIST_FILES}
do
printf %s " Comparing " $i "....." >> ${REGRESSIONTEST_LOG}
printf %s " Comparing " $i "....."

if [ ! -f ${RUNDIR}/$i ] ; then

  echo ".......MISSING file" >> ${REGRESSIONTEST_LOG}
  echo ".......MISSING file"
  test_status='FAIL'

elif [ ! -f ${RTPWD}/${CNTL_DIR}/$i ] ; then

  echo ".......MISSING baseline" >> ${REGRESSIONTEST_LOG}
  echo ".......MISSING baseline"
  test_status='FAIL'

else

  d=`cmp ${RTPWD}/${CNTL_DIR}/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
#  (echo " ......***NOT OK***" ; echo ; echo "   $i differ!   ")>> ${REGRESSIONTEST_LOG}
#   echo " ......***NOT OK***" ; echo ; echo "   $i differ!   " ; exit 2
    echo ".......***NOT OK***" >> ${REGRESSIONTEST_LOG} 
    echo ".......***NOT OK***" 
    test_status='FAIL' 
    if [ ${BAIL_CONDITION}"x" = FILE"x" ]; then 
      echo "BAIL_CONDITION=FILE, Abort testing on failure" 
      exit 2 
    fi
     
  else

    echo "....OK" >> ${REGRESSIONTEST_LOG}
    echo "....OK"
  fi

fi

done

if [ $test_status = 'FAIL' ]; then echo $TEST_NAME >> fail_test ; fi

#
     else
#
# --- create baselines
#

 echo;echo;echo "Moving set ${TEST_NR} files ...."

for i in ${LIST_FILES}
do
  printf %s " Moving " $i "....."
  if [ -f ${RUNDIR}/$i ] ; then
    cp ${RUNDIR}/${i} ${STMP}/${USER}/REGRESSION_TEST/${CNTL_DIR}/${i}
  else
    echo "Missing " ${RUNDIR}/$i " output file"
    echo;echo " Set ${TEST_NR} failed "
    exit 2
  fi
done

# ---
     fi
# ---

echo "Test ${TEST_NR} ${test_status} " >> ${REGRESSIONTEST_LOG}
(echo;echo;echo)                       >> ${REGRESSIONTEST_LOG}
echo "Test ${TEST_NR} ${test_status} "
(echo;echo;echo)

sleep 4
echo;echo

####################################################################################################
# End test
####################################################################################################

exit 0