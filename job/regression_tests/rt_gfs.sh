#!/bin/ksh
set -ua

####################################################################################################
# Make configure and run files
####################################################################################################

JBNME=NEMS_RT_${TEST_NR}_$$

cat gfs_ll.IN       | sed s:_JBNME_:${JBNME}:g   \
                    | sed s:_CLASS_:${CLASS}:g   \
                    | sed s:_ACCNR_:${ACCNR}:g   \
                    | sed s:_TASKS_:${TASKS}:g   \
                    | sed s:_THRDS_:${THRD}:g    >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:${TASKS}:g   \
                    | sed s:_PE1_:${PE1}:g       \
                    | sed s:_WPG_:${WTPG}:g      \
                    | sed s:_WRTGP_:${WRTGP}:g   \
                    | sed s:_THRDS_:${THRD}:g    \
                    | sed s:_NSOUT_:${NSOUT}:g   \
                    | sed s:_QUILT_:${QUILT}:g   \
                    | sed s:_CP2_:${CP2}:g       \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_FDFI_:${FDFI}:g \
                    | sed s:_REDUCEDGRID_:${REDUCEDGRID}:g \
                    | sed s:_ADIAB_:${ADIAB}:g \
                    | sed s:_NUMFILE_:${NUMFILE}:g \
                    | sed s:_NDAYS_:${NDAYS}:g   >  gfs_fcst_run

####################################################################################################
# Copy init files
####################################################################################################

mkdir -p ${RUNDIR}
cp Chem_Registry.rc ${RUNDIR}/Chem_Registry.rc
cp MAPL.rc ${RUNDIR}/MAPL.rc
cp ${RTPWD}/GFS_NODFI/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS_NODFI/sfcanl.2009072400 ${RUNDIR}/.

####################################################################################################
# Submit test
####################################################################################################

llsubmit gfs_ll 2>&1 | grep submitted > /dev/null

echo "Test ${TEST_NR}" >> RegressionTests.log
echo "Test ${TEST_NR}"
echo ${TEST_DESCR} >> RegressionTests.log
echo ${TEST_DESCR}
(echo "GFS, ${TASKS} proc, ${THRD} thread";echo;echo)>> RegressionTests.log
 echo "GFS, ${TASKS} proc, ${THRD} thread";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq -u ${LOGIN} -f %st %jn | grep ${JBNME} | awk '{ print $1}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST ${TEST_NR} is waiting in a queue, Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST ${TEST_NR} is running,            Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST ${TEST_NR} is ready to run,       Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status
else                          echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status
fi

job_running=`llq -u ${LOGIN} -f %st %jn | grep ${JBNME} | wc -l`
  (( n=n+1 ))
done

####################################################################################################
# Check results
####################################################################################################

(echo;echo;echo "Checking test ${TEST_NR} results ....")>> RegressionTests.log
 echo;echo;echo "Checking test ${TEST_NR} results ...."

for i in ${LIST_FILES}

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${CNTL_PTH}/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit 2
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test ${TEST_NR} failed ")>> RegressionTests.log
  echo;echo " Test ${TEST_NR} failed "
  exit 2

fi

done

echo " Test ${TEST_NR} passed " >> RegressionTests.log
echo " Test ${TEST_NR} passed "

sleep 4
clear;echo;echo

####################################################################################################
# End test
####################################################################################################

exit 0
