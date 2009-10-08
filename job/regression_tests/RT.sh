#!/bin/ksh

date > Compile.log
date > RegressionTests.log
echo "Start Regression test" >> RegressionTests.log

#########################################
# RTPWD - path with previous stored data
#########################################

export RTPWD=/meso/noscrub/wx20rv/REGRESSION_TEST
#export RTPWD_GFS=/climate/noscrub/wx20wa/esmf/nems/REGRESSION_TEST

#########################################
# PATHRT - regression test path
#########################################
export PATHRT=`pwd`
cd ../../
#########################################
# PATHTR - NEMS path
#########################################
export PATHTR=`pwd`

RUNDIR_ROOT=/ptmp/${LOGIN}/RT_$$
rm -rf ${RUNDIR_ROOT}

clear;echo;echo

####################################################################################################
# Clean and compile NMMB core
####################################################################################################
echo "Preparing NMMB core for regression tests"
printf %s "Compiling NMMB core (this will take ~10 minutes)......."
cd ush

./clean_stub.sh                        >> ${PATHRT}/Compile.log 2>&1
./clean.sh                             >> ${PATHRT}/Compile.log 2>&1
./compile_configure.sh nmm gfs_physics >> ${PATHRT}/Compile.log 2>&1
./compile.sh                           >> ${PATHRT}/Compile.log 2>&1

if [ -f ../exe/NMM_NEMS.x ] ; then
  echo "   NMMB core Compiled";echo;echo
else
  echo "   NMMB core is NOT compiled" >> ${PATHRT}/RegressionTests.log
  echo "   NMMB core is NOT compiled"
  exit
fi

cd $PATHRT

####################################################################################################
# Submit test 1
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NMM_CNTRL

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_glob_ll.IN   | sed s:_TPN_:32:g         \
                     | sed s:_THRD_:1:g         \
                     | sed s:_GS_:#:g           \
                     | sed s:_RTPWD_:${RTPWD}:g \
                     | sed s:_SRCD_:${PATHTR}:g \
                     | sed s:_RUND_:${RUNDIR}:g  >  nmm_glob_ll

cat nmm_glob_conf.IN | sed s:_INPES_:06:g       \
                     | sed s:_WTPG_:2:g         \
                     | sed s:_FCSTL_:48:g       \
                     | sed s:_NEMSI_:false:g    \
                     | sed s:_RSTRT_:false:g    \
                     | sed s:_gfsP_:false:g     \
                     | sed s:_RGS_:false:g      \
                     | sed s:_WGS_:true:g       >  configure_file

job_id=`llsubmit nmm_glob_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 1" >> RegressionTests.log
echo "Test 1"
echo "Compare results with previous trunk version" >> RegressionTests.log
echo "Compare results with previous trunk version"
(echo "NMMB, global case";echo;echo)>> RegressionTests.log
 echo "NMMB, global case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 1 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 1 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 1 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 1 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 1 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 1 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 1 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
         nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_glob/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 1 failed ")>> RegressionTests.log
  echo;echo " Test 1 failed "
  exit

fi

done

echo " Test 1 passed " >> RegressionTests.log
echo " Test 1 passed "
export timing1=`grep total_tim $PATHRT/err | tail -1 | awk '{ print $5 }'`
export timingc=`cat ${RTPWD}/NMMB_glob/timing.txt`
(echo;echo " Original timing: " $timingc " , test1 timing: " $timing1)>> RegressionTests.log
 echo;echo " Original timing: " $timingc " , test1 timing: " $timing1

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 2
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NMM_NEMSIO

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_glob_ll.IN   | sed s:_TPN_:32:g         \
                     | sed s:_THRD_:1:g         \
                     | sed s:_GS_:#:g           \
                     | sed s:_RTPWD_:${RTPWD}:g \
                     | sed s:_SRCD_:${PATHTR}:g \
                     | sed s:_RUND_:${RUNDIR}:g >  nmm_glob_ll

cat nmm_glob_conf.IN | sed s:_INPES_:06:g       \
                     | sed s:_WTPG_:2:g         \
                     | sed s:_FCSTL_:48:g       \
                     | sed s:_NEMSI_:true:g     \
                     | sed s:_RSTRT_:false:g    \
                     | sed s:_gfsP_:false:g     \
                     | sed s:_RGS_:false:g      \
                     | sed s:_WGS_:false:g      >  configure_file

job_id=`llsubmit nmm_glob_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 2" >> RegressionTests.log
echo "Test 2"
echo "Test NEMSIO as input file" >> RegressionTests.log
echo "Test NEMSIO as input file"
(echo "NMMB, global case";echo;echo)>> RegressionTests.log
 echo "NMMB, global case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 2 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 2 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 2 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 2 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 2 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 2 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 2 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
         nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_glob/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 2 failed ")>> RegressionTests.log
  echo;echo " Test 2 failed "
  exit

fi

done

echo " Test 2 passed " >> RegressionTests.log
echo " Test 2 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 3
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NMM_REST

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_glob_ll.IN   | sed s:_TPN_:32:g         \
                     | sed s:_THRD_:1:g         \
                     | sed s:_GS_:#:g           \
                     | sed s:_RTPWD_:${RTPWD}:g \
                     | sed s:_SRCD_:${PATHTR}:g \
                     | sed s:_RUND_:${RUNDIR}:g >  nmm_glob_ll

cat nmm_glob_conf.IN | sed s:_INPES_:06:g       \
                     | sed s:_WTPG_:2:g         \
                     | sed s:_FCSTL_:48:g       \
                     | sed s:_NEMSI_:false:g    \
                     | sed s:_RSTRT_:true:g     \
                     | sed s:_gfsP_:false:g     \
                     | sed s:_RGS_:false:g      \
                     | sed s:_WGS_:false:g      >  configure_file

job_id=`llsubmit nmm_glob_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 3" >> RegressionTests.log
echo "Test 3"
echo "Test restart run" >> RegressionTests.log
echo "Test restart run"
(echo "NMMB, global case";echo;echo)>> RegressionTests.log
 echo "NMMB, global case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 3 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 3 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 3 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 3 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 3 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 3 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 3 results ...."

#
# NOTE! first hitory files differ because of PSFC, others should be identical:
#
for i in nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_glob/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 3 failed ")>> RegressionTests.log
  echo;echo " Test 3 failed "
  exit

fi

done

echo " Test 3 passed " >> RegressionTests.log
echo " Test 3 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 4
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NMM_REST_NIO

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_glob_ll.IN   | sed s:_TPN_:32:g         \
                     | sed s:_THRD_:1:g         \
                     | sed s:_GS_:#:g           \
                     | sed s:_RTPWD_:${RTPWD}:g \
                     | sed s:_SRCD_:${PATHTR}:g \
                     | sed s:_RUND_:${RUNDIR}:g >  nmm_glob_ll

cat nmm_glob_conf.IN | sed s:_INPES_:06:g       \
                     | sed s:_WTPG_:2:g         \
                     | sed s:_FCSTL_:48:g       \
                     | sed s:_NEMSI_:true:g     \
                     | sed s:_RSTRT_:true:g     \
                     | sed s:_gfsP_:false:g     \
                     | sed s:_RGS_:false:g      \
                     | sed s:_WGS_:false:g      >  configure_file

job_id=`llsubmit nmm_glob_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 4" >> RegressionTests.log
echo "Test 4"
echo "Test restart run from NEMSIO file" >> RegressionTests.log
echo "Test restart run from NEMSIO file"
(echo "NMMB, global case";echo;echo)     >> RegressionTests.log
 echo "NMMB, global case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 4 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 4 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 4 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 4 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 4 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 4 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 4 results ...."

#
# NOTE! first hitory files differ because of PSFC, others should be identical:
#
for i in nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_glob/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 4 failed ")>> RegressionTests.log
  echo;echo " Test 4 failed "
  exit

fi

done

echo " Test 4 passed " >> RegressionTests.log
echo " Test 4 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 5
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NMM_DECOMP

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_glob_ll.IN   | sed s:_TPN_:16:g         \
                     | sed s:_THRD_:1:g         \
                     | sed s:_GS_::g            \
                     | sed s:_RTPWD_:${RTPWD}:g \
                     | sed s:_SRCD_:${PATHTR}:g \
                     | sed s:_RUND_:${RUNDIR}:g >  nmm_glob_ll

cat nmm_glob_conf.IN | sed s:_INPES_:03:g       \
                     | sed s:_WTPG_:1:g         \
                     | sed s:_FCSTL_:24:g       \
                     | sed s:_NEMSI_:false:g    \
                     | sed s:_RSTRT_:false:g    \
                     | sed s:_gfsP_:false:g     \
                     | sed s:_RGS_:true:g       \
                     | sed s:_WGS_:false:g      >  configure_file

job_id=`llsubmit nmm_glob_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 5" >> RegressionTests.log
echo "Test 5"
echo "Test different decomposition" >> RegressionTests.log
echo "Test different decomposition"
(echo "NMMB, global case";echo;echo)>> RegressionTests.log
 echo "NMMB, global case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 5 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 5 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 5 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 5 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 5 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 5 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 5 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 \
         nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_glob/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 5 failed ")>> RegressionTests.log
  echo;echo " Test 5 failed "
  exit

fi

done

echo " Test 5 passed " >> RegressionTests.log
echo " Test 5 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 6
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NMM_THREAD

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_glob_ll.IN   | sed s:_TPN_:32:g         \
                     | sed s:_THRD_:2:g         \
                     | sed s:_GS_::g            \
                     | sed s:_RTPWD_:${RTPWD}:g \
                     | sed s:_SRCD_:${PATHTR}:g \
                     | sed s:_RUND_:${RUNDIR}:g >  nmm_glob_ll

cat nmm_glob_conf.IN | sed s:_INPES_:06:g       \
                     | sed s:_WTPG_:2:g         \
                     | sed s:_FCSTL_:48:g       \
                     | sed s:_NEMSI_:false:g    \
                     | sed s:_RSTRT_:false:g    \
                     | sed s:_gfsP_:false:g     \
                     | sed s:_RGS_:false:g      \
                     | sed s:_WGS_:false:g      >  configure_file

job_id=`llsubmit nmm_glob_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 6" >> RegressionTests.log
echo "Test 6"
echo "Test threading " >> RegressionTests.log
echo "Test threading "
(echo "NMMB, global case";echo;echo)>> RegressionTests.log
 echo "NMMB, global case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 6 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 6 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 6 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 6 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 6 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 6 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 6 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
         nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_glob/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 6 failed ")>> RegressionTests.log
  echo;echo " Test 6 failed "
  exit

fi

done

echo " Test 6 passed " >> RegressionTests.log
echo " Test 6 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 7
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NMM_gfsP

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_glob_ll.IN   | sed s:_TPN_:32:g         \
                     | sed s:_THRD_:1:g         \
                     | sed s:_GS_:#:g           \
                     | sed s:_RTPWD_:${RTPWD}:g \
                     | sed s:_SRCD_:${PATHTR}:g \
                     | sed s:_RUND_:${RUNDIR}:g  >  nmm_glob_ll

cat nmm_glob_conf.IN | sed s:_INPES_:06:g       \
                     | sed s:_WTPG_:2:g         \
                     | sed s:_FCSTL_:24:g       \
                     | sed s:_NEMSI_:false:g    \
                     | sed s:_RSTRT_:false:g    \
                     | sed s:_gfsP_:true:g      \
                     | sed s:_RGS_:false:g      \
                     | sed s:_WGS_:false:g      >  configure_file

job_id=`llsubmit nmm_glob_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 7" >> RegressionTests.log
echo "Test 7"
echo "Test GFS physics package in NMMB" >> RegressionTests.log
echo "Test GFS physics package in NMMB"
(echo "NMMB, global case";echo;echo)>> RegressionTests.log
 echo "NMMB, global case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 7 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 7 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 7 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 7 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 7 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 7 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 7 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024      \
         nmm_b_restart.012 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006   \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_restart_nemsio.012
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_gfsP_glob/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 7 failed ")>> RegressionTests.log
  echo;echo " Test 7 failed "
  exit

fi

done

echo " Test 7 passed " >> RegressionTests.log
echo " Test 7 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 8
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NMM_REG_CTL

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_reg_ll.IN   | sed s:_TPN_:32:g         \
                    | sed s:_THRD_:1:g         \
                    | sed s:_GS_:#:g           \
                    | sed s:_RTPWD_:${RTPWD}:g \
                    | sed s:_SRCD_:${PATHTR}:g \
                    | sed s:_RUND_:${RUNDIR}:g >  nmm_reg_ll

cat nmm_reg_conf.IN | sed s:_INPES_:06:g       \
                    | sed s:_WTPG_:2:g         \
                    | sed s:_FCSTL_:48:g       \
                    | sed s:_NEMSI_:false:g    \
                    | sed s:_RSTRT_:false:g    \
                    | sed s:_gfsP_:false:g     \
                    | sed s:_RGS_:false:g      \
                    | sed s:_WGS_:true:g       >  configure_file

job_id=`llsubmit nmm_reg_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 8" >> RegressionTests.log
echo "Test 8"
echo "Compare results with previous trunk version" >> RegressionTests.log
echo "Compare results with previous trunk version"
(echo "NMMB, regional case";echo;echo)>> RegressionTests.log
 echo "NMMB, regional case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 8 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 8 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 8 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 8 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 8 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 8 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 8 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
         nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_reg/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 8 failed ")>> RegressionTests.log
  echo;echo " Test 8 failed "
  exit

fi

done

echo " Test 8 passed " >> RegressionTests.log
echo " Test 8 passed "
export timing1=`grep total_tim $PATHRT/err | tail -1 | awk '{ print $5 }'`
export timingc=`cat ${RTPWD}/NMMB_reg/timing.txt`
(echo;echo " Original timing: " $timingc " , test1 timing: " $timing1)>> RegressionTests.log
 echo;echo " Original timing: " $timingc " , test1 timing: " $timing1

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 9
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NEM_REG_NEMSIO

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_reg_ll.IN   | sed s:_TPN_:32:g         \
                    | sed s:_THRD_:1:g         \
                    | sed s:_GS_:#:g           \
                    | sed s:_RTPWD_:${RTPWD}:g \
                    | sed s:_SRCD_:${PATHTR}:g \
                    | sed s:_RUND_:${RUNDIR}:g >  nmm_reg_ll

cat nmm_reg_conf.IN | sed s:_INPES_:06:g       \
                    | sed s:_WTPG_:2:g         \
                    | sed s:_FCSTL_:12:g       \
                    | sed s:_NEMSI_:true:g     \
                    | sed s:_RSTRT_:false:g    \
                    | sed s:_gfsP_:false:g     \
                    | sed s:_RGS_:false:g      \
                    | sed s:_WGS_:false:g      >  configure_file

job_id=`llsubmit nmm_reg_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 9" >> RegressionTests.log
echo "Test 9"
echo "Test NEMSIO as input file" >> RegressionTests.log
echo "Test NEMSIO as input file"
(echo "NMMB, regional case";echo;echo)>> RegressionTests.log
 echo "NMMB, regional case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 9 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 9 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 9 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 9 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 9 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 9 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 9 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.009 nmm_b_history.012    \
         nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                   \
         nmm_b_history_nemsio.009 nmm_b_history_nemsio.012
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_reg/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 9 failed ")>> RegressionTests.log
  echo;echo " Test 9 failed "
  exit

fi

done

echo " Test 9 passed " >> RegressionTests.log
echo " Test 9 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 10
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NMM_REG_RST

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_reg_ll.IN   | sed s:_TPN_:32:g         \
                    | sed s:_THRD_:1:g         \
                    | sed s:_GS_:#:g           \
                    | sed s:_RTPWD_:${RTPWD}:g \
                    | sed s:_SRCD_:${PATHTR}:g \
                    | sed s:_RUND_:${RUNDIR}:g >  nmm_reg_ll

cat nmm_reg_conf.IN | sed s:_INPES_:06:g       \
                    | sed s:_WTPG_:2:g         \
                    | sed s:_FCSTL_:48:g       \
                    | sed s:_NEMSI_:false:g    \
                    | sed s:_RSTRT_:true:g     \
                    | sed s:_gfsP_:false:g     \
                    | sed s:_RGS_:false:g      \
                    | sed s:_WGS_:false:g      >  configure_file

job_id=`llsubmit nmm_reg_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 10" >> RegressionTests.log
echo "Test 10"
echo "Test restart run" >> RegressionTests.log
echo "Test restart run"
(echo "NMMB, regional case";echo;echo)>> RegressionTests.log
 echo "NMMB, regional case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 10 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 10 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 10 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 10 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 10 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 10 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 10 results ...."

#
# NOTE! first hitory files differ because of PSFC, others should be identical:
#
for i in nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_reg/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 10 failed ")>> RegressionTests.log
  echo;echo " Test 10 failed "
  exit

fi

done

echo " Test 10 passed " >> RegressionTests.log
echo " Test 10 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 11
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NMM_REG_RST_NIO

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_reg_ll.IN   | sed s:_TPN_:32:g         \
                    | sed s:_THRD_:1:g         \
                    | sed s:_GS_:#:g           \
                    | sed s:_RTPWD_:${RTPWD}:g \
                    | sed s:_SRCD_:${PATHTR}:g \
                    | sed s:_RUND_:${RUNDIR}:g >  nmm_reg_ll

cat nmm_reg_conf.IN | sed s:_INPES_:06:g       \
                    | sed s:_WTPG_:2:g         \
                    | sed s:_FCSTL_:48:g       \
                    | sed s:_NEMSI_:true:g     \
                    | sed s:_RSTRT_:true:g     \
                    | sed s:_gfsP_:false:g     \
                    | sed s:_RGS_:false:g      \
                    | sed s:_WGS_:false:g      >  configure_file

job_id=`llsubmit nmm_reg_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 11" >> RegressionTests.log
echo "Test 11"
echo "Test restart run with NEMSIO file " >> RegressionTests.log
echo "Test restart run with NEMSIO file "
(echo "NMMB, regional case";echo;echo)    >> RegressionTests.log
 echo "NMMB, regional case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 11 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 11 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 11 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 11 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 11 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 11 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 11 results ...."

#
# NOTE! first hitory files differ because of PSFC, others should be identical:
#
for i in nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_reg/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 11 failed ")>> RegressionTests.log
  echo;echo " Test 11 failed "
  exit

fi

done

echo " Test 11 passed " >> RegressionTests.log
echo " Test 11 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 12
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NMM_REG_DECOMP

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_reg_ll.IN   | sed s:_TPN_:16:g         \
                    | sed s:_THRD_:1:g         \
                    | sed s:_GS_::g            \
                    | sed s:_RTPWD_:${RTPWD}:g \
                    | sed s:_SRCD_:${PATHTR}:g \
                    | sed s:_RUND_:${RUNDIR}:g >  nmm_reg_ll

cat nmm_reg_conf.IN | sed s:_INPES_:03:g       \
                    | sed s:_WTPG_:1:g         \
                    | sed s:_FCSTL_:12:g       \
                    | sed s:_NEMSI_:false:g    \
                    | sed s:_RSTRT_:false:g    \
                    | sed s:_gfsP_:false:g     \
                    | sed s:_RGS_:true:g       \
                    | sed s:_WGS_:false:g      >  configure_file

job_id=`llsubmit nmm_reg_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 12" >> RegressionTests.log
echo "Test 12"
echo "Test different decomposition" >> RegressionTests.log
echo "Test different decomposition"
(echo "NMMB, regional case";echo;echo)>> RegressionTests.log
 echo "NMMB, regional case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 12 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 12 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 12 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 12 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 12 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 12 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 12 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.009 nmm_b_history.012

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_reg/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 12 failed ")>> RegressionTests.log
  echo;echo " Test 12 failed "
  exit

fi

done

echo " Test 12 passed " >> RegressionTests.log
echo " Test 12 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 13
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NMM_REG_THREAD

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_reg_ll.IN   | sed s:_TPN_:32:g         \
                    | sed s:_THRD_:2:g         \
                    | sed s:_GS_::g            \
                    | sed s:_RTPWD_:${RTPWD}:g \
                    | sed s:_SRCD_:${PATHTR}:g \
                    | sed s:_RUND_:${RUNDIR}:g >  nmm_reg_ll

cat nmm_reg_conf.IN | sed s:_INPES_:06:g       \
                    | sed s:_WTPG_:2:g         \
                    | sed s:_FCSTL_:48:g       \
                    | sed s:_NEMSI_:false:g    \
                    | sed s:_RSTRT_:false:g    \
                    | sed s:_gfsP_:false:g     \
                    | sed s:_RGS_:false:g      \
                    | sed s:_WGS_:false:g      >  configure_file

job_id=`llsubmit nmm_reg_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 13" >> RegressionTests.log
echo "Test 13"
echo "Test threading " >> RegressionTests.log
echo "Test threading "
(echo "NMMB, regional case";echo;echo)>> RegressionTests.log
 echo "NMMB, regional case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 13 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 13 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 13 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 13 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 13 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 13 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 13 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
         nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_reg/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 13 failed ")>> RegressionTests.log
  echo;echo " Test 13 failed "
  exit

fi

done

echo " Test 13 passed " >> RegressionTests.log
echo " Test 13 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 14
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/NMM_REG_gfsP

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_reg_ll.IN   | sed s:_TPN_:32:g         \
                    | sed s:_THRD_:1:g         \
                    | sed s:_GS_:#:g           \
                    | sed s:_RTPWD_:${RTPWD}:g \
                    | sed s:_SRCD_:${PATHTR}:g \
                    | sed s:_RUND_:${RUNDIR}:g >  nmm_reg_ll

cat nmm_reg_conf.IN | sed s:_INPES_:06:g       \
                    | sed s:_WTPG_:2:g         \
                    | sed s:_FCSTL_:24:g       \
                    | sed s:_NEMSI_:false:g    \
                    | sed s:_RSTRT_:false:g    \
                    | sed s:_gfsP_:true:g      \
                    | sed s:_RGS_:false:g      \
                    | sed s:_WGS_:false:g      >  configure_file

job_id=`llsubmit nmm_reg_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 14" >> RegressionTests.log
echo "Test 14"
echo "Test GFS physics package in NMMB" >> RegressionTests.log
echo "Test GFS physics package in NMMB"
(echo "NMMB, regional case";echo;echo)>> RegressionTests.log
 echo "NMMB, regional case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 14 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 14 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 14 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 14 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 14 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 14 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 14 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024     \
         nmm_b_restart.012 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006  \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_restart_nemsio.012
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_gfsP_reg/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 14 failed ")>> RegressionTests.log
  echo;echo " Test 14 failed "
  exit

fi

done

echo " Test 14 passed " >> RegressionTests.log
echo " Test 14 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Clean and compile GFS core
####################################################################################################
echo "Preparing GFS core for regression tests" >> RegressionTests.log
echo "Preparing GFS core for regression tests"
printf %s "Compiling GFS core (this will take ~10 minutes)......."
cd ${PATHTR}/ush

./clean_stub.sh            >> ${PATHRT}/Compile.log 2>&1
./clean.sh                 >> ${PATHRT}/Compile.log 2>&1
./compile_configure.sh gfs >> ${PATHRT}/Compile.log 2>&1
./compile.sh               >> ${PATHRT}/Compile.log 2>&1

if [ -f ../exe/GFS_NEMS.x ] ; then
  echo "   GFS core Compiled";echo;echo
else
  echo "   GFS core is NOT compiled" >> ${PATHRT}/RegressionTests.log
  echo "   GFS core is NOT compiled"
  exit
fi

cd $PATHRT

####################################################################################################
# Submit test 15
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/GFS_32
mkdir -p ${RUNDIR}
cp gfs_configfile ${RUNDIR}/configure_file
cp ${RTPWD}/GFS/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS/sfcanl.2009072400 ${RUNDIR}/.
rm -f err out configure_file nmm_glob_ll nmm_reg_ll gfs_fcst_run  gfs_ll

cat gfs_ll.IN       | sed s:_TASKS_:32:g         \
                    | sed s:_THRDS_:1:g          >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:32:g         \
                    | sed s:_PE1_:30:g           \
                    | sed s:_WPG_:2:g            \
                    | sed s:_THRDS_:1:g          \
                    | sed s:_NSOUT_:0:g          \
                    | sed s:_QUILTING_:.true.:g  \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_NDAYS_:2:g          >  gfs_fcst_run

job_id=`llsubmit gfs_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 15" >> RegressionTests.log
echo "Test 15"
echo "Compare results with previous trunk version" >> RegressionTests.log
echo "Compare results with previous trunk version"
(echo "GFS, 32 proc, 1 thread";echo;echo)>> RegressionTests.log
 echo "GFS, 32 proc, 1 thread";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST 15 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 15 is running,            ID: " $job_id " Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST 15 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 15 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 15 is finished,           ID: " $job_id " Status: " $status
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 15 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 15 results ...."

for i in sigf03 sigf06 sigf09 sigf12 sigf24 sigf48 \
         sfcf03 sfcf06 sigf09 sfcf12 sfcf24 sfcf48 \
         flxf03 flxf06 sigf09 flxf12 flxf24 flxf48

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/GFS/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 15 failed ")>> RegressionTests.log
  echo;echo " Test 15 failed "
  exit

fi

done

echo " Test 15 passed " >> RegressionTests.log
echo " Test 15 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 16
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/GFS_60
mkdir -p ${RUNDIR}
cp gfs_configfile ${RUNDIR}/configure_file
cp ${RTPWD}/GFS/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS/sfcanl.2009072400 ${RUNDIR}/.
rm -f err out configure_file nmm_glob_ll nmm_reg_ll gfs_fcst_run  gfs_ll

cat gfs_ll.IN       | sed s:_TASKS_:60:g         \
                    | sed s:_THRDS_:1:g          >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:60:g         \
                    | sed s:_PE1_:58:g           \
                    | sed s:_WPG_:2:g            \
                    | sed s:_THRDS_:1:g          \
                    | sed s:_NSOUT_:0:g          \
                    | sed s:_QUILTING_:.true.:g  \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_NDAYS_:2:g          >  gfs_fcst_run

job_id=`llsubmit gfs_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 16" >> RegressionTests.log
echo "Test 16"
echo "Test different decomposition" >> RegressionTests.log
echo "Test different decomposition"
(echo "GFS, 60 proc, 1 thread";echo;echo)>> RegressionTests.log
 echo "GFS, 60 proc, 1 thread";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST 16 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 16 is running,            ID: " $job_id " Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST 16 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 16 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 16 is finished,           ID: " $job_id " Status: " $status
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 16 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 16 results ...."

for i in sigf03 sigf06 sigf09 sigf12 sigf24 sigf48 \
         sfcf03 sfcf06 sigf09 sfcf12 sfcf24 sfcf48 \
         flxf03 flxf06 sigf09 flxf12 flxf24 flxf48

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/GFS/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 16 failed ")>> RegressionTests.log
  echo;echo " Test 16 failed "
  exit

fi

done

echo " Test 16 passed " >> RegressionTests.log
echo " Test 16 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 17
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/GFS_16
mkdir -p ${RUNDIR}
cp gfs_configfile ${RUNDIR}/configure_file
cp ${RTPWD}/GFS/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS/sfcanl.2009072400 ${RUNDIR}/.
rm -f err out configure_file nmm_glob_ll nmm_reg_ll gfs_fcst_run  gfs_ll

cat gfs_ll.IN       | sed s:_TASKS_:16:g         \
                    | sed s:_THRDS_:2:g          >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:16:g         \
                    | sed s:_PE1_:14:g           \
                    | sed s:_WPG_:2:g            \
                    | sed s:_THRDS_:2:g          \
                    | sed s:_NSOUT_:0:g          \
                    | sed s:_QUILTING_:.true.:g  \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_NDAYS_:2:g          >  gfs_fcst_run

job_id=`llsubmit gfs_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 17" >> RegressionTests.log
echo "Test 17"
echo "Test threads" >> RegressionTests.log
echo "Test threads"
(echo "GFS, 16 proc, 2 threads";echo;echo)>> RegressionTests.log
 echo "GFS, 16 proc, 2 threads";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST 17 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 17 is running,            ID: " $job_id " Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST 17 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 17 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 17 is finished,           ID: " $job_id " Status: " $status
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 17 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 17 results ...."

for i in sigf03 sigf06 sigf09 sigf12 sigf24 sigf48 \
         sfcf03 sfcf06 sigf09 sfcf12 sfcf24 sfcf48 \
         flxf03 flxf06 sigf09 flxf12 flxf24 flxf48

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/GFS/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 17 failed ")>> RegressionTests.log
  echo;echo " Test 17 failed "
  exit

fi

done

echo " Test 17 passed " >> RegressionTests.log
echo " Test 17 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 18
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/GFS_01
mkdir -p ${RUNDIR}
cp gfs_configfile ${RUNDIR}/configure_file
cp ${RTPWD}/GFS/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS/sfcanl.2009072400 ${RUNDIR}/.
rm -f err out configure_file nmm_glob_ll nmm_reg_ll gfs_fcst_run  gfs_ll

cat gfs_ll.IN       | sed s:_TASKS_:1:g          \
                    | sed s:_THRDS_:1:g          >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:1:g          \
                    | sed s:_PE1_:1:g            \
                    | sed s:_WPG_:1:g            \
                    | sed s:_THRDS_:1:g          \
                    | sed s:_NSOUT_:0:g          \
                    | sed s:_QUILTING_:.false.:g  \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_NDAYS_:1:g          >  gfs_fcst_run

job_id=`llsubmit gfs_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 18" >> RegressionTests.log
echo "Test 18"
echo "Test single processor" >> RegressionTests.log
echo "Test single processor"
(echo "GFS, 1 proc, 1 threads, no quilting";echo;echo)>> RegressionTests.log
 echo "GFS, 1 proc, 1 threads, no quilting";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST 18 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 18 is running,            ID: " $job_id " Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST 18 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 18 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 18 is finished,           ID: " $job_id " Status: " $status
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 18 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 18 results ...."

for i in sigf03 sigf06 sigf09 sigf12 sigf24 \
         sfcf03 sfcf06 sigf09 sfcf12 sfcf24 \
         flxf03 flxf06 sigf09 flxf12 flxf24

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/GFS/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 18 failed ")>> RegressionTests.log
  echo;echo " Test 18 failed "
  exit

fi

done

echo " Test 18 passed " >> RegressionTests.log
echo " Test 18 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 19
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/GFS_01_NSOUT
mkdir -p ${RUNDIR}
cp gfs_configfile ${RUNDIR}/configure_file
cp ${RTPWD}/GFS/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS/sfcanl.2009072400 ${RUNDIR}/.
rm -f err out configure_file nmm_glob_ll nmm_reg_ll gfs_fcst_run  gfs_ll

cat gfs_ll.IN       | sed s:_TASKS_:1:g          \
                    | sed s:_THRDS_:1:g          >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:1:g          \
                    | sed s:_PE1_:1:g            \
                    | sed s:_WPG_:1:g            \
                    | sed s:_THRDS_:1:g          \
                    | sed s:_NSOUT_:1:g          \
                    | sed s:_QUILTING_:.false.:g  \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_NDAYS_:1:g          >  gfs_fcst_run

job_id=`llsubmit gfs_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 19" >> RegressionTests.log
echo "Test 19"
echo "Test single processor" >> RegressionTests.log
echo "Test single processor"
(echo "GFS, 1 proc, 1 threads,no quilting,nsout=1";echo;echo)>> RegressionTests.log
 echo "GFS, 1 proc, 1 threads,no quilting,nsout=1";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST 19 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 19 is running,            ID: " $job_id " Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST 19 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 19 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 19 is finished,           ID: " $job_id " Status: " $status
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 19 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 19 results ...."

for i in sigf03 sigf06 sigf09 sigf12 sigf24 \
         sfcf03 sfcf06 sigf09 sfcf12 sfcf24 \
         flxf03 flxf06 sigf09 flxf12 flxf24

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/GFS/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 19 failed ")>> RegressionTests.log
  echo;echo " Test 19 failed "
  exit

fi

done

echo " Test 19 passed " >> RegressionTests.log
echo " Test 19 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 20
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/GFS_16_NOQUILT_NSOUT
mkdir -p ${RUNDIR}
cp gfs_configfile ${RUNDIR}/configure_file
cp ${RTPWD}/GFS/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS/sfcanl.2009072400 ${RUNDIR}/.
rm -f err out configure_file nmm_glob_ll nmm_reg_ll gfs_fcst_run  gfs_ll

cat gfs_ll.IN       | sed s:_TASKS_:16:g         \
                    | sed s:_THRDS_:2:g          >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:16:g         \
                    | sed s:_PE1_:16:g           \
                    | sed s:_WPG_:1:g            \
                    | sed s:_THRDS_:2:g          \
                    | sed s:_NSOUT_:2:g          \
                    | sed s:_QUILTING_:.false.:g  \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_NDAYS_:2:g          >  gfs_fcst_run

job_id=`llsubmit gfs_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 20" >> RegressionTests.log
echo "Test 20"
echo "Test threads" >> RegressionTests.log
echo "Test threads"
(echo "GFS, 16 proc, 2 threads,no quilt, output every 2 time steps";echo;echo)>> RegressionTests.log
 echo "GFS, 16 proc, 2 threads,no quilt, output every 2 time steps";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST 20 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 20 is running,            ID: " $job_id " Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST 20 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 20 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 20 is finished,           ID: " $job_id " Status: " $status
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 20 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 20 results ...."

for i in sigf03 sigf06 sigf09 sigf12 sigf24 sigf48 \
         sfcf03 sfcf06 sigf09 sfcf12 sfcf24 sfcf48 \
         flxf03 flxf06 sigf09 flxf12 flxf24 flxf48

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/GFS/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 20 failed ")>> RegressionTests.log
  echo;echo " Test 20 failed "
  exit

fi

done

echo " Test 20 passed " >> RegressionTests.log
echo " Test 20 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 21
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/GFS_60_NOQUILT
mkdir -p ${RUNDIR}
cp gfs_configfile ${RUNDIR}/configure_file
cp ${RTPWD}/GFS/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS/sfcanl.2009072400 ${RUNDIR}/.
rm -f err out configure_file nmm_glob_ll nmm_reg_ll gfs_fcst_run  gfs_ll

cat gfs_ll.IN       | sed s:_TASKS_:60:g         \
                    | sed s:_THRDS_:1:g          >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:60:g         \
                    | sed s:_PE1_:60:g           \
                    | sed s:_WPG_:1:g            \
                    | sed s:_THRDS_:1:g          \
                    | sed s:_NSOUT_:0:g          \
                    | sed s:_QUILTING_:.false.:g  \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_NDAYS_:2:g          >  gfs_fcst_run

job_id=`llsubmit gfs_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 21" >> RegressionTests.log
echo "Test 21"
echo "Test different decomposition" >> RegressionTests.log
echo "Test different decomposition"
(echo "GFS, 60 proc, 1 thread, no quilt";echo;echo)>> RegressionTests.log
 echo "GFS, 60 proc, 1 thread, no quilt";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST 21 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 21 is running,            ID: " $job_id " Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST 21 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 21 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 21 is finished,           ID: " $job_id " Status: " $status
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 21 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 21 results ...."

for i in sigf03 sigf06 sigf09 sigf12 sigf24 sigf48 \
         sfcf03 sfcf06 sigf09 sfcf12 sfcf24 sfcf48 \
         flxf03 flxf06 sigf09 flxf12 flxf24 flxf48

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/GFS/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 21 failed ")>> RegressionTests.log
  echo;echo " Test 21 failed "
  exit

fi

done

echo " Test 21 passed " >> RegressionTests.log
echo " Test 21 passed "

sleep 4
clear;echo;echo
#
####################################################################################################
# Submit test 22
####################################################################################################

export RUNDIR=${RUNDIR_ROOT}/GFS_32_NOQUILT
mkdir -p ${RUNDIR}
cp gfs_configfile ${RUNDIR}/configure_file
cp ${RTPWD}/GFS/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS/sfcanl.2009072400 ${RUNDIR}/.
rm -f err out configure_file nmm_glob_ll nmm_reg_ll gfs_fcst_run  gfs_ll

cat gfs_ll.IN       | sed s:_TASKS_:32:g         \
                    | sed s:_THRDS_:1:g          >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:32:g         \
                    | sed s:_PE1_:32:g           \
                    | sed s:_WPG_:1:g            \
                    | sed s:_THRDS_:1:g          \
                    | sed s:_NSOUT_:4:g          \
                    | sed s:_QUILTING_:.false.:g  \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_NDAYS_:2:g          >  gfs_fcst_run

job_id=`llsubmit gfs_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 22" >> RegressionTests.log
echo "Test 22"
echo "Compare results with previous trunk version" >> RegressionTests.log
echo "Compare results with previous trunk version"
(echo "GFS, 32 proc, 1 thread, no quilt, output every 4 timestep";echo;echo)>> RegressionTests.log
 echo "GFS, 32 proc, 1 thread, no quilt, output every 4 timestep";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST 22 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 22 is running,            ID: " $job_id " Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST 22 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 22 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 22 is finished,           ID: " $job_id " Status: " $status
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 22 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 22 results ...."

for i in sigf03 sigf06 sigf09 sigf12 sigf24 sigf48 \
         sfcf03 sfcf06 sigf09 sfcf12 sfcf24 sfcf48 \
         flxf03 flxf06 sigf09 flxf12 flxf24 flxf48

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/GFS/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 22 failed ")>> RegressionTests.log
  echo;echo " Test 22 failed "
  exit

fi

done

echo " Test 22 passed " >> RegressionTests.log
echo " Test 22 passed "

sleep 4
clear;echo;echo
####################################################################################################
####################################################################################################
####################################################################################################


rm -f err out configure_file nmm_glob_ll nmm_reg_ll gfs_fcst_run  gfs_ll
cd ${PATHTR}/ush
./clean_stub.sh      > /dev/null 2>&1
./clean.sh > /dev/null 2>&1

date >> ${PATHRT}/RegressionTests.log

rm -rf ${RUNDIR_ROOT}
exit
