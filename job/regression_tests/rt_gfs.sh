#!/bin/ksh
set -ua

export GEFS_ENSEMBLE=${GEFS_ENSEMBLE:-0}
echo "GEFS_ENSEMBLE=" $GEFS_ENSEMBLE

mkdir -p ${RUNDIR}

if [ $GEFS_ENSEMBLE = 0 ] ; then

####################################################################################################
# For the stand alone GFS regression tests.
####################################################################################################

####################################################################################################
# Make configure and run files
####################################################################################################

## determine GOCART and TRACER from gocart_aerosol and passive_tracer 
export gocart_aerosol=${gocart_aerosol:-NO}
export passive_tracer=${passive_tracer:-NO}
if [ $gocart_aerosol = 'YES' ]; then
 export GOCART=1 
else
 export GOCART=0 
fi
if  [ $passive_tracer = 'YES' ]; then
 export TRACER=.true.
else
 export TRACER=.false.
fi
##


cat gfs_fcst_run.IN | sed s:_TASKS_:${TASKS}:g   \
                    | sed s:_PE1_:${PE1}:g       \
                    | sed s:_WPG_:${WTPG}:g      \
                    | sed s:_WRTGP_:${WRTGP}:g   \
                    | sed s:_THRDS_:${THRD}:g    \
                    | sed s:_NSOUT_:${NSOUT}:g   \
                    | sed s:_QUILT_:${QUILT}:g   \
                    | sed s:_IAER_:${IAER}:g       \
                    | sed s:_wavecoef_:${wavecoef}:g     \
                    | sed s:_wavegrid_:${wavegrid}:g     \
                    | sed s:_lm_:${lm}:g         \
                    | sed s:_lsoil_:${lsoil}:g   \
                    | sed s:_MEMBER_NAMES_:${MEMBER_NAMES}:g   \
                    | sed s:_CP2_:${CP2}:g       \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_FDFI_:${FDFI}:g \
                    | sed s:_FHRES_:${FHRES}:g \
                    | sed s:_REDUCEDGRID_:${REDUCEDGRID}:g \
                    | sed s:_ADIAB_:${ADIAB}:g \
                    | sed s:_NSTFCST_:${NST_FCST}:g \
                    | sed s:_GOCART_:${GOCART}:g \
                    | sed s:_TRACER_:${TRACER}:g \
                    | sed s:_SFCPRESSID_:${SFCPRESS_ID}:g \
                    | sed s:_THERMODYNID_:${THERMODYN_ID}:g \
                    | sed s:_IDVC_:${IDVC}:g \
                    | sed s:_NDSLFV_:${NDSLFV}:g \
                    | sed s:_SPECTRALLOOP_:${SPECTRALLOOP}:g \
                    | sed s:_NDAYS_:${NDAYS}:g   >  gfs_fcst_run

####################################################################################################
# Copy init files
####################################################################################################

cp atmos.configure_gfs ${RUNDIR}/atmos.configure
cp MAPL.rc ${RUNDIR}/MAPL.rc
cp Chem_Registry.rc ${RUNDIR}/Chem_Registry.rc

if [ $GOCART = 1 ] ; then
 export EXTDIR=/global/save/wx23lu/NEMS/fix
 export RCSDIR=/global/save/wx23lu/NEMS/Chem_Registry
 cp ${RCSDIR}/*.rc ${RUNDIR}/.
 cp -r  ${EXTDIR}/ExtData ${RUNDIR}/.
fi

if [ $IDVC = 2 ] ; then
  cp ${RTPWD}/GFS_DFI_REDUCEDGRID_HYB/gfsanl.2010010100 ${RUNDIR}/.
  cp ${RTPWD}/GFS_DFI_REDUCEDGRID_HYB/sfcanl.2010010100 ${RUNDIR}/.

#  cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/t62hyb/gfsanl.2010010100 ${RUNDIR}/.
#  cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/t62hyb/sfnanl.2010010100 ${RUNDIR}/sfcanl.2010010100


elif [ $IDVC = 3 ] ; then
  cp ${RTPWD}/GFS_NODFI/gfsanl.2010010100 ${RUNDIR}/.
  cp ${RTPWD}/GFS_NODFI/sfcanl.2010010100 ${RUNDIR}/.

#  cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/t62/gfsanl.2010010100 ${RUNDIR}/.
#  cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/t62/sfnanl.2010010100 ${RUNDIR}/sfcanl.2010010100
fi

else


####################################################################################################
# For the concurrency ensemble GEFS regression test.
####################################################################################################

cd $PATHRT

echo 'RUNDIR=' $RUNDIR

 cp ${RTPWD}/GEFS_data_2008082500/* $RUNDIR
# cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/GEFS_data_2008082500/gfsanl* $RUNDIR
# cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/GEFS_data_2008082500/sfcanl* $RUNDIR
cp $PATHRT/gfs_configfile_190 $RUNDIR/configure_file


cat gfs_fcst_run_GEFS.IN \
                    | sed s:_SRCDIR_:${PATHTR}:g \
                    | sed s:_NDSLFV_:${NDSLFV}:g \
                    | sed s:_RUNDIR_:${RUNDIR}:g > gfs_fcst_run

cp Chem_Registry.rc ${RUNDIR}/Chem_Registry.rc
cp atmos.configure_gfs ${RUNDIR}/atmos.configure

fi

####################################################################################################
# Submit test
####################################################################################################

JBNME=NEMS_RT_${TEST_NR}_$$

cat gfs_ll.IN       | sed s:_JBNME_:${JBNME}:g   \
                    | sed s:_CLASS_:${CLASS}:g   \
                    | sed s:_GROUP_:${GROUP}:g   \
                    | sed s:_ACCNR_:${ACCNR}:g   \
                    | sed s:_WLCLK_:${WLCLK}:g   \
                    | sed s:_TASKS_:${TASKS}:g   \
                    | sed s:_THRDS_:${THRD}:g    >  gfs_ll

llsubmit gfs_ll 2>&1 | grep submitted > /dev/null

echo "Test ${TEST_NR}" >> RegressionTests.log
echo "Test ${TEST_NR}"
echo ${TEST_DESCR} >> RegressionTests.log
echo ${TEST_DESCR}
(echo "GFS, ${TASKS} proc, ${THRD} thread";echo;echo)>> RegressionTests.log
 echo "GFS, ${TASKS} proc, ${THRD} thread";echo;echo

# wait for the job to enter the queue
job_running=0
until [ $job_running -eq 1 ]
do
echo "TEST is waiting to enter the queue"
job_running=`llq -u ${LOGIN} -f %st %jn | grep ${JBNME} | wc -l`;sleep 5
done

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

export status=`llq -u ${LOGIN} -f %st %jn | grep ${JBNME} | awk '{ print $1}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST ${TEST_NR} is waiting in a queue, Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST ${TEST_NR} is running,            Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST ${TEST_NR} is ready to run,       Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status
else                          echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status
fi

sleep 60

job_running=`llq -u ${LOGIN} -f %st %jn | grep ${JBNME} | wc -l`
  (( n=n+1 ))
done

####################################################################################################
# Check results
####################################################################################################

(echo;echo;echo "Checking test ${TEST_NR} results ....")>> RegressionTests.log
 echo;echo;echo "Checking test ${TEST_NR} results ...."

#
     if [ ${CREATE_BASELINE} = false ]; then
#
# --- regression test comparison ----
#

for i in ${LIST_FILES}

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/${CNTL_DIR}/$i ${RUNDIR}/$i | wc -l`

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
    cp ${RUNDIR}/${i} /stmp/${LOGIN}/REGRESSION_TEST/${CNTL_DIR}/${i}
  else
    echo "Missing " ${RUNDIR}/$i " output file"
    echo;echo " Set ${TEST_NR} failed "
    exit 2
  fi
done

# ---
     fi
# ---

echo " Test ${TEST_NR} passed " >> RegressionTests.log
echo " Test ${TEST_NR} passed "

sleep 4
clear;echo;echo

####################################################################################################
# End test
####################################################################################################

exit 0
