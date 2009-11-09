#!/bin/ksh

################################################################################
# Script name:         prepare
# Script description:  Sets up Load Leveler commands for serial batch jobs
#                      in /scripts needed to run the WRF-NMM launcher
#
# Authors:    Ed Colon
################################################################################
#
#  Used to set the following Load Leveler commands in serial batch jobs in the
#  /scripts subdirectory
#
#
#-- *** IMPORTANT *** IMPORTANT *** IMPORTANT *** IMPORTANT ***
#
#   Make sure that the settings used here are consistent with those set in prelaunch!!!
#
# ------ NCEP IBM SP ----------------
#=======================================================================

    cd ../src/main && gmake -f Makefile_main clean
    rm -f Makefile_main
    rm -f Makefile_stub
    rm -f atmos/makefile
    rm -f atmos/makefile_stub 
echo 'NEMS cleaned'

exit
