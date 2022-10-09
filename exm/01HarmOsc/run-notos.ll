#!/bin/bash
# @ job_name = rschr-01
# @ class = powiew
# @ account_no = G15-9
# @ error = a.err
# @ output = a.out
# @ environment = COPY_ALL
# @ wall_clock_limit = 40:00:00
# @ notification = never
# @ notify_user = $(user)@icm.edu.pl
# @ job_type = bluegene
# @ bg_size = 32
# @ bg_connection = PREFER_TORUS
# @ queue

# bg_connection = TORUS
# bg_connection = MESH




/bgsys/drivers/ppcfloor/bin/mpirun -np 32 \
-exe ${PWD}/../../bin/rschr.x             \
-mode SMP                                  \
-mapfile XYZT                             \
-cwd ${PWD}                               \
-args "-a ./param.inp  -b ./pot.def" > rschr.out 


#/bgsys/drivers/ppcfloor/bin/mpirun -np 128 \
#-exe ${PWD}/../../bin/rschr.x             \
#-mode VN                                  \
#-mapfile TXYZ                             \
#-cwd ${PWD}                               \
#-args "-a ./param.inp  -b ./pot.def" > rschr.out 


#
# Submit the job to the queuing system
# llsubmit run-notos.ll
#

