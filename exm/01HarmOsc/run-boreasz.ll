#!/bin/bash
# @ job_name = rschr-01
# @ output = $(job_name)_$(jobid).out
# @ error =  $(job_name)_$(jobid).err
# @ account_no = G15-9
# @ class = powiew
# @ node = 1
# @ tasks_per_node = 4
# @ wall_clock_limit = 08:00:00
# @ network.MPI = sn_all,not_shared,US,HIGH
# @ notification = never
# @ environment = COPY_ALL
# @ job_type = parallel
# @ queue

#
# Influential evironment parameters
#
# export MP_BUFFER_MEM=64M
# export MP_EAGER_LIMIT=2K
# export MP_CSS_INTERRUPT=no
# export MP_STATISTICS=print
# export MP_PRINTENV=yes 
export MP_CLOCK_SOURCE=ntp 



mpiexec ../../bin/rschr.x -vecscatter_alltoall -a ./param.inp -b ./pot.def > rschr.out

#
# Use MPI_Win_ interface. It works fine.
#
# mpiexec ../../bin/rschr.x -vecscatter_window -a ./param.inp -b ./pot.def > rschr.out




# 
# Following does not work on boreasz
# 
# mpiexec ../../bin/rschr.x -vecscatter_ssend -a ./param.inp -b ./pot.def > rschr.out
# mpiexec ../../bin/rschr.x -vecscatter_rsend -a ./param.inp -b ./pot.def > rschr.out
# mpiexec ../../bin/rschr.x -a ./param.inp -b ./pot.def > rschr.out


#
# Submit the job to the queuing system:
# llsubmit run-boreasz.ll
#
