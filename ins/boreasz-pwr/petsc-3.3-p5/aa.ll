#@ job_name = petsc-configure
#@ output = petsc-configure.out
#@ error =  petsc-configure.err
#@ account_no = G15-9
#@ class = powiew
#@ node =  1
#@ tasks_per_node = 1
#@ wall_clock_limit = 01:00:00
#@ network.MPI = sn_all,not_shared,US,HIGH
#@ notification = never
#@ environment = COPY_ALL
#@ job_type = parallel
#@ queue


mpiexec ./conftest-aa-aa


