# @ job_name = config-petsc
# @ class = powiew
# @ account_no = G15-9
# @ error = config-petsc.err
# @ output = config-petsc.out
# @ environment = COPY_ALL
# @ wall_clock_limit = 01:00:00
# @ notification = error
# @ notify_user = $(user)@icm.edu.pl
# @ job_type = bluegene
# @ bg_size = 1 
# @ bg_connection = MESH
# @ queue

/bgsys/drivers/ppcfloor/bin/mpirun -np 1 ./conftest-aa-aa




