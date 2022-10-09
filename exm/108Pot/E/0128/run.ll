# @ job_name = 108-E-128
# @ class = powiew
# @ account_no = G15-9
# @ error = a.err
# @ output = a.out
# @ environment = COPY_ALL
# @ wall_clock_limit = 40:00:00
# @ notification = never
# @ job_type = bluegene
# @ bg_size = 128
# @ bg_connection = PREFER_TORUS
# @ queue





/bgsys/drivers/ppcfloor/bin/mpirun -np 128                \
-exe /bgdata/home/users/romz/rschr-5.4/bin/rschr.x       \
-mode SMP                                                \
-mapfile TXYZ                                           \
-cwd /bgdata/home/users/romz/rschr-5.4/exm/108Pot/E/0128    \
-args "../param.inp  ../pot.def" > rschr.out 


