set -xe

mpicc -I../include -o allpair allpair.c ../util/mtest.c ../util/mtest_common.c ../util/mtest_single.c

export MPIR_CVAR_ODD_EVEN_CLIQUES=1
export MPIEXEC_TIMEOUT=30
# export MPITEST_VERBOSE=1

export MPIR_CVAR_CH4_NUM_VCIS=4

for i in {1..1} ; do
    echo ==== $i ====
    time mpirun -l -n 2 ./allpair || exit
done
