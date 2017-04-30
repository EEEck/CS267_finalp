$ cat scaling.sh
#!/bin/bash
#SBATCH -N 1
#SBATCH -p regular
#SBATCH -C knl,quad,cache
#SBATCH -J first_test
#SBATCH -t 01:59:20

#this is the KNL script
echo ncores,nht,arch,time > threadscale_knl.csv
for nc in 1 2 4 8 16 32 64; do
for nht in 1 2 4; do

export OMP_NUM_THREADS=$(( ${nc} * ${nht} ))
export OMP_PLACES=cores"(${nc})"
export OMP_PROC_BIND=spread
#module load memkind 
echo $nc $nht >> output
srun -n 1 -c 272 --cpu_bind=cores ./fast >> output
done
done

