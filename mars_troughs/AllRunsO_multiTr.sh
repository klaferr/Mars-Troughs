#!/usr/bin/bash
#SBATCH -A bramsona
#SBATCH -t 0-72:00:00
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH -n 64

for ia in $(seq 1 4)
do
for ir in $(seq 2 5)
do
 srun --exclusive -c4 -n1  python ../scripts/tmp_inversion_multiTr.py -acc $ia -retr $ir -steps 500000 -thin_by 1 -dir '/home/klaferri/Desktop/Research/TMP3/Slice6/' -tmp 3 -trslice 6 -retreat "fit" -data 'Obliquity' &
done

done

wait
