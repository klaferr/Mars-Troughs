#!/usr/bin/bash
#SBATCH -A bramsona
#SBATCH -t 0-48:00:00
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH --mem-per-cpu=10G

python /home/klaferri/Desktop/Research/Mars-Troughs/scripts/figures_mcmc.py -objpath '/home/klaferri/Desktop/Research/TMP3/Slice6/Obliquity/TMP3/obj/' -plotdir '/home/klaferri/Desktop/Research/TMP3/Slice6/' -nmodels 50 -stepEnsemble 1


