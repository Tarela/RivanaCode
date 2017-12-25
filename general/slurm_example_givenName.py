SRRs=['SRR1533850','SRR1533851','SRR1533852','SRR2085919','SRR2096437']
SPEs=['mm10','mm10','mm10','hg38','hg38']

for i in range(5):
    SRR = SRRs[i]
    species = SPEs[i]
    model="""#!/bin/bash
#SBATCH -n 8
#SBATCH --mem=100000
#SBATCH -t 18:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o %s.out
#SBATCH -e %s.err
#Run program

./downmap %s %s
"""%(SRR,SRR,SRR,species)
    outf=  open('%s.slurm'%(SRR),'w')
    outf.write(model)
    outf.close()

