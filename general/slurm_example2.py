inf = open('meta_cmd.sh')
for line in inf:
    SRR = line.split()[1]
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
"""%(SRR,SRR,SRR,'hg38')
    outf=  open('%s.slurm'%(SRR),'w')
    outf.write(model)
    outf.close()

