jobname = raw_input('Enter Job name:')
partition = raw_input('Enter partition (example cpu,knl,gpu) comma seperated: ')
memreq = raw_input('Enter Memory Required for your Program(in MB): ')
nodes = raw_input('Enter nodes required by your program: ')
tasks = raw_input('Enter tasks to be created( N tasks indicates your program will be executeed N times in parallel) : ')
taskspn = raw_input('Enter number of tasks to be allocated per node: ')
slurmScript = raw_input('Enter SLURM Batch Script output filename: ')

with open(slurmScript,"w") as f:
	f.write('#!/bin/bash\n\n')
	f.write('#SBATCH -p '+partition+'\n')
	f.write('#SBATCH -t 00:20:00\n')
	f.write('#SBATCH -m='+memreq+'\n')
	f.write('#SBATCH --nodes='+nodes+'\n')
	f.write('#SBATCH --tasks='+tasks+'\n')
	f.write('#SBATCH --ntasks-per-node='+taskspn+'\n')
	f.write('#SBATCH --job-name='+jobname+'\n')
	f.write('#SBATCH -o '+jobname+'_%j.o\n')
	f.write('\n\n')
	f.write('srun hostname')

f.close()


print('SLURM script created....Enter your code at srun line in '+slurmScript +' file\n\n')
print('Use ' + '\x1b[6;30;42m' + 'sbatch '+slurmScript + '\x1b[0m'+ ' to run in OpenHPC Cluster. Output files will be prepended with jobname')
