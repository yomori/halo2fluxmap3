#!/bin/bash

for i in {1..10}
do
echo '#!/bin/bash' > submit_job
echo "#SBATCH -t 01:00:00">> submit_job
echo "#SBATCH --partition=chihway" >> submit_job
echo "#SBATCH --account=pi-chihway" >> submit_job
echo "#SBATCH --exclusive" >> submit_job
echo "#SBATCH --nodes=1" >> submit_job
echo "#SBATCH --job-name=MDPL_CIB${i}" >> submit_job
echo "#SBATCH --exclude=midway2-0282,midway2-0283,midway2-0284" >> submit_job
#echo "#SBATCH --mem 58000" >> submit_job
echo "cd $PWD" >> submit_job
echo "source /project2/chihway/setup/setup_intel.sh" >> submit_job
echo python build_cibmap_universemachine.py 4096 $((i)) 50 >> submit_job
sbatch submit_job
done

