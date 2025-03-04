### simulate data 
source("data_processing_BMsim11Y.R")

### write stan files
source("BMsim11Y_stan.R")

### estimate models L = 1, 2, and 3
# cd ./AD
# sbatch bash2.sh
# squeue -u yxu2
if(0){
  #!/bin/bash
  #SBATCH --time=2-00:00:00
  #SBATCH --mem-per-cpu=5G
  #SBATCH --job-name=simBM1
  #SBATCH --array=1-300
  #SBATCH --output=./out_err/slurm_%a.out
  #SBATCH --error=./out_err/slurm_%a.err
  
  module load conda_R/4.3.x
  Rscript BMsim11Y.R $SLURM_ARRAY_TASK_ID 

  
  # cd ./AD
  # qsub bash1.sh
  
  #!/bin/bash
  #$ -N sim2 
  # job name
  #$ -cwd 
  # execute in the current directory
  #$ -l mem_free=5G,h_vmem=5G 
  #$ -l h_fsize=6G
  # 200 cores 768G RAM
  # node has >= mem_free RAM, kill job when uses > h_vmem
  # h_fsize file size limit
  #$ -R y
  # job reservation - yes; -r y: rerun if aborted without consistent exit state
  #$ -t 1-300
  # array job specify the range of tasks
  # assign to variable $SGE_TASK_ID in script
  #$ -M yxu143@jhu.edu
  # list mail address
  
  module load conda_R/4.0
  Rscript BMsim11Y.R
  # http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html
}

### create posterior summaries postsumm1,2,3.RData
# By p9.pdf, we see that L=2 has optimal and suboptimal solutions, with 
# average posteriors  -30792.54 and -30920.29, respectively. The suboptimal
# solution has alpha0 < 0.5, as in p10.pdf.
# Hence, the criteria of including a chain from a simulation replicate are
# 1. parameters converged, e.g. sigma^2 for all outcomes < 100
# 2. if L = 2, exclude suboptimal solution, e.g. alpha0 < 0.5
# 3. if both 1 and 2 satisfied, choose the chain (out of the two chains) with the larger posterior

source("BMsim_table_final.R")

### create summary plots
source("BMsim_plot.R")

