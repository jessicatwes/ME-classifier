#!/bin/bash
#SBATCH --job-name=meclassifier # Job name
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jewe1055@colorado.edu # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=8 # Number of CPU
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --mem=8gb # Memory limit
#SBATCH --output=/Users/jewe1055/experiments/exp76/eofiles/%x_%j.out
#SBATCH --error=/Users/jewe1055/experiments/exp76/eofiles/%x_%j.err

### Displays the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Script directory is `pwd`
echo Input directory is ${INDIR}
echo Output SAM files directory is ${OUTDIR}
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes


### Assigns path variables
# INDIR directory for motif and tfit
INDIR_MOTIF='/Users/jewe1055/experiments/exp76/ME-classifier/datasample/hg38_v11_p1e-06'
INDIR_MUS='/Users/jewe1055/experiments/exp76/ME-classifier/datasample/tfit'
BASE_DIR='/Users/jewe1055/experiments/exp76/ME-classifier/output/'
SCRIPTS='/Users/jewe1055/experiments/exp76/ME-classifier'
ENV_DIR='/Users/jewe1055/jhub_venv/bin/activate'


######### Inputs tfit and motif call and outputs distance table ##########
mkdir -p $BASE_DIR
i=1
#for d in $NASCENT_DIR/*; do
#    if [[ -d $d ]]; then
#        if [[ -d $d/tfit/ ]]; then
#            for SRR_FILE in $(find $d/tfit/ -name '*.sorted_split_bidir_predictions.bed'); do
#                if [[ -d $d/tfit/ ]]; then  
#                    for MOTIF_FILE in `ls -1 $INDIR_MOTIF/*.bed`; do
#                        SRR="$(basename -- $SRR_FILE)"
#                        SRR="${SRR%.*}"
#                        MOTIF="$(basename -- $MOTIF_FILE)"
#                        MOTIF="${MOTIF%.*}"
#                        OUTDIR="$BASE_OUTDIR/$SRR/$MOTIF/"
#                        if ! test -f "$OUTDIR/raw_barcode_vals.csv"; then # folder not processed or not completed      
#                            rm "$OUTDIR/*"                        
#                            echo "================ calculating distances ================"                            
#                            echo SRR_FILE is $SRR_FILE
#                            echo MOTIF_FILE is $MOTIF_FILE
#                            echo Outdir is $OUTDIR
#                            mkdir -p $OUTDIR
#                            Rscript $SCRIPTS/dist_table_1motif.R $OUTDIR $MOTIF_FILE $SRR_FILE             
#                        fi    
#                    done
#                fi
#            done
#        fi
#    fi
#    ((i=i+1))
#done    



######### Org code for single mus folder ##########
for SRR_FILE in `ls -1 $INDIR_MUS/*-1_prelim_bidir_hits.bed`; do
  for MOTIF_FILE in `ls -1 $INDIR_MOTIF/*.bed`; do
    SRR="$(basename -- $SRR_FILE)"
    SRR="${SRR%.*}"
    MOTIF="$(basename -- $MOTIF_FILE)"
    MOTIF="${MOTIF%.*}"
    OUTDIR="$BASE_DIR$SRR/$MOTIF/"
    echo $OUTDIR
    if ! test -f "$OUTDIR/raw_barcode_vals.csv"; then # folder not processed or not completed
      rm "$OUTDIR/*"
      echo "================ calcuating distances ================"
      echo $SRR
      echo $MOTIF
      mkdir -p $OUTDIR
      Rscript $SCRIPTS/1_calculate_distances.R $OUTDIR $MOTIF_FILE $SRR_FILE
    fi
    ((i=i+1))
  done
done


######### Input distance table and output counts per distance ##########
source $ENV_DIR
python3 ${SCRIPTS}/2_count_distances.py $BASE_DIR $BASE_DIR



######### Inputs MD count table. Convert to z-score and call peak finder on both MD and z-score #########
python3 ${SCRIPTS}/3_zscore_pf_calc.py $BASE_DIR $BASE_DIR

