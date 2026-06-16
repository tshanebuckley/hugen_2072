module load plink/1.90b7.7

# edit this to point to your directory 
OUT_FOLDER=$HOME/p8

CLASS_FOLDER=/ix1/hugen2072-2026s/p8/

plink --bfile $CLASS_FOLDER/p8_study1 \
--score  $CLASS_FOLDER/p8_pgs002975_scores.txt header sum \
--out $OUT_FOLDER/p8_study1_scored 