# ME-classifier
## Calculating the co-occurence and distance between TF DNA binding motifs and enhancer RNA (eRNA) calls

This github repository has 3 main components:
1) 1_calculate_distances.R: This R script uses GRange to calculate the distance between two files; TF motif calls (eg HOCOMOCO) and bidir calls (from <a href="https://github.com/Dowell-Lab/Tfit">Tfit</a>). The output will include a metadata table, a motif GRange bed, a eRNA GRange bed, and a distance table for a single motif and eRNA sample.
2) 2_count_distances.py: This python script takes the distance table output from step 1 and compiles them together into a single motif frequency table.
3) 3_zscore_pf_calc.py: This python script will take the motif frequency table (md_table) and converts the table into a secondary table that uses z-score to do a false discovery rate correction to reduce background. The motif frequncy and z-score frequency table will be process through a peak finding algorithm to find regions of enriched signal and then classify each TF based on number of peaks.

To run all the script, a run.sh is a bash script that will run the 3 scripts above to produce frequency table and classification output.

