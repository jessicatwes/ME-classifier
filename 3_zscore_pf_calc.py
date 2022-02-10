import sys
import csv
from pf_lib import main
from pf_lib import input_parser
from pf_lib import find_peaks_and_or_valleys
from pf_lib import plot_raw_data
from pf_lib import plot_peaks_and_valleys
from pf_lib import select_method
from pf_lib import run_pf_plot
from pf_lib import draft_categorize_msd_file_into_tf_regulation
from pf_lib import zscore_conv

# set pathway
md_table_file = 'md_table.csv'
md_tfpf_file = 'md_tfpf.csv'
zscore_table_file = 'zscore_table.csv'
zscore_tfpf_file = 'zscore_tfpf.csv'

BASE_DIR = sys.argv[1]

# initialize zscore
all_zscores = []

# read csv file
with open(BASE_DIR+md_table_file) as csvfile:
    read_csv = csv.reader(csvfile, delimiter=',')
    header = header = next(read_csv)
    i = 0
    for md_row in read_csv:
        # get the data portion of the csv row and convert to integers
        md_data = [int(i) for i in md_row[2:]]
        tf_data = md_row[0]
        srr_data = md_row[1]
        zscores = zscore_conv(md_data)
        if zscores: # toss rows with 2 or less count values
            zscores = [tf_data, srr_data] + zscores #concatnate 2 lists
            all_zscores.append(zscores)

with open(BASE_DIR+zscore_table_file, 'w', encoding='UTF8', newline='') as f:
    write_csv = csv.writer(f)
    write_csv.writerow(header)
    for row in all_zscores:
        write_csv.writerow(row)

#### Run peak finder and tf regulation categories
# for md
draft_categorize_msd_file_into_tf_regulation(BASE_DIR+md_table_file, BASE_DIR+md_tfpf_file)
# for z-score
draft_categorize_msd_file_into_tf_regulation(BASE_DIR+zscore_table_file, BASE_DIR+zscore_tfpf_file)
