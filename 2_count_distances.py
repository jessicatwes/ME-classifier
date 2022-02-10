""" reads the distance data output from dist_table_1motif.R 
and counts the calls ate each distance """
import os
import sys
import csv


def count_distances(input_file_name):
    input_file_name = open(input_file_name)
    input_csv = csv.reader(input_file_name, delimiter=',')
    line = None
    counts = [0 for _ in range(2999)]
    next(input_csv) # remove header
    for row in input_csv:
        counts[1499+int(row[1])] += int(row[2])
    return counts
    print(counts)
    


def write_distance_counts(output_path, all_counts):
    output_file_name =  os.path.join(output_path, 'md_table.csv')
    if os.path.exists(output_file_name):
        os.remove(output_file_name)
    output_file = open(output_file_name, 'w')
    header = ['ID', 'SRZ'] + [str(i) for i in range(-1499,1499)]
    output_file.write(','.join(header) + '\n')
    for counts in all_counts:
        print(counts)
        output_file.write(','.join(counts) + '\n')
    output_file.close()


def count_all_distances(input_path):
    all_counts = []
    for srr_folder in os.listdir(input_path):
        srr_path = os.path.join(input_path, srr_folder)
        if not os.path.isdir(srr_path):
            continue
        for motif_folder in os.listdir(srr_path):
            motif_path = os.path.join(srr_path, motif_folder)
            if not os.path.isdir(motif_path):
                continue
            print('counting calls for', srr_folder, 'and', motif_folder)
            input_file_name = os.path.join(motif_path, 'raw_barcode_vals.csv')
            if not os.path.isfile(input_file_name):
                continue
            counts = count_distances(input_file_name)
            all_counts.append([motif_folder, srr_folder] + [str(i) for i in counts])
    return all_counts


if __name__ == '__main__':
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    all_counts = count_all_distances(input_path)
    write_distance_counts(output_path, all_counts)
