# Import libraries/packages
import os
import argparse
import csv
import math
import statistics
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from statistics import *
from itertools import groupby


def main(file_name, selected_methods):
    '''
    Main function call select_method() that contains input of
    argparse and run the peak finder algorithm.

    Input parameters
    ---------------
    filename           : str
                        Path to motif distribution csv file
                        that will be our raw input data
    selected_methods   : list
                        list of methods that user wants for peak finding.

    Output
    ---------------
    return None
    '''

    # get the input from argparse
    my_args = input_parser(selected_methods)
    motif_distribution = my_args.motif_distribution
    exponential = my_args.exponential
    self_power = my_args.self_power
    threshold = my_args.threshold
    fast_fourier = my_args.fast_fourier

    # select_method() function trigger to run selected methods from
    # above parser
    select_method(file_name, motif_distribution, exponential, self_power,threshold, fast_fourier)
    return None


def select_method(file_name, motif_distribution, exponential, self_power, threshold, fast_fourier):
    '''
    This function can let the user choose which method that want to
    use for finding peak

    Input parameters
    ---------------
    filename           : str
                         This parameter contain the naming of input file
    motif_distribution : str
                         This parameter contain whether use
                         motif_distribution method by argparse
    exponential        : str
                         This parameter contain whether use exponential
                         method by argparse
    self_power         : str
                         This parameter contain whether use self_power
                         method by argparse
    threshold          : str
                         This parameter contain whether use threshold method
                         by argparse
    fast_fourier       : str
                         This parameter contain whether to use fast fourier transform method
                         by argparse

    Output returns
    ---------------
    raw_data           : list
                         Store the motif distribution raw data as a list
    '''
    reduce_size = 10
    smooth_size = 10

    # Convert the original format to a 2D list
    raw_data = _preprocess_mds_file(file_name)
    raw_data = list(map(list, zip(*raw_data)))

    for i in range(len(raw_data)):
        count = 0
        for raw_row in raw_data[i]:
            # Print which method is based on parser
            if count == 0 and motif_distribution == "True":
                print("Original Method")
                run_pf_plot(raw_row, reduce_size, smooth_size)
            elif count == 1 and exponential == "True":
                print("exp method")
                run_pf_plot(raw_row, reduce_size, smooth_size)
            elif count == 2 and self_power == "True":
                print("self-power method")
                run_pf_plot(raw_row, reduce_size, smooth_size)
            elif count == 3 and threshold == "True":
                print("threshold method")
                run_pf_plot(raw_row, reduce_size, smooth_size)
            elif count == 4 and fast_fourier == "True":
                print("fast fourier transform method")
                run_pf_plot(raw_row, reduce_size, smooth_size)
            count = count + 1
    return raw_data


def run_pf_plot(raw_row, reduce_size, smooth_size):
    '''
    This function calls the find_peaks_and_or_valleys() to find the
    peak/valley based on statistical threshold of 1.5 standard deviation
    from the mean. It prints the output data after running the function.

    Parameters
    ---------------
    raw_row            : list
                         This parameter contain the raw data
    reduce_size        : int
                         This parameter contain the size of reducing raw data
    smooth_size        : int
                         This parameter contain the size of smoothing raw
                         data

    Returns
    ---------------
    Prints the data in the jupyter notebook.
    '''
    reduced_row, peaks, valleys, peak_properties, valley_properties = find_peaks_and_or_valleys(
        raw_row, reduce_size, smooth_size
    )

    peaks_bp = []
    for old_loc in peaks:
        new_loc = (old_loc + int(smooth_size/2))*reduce_size
        peaks_bp.append(new_loc)

    valleys_bp = []
    for old_loc in valleys:
        new_loc = (old_loc + int(smooth_size/2))*reduce_size
        valleys_bp.append(new_loc)

    plot_raw_data(raw_row, peaks_bp, valleys_bp)
    plot_peaks_and_valleys(reduced_row, peaks, valleys)
    print("---------------------------------------------------")
    return raw_row, reduced_row, peaks, valleys, peaks_bp, valleys_bp


def find_peaks_and_or_valleys(raw_row, reduce_size, smooth_size):
    '''
    Take raw data and call _determine_scans() to scan for peaks and valley
    based on detection of continuous signal. Calculate the statistic or
    reduced data to define outliers

    Input parameters
    ---------------
    raw_row           : list
                        motif distribution data as a list
    reduce_size       : int
                        size for bucket in reduction
    smooth_size       : int
                        size for smoothing function

    Output return
    ---------------
    reduced_row       : list
                        processed list after _reduce_data_for_peak_detection()
    peaks             : list
                        list of called peaks
    valleys           : list
                        list of called valleys
    peak_properties   : list
                        list of peak properties including the number
    valley_properties : list
    '''
    reduced_row = _reduce_data_for_peak_detection(raw_row, reduce_size, smooth_size)
    peak_scan, valley_scan = _determine_scans(reduced_row)

    peaks = []
    valleys = []
    peak_properties = []
    valley_properties = []
    ymax = max(reduced_row)
    mean_y = mean(reduced_row)
    median_y = median(reduced_row)
    stdev_y = stdev(reduced_row)
    stdev_15 = stdev_y * 1.5
    stdev_20 = stdev_y * 2
    y = np.array(reduced_row)
    # print("StdDev: ", round(stdev_y,3))
    # invert data before finding peaks
    if peak_scan == True:
        '''
        Use Scipy find_peaks library with the set parameters (these
        can be alter for optimization). Returns indexes of y that are
        detected as peaks and the corresponding properties.
        '''
        peaks, peak_properties = find_peaks(
            y, height=mean_y+stdev_15, distance=20,
            prominence=4, width=3
        )
        # print("Peak properties: ", peak_properties)

    if valley_scan == True:
        '''
        Use Scipy find_peaks library with the set parameters
        (these can be alter for optimization). Returns indexes
        of y that are detected as valleys and the corresponding
        properties.
        '''
        y_inverse = ymax - y
        mean_y_inverse = ymax - mean_y
        valleys, valley_properties = find_peaks(
            y_inverse, height=mean_y_inverse+stdev_15, distance=20,
            prominence=4, width=3)
        # print("Valley properties: ", valley_properties)

    # TODO calculate peak and valley positions as diatance from center not
    # the reduced data location
    return reduced_row, peaks, valleys, peak_properties, valley_properties


def input_parser(selected_methods):
    '''
    This function parse user's input to select which method to implement;
    motif distribution without manipulation or motif distribution after
    exponential, self-power, or threshold. Function is called in main().

    Input parameters
    ---------------
    None

    Output
    ---------------
    method_inputs       : arguments to determine which methods to call
    '''
    parser = argparse.ArgumentParser()

    # Motif distribution (MD) method
    parser.add_argument(
        '-m',
        '--motif_distribution',
        help="This option takes the original motif distribution to process through the peak finder algorithm.",
        required=False
    )

    # Exponential method
    parser.add_argument(
        '-e',
        '--exponential',
        help="This option implement the exponential method on the motif distribution prior to processing it through the peak finder algorithm.",
        required=False
    )

    # Self power method
    parser.add_argument(
        '-s',
        '--self_power',
        help="This option implement the  self-power method on the motif distribution prior to processing it through the peak finder algorithm.",
        required=False
    )

    # threshold method
    parser.add_argument(
        '-t',
        '--threshold',
        help="This option implement a threshold cutoff on the motif distribution prior to processing it through the peak finder algorithm.",
        required=False
    )

    # fast fourier transform method
    parser.add_argument(
        '-f',
        '--fast_fourier',
        help="This option implement a fast fourier transform on the motif distribution prior to processing it through the peak finder algorithm.",
        required=False
    )
    
    # add parser by list because of the limitation of Jupyter Notebook
    # to use argparse
    method_inputs = parser.parse_args(args=selected_methods)

    return method_inputs


def _preprocess_mds_file(file_name):
    '''
    This function will read a motif distribution (MD) csv file. After
    reading the file, the function preprocess the data into the format
    required for the different methods

    Input parameters
    ---------------
    filename           : str
                         This parameter contain the name of input file
                         as a string

    Output returns
    ---------------
    raw_data           : list
                         Return the raw motif distribution data as a list
    exp_data           : list
                         Return the exponential data as a list
    selfpower_data     : list
                         Return the self-power data as a list
    threshold_data     : list
                         Return the threshold data as a list
    '''
    
    ids = []
    srrs = []
    raw_data = []
    exp_data = []
    selfpower_data = []
    threshold_data = []
    fft_data = []
    with open(file_name, newline='') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        next(csv_reader)
        for row in csv_reader:
#             raw_row = [int(value) for value in row[1:]]
            raw_row = [float(value) for value in row[2:]] # for md_table.csv
            # the above list comprehension is equivalent to:
    #         raw_row = []
    #         for value in row[1:]:
    #             raw_row.append(int(value))
            raw_data.append(raw_row)
            ids.append(row[0]) # Code to split ID string.split('_')[0]+'.'+row[0].split('.')[3])
            srrs.append(row[1].split('.')[0])

            # exp method
            exp_row = _exp_method(raw_row)
            exp_data.append(exp_row)

            # power method
            power_row = _selfpower_method(raw_row)
            selfpower_data.append(power_row)

            # threshold method - there are some bugs here
            threshold_row = _threshold_method(raw_row)
            threshold_data.append(threshold_row)
            
            # fast fourier transform method
            fft_row = _fft_method(raw_row)
            fft_data.append(fft_row)

    return ids, srrs, raw_data, exp_data, selfpower_data, threshold_data, fft_data


def _exp_method(raw_row):
    '''
    This function takes in the motif distribution as a list and
    process it using numpy exponential to get a list of the
    exponential value. The function returns a list of values after
    taking the exponential.

    Input parameters
    ---------------
    raw_row            : list
                         This parameter contain raw motif distribution
                         as a list

    Output returns
    ---------------
    exp_data           : list
                         This parameter is a list of values after the
                         exponential of the motif distribution is taken
    '''
    e = math.exp(1)  # returns E raised to the power of x
    exp_list = np.linspace(e, e, len(raw_row))
    exp_data = np.multiply(exp_list, raw_row)
    return exp_data


def _selfpower_method(raw_row):
    '''
    This function takes in the motif distribution as a list, takes the
    self-power of each value and returns a list.

    Input parameters
    ---------------
    raw_row            : list
                         This parameter contain raw raw motif distribution
                         as a list

    Output returns
    ---------------
    selfpower_data     : list
                         This parameter is a list of values after
                         taking self power
    '''
#     e = math.exp(1)
    selfpower_data = np.multiply(raw_row, raw_row)
    selfpower_data = np.array(selfpower_data, dtype=np.float64)
    return selfpower_data


def _threshold_method(raw_row):
    '''
    This function takes in the motif distribution as a list. The
    function takes each row and finds the max value to set the
    threshold as 50% of the max value. If the value is above the
    threshold, it is appended to the new list. If it is less than
    threshold, a value of 0 will be append to the list.

    Input Parameters
    ---------------
    raw_row            : list
                         This parameter contain raw raw motif distribution
                         as a list

    Output returns
    ---------------
    threshold_data     : list
                         This parameter filters the list to keep only
                         values that are above the threshold set at 50%
                         of the max.
    '''
    count = 0
    threshold_raw = []
    for i in raw_row:
        # print(raw_row)
        if (float(i) > max(raw_row) * 0.5):
            threshold_raw.append(raw_row[count])
        else:
            threshold_raw.append(0)
        count = count + 1
    return threshold_raw


def _fft_method(raw_row):
    '''
    This function takes in the motif distribution as a list, and computes a 1D n-point discrete
    Fourier Transform (DFT) with the efficient Fast Fourier Transform (FFT) algorithm.

    Input parameters
    ---------------
    raw_row            : list
                         This parameter contain raw raw motif distribution
                         as a list

    Output returns
    ---------------
    ddt_data          : list
                         This parameter is a list of values after taking fast fourier transform. 
                         Use input length as n
    '''
    fft_data = np.fft.fft(raw_row, n=None, axis=-1, norm=None)
    return fft_data


def _reduce_data_for_peak_detection(raw_row, reduce_size, smooth_size):
    '''
    Preprocess the raw MD distribution data; 1) Reduce the MD data to
    accentuate the signal, 2) Smoothing function to improve peak detection

    Input parameters
    --------------------
    raw_row            : list
                        input data is a .csv file of motif distribution.
                        Each row is a single data point.
    reduce_size        : int
                        size to reduce list into n-sized bin
    smooth_size        : int
                        size to smooth bucket from bucket_reduction()

    Output returns
    --------------------
    smoothed_row       : list
                        Returns a list with reduced and smoothed data
    '''

    def bucket_reduction(lst, n):
        '''
        Yield successive n-sized chunks from lst(list)

        Input parameters
        --------------------
        lst            : list
                        input data as a list
        n              : int
                        define the size for the buckets
        '''
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    def bucket_smooth(lst, n):
        '''
        Smoothing function through rolling bins

        Input parameters
        --------------------
        lst            : list
                        input data as a list
        n              : int
                        define the size for the buckets
        '''
        for i in range(0, len(lst) - n, 1):
            yield lst[i:i + n]

    ''' Data reduction of raw data into buckets of num_bucket units each
    Data divided into buckets of num_bucket units
    eg., [[0,5], [6, 10], [11,15], ...]
    Sum each bucket and generate new dataset with bucketed data.
    This allows us to accentuate peak data and create larger gaps between
    signal and noise
    '''
    bucketed_row = []
    for i in bucket_reduction(raw_row, reduce_size):
        bucketed_row.append(sum(i))

    ''' Smooth out the histogram created previously
    Generating rolling buckets
    eg., [[0,n], [1,n+1], [2, n+2], ...]
    Smooth out to make it easier to detect real peaks and filtering out
    1 unit wide 'peaks'.
    Smooths out noisy peak and improves detection rates.
    '''
    smoothed_row = []
    for i in bucket_smooth(bucketed_row, smooth_size):
        smoothed_row.append(mean(i))

    return smoothed_row


def draft_categorize_msd_file_into_tf_regulation(mds_file_path, tf_regulation_file):
    '''
    This function is not working properly at the moment. Read an
    mds file containing the motfif distribution and use peak
    detection to categorize tf regulation

    Input parameters
    ---------------
    mds_file_path      : string
                         Contains the full file path of the MDs file
    tf_regulation_file : file object
                         Object for writing tf regulation output
    smooth_size        : int
                         This parameter contain the size of smoothing
                         raw data

    Output returns
    ---------------
    None
    outputs to the file specified by the tf_regulation_file file object
    '''
    reduce_size = 10
    smooth_size = 10

    #sra = os.path.basename(mds_file_path).split('_')[0]
    ids, srrs, raw_data, exp_data, selfpower_data, threshold_data, fft_data = _preprocess_mds_file(mds_file_path)
    f = open(tf_regulation_file, "w")

    rows_by_counts = {}
    for i in range(0, len(raw_data)):
        raw_data_row = raw_data[i]
        # plot_raw_data(raw_data_row)
        reduced_row, peaks, valleys, peak_properties, valley_properties = find_peaks_and_or_valleys(
            raw_data_row, reduce_size, smooth_size
        )
        peak_count = len(peaks)
        valley_count = len(valleys)
        if (peak_count, valley_count) not in rows_by_counts:
            rows_by_counts[(peak_count, valley_count)] = []
        rows_by_counts[(peak_count, valley_count)].append(i)
        id = ids[i] #.replace('HO_', '').replace('_HUMAN.H10MO', '')
        srr = srrs[i]
        if peak_count == 0 and valley_count == 0:
            category = 'inactive'
        elif peak_count == 1 and valley_count == 0:
            category = 'active'
        elif peak_count == 0 and valley_count == 1:
            category = 'repressed'
        elif (peak_count == 2 or peak_count == 4) and valley_count == 0:
            category = 'offset'
        else:
            category = 'uncategorized'

        tf_regulation_line = '{},{},{},{},{}'.format(
            peak_count, valley_count, id, srr, category
        )

        peaks_bp = []
        for old_loc in peaks:
            new_loc = (old_loc + int(smooth_size/2))*reduce_size - 1500
            peaks_bp.append(new_loc)
        valleys_bp = []
        for old_loc in valleys:
            new_loc = (old_loc + int(smooth_size/2))*reduce_size - 1500
            valleys_bp.append(new_loc)

        for i in range(4):
            if len(peaks_bp) > i:
                tf_regulation_line += ',{}'.format(peaks_bp[i])
            else:
                tf_regulation_line += ','
        if len(valleys_bp) > 0:
            tf_regulation_line += ',{}'.format(valleys_bp[0])
        else:
            tf_regulation_line += ','
        tf_regulation_line += '\n'
        f.write(tf_regulation_line)
    f.close()
    return None


def _determine_scans(reduced_row):
    '''
    Determine if we should scan for peaks, valleys or both
    Uses 1.5 Std Dev outlier detection to determine if there is a
    substantial negative peak. If n points are less than
    mean(y) - stddev(y), then the "peak" will be treated as a valley

    Input parameters
    ---------------
    reduced_row        : list
                         reduced data to be used in peak detection

    Output returns
    ---------------
    peak_scan          : boolean value
                        Value to specify if we should look for peaks
    valley_scan        : boolean value
                        Value to specify if we should look for valleys
    '''
    mean_y = mean(reduced_row)
    stdev_y = stdev(reduced_row)
    stdev_15 = stdev_y * 1.5
    stdev_20 = stdev_y * 2

    peak_scan = False
    valley_scan = False
    scan_threshold = 5  # need n continuous signal
    signals = []
    for i in reduced_row:
        ''' generate signal array of peaks and valleys based on
        outlier calculation via Std Dev
        '''
        if i > (mean_y + stdev_15):
            signals.append(1)
        elif i < (mean_y - stdev_15):
            signals.append(-1)
        else:
            signals.append(0)

    '''
    Group by signal level to get the length of continuous signals
    [0, 0, 1, 1, 1, 0, 1, 1] => [(0,2),(1,3),(0,1),(1,2)]
    If there is a long continuous signal greater than scan_threshold,
    then indicate peak/valley
    '''
    grouped_signal_counts = [
        (k, sum(1 for i in g)) for k, g in groupby(signals)
    ]
    for (signal_value, count) in grouped_signal_counts:
        '''
        scan for peaks and valleys based on continuous
        signals > scan_threshold
        '''
        if signal_value == 1 and count > scan_threshold:
            peak_scan = True
        if signal_value == -1 and count > scan_threshold:
            valley_scan = True
    return peak_scan, valley_scan


def plot_peaks_and_valleys(reduced_row, peaks, valleys):
    '''
    Plot the peaks and valleys found using scipy find_peaks library

    Input parameters
    ---------------
    reduced_row        : list
                        reduced data to be used in peak detection
    peaks              : list
                        list of called peaks from find_peaks_and_or_valleys()
    valleys            : list
                        list of called valleys from find_peaks_and_or_valleys()

    Output returns
    ---------------
    None
    '''
    y = np.array(reduced_row)
    stdev_y = stdev(reduced_row)
    stdev_15 = stdev_y * 1.5
    mean_y = mean(reduced_row)
    median_y = median(reduced_row)
    meanArr = np.array([mean_y] * len(y))
    medianArr = np.array([median_y] * len(y))

    plt.rcParams["figure.figsize"] = (10, 3)
    plt.plot(y)
    plt.plot(meanArr, color="red")
    plt.plot(medianArr, color="red", linestyle='dashed')
    plt.plot(meanArr + stdev_15, color="orange")
    plt.plot(meanArr - stdev_15, color="orange")
    plt.plot(peaks, y[peaks], "xm")
    plt.plot(valleys, y[valleys], "xm")
    plt.show()
    return None


def plot_raw_data(raw_row, peaks_bp, valleys_bp):
    """ Plot the raw MD data """
    plt.rcParams["figure.figsize"] = (10, 3)
    y = np.array(raw_row)
    plt.plot(y)
    plt.plot(peaks_bp, y[peaks_bp], "xm")
    plt.plot(valleys_bp, y[valleys_bp], "xm")
    plt.show()


def zscore_conv(tf_data):
    """ input = a list of integer motif counts
        output = a list of float zscores """
    max_motif_count = max(tf_data) # get the highest count
    if max_motif_count <= 3: # ignore lines with too low of a count
        return False
    num_bins = max_motif_count+1
    num_25_percentile = max([int(math.floor(num_bins/4)),1]) # how many bins to ignore on 
    num_50_percentile = num_bins - num_25_percentile*2
    # calculate the center bins to use on sorted data to get the iqr
    iqr_bins = [i for i in range(num_25_percentile, num_25_percentile+num_50_percentile)]
    # motif_freqencies is a list of [bin_num, freqency]
    motif_freqencies = [[i,0] for i in range(num_bins)]
    for motif_count in tf_data:
        motif_freqencies[motif_count][1] += 1
    # sort the count bins by frequency
    motif_freqencies.sort(key=lambda bin: bin[1])
    # extract the value for the bins inside the iqr
    iqr_values = [motif_freqencies[bin][0] for bin in iqr_bins]
    iqr_mean = statistics.mean(iqr_values)
    iqr_stdev = statistics.stdev(iqr_values)
    zscores = [(motif_count - iqr_mean) / iqr_stdev for motif_count in tf_data]
    return zscores


if __name__ == "__main__":
    main()
