
import numpy as np
import pandas as pd
import random
import multiprocessing
import sys

## bash line to get the condensed sam file

## bash line to get the names
# bioawk -c'fastx' '{print $name" "length($seq)}' discovar_5th_contigs_mt9k.fa > contig_names_mt9kb.txt

def load_input_data(filename):
    ## put here as a comment the bash line needed to create this
    data = pd.read_csv(filename, delimiter=" ", names=["name", "p1", "p2", "cuenta"])
    return data

def filter_data_length(data, contig_names, min_contig_length):
    length_filtered_contig_names = {name: largo for name, largo in contig_names.iteritems() if largo>min_contig_length}
    filtered_data = data[data['name'].isin(length_filtered_contig_names.keys())]
    return filtered_data

def shuffle_data(data, n_sample, contig_names):
    random_names = random.sample(list(data['name'].unique()), n_sample)
    filtered_data = data[data['name'].isin(random_names)]
    random_contig_dict = {name: contig_names[name] for name in random_names}
    return filtered_data, random_contig_dict

def link_coverage_pileup(data, contig_names):
    bulk_data=[]
    contig_link_coverage={}

    for contig, largo in contig_names.iteritems():
        fdata = data[data["name"]==contig]

        if fdata.size > 0:

            bin_number = int(largo/BIN_SIZE)+1
            mdata = [0]*bin_number

            for k, v in fdata.iterrows():
                for i in xrange(int(v.p1/BIN_SIZE), int(v.p2/BIN_SIZE)-1):
                    mdata[i]+=v.cuenta

            BIN_BUFFER=int(LIBRARY_SIZE/BIN_SIZE) # Bins to ignore in each end because of the border effect
            for p in xrange(int(BIN_BUFFER), int(bin_number-BIN_BUFFER)):
                bulk_data.append(mdata[p])

            contig_link_coverage[contig] = mdata

    return np.array(bulk_data), contig_link_coverage

def calculate_threshold(data, contig_names, percentile):
    bulk_data, link_pileup=link_coverage_pileup(data, contig_names)
    threshold = np.percentile(bulk_data, percentile)
    return threshold, bulk_data

def get_threshold((data, contig_names, percentile)):
    ## select a random sample from the data
    random_data, random_contig_dict = shuffle_data(data, 100, contig_names)
    ## calculate threshold
    threshold, bdata = calculate_threshold(data, random_contig_dict, percentile)
    print ".",
    return threshold

def compute_threshold(data, contig_names, reps, thredas, percentile):
    """ Get filtered data and compute the threshold multiple times """
    input_data= tuple([(data, contig_names, percentile) for x in xrange(reps)])
    p=multiprocessing.Pool(thredas)
    threshold_list = p.map(get_threshold, input_data)
    threshold = np.mean(threshold_list)
    return threshold, threshold_list

def process_contig(coverage_vector, threshold):

    x = np.array(coverage_vector)
    fx = x[BIN_BUFFER:len(x)-BIN_BUFFER]
    masked_x=fx>threshold

    true_count = masked_x.sum()
    false_count = len(fx)-masked_x.sum()

    return fx, masked_x, true_count, false_count

def split_in_chunks(fx, masked):
    """ Split the contigs into the new broken chunks """
    current = masked[0]
    changes = []
    prev_change = 0
    good_intervals=[]
    bad_intervals=[]
    for x in xrange(len(masked)):
        if masked[x] != current:
            changes.append(x)
            if current == True:
                good_intervals.append([prev_change, x])
            current = masked[x]
            prev_change=x
    if current == True:
        good_intervals.append([prev_change, x])
    return good_intervals

def contig_bins2bases(splited_coords):
    if len(splited_coords) == 1:
        return [(splited_coords[0]*BIN_SIZE)+(BIN_BUFFER*BIN_SIZE*2), ]
    else:
        part_sizes=[(splited_coords[0]*BIN_SIZE)+(BIN_BUFFER*BIN_SIZE),]
        for p in xrange(1, len(splited_coords)-1):
            part_sizes.append(splited_coords[p]*BIN_SIZE)
        part_sizes.append((splited_coords[-1]*BIN_SIZE)+(BIN_BUFFER*BIN_SIZE))
    return part_sizes

def measure_contig((data, name, contig_length, threshold)):

    largo = contig_length
    fdata = data[data["name"] == name]
    # if fdata.size[0] == 0:
    #     return {}

    bin_number = int(largo / BIN_SIZE) + 1
    mdata = [0] * bin_number

    for k, v in fdata.iterrows():
        for i in xrange(int(v.p1 / BIN_SIZE), int(v.p2 / BIN_SIZE) - 1):
            mdata[i] += v.cuenta

    contig_metrics={}
    splited_contigs_length=[]
    accumulated_supported_content=[] ## TODO hacer que lo reyene
    breakpoints_per_contig=[]
    
    fx, masked_x, tc, fc = process_contig(mdata, threshold)
    if tc == 0:
        return {"breakpoints_count":0, 
                "number_of_parts": 1, 
                "size_of_contigs": contig_length,
                "original_contig_size": contig_length,
                "breakpoints_per_kbp": 0}

    sections=split_in_chunks(fx, masked_x)
    for gi in sections:
        splited_contigs_length.append(gi[1] - gi[0])
    breakpoints_per_contig.append(len(sections)-1)

    contig_metrics['breakpoints_count']=len(sections)-1
    contig_metrics['number_of_parts']=len(splited_contigs_length)
    contig_metrics['size_of_contigs']=contig_bins2bases(splited_contigs_length)
    contig_metrics['original_contig_size']=contig_length
    contig_metrics['breakpoints_per_kbp']=(len(sections)-1)/(contig_length/1000.0)

    return contig_metrics

def compute_contig_measurements(data, contig_names, threshold, threads, MIN_CONTIG_LENGTH):
    """ `Compute the measurement for all the contigs present in the dataset """
    print "Number of uinque scaffolds present: %s" %(len(data['name'].unique()), )
    # bulkdata, linkdata = link_coverage_pileup(data, contig_names)

    print "Done pileup collectiong args to run the multip" 
    args=[]
    for k, v in contig_names.iteritems():
        if v > MIN_CONTIG_LENGTH:
            args.append((data, k, v, threshold))
    print len(args)
    print "Starting pool"
    p=multiprocessing.Pool(threads)
    contig_measurements = p.map(measure_contig, args)
    print "Done with pool"
    return contig_measurements   


if __name__ == "__main__":

  ## Sam processing line
  # cat <SAMFILE>.sam |
  # grep - v '@SQ' |
  # grep - v '@PG' |
  # awk '{if($9>=60 && $7=="=")print $0}' |
  # awk '{print $3" "int($4/10)*10 " "int($8/10)*10}' |
  # sort - k1, 2 |
  # uniq - c |
  # awk '{print $2" "$3" "$4" "$1}' > link_coverage.txt

  filename = sys.argv[1]
  contigs_name_file = sys.argv[2]
  out_filename = sys.argv[3]
  percentile = float(sys.argv[4])
  threads = int(sys.argv[5])

  LIBRARY_SIZE=11000
  BIN_SIZE=1000
  BIN_BUFFER=int(LIBRARY_SIZE/BIN_SIZE)
  MIN_CONTIG_LENGTH=LIBRARY_SIZE*3

  print "Loading contig names..."
  ## Load contig name_length dict
  contig_names = {x.split(" ")[0]: int(x.strip().split(" ")[1]) for x in open(contigs_name_file)}

  print "Loading data..."
  ## load data form file into data
  data = load_input_data(filename)
  print "All table size: %s" %(data.shape, )

  print "filtering data..."
  ## filter data by size
  data = filter_data_length(data, contig_names, MIN_CONTIG_LENGTH)
  print "Length filtered size: %s" %(data.shape, )

  m_threshold, threshold_reps=compute_threshold(data, contig_names, 100, threads, percentile)
  contig_measurements = compute_contig_measurements(data, contig_names, m_threshold, threads, MIN_CONTIG_LENGTH)

  print "Writing the files..."
  outfile = open(out_filename, 'w')
  outfile.write("#PARAMETERS:\n")
  outfile.write("#FILENAME: %s\n" %(filename, ))
  outfile.write("#CONTIG FILENAME: %s\n" %(contigs_name_file, ))
  outfile.write("#LIBRARY SIZE: %s\n" %(LIBRARY_SIZE, ))
  outfile.write("#BIN SIZE: %s\n" %(BIN_SIZE,))
  outfile.write("#MIN CONTIG SIZE: %s\n" %(MIN_CONTIG_LENGTH,))
  outfile.write("#PERCENTILE: %s\n" %(percentile, ))
  outfile.write("#THREADS: %s\n" %(threads,))
  outfile.write("#THRESHOLD_MEAN: %s\n" % (m_threshold,))
  outfile.write("#THRESHOLD_REPS: %s\n" % (threshold_reps,))
  for x in contig_measurements:
    outfile.write("%s\n" %(x, ))
