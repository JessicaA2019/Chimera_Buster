## Chimera_buster
## By: Jessica L Albert
## Last edited : 3/8/24

import csv
import math
import time
import edlib


#sample_file = "barcode5_clusters_consensus.fasta"
#size_file = "barcode5_pre_clusters_consensus.fasta"
#mismatch_tolerance = 2
#output_name = "barcode5_2mm_pre_post_ds2_"

# Take cluster concensus fasta and create list for each entry
## [name, fwd umi, rev umi, size]

def readFastq(filename):
    #Reads FASTQ file and remove the special characters!
    samples = []
    with open(filename) as fh:
        while True:
            name = fh.readline().rstrip() # read name
            seq = fh.readline().rstrip() # read base sequence
            if len(seq) == 0:
                break
            samples.append([name,seq])
    return samples

#changes a list of lists into a sigle list

def flatten_list(data):
    flat_list = []
    for element in data:
        if type(element) == list:
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list

def check_mismatch(ref_seq, test_seq, mismatch_tolerance):
    alignment = edlib.align(ref_seq,test_seq,k=(mismatch_tolerance+1))
    mismatching = alignment['editDistance']
    '''if len(test_seq) == len(ref_seq):
        pos_base = 0
        mismatching = 0
        for letter in test_seq:
            if letter != ref_seq[pos_base]:
                mismatching = mismatching + 1
                pos_base = pos_base + 1
            else:
                pos_base = pos_base + 1
    else:
        test_base = 0
        ref_base = 0
        mismatching = 0
        while test_base < len(test_seq) and ref_base < len(ref_seq):
            if test_seq[test_base] != ref_seq[ref_base]:
                #check if deletion in test_seq
                if ref_seq[ref_base+1:ref_base+4] == test_seq[test_base:test_base+3]:
                    mismatching = mismatching + 1
                    test_base = test_base + 3
                    ref_base = ref_base + 4
                #check if insertion in test_seq
                elif ref_seq[ref_base:ref_base+3] == test_seq[test_base+1:test_base+4]:
                    mismatching = mismatching + 1
                    test_base = test_base + 4
                    ref_base = ref_base + 3
                #just a mismatch/multiple mismatches so just go to next nucleotide
                else:
                    mismatching = mismatching + 1
                    test_base = test_base + 1
                    ref_base = ref_base + 1
            else:
                test_base = test_base + 1
                ref_base = ref_base + 1
        '''
    return mismatching

def chimera_buster(sample_file, size_file, output_name, mismatch_tolerance):
    start_time = time.time()

    print ("Loading and preparing samples...")
    sample_list = readFastq(sample_file)
    sample_size_list = readFastq(size_file)

    entry_num = 0

    for entry in sample_list:
        sample_list[entry_num] = [entry[0].split(";"), entry[1]]
        sample_list[entry_num] = flatten_list(sample_list[entry_num])
        entry_num = entry_num + 1
        
    entry_num = 0

    for entry in sample_size_list:
        sample_size_list[entry_num] = entry[0].split(";")
        sample_size_list[entry_num] = flatten_list(sample_size_list[entry_num])
        entry_num = entry_num + 1
        
    size_dict = {}

    for sublist in sample_size_list:
        size_dict[sublist[-1][10:]] = sublist[-2][5:]
        
    working_list = []

    for sublist in sample_list:
        sample = []
        #add cluster name
        name = sublist[0][10:-9]
        sample.append(name)
        #add fwd and rev umi seqs from header
        umi_fwd = sublist[4][12:]
        umi_rev = sublist[5][12:]
        sample.append(umi_fwd)
        sample.append(umi_rev)
        # add final cluster size
        size = sublist[-3][5:]
        sample.append(int(size))
        # add earlier cluster size
        pre_size = size_dict[name]
        sample.append(int(pre_size))
        if len(sample) != 5:
            print ("Error")
        working_list.append(sample)

    
    working_list.sort(reverse=True, key=lambda working_list: (working_list[3], working_list[4]))

    lap_time = time.time()
    print ("Beginning filtering...")

    sample_num = 0

    fwd_umi = []
    rev_umi = []
    chimeras = []
    nonchimeras = []

    for sample in working_list:
        chimera_stat = False
        for seq in fwd_umi:
            if seq == sample[2]:
                chimera_stat = True
                break
            else:
                mismatching = check_mismatch(seq, sample[1], mismatch_tolerance)
                if mismatching <= mismatch_tolerance and mismatching >= 0:
                    chimera_stat = True
                    break
        for seq in rev_umi:
            if seq == sample[2]:
                chimera_stat = True
                break
            else:
                mismatching = check_mismatch(seq, sample[2], mismatch_tolerance)
                if mismatching <= mismatch_tolerance and mismatching >= 0:
                    chimera_stat = True
                    break
        if sample[1] not in fwd_umi:
            fwd_umi.append(sample[1])
        if sample[2] not in rev_umi:
            rev_umi.append(sample[2])
        if (len(chimeras)+len(nonchimeras)) % 1000 == 0:
            print(str(len(chimeras)+len(nonchimeras)) + " of " + str(len(working_list)) + " samples processed. (%s seconds)             " % (time.time() - lap_time), end = '\r' )
        if chimera_stat == False:
            nonchimeras.append(sample)
        elif chimera_stat == True:
            chimeras.append(sample)
        
                   
    print(str(len(chimeras)+len(nonchimeras)) + " of " + str(len(working_list)) + " samples processed. (%s seconds)              " % (time.time() - lap_time))
    lap_time = time.time()
    print("Writing outputs...")
    
    with open(output_name+"chimera_v_nonchimera.tsv",'w') as tfile:
        tfile.write("Sample \tUMI_fwd \tUMI_rev \tFinal cluster # \t Prelim Cluster # \t Chimera Status\n")
        for sample in nonchimeras:
            for item in sample:
                tfile.write(str(item) + "\t")
            tfile.write("N")
            tfile.write("\n")
        for sample in chimeras:
            for item in sample:
                tfile.write(str(item) + "\t")
            tfile.write("Y")
            tfile.write("\n")
            
    with open(output_name+"chimera_list.txt",'w') as tfile:
        for sample in chimeras:
            tfile.write(str(sample[0]) + "\n")
            
    with open(output_name+"nonchimera_list.txt",'w') as tfile:
        for sample in nonchimeras:
            tfile.write(str(sample[0]) + "\n")
            
    print("# of Chimeras: " + str(len(chimeras)))
    print("# of Non-Chimeras: " + str(len(nonchimeras)))
    print ("Complete.\n--- %s seconds total ---" % (time.time() - start_time))

