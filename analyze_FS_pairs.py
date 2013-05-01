#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import copy
import numpy
import cPickle
import ConfigParser


E = os.path.exists
J = os.path.join


import IlluminaUtils.lib.fastqlib as u
from IlluminaUtils.utils.runconfiguration import RunConfiguration
from IlluminaUtils.utils.helperfunctions import ReadIDTracker 
from IlluminaUtils.utils.helperfunctions import QualityScoresHandler
from IlluminaUtils.utils.helperfunctions import visualize_qual_stats_dict
from IlluminaUtils.utils.helperfunctions import store_cPickle_obj
from IlluminaUtils.utils.helperfunctions import load_cPickle_obj



# reverse and forward primers for V6
reverse_primer = re.compile('AGGTG.TGCATGG[C,T][C,T]GTCG') 

# forward primers:
# 
# 967F-AQ    C    T A A    C    CG A    .    G    AACCT [CT] ACC
# 967F-UC3   A    T A C    G    CG A    [AG] G    AACCT T    ACC
# 967F-PP    C    . A C    G    CG A    A    G    AACCT T    A.C
# 967F-UC12* C    A A C    G    CG [AC] A    [AG] AACCT T    ACC
#  COMBINED: [CA] . A [AC] [CG] CG [AC] .    [AG] AACCT [CT] A.C
forward_primer = re.compile('[CA].A[AC][CG]CG[AC].[AG]AACCT[CT]A.C')
 


class ConfigError(Exception):
    pass


def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False


def get_rp_match(seq, reverse_primer):
    offset = 70
    rp_match = None

    while 1:
        m = reverse_primer.search(seq, offset)
                
        if not m:
            break
        else:
            rp_match = m
            offset = m.start() + 1

    return rp_match

def get_percent(x, y):
    if y == 0:
        return 0.0
    else:
        return x * 100.0 / y

def reverse_complement(seq):
    conv_dict = {'A': 'T',
                 'T': 'A',
                 'C': 'G',
                 'G': 'C',
                 'N': 'N'}

    return ''.join(reversed([conv_dict[n] for n in seq]))


def main(config, visualize_quality_curves = False):    

    #####################################################################################
    # dealing with output file pointers..
    #####################################################################################
    
    GetFilePath = lambda p: os.path.join(config.output_directory, config.project_name + '-' + p)
    
    errors_fp                 = open(GetFilePath('STATS.txt'), 'w')
    
    missing_sequencing_oligo_R1  = u.FastQOutput(GetFilePath('MISSING_SEQUENCING_OLIGO_R1.fastq'))
    missing_sequencing_oligo_R2  = u.FastQOutput(GetFilePath('MISSING_SEQUENCING_OLIGO_R2.fastq'))
   
    sequencing_oligo_anomaly_R1  = u.FastQOutput(GetFilePath('SEQUENCING_OLIGO_ANOMALY_R1.fastq'))
    sequencing_oligo_anomaly_R2  = u.FastQOutput(GetFilePath('SEQUENCING_OLIGO_ANOMALY_R2.fastq'))

    has_sequencing_error_R1  = u.FastQOutput(GetFilePath('HAS_SEQUENCING_ERROR_R1.fastq'))
    has_sequencing_error_R2  = u.FastQOutput(GetFilePath('HAS_SEQUENCING_ERROR_R2.fastq'))

    tag_length_anomaly_R1  = u.FastQOutput(GetFilePath('TAG_LENGTH_ANOMALY_R1.fastq'))
    tag_length_anomaly_R2  = u.FastQOutput(GetFilePath('TAG_LENGTH_ANOMALY_R2.fastq'))
    
    missing_primer_R1 = u.FastQOutput(GetFilePath('MISSING_V6_PRIMER_R1.fastq'))
    missing_primer_R2 = u.FastQOutput(GetFilePath('MISSING_V6_PRIMER_R2.fastq'))
    
    passed_R1  = u.FastQOutput(GetFilePath('PASSED_R1.fastq'))
    passed_R2  = u.FastQOutput(GetFilePath('PASSED_R2.fastq'))
    
    tag_length_distribution_output_path   = open(GetFilePath('TAG_LENGTH_DIST.txt'), 'w')
    
    perfect_reads_fasta_fp    = open(GetFilePath('PERFECT_reads.fa'), 'w')
    id_tracker_output_path    = GetFilePath('READ_IDs.cPickle.z')
    tags_dict_output_path    = GetFilePath('TAGS_DICT.cPickle')
    
    if visualize_quality_curves:
        qual_dict_output_path = GetFilePath('Q_DICT.cPickle.z')

    #####################################################################################
    # some useful variables before we begin..
    #####################################################################################
   
    number_of_pairs = 0
    number_of_pairs_passed = 0
    missing_sequencing_oligo = 0
    sequencing_oligo_anomaly = 0
    has_sequencing_error = 0
    perfect_compliance = 0
    tag_length_anomaly = 0
    rp_failed = 0
    fp_failed = 0

    if visualize_quality_curves:
        qual_dict = QualityScoresHandler()
        
    read_id_tracker = ReadIDTracker()

    tag_length_distribution = {}
    
    tags_dict = {}


    #####################################################################################
    # main loop per file listed in config:
    #####################################################################################
    COMPRESSED = lambda x: os.path.basename(config.pair_1[index]).endswith('.gz')
    
    for index in range(0, len(config.pair_1)):
        try:
            pair_1 = u.FastQSource(config.pair_1[index], compressed = True if COMPRESSED(config.pair_1[index]) else False)
            pair_2 = u.FastQSource(config.pair_2[index], compressed = True if COMPRESSED(config.pair_2[index]) else False)

        except u.FastQLibError, e:
            print "FastQLib is not happy.\n\n\t", e, "\n"
            sys.exit()


        #####################################################################################
        # main loop per read:
        #####################################################################################
       
        while pair_1.next() and pair_2.next():
            number_of_pairs += 1
            h1, s1, q1 = pair_1.entry.header_line, pair_1.entry.sequence, pair_1.entry.qual_scores
            h2, s2, q2 = pair_2.entry.header_line, reverse_complement(pair_2.entry.sequence), ''.join(reversed(pair_2.entry.qual_scores))

            if pair_1.p_available:
                pair_1.print_percentage()


            #####################################################################################
            # PHASE 1: getting rid of sequencing primers.
            #
            #                   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx     <- s1
            #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                     <- rc'd s2
            #       ^----------^                                       ^--------------^
            #  sequencing_oligo_in_s2                               sequencing_oligo_in_s1
            #
            #
            #####################################################################################
            
            # trim from the end to get more matches:
            sequencing_oligo_in_s1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

            # trim from the beginning to get more matches:
            sequencing_oligo_in_s2 = 'TCTTTCCCTACACGACGCTCTTCCGATCT'


            sequencing_oligo_in_s1_pos = s1.find(sequencing_oligo_in_s1[:15])
            sequencing_oligo_in_s2_pos = s2.find(sequencing_oligo_in_s2[15:])
            
            if sequencing_oligo_in_s1_pos == -1 or sequencing_oligo_in_s2_pos == -1:
                # missing sequencing oligo
                missing_sequencing_oligo_R1.store_entry(pair_1.entry)
                missing_sequencing_oligo_R2.store_entry(pair_2.entry)
                missing_sequencing_oligo += 1
                
                read_id_tracker.update(pair_1, pair_2, 'MISS_SEQ_OLIGO')
                if visualize_quality_curves:
                    qual_dict.update(pair_1, pair_2, 'MISS_SEQ_OLIGO')
                    
                continue
            
            if sequencing_oligo_in_s1_pos < sequencing_oligo_in_s2_pos:
                # sequencing oligo anomaly. we gotta figure out why this happens. 
                sequencing_oligo_anomaly_R1.store_entry(pair_1.entry)
                sequencing_oligo_anomaly_R2.store_entry(pair_2.entry)
                sequencing_oligo_anomaly += 1
                
                read_id_tracker.update(pair_1, pair_2, 'SEQ_OLIGO_ANOMALY')
                if visualize_quality_curves:
                    qual_dict.update(pair_1, pair_2, 'SEQ_OLIGO_ANOMALY')
    
                continue
           
            # update s1 and s2: 
            s1 = s1[:sequencing_oligo_in_s1_pos]
            s2 = s2[sequencing_oligo_in_s2_pos + len(sequencing_oligo_in_s2[15:]):]
            
            if s1 != s2:
                has_sequencing_error += 1
                has_sequencing_error_R1.store_entry(pair_1.entry)
                has_sequencing_error_R2.store_entry(pair_2.entry)

                read_id_tracker.update(pair_1, pair_2, 'HAS_SEQ_ERR')
                if visualize_quality_curves:
                    qual_dict.update(pair_1, pair_2, 'HAS_SEQ_ERR')

                continue
            
            perfect_compliance += 1

            #####################################################################################
            # PHASE 2: finding primers and isolating tags at both ends. now we can work only with
            #          s1 since both s1 and s2 are identical to each other.
            #
            #     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   <- s1
            #     ^------------^                                               ^--------------^
            #      random tag1  ^-------^                              ^------^  random tag2
            #                      FP    ^----------------------------^   RP
            #                                          V6
            #
            #
            #####################################################################################


            #####################################################################################
            # find reverse primer, retain 'random tag2'
            #####################################################################################

            # find all instances of matching
            rp_obj = get_rp_match(s1, reverse_primer)

            if (not rp_obj):
                # reverse primer wasn't there
                rp_failed += 1

                read_id_tracker.update(pair_1, pair_2, 'FAILED_RP')
                if visualize_quality_curves:
                    qual_dict.update(pair_1, pair_2, 'FAILED_RP')

                continue

            tag_2 = s1[rp_obj.end():]
            rp = s1[rp_obj.start():rp_obj.end()]
            s1 = s1[:rp_obj.start()]
           

            #####################################################################################
            # find forward primer, retain 'random tag1'
            #####################################################################################

            fp_obj = forward_primer.search(s1)
            
            if (not fp_obj):
                # forward primer wasn't there
                fp_failed += 1

                read_id_tracker.update(pair_1, pair_2, 'FAILED_FP')
                if visualize_quality_curves:
                    qual_dict.update(pair_1, pair_2, 'FAILED_FP')
 
                continue

            tag_1 = s1[:fp_obj.start()]
            fp = s1[fp_obj.start():fp_obj.end()]
            s1 = s1[fp_obj.end():]

            tag = tag_1 + '-' + tag_2
            tag_len = len(tag_1) + len(tag_2)
            
            if tag_length_distribution.has_key(tag_len):
                tag_length_distribution[tag_len] += 1
            else:
                tag_length_distribution[tag_len] = 1
            
            # we expect random tags to be between 13 and 23 in length.
            if tag_len < 13 or tag_len > 23:
                tag_length_anomaly += 1
                tag_length_anomaly_R1.store_entry(pair_1.entry)
                tag_length_anomaly_R2.store_entry(pair_2.entry)

                read_id_tracker.update(pair_1, pair_2, 'TAG_LEN_ANOMALY')
                if visualize_quality_curves:
                    qual_dict.update(pair_1, pair_2, 'TAG_LEN_ANOMALY')

                continue

            #####################################################################################
            # PHASE 3: we have everything we need:
            # 
            # tag_1 = random tag 1
            # tag_2 = random tag 2
            # tag = tag_1 + '-' + tag_2
            # fp = forward primer
            # rp = distal primer
            # s1 = v6 read
            #
            #####################################################################################          
            
            number_of_pairs_passed += 1
            passed_R1.store_entry(pair_1.entry)
            passed_R2.store_entry(pair_2.entry)
            
            read_id_tracker.update(pair_1, pair_2, 'PASSED')
            if visualize_quality_curves:
                qual_dict.update(pair_1, pair_2, 'PASSED')
            
            #read = fp + s1 + rp
            read = s1
            
            if tags_dict.has_key(tag):
                if tags_dict[tag].has_key(read):
                    tags_dict[tag][read] += 1
                else:
                    tags_dict[tag][read] = 1
            else:
                tags_dict[tag] = {read: 1}
            

            perfect_reads_fasta_fp.write('>%s\n%s\n' % (pair_1.entry.header_line, s1))

    # DONE.


    total_pairs_failed = number_of_pairs - number_of_pairs_passed
    errors_fp.write('number of pairs                 : %d\n' % (number_of_pairs))
    errors_fp.write('perfect compliance observed     : %d (%%%.2f of all pairs)\n' % (perfect_compliance, get_percent(perfect_compliance, number_of_pairs)))
    errors_fp.write('total pairs passed all phases   : %d (%%%.2f of all pairs)\n' % (number_of_pairs_passed, get_percent(number_of_pairs_passed, number_of_pairs)))
    errors_fp.write('total pairs failed              : %d (%%%.2f of all pairs)\n' % (total_pairs_failed, get_percent(total_pairs_failed, number_of_pairs)))
    errors_fp.write('  missing sequencing oligo      : %d (%%%.2f of all failed pairs)\n' % (missing_sequencing_oligo, get_percent(missing_sequencing_oligo, total_pairs_failed)))
    errors_fp.write('  sequencing oligo anomaly      : %d (%%%.2f of all failed pairs)\n' % (sequencing_oligo_anomaly, get_percent(sequencing_oligo_anomaly, total_pairs_failed)))
    errors_fp.write('  has sequencing error          : %d (%%%.2f of all failed pairs)\n' % (has_sequencing_error, get_percent(has_sequencing_error, total_pairs_failed)))
    errors_fp.write('  forward primer was missing    : %d (%%%.2f of all failed pairs)\n' % (fp_failed, get_percent(fp_failed, total_pairs_failed)))
    errors_fp.write('  reverse primer was missing    : %d (%%%.2f of all failed pairs)\n' % (rp_failed, get_percent(rp_failed, total_pairs_failed)))
    
    for fate in [f for f in read_id_tracker.fates if f.startswith('FAILED')]:
        errors_fp.write('  %s%s: %d (%%%.2f of all failed pairs)\n' % (fate,
                                                                       (' ' * (30 - len(fate))),
                                                                       len(read_id_tracker.ids[fate]),
                                                                       get_percent(len(read_id_tracker.ids[fate]), total_pairs_failed)))

    print

   
    sys.stderr.write('\rStoring tag length distribution ...')
    for i in sorted(tag_length_distribution.keys()):
        tag_length_distribution_output_path.write('%s\t%s\n' % (i, tag_length_distribution[i]))
    tag_length_distribution_output_path.close()
    sys.stderr.write('\n')
    
    sys.stderr.write('\rRead ID tracker dict is being stored ...')
    sys.stderr.flush()
    read_id_tracker.store(id_tracker_output_path)
    sys.stderr.write('\n')

    sys.stderr.write('\rTags dict is being stored ...')
    cPickle.dump(tags_dict, open(tags_dict_output_path, 'wb'))
    sys.stderr.write('\n')

    if visualize_quality_curves:
        sys.stderr.write('\rQuality scores dict is being stored ...')
        sys.stderr.flush()
        qual_dict.store_dict(qual_dict_output_path)

        for entry_type in qual_dict.entry_types:
            sys.stderr.write('\rQuality scores visualization in progress: %s%s' % (entry_type, ' ' * 20))
            sys.stderr.flush()
            visualize_qual_stats_dict(qual_dict.data[entry_type], GetFilePath(entry_type),\
                            title = 'Mean PHRED scores for pairs tagged as "%s"' % entry_type)
        sys.stderr.write('\n')

    errors_fp.close()
    perfect_reads_fasta_fp.close()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Dealing with V6 perfect overlaps')
    parser.add_argument('user_config', metavar = 'CONFIG_FILE',
                                        help = 'User configuration to run. See the source code to\
                                                see an example.')
    parser.add_argument('--visualize-quality-curves', action = 'store_true', default = False,
                                        help = 'When set, mean quality score for individual bases will be\
                                                stored and visualized for each group of reads.')


    args = parser.parse_args()
    user_config = ConfigParser.ConfigParser()
    user_config.read(args.user_config)

    try: 
        config = RunConfiguration(user_config)
    except ConfigError, e:
        print "There is something wrong with the config file. This is what we know: \n\n", e
        print
        sys.exit()

    sys.exit(main(config,
                  visualize_quality_curves = args.visualize_quality_curves))
