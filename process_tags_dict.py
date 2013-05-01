import sys
import numpy
import cPickle
import operator


class Errors:
    def __init__(self, downstream = 1, upstream = 1):
        self.errors = {'upstream': {}, 'downstream': {}}
        self.error_locations = {}
        self.downstream = downstream
        self.upstream = upstream
    
    def add_error(self, seq, pos):
        if self.error_locations.has_key(pos):
            self.error_locations[pos] += 1
        else:
            self.error_locations[pos] = 1
            
        if not pos < self.downstream:
            downstream_seq = seq[pos - self.downstream:pos]
            if self.errors['downstream'].has_key(downstream_seq):
                self.errors['downstream'][downstream_seq] += 1
            else:
                self.errors['downstream'][downstream_seq] = 1

        if not len(seq) < pos + self.upstream + 1:
            upstream_seq = seq[pos + 1:pos + 1 + self.upstream]
            if self.errors['upstream'].has_key(upstream_seq):
                self.errors['upstream'][upstream_seq] += 1
            else:
                self.errors['upstream'][upstream_seq] = 1
    

    def show_error_distribution(self):
        err_dist = []
        for i in range(0, max(self.error_locations.keys())):
            if not self.error_locations.has_key(i):
                err_dist.append(0)
            else:
                err_dist.append(self.error_locations[i])
        
        print
        print err_dist
        
        import matplotlib.pyplot as plt
        plt.plot(err_dist)
        plt.xlim(xmin = 0, xmax = len(err_dist))
        plt.show()
    
    def show_error_downstream(self):
        print sorted(self.errors['downstream'].iteritems(), key=operator.itemgetter(1), reverse = True)
        
    def show_error_upstream(self):
        print sorted(self.errors['upstream'].iteritems(), key=operator.itemgetter(1), reverse = True)
        
        
def pp(n):
    """Pretty print function for very big integers"""
    if type(n) != int:
        return n

    ret = []
    n = str(n)
    for i in range(len(n) - 1, -1, -1):
        ret.append(n[i])
        if (len(n) - i) % 3 == 0:
            ret.append(',')
    ret.reverse()
    return ''.join(ret[1:]) if ret[0] == ',' else ''.join(ret)

def main(tags_dict_path, min_tag_abundance = 10, dry_run = False, save_individual_fasta = False, print_pcr_errors = False, downstream = 1, upstream = 1):
    tags_dict = cPickle.load(open(tags_dict_path))
    
    num_tags = len(tags_dict)
    tags_sorted = []
    abundances = []
    for tag in tags_dict:
        abundance = sum(tags_dict[tag].values())
        tags_sorted.append((abundance, tag), )
        abundances.append(abundance)
    tags_sorted.sort(reverse = True)
    
    total_number_of_reads = sum([x[0] for x in tags_sorted])
    tags_final = [x for x in tags_sorted if x[0] > min_tag_abundance]
    tags_discarded = [x for x in tags_sorted if x[0] <= min_tag_abundance]
    reads_represented = sum([x[0] for x in tags_final])
    reads_discarded = sum([x[0] for x in tags_discarded])
    
    print 'Number of tags ...........................: %s' % pp(num_tags)
    print 'Number of reads ..........................: %s' % pp(total_number_of_reads)
    print 'Mean tag abundance .......................: %.2f' % numpy.mean(abundances) 
    print 'Standard deviation .......................: %.2f' % numpy.std(abundances) 
    print 'Top 5 most abundant tags .................: %s' % ', '.join(['%s (%d)' % (x[1], x[0]) for x in tags_sorted[:5]])
    print 'Minimum tag abundance ....................: %s' % pp(min_tag_abundance)
    print 'Number of tags above --min-tag-abundance .: %s' % pp(len(tags_final))
    print '  Number of reads represented ............: %s (%.2f%% of all reads)' % (pp(reads_represented), reads_represented * 100.0 / total_number_of_reads)
    print 'Number of tags below --min-tag-abundance .: %s' % pp(len(tags_discarded))
    print '  Number of reads discarded ..............: %s (%.2f%% of all reads)' % (pp(reads_discarded), reads_discarded * 100.0 / total_number_of_reads)

    if dry_run:
        sys.exit()
    
    errors = Errors(upstream = upstream, downstream = downstream)
    
    tags = [t[1] for t in tags_final]
    
    num_indels_ignored = 0
    
    for i in range(0, len(tags)):
        tag = tags[i]
        
        if save_individual_fasta:
            f = open('TAG_%.5d_%s.fa' % (i, tag), 'w')
            counter = 0
            for seq in tags_dict[tag]:
                counter += 1
                for i in range(0, tags_dict[tag][seq]):
                    f.write('>%s | seq: %d | count: %d\n%s\n' % (tag, counter, i, seq))
                f.close()
        
        seqs_sorted_by_abundance = sorted(tags_dict[tag].iteritems(), key=operator.itemgetter(1), reverse = True)
        most_abundant_seq = seqs_sorted_by_abundance[0][0]
        
        for tpl in seqs_sorted_by_abundance[1:]:
            current_seq = tpl[0]
            if len(most_abundant_seq) != len(current_seq):
                # WE ARE IGNOREING IN/DELS
                num_indels_ignored += 1
                continue

            if print_pcr_errors:
                print '--'
                print most_abundant_seq
                print current_seq
                comp_line = ''
                num_line = ''
                for pos in range(0, len(current_seq)):
                    if most_abundant_seq[pos] != current_seq[pos]:
                        comp_line += '*'
                        num_line += '%d' % pos
                    else:
                        comp_line += ' '
                        num_line += ' '
                print comp_line
                print num_line
            
            for pos in range(0, len(current_seq)):
                if most_abundant_seq[pos] != current_seq[pos]:
                    errors.add_error(most_abundant_seq, pos)
                

    print 'downstream:' 
    errors.show_error_downstream()

    print 'upstream:' 
    errors.show_error_upstream()
    
    print 'num in/dels ignored:'
    print num_indels_ignored
            
        
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Parse and make sense of the Tags dict')
    parser.add_argument('tags_dict_path', metavar = 'TAGS_DICT',
                                        help = 'Serialized Tags dict to analyze')
    parser.add_argument('-m', '--min-tag-abundance', type=int, default = 10,
                                        help = 'minimum number of reads from one tag')
    parser.add_argument('--dry-run', action = 'store_true', default = False,
                                        help = 'When declared, program would give simple stats and exit')
    parser.add_argument('--save-individual-fasta', action = 'store_true', default = False,
                                        help = 'When declared, program stores sequence information per tag\
                                                that passes --min-abundance-tag parameter')
    parser.add_argument('--print-pcr-errors', action = 'store_true', default = False,
                                        help = 'If declared there will be a lot of output :/')
    parser.add_argument('--downstream', type=int, default = 1, metavar = 'INT',
                                        help = 'How many bases downstream of an error shall be stored')
    parser.add_argument('--upstream', type=int, default = 1, metavar = 'INT',
                                        help = 'How many bases upstream of an error shall be stored')
    


    args = parser.parse_args()

    sys.exit(main(args.tags_dict_path, args.min_tag_abundance, args.dry_run, args.save_individual_fasta, args.print_pcr_errors, args.downstream, args.upstream))
