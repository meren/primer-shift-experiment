import sys
import numpy
import cPickle
import operator
import itertools


def bump_dict(d, key, val = 1):
    if d.has_key(key):
        d[key] += val
    else:
        d[key] = val
    
    return d

def get_sorted(d):
    return sorted(d.iteritems(), key=operator.itemgetter(1), reverse = True)

class Errors:
    def __init__(self, tags_dict, downstream = 1, upstream = 1):
        self.errors = {'upstream': {}, 'downstream': {}, 'linked': {}, 'mutation': {}, 'word': {}}
        self.null_stats = {'upstream': {}, 'downstream': {}, 'linked': {}, 'word':{}, 'bases': {}}
        self.error_locations = {}
        self.downstream = downstream
        self.upstream = upstream
        self.tags_dict = tags_dict
        
        self.bases = ['A', 'T', 'C', 'G']
        
        self.compute_null_stats()
    
    def compute_null_stats(self):
        
        u = [self.bases] * self.upstream
        upstream_combinations = [''.join(x) for x in list(itertools.product(*u))]

        d = [self.bases] * self.downstream
        downstream_combinations = [''.join(x) for x in list(itertools.product(*d))]
        
        l = [self.bases] * (1 + self.downstream + self.upstream)
        linked_combinations = [''.join(x) for x in list(itertools.product(*l))]
        linked_combinations = list(set([x[0:self.upstream] + '|' + x[self.upstream + 1:] for x in linked_combinations]))
        
        
        downstream_combinations = [] if downstream_combinations == [''] else downstream_combinations
        upstream_combinations = [] if upstream_combinations == [''] else upstream_combinations
        
        for tag in self.tags_dict:
            for seq in self.tags_dict[tag]:
                freq = self.tags_dict[tag][seq]
                
                # BASES:
                for base in self.bases:
                    self.null_stats['bases'] = bump_dict(self.null_stats['bases'], base, seq.count(base) * freq)
                
                
                # UPSTREAM AND DOWNSTREAM COMBINATIONS
                for upstream_combination in upstream_combinations:
                    uc = seq.count(upstream_combination) * freq
                    if uc:
                        self.null_stats['upstream'] = bump_dict(self.null_stats['upstream'], upstream_combination, uc)
                        
                for downstream_combination in downstream_combinations:
                    dc = seq.count(downstream_combination) * freq
                    if dc:
                        self.null_stats['downstream'] = bump_dict(self.null_stats['downstream'], downstream_combination, dc)
    
    
        print self.null_stats
    
    
    def add_error(self, abundant_seq, erronous_seq, pos):
        if self.error_locations.has_key(pos):
            self.error_locations[pos] += 1
        else:
            self.error_locations[pos] = 1
        
        downstream_seq, upstream_seq = '', '' 
        
        if self.upstream and not pos < self.upstream:
            upstream_seq = abundant_seq[pos - self.upstream:pos]
            self.errors['upstream'] = bump_dict(self.errors['upstream'], upstream_seq)

        if self.downstream and not len(abundant_seq) < pos + self.downstream + 1:
            downstream_seq = abundant_seq[pos + 1:pos + 1 + self.downstream]
            self.errors['downstream'] = bump_dict(self.errors['downstream'], downstream_seq)
    
        if upstream_seq or downstream_seq:
            linked_seq = downstream_seq + '|' + upstream_seq
            self.errors['linked'] = bump_dict(self.errors['linked'], linked_seq)
            
        if upstream_seq or downstream_seq:
            from_seq = downstream_seq + abundant_seq[pos].lower() + upstream_seq
            to_seq = downstream_seq + erronous_seq[pos].lower() + upstream_seq
            
            self.errors['word'] = bump_dict(self.errors['word'], '%s:%s' % (from_seq, to_seq))
 
        mutation = abundant_seq[pos] + ':' + erronous_seq[pos]
        self.errors['mutation'] = bump_dict(self.errors['mutation'], mutation)
            

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
    
    def show_errors_downstream(self):
        print get_sorted(self.errors['downstream'])
        
    def show_errors_upstream(self):
        print get_sorted(self.errors['upstream'])
        
    def show_errors_linked(self):
        print get_sorted(self.errors['linked'])
        
    def show_mutations(self):
        print get_sorted(self.errors['mutation'])
        
    def show_words(self):
        print get_sorted(self.errors['word'])
        
    def show(self):
        print
        print
        print 'Downstream:' 
        self.show_errors_downstream()

        print
        print
        print 'Upstream:' 
        self.show_errors_upstream()
    
        print
        print
        print 'Linked:' 
        self.show_errors_linked()
    
        print
        print
        print 'Mutations:' 
        self.show_mutations()
        
        print
        print
        print 'Words:'
        self.show_words()
    
        
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
    
    tags = [t[1] for t in tags_final]
    tags_dict_final = {}
    for tag in tags:
        tags_dict_final[tag] = tags_dict[tag] 
        
    errors = Errors(tags_dict = tags_dict_final, upstream = upstream, downstream = downstream)
    
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
                    errors.add_error(most_abundant_seq, current_seq, pos)
                

    print '* num in/dels ignored:', num_indels_ignored
    errors.show()
            
        
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
