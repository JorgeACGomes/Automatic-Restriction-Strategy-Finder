from pydna.common_sub_strings import common_sub_strings as overlaps
# from Bio.Restriction import *
from operator import itemgetter as _itemgetter
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from Bio import SeqIO
from Bio.Restriction import CommOnly as Commercial
# AllEnzymes is empty
# from Bio.Restriction import AllEnzymes


class Restriction_Finder:

    def __init__(self, seqs, min_size=300, r_enzymes=Commercial):
        self.seqs = seqs
        self.contiguous = self.longest_contiguous_sequence()
        self.re = r_enzymes
        self.min_size = min_size
        self.restriction_finder()

    # gets the longest contiguous sequence in a list of sequences
    def longest_contiguous_sequence(self):
        init = overlaps(str(self.seqs[0].seq).lower(), str(self.seqs[1].seq).lower(), 25)
        ref = get_longest(str(self.seqs[1].seq).lower(), init)
        for s in self.seqs[2:]:
            contiguous = get_longest(str(s.seq).lower(), overlaps(ref, str(s.seq).lower(), 25))
            ref = contiguous
        return Dseqrecord(ref)

    def restriction_finder(self):
        # firstly lets find all of the enzymes that cut through the contiguous sequence

        # this dictionary is made up from all of the re that cut through the contiguous sequence
        cuts_contiguous = self.re.search(self.contiguous.seq, linear=False)
        enzymes = []
        # print(self.contiguous.seq)

        # to see if there are enzymes that cut more than once (produce two or more fragments)
        # print(len([x for x in cuts_contiguous.keys() if len(self.contiguous.cut(x)) > 2]))

        for re, sites in cuts_contiguous.items():
            # if it cuts
            if len(sites) > 0:
                # if it produces fragments, why do some not produce?? // now for 2 to restrict to single cutters
                # test for single cutters on contiguous sequence

                if len(self.contiguous.cut(re)) == 2:
                    # test for fragment size
                    if self.has_size(self.contiguous.cut(re)):
                        enzymes.append(re)

        print(enzymes)
        print(len(enzymes))
        # until here we have all the enzymes that execute single cuts through the contiguous sequence

        # now enzymes that cut 2 or more times through all sequences
        last_checked = ''
        cuts_all = self.re.search(self.seqs[0].seq, linear=False)
        for s in self.seqs[1:]:
            if s != last_checked:  # gain efficiency
                cuts_s = self.re.search(s.seq, linear=False)
                last_checked = s
                # common with the ones from last different
                cuts_all = {x: cuts_all[x] for x in cuts_all if x in cuts_s}

        # only double_cutters or more
        enzymes2 = [x for x in cuts_all if (len(self.seqs[0].cut(x)) > 2)]
        to_remove = []

        # fragment size verification
        for s in self.seqs:
            for e in range(len(enzymes2)):
                for f in s.cut(enzymes2[e]):
                    if len(str(f.seq)) < self.min_size:
                        to_remove.append(e)
        enzymes2 = [enzymes2[x] for x in range(len(enzymes2)) if x not in to_remove]

        print(len(enzymes2))

        # if everything went well we now have all the enzymes that cut on the contiguous sequence just one time
        # and the enzymes that cut through all the sequences more than just once
        # now, how should i pick the enzymes for the selective digestion?
        # hipothesis: check all permutations of enzymes and enzymes2 and test for fragment size, with subsequent digest?
        # hipothesis: get a pair of enzymes that produce different fragments for different seqs
        # should i maximize the difference(in bp) between the fragments produced?

        possible_results = {}
        for e1 in enzymes:
            for e2 in enzymes2:
                frags1 = self.seqs[0].cut([e1, e2])
                frags2 = self.seqs[3].cut([e1, e2])
                if len(frags1) != len(frags2):
                    print(e1,e2)



        return

    def has_size(self, frags):
        for f in frags:
            if len(str(f.seq)) < self.min_size:
                return False
        return True


# auxiliary functions
# gets the tuple with the maximum value for the 2nd index
# and slices the sequence with those given parameters
def get_longest(seq, out_list):
    coords = max(out_list, key=_itemgetter(2))  # output do overlaps
    start = coords[1]
    end = coords[1] + coords[2]
    return seq[start:end]




# transforms string sequences into Dseqrecord Objects


# testing        
if __name__ == '__main__':
    testseqs = ['aactgatcgtcgaactgaaactgatcg',
                 'aactgacaaaatcgaactgactgacaa', 
                 'aacttagactcgaactgaaacttagac',
                 'aactcgaacaatcgaactgactcgaac', 
                 'aactgatacatcgaactgaactgatac']
    # as i will be testing with plasmids with over 2k bp might as well use a fasta file as input
    fileseqs = [Dseqrecord(seq) for seq in SeqIO.parse('seqs.txt', 'fasta')]
    testseqrecs = [Dseqrecord(x, circular=True) for x in testseqs]
    
    rf = Restriction_Finder(fileseqs, 300, Commercial)
