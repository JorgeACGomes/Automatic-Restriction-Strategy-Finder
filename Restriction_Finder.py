from pydna.common_sub_strings import common_sub_strings as overlaps
# from Bio.Restriction import *
from operator import itemgetter as _itemgetter
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq

# from warnings import warn


class Restriction_Finder:

    def __init__(self, seqs, min_size=300, r_enzymes='Commercial'):
        self.seqs = seqs
        self.Dseqs = toDseqs(seqs)
        # teste de eficicencia
        self.conti = 'aact'
        # manter este
        self.contiguous = self.longest_contiguous_sequence()
        self.re = self.getRenzymes(r_enzymes)
        self.min_size = min_size
        self.restriction_finder()

    # given the user input it decides which set of enzymes to use, either All of the know enzymes in  REBase or
    # only those which are commercially available

    def getRenzymes(self, re):
        if re.lower() == 'commercial' or re.lower() == 'comm':
            from Bio.Restriction import CommOnly as Commercial
            return Commercial
        elif re.lower() == 'all':
            from Bio.Restriction import AllEnzymes as All
            return All
        else:
            from Bio.Restriction import CommOnly as Commercial
            # warning
            print('WARNING: Given restriction enzyme group is not available!\n'
                  'Enzyme list set to: Commercialy available\n')
            return Commercial

    # gets the longest contiguous sequence in a list of sequences
    def longest_contiguous_sequence(self):
        # see whats best to set as the limit because of efficiency
        init = overlaps(self.seqs[0], self.seqs[1], 4)
        ref = get_longest(self.seqs[1], init)
        for sequence in self.seqs[2:]:
            contiguous = get_longest(sequence, overlaps(ref, sequence, 2))
            ref = contiguous
        # return an Dseq object so it can be manipulated by restriction package from BioPython
        return toDseqs(ref)

    def restriction_finder(self):
        # very early yet, i was just messing around with it to try and understand how does Bio.Restirctio package works,
        # but i'm getting an error when i try to use the catalyse() function with the enzymmes that cut the contiguous
        # sequence

        # this dictionary is made up from all of the re that cut through the contiguous sequence
        dicti = self.re.search(self.contiguous, linear=False)
        for x, y in dicti.items():
            if len(y) != 0:
                # print(x)
                # print(dicti)

                # This was supposed to retrieve 2 Dseqs objects, but i get an error when i try to use it
                print(x.catalyse(self.contiguous))
                # >>AttributeError: 'IUPACAmbiguousDNA' object has no attribute 'split'

                # print(type(x))


# auxiliary functions

# gets the tuple with the maximum value for the 2nd index
# and slices the sequence with those given parameters
def get_longest(seq, out_list):
    coords = max(out_list, key=_itemgetter(2))  # output do overlaps
    start = coords[1]
    end = coords[1] + coords[2]
    return seq[start:end]


# transforms string sequences into Dseqrecord Objects
def toDseqs(seqs):
    if type(seqs) == list:
        return [Dseqrecord(x) for x in seqs]
    else:
        return Dseq(seqs)


# testing        
if __name__ == '__main__':
    testseqs = ['aactgatcgtcgaactgaaactgatcg', 'aactgacaaaatcgaactgactgacaa', 'aacttagactcgaactgaaacttagac',
                'aactcgaacaatcgaactgactcgaac', 'aactgatacatcgaactgaactgatac']
    rf = Restriction_Finder(testseqs, 300, 'comm')

# z=dseqrecord('dasdasd')
# z.cut(bamHI)
# overlaps só aceita strings
# limit = 25
# overlaps é case sensitive