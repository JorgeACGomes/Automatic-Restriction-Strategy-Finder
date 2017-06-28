from pydna.common_sub_strings import common_sub_strings as overlaps
from operator import itemgetter as _itemgetter
from pydna.dseqrecord import Dseqrecord
from Bio import SeqIO
from Bio.Restriction import CommOnly as Commercial
# AllEnzymes is empty
# from Bio.Restriction import AllEnzymes


class Plasmid:
    # total is the complete sequence, contiguous is the contiguous part, and insert is the non-contiguous
    def __init__(self, Dseqrecord, contiguous):
        self.total = Dseqrecord
        self.contiguous = contiguous
        self.insert = self.get_insert()
        self.cutters = None

    def get_insert(self):
        start = str(self.total.seq).lower().find(str(self.contiguous.seq).lower())
        end = start + len(self.contiguous.seq) - 1
        insert = str(self.total.seq)[end + 1::] + str(self.total.seq)[0:start]
        return Dseqrecord(insert)

    def set_cutters(self,cutters):
        self.cutters = cutters


class Restriction_Finder:

    def __init__(self, seqs, min_size=300, r_enzymes=Commercial):
        self.seqs = seqs
        self.contiguous = self.longest_contiguous_sequence()
        self.re = r_enzymes
        self.min_size = min_size
        self.linear = self.seqs[0].linear
        self.plasmids = [Plasmid(x, self.contiguous) for x in self.seqs]
        # doubt here
        if self.linear:
            self.restriction_finder_linear()
        #else:
        #    self.restriction_finder_circular()
        # self.restriction_finder()

    # gets the longest contiguous sequence in a list of sequences
    # pairwise alignment
    def longest_contiguous_sequence(self):
        init = overlaps(str(self.seqs[0].seq).lower(), str(self.seqs[1].seq).lower(), 25)
        ref = get_longest(str(self.seqs[1].seq).lower(), init)
        for s in self.seqs[2:]:
            contiguous = get_longest(str(s.seq).lower(), overlaps(ref, str(s.seq).lower(), 25))
            ref = contiguous
        return Dseqrecord(ref)

    def restriction_finder_linear(self):
        # firstly lets find all of the enzymes that cut through the contiguous sequence

        # this dictionary is made up from all of the re that cut through the contiguous sequence
        cuts_contiguous = self.re.search(self.contiguous.seq, linear=True)
        enzymes = []

        # to see if there are enzymes that cut more than once (produce two or more fragments)

        for re, sites in cuts_contiguous.items():
            # if it cuts
            if len(sites) > 0:
                # now for single cutters
                # test for single cutters on contiguous sequence
                if len(self.contiguous.cut(re)) == 2:
                    # test for fragment size
                    if self.has_size(self.contiguous.cut(re)):
                        enzymes.append(re)

        # until here we have all the enzymes that execute single cuts through the contiguous sequence
        # stored in enzymes

        # now enzymes that cut 2 or more times through all sequences' non-contiguous site
        # which can be acessed by the plasmid.insert attribute
        for s in self.plasmids:
                cuts_s = self.re.search(s.insert.seq, linear=True)
                # double cutters, and size verification
                cuts_s = [x for x in cuts_s if ((len(s.insert.cut(x)) > 2) and self.has_size(s.insert.cut(x), 3))]
                s.set_cutters(set(cuts_s))

        # see if there are enzymes that cuts through all
        # which would be the best case scenario as 1 enzyme would be enough, and allow to save money
        # as the plasmids are objects with no id i might as well store them as their index on self.plasmids
        enz_all = [x.cutters for x in self.plasmids]
        if len(enz_all[0].intersection(*enz_all[1:])) != 0:
            enzymes2 = enz_all[0].intersection(*enz_all[1:])

        # cause if there aren't:
        # we must search for a group of enzymes (up to the number of inputed sequences)
        # that are able to cut through all of the inserts
        else:
            # enzyme cuts mapping (indexes of plasmid objects in self.plasmids)
            mapping = {}
            for e in list(set().union(*enz_all)):
                mapping[e] = []
                for p in range(len(self.plasmids)):
                    if e in self.plasmids[p].cutters:
                        mapping[e].append(p)

            # now selecting the enzymes
            # how many sequences does the most versatile enzyme cut?
            mx = len(mapping[max(mapping, key=lambda x: len(mapping[x]))])
            # get the most versatiles
            most_versatile = [x for x in mapping.keys() if len(mapping[x]) == mx]
            # get the possible combinations of enzymes using the most versatile ones (save money)
            combos = []
            for e in most_versatile:
                # gets the plasmids which insert is not cut by the most versatile on each iteration
                missing = set(range(len(self.plasmids))) - set((mapping[e]))





        # if everything went well we now have all the enzymes that cut on the contiguous sequence just one time
        # and the enzymes that cut through all the sequences more than just once
        # now, how should i pick the enzymes for the selective digestion?
        # hipothesis: check all permutations of enzymes and enzymes2 and test for fragment size, with subsequent digest?
        # hipothesis: get a pair of enzymes that produce different fragments for different seqs
        # should i maximize the minimum difference(in bp) between the fragments produced?

        #possible_results = {}
        #for e1 in enzymes:
        #    for e2 in enzymes2:
        #        frags1 = self.seqs[0].cut([e1, e2])
        #        frags2 = self.seqs[3].cut([e1, e2])
        #        if len(frags1) != len(frags2):
        #            print(e1,e2)
        # return

    def has_size(self, frags, shrink=1):
        for f in frags:
            if len(str(f.seq)) < self.min_size/shrink:
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

# testing

if __name__ == '__main__':
    # as i will be testing with plasmids with over 2k bp might as well use a fasta file as input
    fileseqs = [Dseqrecord(seq, circular=False) for seq in SeqIO.parse('seqs.txt', 'fasta')]

    rf = Restriction_Finder(fileseqs, 300, Commercial)
    # insert sizes:
    # 446
    # 568
    # 2317
    # 500