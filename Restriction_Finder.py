from pydna.common_sub_strings import common_sub_strings as overlaps
from operator import itemgetter as _itemgetter
from pydna.dseqrecord import Dseqrecord
from Bio import SeqIO
from Bio.Restriction import CommOnly as Commercial
from itertools import permutations as _permutations
# AllEnzymes is empty
# from Bio.Restriction import AllEnzymes


class Plasmid:
    # total is the complete sequence, contiguous is the contiguous part, and insert is the non-contiguous
    # cutters are a set of enzymes that cut through the plasmids' insert
    def __init__(self, Dseqrecord, contiguous):
        self.total = Dseqrecord
        self.contiguous = contiguous
        self.insert = self.get_insert()
        self.cutters = None

    def get_insert(self):
        # need to check if this is correct
        start = str(self.total.seq).lower().find(str(self.contiguous.seq).lower())
        end = start + len(self.contiguous.seq) - 1
        insert = str(self.total.seq)[end + 1::] + str(self.total.seq)[0:start]
        return Dseqrecord(insert)

    def set_cutters(self,cutters):
        self.cutters = cutters


class Restriction_Finder:

    def __init__(self, seqs, min_size=300, r_enzymes=Commercial, budget = True):
        self.seqs = seqs
        self.contiguous = self.longest_contiguous_sequence()
        self.re = r_enzymes
        self.min_size = min_size
        self.budget = budget
        self.linear = self.seqs[0].linear
        self.plasmids = [Plasmid(x, self.contiguous) for x in self.seqs]
        # doubt here
        if self.linear:
            self.restriction_finder()

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

    def warning(self):
        if self.budget:
            print('It was not possible to find a suitable combination of enzymes to perform a diagnostic digest'
                  'Try setting budget = False')
        else:
            print('It was not possible to find a suitable combination of enzymes to perform a diagnostic digest')


    def restriction_finder(self):
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

        # !!!!!!
        # could i include a parameter here that could select the else approach if the results aren't good for
        # sequence differentiation?
        # !!!!!!

        enz_all = [x.cutters for x in self.plasmids]

        # budget option --> tries to obtain better results using more enzymes

        if len(enz_all[0].intersection(*enz_all[1:])) != 0 and self.budget:
            enzymes2 = enz_all[0].intersection(*enz_all[1:])

        # cause if there aren't:
        # we must search for a group of enzymes (up to the number of inputed sequences)
        # that are able to cut through all of the inserts
        else:
            # enzyme cuts mapping (indexes of plasmid objects in self.plasmids)
            # dictionary
            mapping = {}
            for e in list(set().union(*enz_all)):
                mapping[e] = []
                for p in range(len(self.plasmids)):
                    if e in self.plasmids[p].cutters:
                        mapping[e].append(p)
            # now selecting the enzymes

            # how many sequences does the most versatile enzyme cut?
            mx = len(mapping[max(mapping, key=lambda x: len(mapping[x]))])

            # get the most versatile ones // the lesser versatile ones
            most_versatile = [x for x in mapping.keys() if len(mapping[x]) == mx]
            # if mx = 1 the most versatile enzymes will be the same as the less versatile
            if mx != 1:
                less_versatile = [x for x in mapping.keys() if x not in most_versatile]

            # get the possible combinations of enzymes using the most versatile ones (save money)
            # how many enzymes do I need in the worst case scenario?
            # will be increasing as we try to cut all inserts with the minimum number of enzymes
            if mx != 1:
                e_count = 1
                wcs = len(self.plasmids) - len(mapping[most_versatile[0]])
            else:
                # this would be the worst case scenario and in that case we won't need to loop in an ascending
                # number of enzymes, so we can jump right to the case where we need 1 enzyme to each plasmid
                e_count = len(self.plasmids)
                wcs = len(self.plasmids)

            # break while
            found = False
            # will store the combinations of enzymes that provide the wanted results , for later analysis
            # also doesn't need to be redefined as it will break the cycle the first time it is filled
            combos = []

            while e_count <= wcs and not found:
                # for the cases where there are no enzymes cutting more than one plasmid
                if mx != 1:
                    # all possible combination of e_count enzymes
                    if e_count == 1:
                        # for the 2 enzyme scenario we don't need permutations
                        # so there will be a couple of steps that will be different and beacuse of that
                        # will be preceded by the same if clause used here
                        # if e_count is 1 we don't need permutations
                        perm = less_versatile
                    else:
                        perm = _permutations(less_versatile, e_count)

                    # we will be picking one of the most versatile and try to pair it with another group
                    # to see if we can cut through all plasmids' insert
                    for e in most_versatile:
                        missing = list(set(range(len(self.plasmids))) - set((mapping[e])))
                        for p in perm:
                            # for the if clause p is an enzyme
                            if e_count == 1:
                                # get the enzymes that cut through the missing plasmids
                                fill_gap = [x for x in mapping if sorted(mapping[x]) == sorted(missing)]
                                if len(fill_gap) > 0:
                                    for e2 in fill_gap:
                                        combi = (e, e2)
                                        if combi not in combos:
                                            combos.append(combi)

                            else:
                                # for the else clause p is a combinationn of enzymes
                                p_cuts = []
                                # p_cuts will contain all of the plasmids cut by the combination of enzymes
                                for e2 in p:
                                    p_cuts.extend(mapping[e2])
                                if sorted(p_cuts) == sorted(missing):
                                    combi = tuple([e] + [x for x in p])
                                    if combi not in combos:
                                        combos.append(combi)

                    if len(combos) != 0:
                        # if combos has enzymes then we don't need to search further
                        found = True

                    else:
                        e_count += 1

                # for the odd case where we need one enzyme for each
                else:
                    perm = _permutations(most_versatile, e_count)
                    for p in perm:
                        p_cuts = []
                        for e2 in p:
                            p_cuts.extend(mapping[e2])
                        # as in this case all of the plasmids are missing cuts, we need to compare with all plasmids
                        if sorted(p_cuts) == list(range(len(self.plasmids))):
                            combi = [x for x in p]
                            if combi not in combos:
                                combos.append(combi)
                    if len(combos) != 0:
                        # if combos has enzymes then we don't need to search further
                        found = True
                    else:
                        # ver melhor isto
                        return self.warning()

            enzymes2 = combos

        # if everything went well we now have all the enzymes that cut on the contiguous sequence just once (enzymes)
        # and the enzymes that cut through all the sequences more than just once (enzymes2)
        # now, how should i pick the enzymes for the selective digestion?
        self.best_set(enzymes, enzymes2)


    def best_set(self, enzymes1, enzymes2):
        # this function should retrieve the best combination of enzymes for gel analysis
        # hipothesis: check all permutations of enzymes and enzymes2 and test for fragment size, with subsequent digest?
        # hipothesis: get a pair of enzymes that produce different fragments for different seqs
        # should i maximize the minimum difference(in bp) between the fragments produced?

        # i guess we should only select re from enzymes that do not feature in enzymes2
        enzymes_f = [x for x in enzymes1 if x not in enzymes2]
        pass

    def has_size(self, frags, shrink=1):
        for f in frags:
            if len(str(f.seq)) < int(self.min_size/shrink):
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

    rf = Restriction_Finder(fileseqs, 300, Commercial, budget=False)
    # insert sizes:
    # 446
    # 568
    # 2317
    # 500
    # amo-te sabias?