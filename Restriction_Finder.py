from pydna.common_sub_strings import common_sub_strings as overlaps
from operator import itemgetter as _itemgetter
from itertools import permutations as _permutations
from itertools import chain as _chain
from collections import Counter
from pydna.dseqrecord import Dseqrecord
from pydna.gel import Gel
from pydna.gel import weight_standard_sample
from Bio import SeqIO
from Bio.Restriction import CommOnly as Commercial
from Bio.Restriction import NonComm as NonComm
from random import choice,sample

# this is a workaround for the AllEnzymes bug , which is recognized as a set in Python 3 instead
# of a restrictionBatch object
All = NonComm
All.update(Commercial)


class Plasmid:
    
    """Class to store information on each sequence considering its plasmid
    configuration namely the presence of a contiguous sequence shared with 
    other sequences but also a different region commonly named as insert.
    Cutters attribute stores which enzymes do cut on insert sequence
    Further documenting will not be made since it is not supposed to be used
    isolated from Restriction_Finder, and as a dependency of it"""
    
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

    def __init__(self, seqs, min_size=150, r_enzymes=Commercial, optim=True, gel=True, iso=True):
        self.seqs = seqs
        
        # verifications for sequences
        if len(self.seqs) < 2:
            warning('seq_num')
            return
        for seq in self.seqs:
            if type(seq) != Dseqrecord:
                warning('Dseqrecord')
                return
        lens = [len(x.seq) for x in self.seqs]
        if min(lens) < 750:
            warning('seq_size')
            return
        
        self.contiguous = self.longest_contiguous_sequence()
        self.results = None
        self.best = None

        if r_enzymes not in [All, Commercial]:
            self.re = Commercial
            warning('enzymes')
        else:
            if r_enzymes == All:
                warning('enzymes2')
            self.re = r_enzymes
            
        if type(min_size) not in [int, float]:
            self.min_size = 150
            warning('min_size')
        else:
            self.min_size = min_size
            
        if type(optim) is not bool:
            self.optim = True
            warning('optim')
        else:
            if optim is False:
                warning('optim2')
            self.optim = optim
            
        if type(gel) is not bool:
            self.gel = True
            warning('gel')
        else:    
            self.gel = gel
            
        if type(iso) is not bool:
            self.iso = True
            warning('iso')
        else:
            if iso == False:
                warning('iso2')
            self.iso = iso
                
        self.plasmids = [Plasmid(x, self.contiguous) for x in self.seqs]
        self.restriction_finder()


    def longest_contiguous_sequence(self):
        
        """gets the longest contiguous sequence from the list of sequences in
        the Restriction_Finder object. Contiguous sequence was found using 
        pairwise alignment"""
        
        if len(self.seqs[0].seq)>50 :
            lim = 25
        else:
            lim = int(len(self.seqs[0].seq)/2)
        try:
            init = overlaps(str(self.seqs[0].seq).lower(), str(self.seqs[1].seq).lower(), lim)
            ref = get_longest(str(self.seqs[1].seq).lower(), init)
            for s in self.seqs[2:]:
                contiguous = get_longest(str(s.seq).lower(), overlaps(ref, str(s.seq).lower(), lim))
                ref = contiguous
            return Dseqrecord(ref)
        except ValueError:
            return warning('similarity')

    def restriction_finder(self):
        
        """Initial part of the algorithm that searches 1st for single cutters 
        on the contiguous area, and 2nd enzymes that cut 2 or more times on the non
        shared sequence. Most of the function complies the second search since
        it must search for every combination of possible enzymes, covering every
        possible case. If it finds suitable enzymes they're sent to best_set.
        This function is heavily commented to make up for its extension and high
        possibility of getting lost while trying to understand it"""
        
        # firstly lets find all of the enzymes that cut through the contiguous sequence

        # this dictionary is made up from all of the re that cut through the contiguous sequence
        cuts_contiguous = self.re.search(self.contiguous.seq, linear=True)
        enzymes1 = []

        # to see if there are enzymes that cut more than once (produce two or more fragments)

        for re, sites in cuts_contiguous.items():
            # if it cuts
            try:
                if len(sites) > 0:
                    # now for single cutters
                    # test for single cutters on contiguous sequence
                    # as it is set as two fragments, will work on both linear and circular
                    if len(self.contiguous.cut(re)) == 2:
                        # test for fragment size
                        if self.has_size(self.contiguous.cut(re)):
                            enzymes1.append(re)
            except TypeError:
                # this block was made to avoid some enzymes that cause issues when using cut
                pass

        # until here we have all the enzymes that execute single cuts through the contiguous sequence
        # stored in enzymes

        # now enzymes that cut 2 or more times through all sequences' non-contiguous site
        # which can be acessed by the plasmid.insert attribute
        for s in self.plasmids:
                cuts_s = self.re.search(s.insert.seq, linear=s.total.linear)
                # double cutters or more
                cuts_s = [x for x, y in cuts_s.items() if len(y) >= 2]
                s.set_cutters(set(cuts_s))

        # see if there are enzymes that cuts through all
        # which would be the best case scenario as 1 enzyme would be enough, and allow to save money
        # as the plasmids are objects with no id i might as well store them as their index on self.plasmids

        enz_all = [x.cutters for x in self.plasmids]

        if len(enz_all[0].intersection(*enz_all[1:])) != 0:
            enzymes2 = list(enz_all[0].intersection(*enz_all[1:]))

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

            # if mx = 1 then the most versatile enzymes will be the same as the less versatile
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
                        return warning('final')

            enzymes2 = combos

        # if everything went well we now have all the enzymes that cut on the contiguous sequence just once (enzymes)
        # and the enzymes that cut through all the sequences more than just once (enzymes2)
        # now, we should pick the enzymes that provide the best results

        # print(enzymes1, enzymes2)
        self.best_set(enzymes1, enzymes2)

    def best_set(self, enzymes1, enzymes2):
        
        """this function should retrieve the best combination of enzymes for 
        gel analysis. Since we are trying to minimize the costs let's see which
        re's from enzymes1 are present in enzymes 2 and if possible try to 
        solve it only with enzymes that are shared between them. As previously 
        done, all of the possible cases are covered up, thus explaining the 
        extension of the code. Enzymes will be tested for fragment size and
        optionally gel bands quality."""

        # the approach will be different if enzymes2 comes without the mapping step
        if type(enzymes2[0]) == tuple:
            no_map = False
        else:
            no_map = True

        if no_map:
            common = [x for x in enzymes1 if x in enzymes2]
            uncommon = [x for x in enzymes1 if x not in common]
        else:
            common = []
            for tup in enzymes2:
                for e in tup:
                    if e in enzymes1:
                        common.append(tup)
            uncommon = [x for x in enzymes1 if x not in common]

        # since we now have the common ones lets see what results we get from using those
        res = []
        for e in common:
            frags = []  # list of lists of fragments
            for p in self.plasmids:
                try:
                    f = p.total.cut(*[e])  # all fragments for each plasmid // e can be a tuple
                    proceed = self.has_size(f)  # for each plasmid check if all frags meet size criteria
                    if proceed:
                        frags.append(f)
                    else:  # if any of the plasmids cut by e enzyme does not meet the size criteria, e is not considered
                        break
                except (TypeError, IndexError) as error:
                    # typeError can be raised when using some enzymes of AllEnzymes
                    # IndexError is raised sometimes by cut method from Dseq, reason Unknown
                    pass

            if len(frags) == len(self.plasmids):
                if bands_differ(frags):
                    res.append(e)  # if the enzyme fills all of the requirements it can be used

        if len(res) != 0:
            # If we already have enzymes in res we don't need to search any further
            self.results = res
            if self.optim:
                return self.optimize_bands(res)
            else:
                return self.to_gel(choice(res))

        else:
            # If it reaches this point we will try to get the best possible results with enzymes
            # that weren't used before so we will be mixing between uncommon and re's from
            # enzymes2
            last_chance = []
            for e in uncommon:
                for etup in enzymes2:
                    if type(etup) == tuple:
                        s = [e]+list(etup)
                    else:
                        s = [e, etup]

                    if sorted(s) not in last_chance:
                        last_chance.append(sorted(s))
            
            # we could only take a sample since this list could get really big
            # and make the program run for days, although it takes out some 
            # of the results reproduceability, it is the only way to make it
            # work in useful time
            try:
                last_chance = sample(last_chance,50)
            except ValueError:
                last_chance = last_chance
            for e in last_chance:
                frags = []
                for p in self.plasmids:
                    try:
                        f = p.total.cut(*[e])
                        proceed = self.has_size(f)
                        if proceed:
                            frags.append(f)
                        else:
                            break
                    except (TypeError, IndexError) as error:
                        pass

        if len(res) != 0:
            self.results = res
            if self.optim:
                return self.optimize_bands(res)
            else:
                # may not produce the best result but saves time
                return self.to_gel(choice(res))
        else:
            return warning('min_size3')

    def optimize_bands(self, res):
        
        """This function receives a list of enzymes and retrieves the gel result
        for the enzyme that produces the maximum of the minimum differences between
        each band, to facilitate the gel analysis
        res : list
        a list of enzymes, or a list of lists(or tuples) each with a set of enzymes"""
        
        if len(res) > 1:
            diffs = {}
            for e in res:  # for each enzyme given in res
                lsfrags = [p.total.cut(*[e]) for p in self.plasmids]
                flat_frags = _chain(*lsfrags)
                dif = 50000  # arbitrary value
                for band in flat_frags:
                    for band2 in flat_frags:
                        if band.seq != band2.seq:
                            l_dif = abs(len(band.seq)-len(band2.seq))
                            if l_dif < dif:
                                dif = l_dif
                diffs[e] = l_dif
            
            best_possible = max(diffs.items(), key=_itemgetter(1))[0]
            if self.iso is False:
                val_counter = Counter(diffs.values())
                best_unique = [(k,v) for k, v in diffs.items() if val_counter[v] == 1]
                if len(best_unique) != 0:
                    return self.to_gel(max(best_unique, key=_itemgetter(1))[0])
                else:
                    warning('enzymes4')

            return self.to_gel(best_possible)
        else:
            return self.to_gel(res[0])

    def to_gel(self, enzymes):
        """Executes the gel, or prints out the solutions 
        for the diagnostic digest"""
        
        self.best = enzymes
        if self.gel:
            lsfrags = [p.total.cut(*[enzymes]) for p in self.plasmids]
            st = weight_standard_sample('1kb+_GeneRuler')
            lsfrags.insert(0, st)
            print('Gel obtained using', enzymes)
            Gel(lsfrags, gel_len=16, wellsep = 8).run()
            # self.solutions()

        else:
            self.solutions()

    def has_size(self, frags):
        """Size verification for all frags"""
        for f in frags:
            if len(str(f.seq)) < self.min_size:
                return False
        return True
    
    def lanes(self):
        """Retrieves a list with sequence IDs, 
        in the same order as they appear on the gel"""
        Ids = []
        counter = 1
        if self.plasmids[0].total.name == 'name?':
            warning('Ids')
        for seq in self.plasmids:
            if seq.total.name == 'name?':
                Ids.append('Seq '+str(counter))
                counter += 1
            else:
                Ids.append(seq.total.name)
                counter += 1
        print('Lane order:\n', Ids)
    
    def solutions(self):
        """Prints the solutions obtained for the diagnostic digest"""
        if self.results is None:
            warning('results')
            return 
        if self.are_isoschizomers():
            iso = True
        else:
            iso = False
        if type(self.best) == list:
            print('Results obtained using the following group of enzymes:')
            for e in self.best:
                print(e, sep=',')
            print('Other possible solutions:')
            for el in self.results:
                if el != self.best:
                    print(el)
        else:
            if iso:
                warning('enzyme3')
            print('Best enzyme:\n', self.best)
            print('Other possible solutions:')
            b = self.results.index(self.best)
            print(self.results[0:b]+ self.results[b+1:])

        return self.results

    def are_isoschizomers(self):
        """Defines if any of the results are isoschizomers"""
        if self.results is None:
            warning('results')
            return
        for e in self.results:
            if self.best in e.isoschizomers():
                return True
        return False


# -------------- auxiliary functions --------------------- #


def get_longest(seq, out_list):
    """gets the tuple with the maximum value for the 2nd index
    and slices the sequence with those given parameters"""
    
    coords = max(out_list, key=_itemgetter(2))
    start = coords[1]
    end = coords[1] + coords[2]
    return seq[start:end]



def bands_differ(lsfrags):
    """given a list of lists of fragments it sees if there's 
    a unique band in each list of fragments"""
    # lsfrags is a list made of fragments for each plasmid, 
    # so len(lsfrags) == len(self.plasmids)
    for p in range(len(lsfrags)):
        others = [len(x.seq) for x in _chain(*(lsfrags[:p] + lsfrags[(p+1):]))]  # list of lists to single list
        count = 0
        for band in lsfrags[p]:
            if len(band.seq) not in others:
                count += 1
        if count == 0:
            return False
    return True


def warning(code):
    """for warnings regarding inputs and results"""
    
    abort = ['Dseqrecord', 'Seq_Size','seq_num']
    d = {'enzyme': 'The introduced enzyme set is not recognized. Using Commercial\n',
         'enzyme2': 'Using all of the known enzymes, some of which may not be available commercially\n, use Commercial instead\n',
         'enzyme3': 'The best enzyme is an isoschizomer so if you run Restriction_Finder() multiple times you will '
                    'get different best enzymes but same digestion result\n',
        'enzymes4': 'All of the possible results are isoschizomers, using the best one\n',
         'min_size': 'Introduced minimun fragment size is not numeric. Using 150bp\n',
         'min_size2': 'Minimum fragment size is too big, results may not exist, try lowering it next time\n',
         'min_size3': 'It was not possible to find a suitable combination of enzymes to perform a diagnostic digest.\n'
                      'A lower minimum fragment size is recommended\n',
         'gel': 'Gel parameter must be either True or False. Using True\n',
         'optim': 'Optim parameter must be either True or False. Using True\n',
         'optim2': 'Band optimization is disabled, results will be faster, although gel readibility may be worse\n',
         'final': 'It was not possible to find a suitable combination of enzymes to perform a diagnostic digest\n',
         'Ids': 'Sequences were not identified so consider the same order of input\n',
         'iso': 'Isoschizomers parameter must be either True or False. Using True\n',
         'iso2': 'Isoschizomers parameter set to False. Unless all of the possible results are isoschizomers, '
                 'results will not be optimal\n',
         'similarity': 'Sequences do not share big similarity, which diverges from the objective of this program\n',
         'seq_size': 'Inputted sequences are too small and thus will not provide visible bands in an agarose gel\n',
         'Dseqrecord': 'At least one of the inputted seqs is not a Dseqrecords\n',
         'seq_num': 'A minimum of 2 sequences is required\n',
         'results': 'There were no results for the given parameters and sequences\n'
         }
    print(d[code])
    if code in abort:
        print('Program will now abort')
        return
    
# testing

if __name__ == '__main__':
    fileseqs = [Dseqrecord(seq, circular=True) for seq in SeqIO.parse('seqs_vegas.txt', 'fasta')]
    import time
    ti = time.time()
    rf = Restriction_Finder(fileseqs, 150, Commercial, optim=False, iso=True)
    #for p in rf.plasmids:
    #   print(len(p.insert.seq))
    #print(len(rf.contiguous.seq))
    #from numpy import mean 
    #print(mean([len(x.total) for x in rf.plasmids]))
    rf.solutions()
    #rf.to_gel(BcoDI)
    #rf.lanes()
    print(rf.are_isoschizomers())
    rf.lanes()
    #sol = rf.solutions()
    #fileseqs[0].cut(*rf.results)
    print(round(time.time()-ti, 3))
