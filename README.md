# Bioinformatics-Project-
Will store the code developed for an ongoing Bioinformatics Project and keep track of the changes made throughout the project development

Update History:

27/06/2017 
Restriction_Finder added to repository

It's a super preliminary approach, just to get around with the way the packages from pydna and biopython work.
Got some compatibility issues with the sequence representation, namely the Dseq format, not being able to use the  enzymes' catalyze method from Bio.Restriction package.
It should be pushed sooner, however I wasn't yet too capable with git.


27/06/2017 (2)
Restriction_Finder updated in repository

Compatibility issues all solved, now using Dseqrecord for all sequence representation, and cut method from Dseq instead of catalyse from bio.restriction package.
Enzyme import and sequence input were changed.
Above changes were made after meeting with professor Bjorn who helped me figure that out.
Algorithm not yet implemented, but it should begin to look better during the week.

27/06/2017 (3)
Restriction_Finder updated in repository

Algorithm already gets all single cutters for the contiguous sequence, that meet the minimum size criteria.
Also, all of the enzymes that are able to cut twice or more on all of the sequences.
Next step should be, choosing the enzyme pair that produces the best fragments for gel analysis.


27/06/2017 (4)
seqs.txt added to repository

File contains fasta sequences from cloning vectors that will be used for testing.