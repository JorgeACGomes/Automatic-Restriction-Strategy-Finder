![alt text](https://upload.wikimedia.org/wikipedia/commons/thumb/9/93/EEUMLOGO.png/200px-EEUMLOGO.png)

# Bioinformatics Project 

###### Project Supervisor: @BjornFJohansson

## Automatic restriction strategy finder for Synthetic Biology Constructs

Will store the code developed for an ongoing Bioinformatics Project and keep track of the changes made throughout the project development.
The main goal of this project is to design an algorithm for the automatic selection of the most effective restriction enzymes for DNA analysis based on user defined criteria. With those enzymes we will be able to perform a diagnostic digest over n user inputted sequences.
The algorithm will then be implemented in Python and embedded within Pydna software package, or make use of some of its resources.

Will also include the preliminary report, which consists of a more theoretical review on Synthetic Biology, DNA cloning and restriction, and the final report
which will add the algorithm explanation, the idea behind the chosen approach, useful explanations, usage examples and also results obtained with chosen testcases.

Pydna repository: https://github.com/BjornFJohansson/pydna



###### Changelog:

###### 27/06/2017 

Restriction_Finder added to repository

It's a super preliminary approach, just to get around with the way the packages from pydna and biopython work.
Got some compatibility issues with the sequence representation, namely the Dseq format, not being able to use the  enzymes' catalyze method from Bio.Restriction package.
It should be pushed sooner, however I wasn't yet too capable with git.


###### 27/06/2017 (2)

Restriction_Finder updated in repository

Compatibility issues all solved, now using Dseqrecord for all sequence representation, and cut method from Dseq instead of catalyse from bio.restriction package.
Enzyme import and sequence input were changed.
Above changes were made after meeting with professor Bjorn who helped me figure that out.
Algorithm not yet implemented, but it should begin to look better during the week.

###### 27/06/2017 (3)

Restriction_Finder updated in repository

Algorithm already gets all single cutters for the contiguous sequence, that meet the minimum size criteria.
Also, all of the enzymes that are able to cut twice or more on all of the sequences.
Next step should be, choosing the enzyme pair that produces the best fragments for gel analysis.


###### 27/06/2017 (4)

seqs.txt added to repository

File contains fasta sequences from cloning vectors found in genbank that will be used for testing.


###### 28/06/2017

Restriction_Finder updated in repository

Plasmid class has been added to make the contiguous/non-contiguous sequences manipulation easier.
Several changes to enzyme search algorithm.

###### 29/06/2017

Restriction_Finder updated in repository

Enzyme selection algorithm already works, now i just need to get the set that produces the best results, for a future gel analysis.

###### 29/06/2017(2)

Preliminary_report.pdf added to repository
seqs_vegas.txt added to repository

As suggested by professor Bjorn I will also include the written reports made for the project.
As I needed a different testing set i included the seqs_vegas.txt which includes sequences from 
yeast plasmid cloning experiments, retrieved from:

Mitchell, Leslie A., et al. "Versatile genetic assembly system (VEGAS) to assemble pathways for expression in S. cerevisiae." Nucleic acids research 43.13 (2015): 6620-6630.

###### 01/07/2017

Example.ipynb added to repository

Restriction_Finder updated in repository

Added a jupyter/IPython notebook on how to use Restriction_Finder.
Restriction_Finder is now working properly, and in its final form, for now.
It's lacking documenting, which will be added in the next days. 
However the code is heavily commented and can be easily understood without proper documentation.

###### 03/07/2017

Restriction_Finder updated in repository
final_report added to repository

Some minor corrections were introduced in Restriction_Finder. More warnings added. Documentation still lacking.
Final report was added. Will be updated in the next days.


###### 04/07/2017

Restriction_Finder updated in repository
final_report updated in repository
Example.ipynb renamed to Restriction_Finder_Guide_.ipynb
Algorithm_Pipeline added to repository

Restriction_Finder is now completely documented, all_best function added, to retrieve all of the best results when the best possible result is an isoschizomer.
final_report underwent some corrections, according to final version delivered to the Professors responsible for the project.
IPython notebook file was updated to match latest changes in the code. Renamed for convenience.
Graphical pipeline of the algorithm was added for better understanding at algorithm's main steps.


### Final report includes:
- *methodologies*
- *algorithm design and implementation*
- *test building*
- *Test results*
- *Results Discussion*
- *Main conclusions on algorithm's strengths and weaknesses*

### Status:

- [x] Finished
- [x] Tested
- [x] Guide Written
- [x] Fully Functional