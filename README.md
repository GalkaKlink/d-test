# d-test
phylogeny-based test for detecting changes in amino acid preferences
References:
1) https://www.biorxiv.org/content/10.1101/589572v2 - this article contains a precise description of the test and its application to Env gene of HIV-1
Brief description:   
The test finds amino acids, substitutions to which are biased towards or away from any chosen focal nodes (external or internal) of the phylogenetic tree. The significance of an amino acid X to be “proximal” or “distal” for a node N is calculated by comparing the mean evolutionary distance (i.e. the way between two points of the phylogenetic tree) between N and substitutions to X with the null-distribution obtained by calculating evolutionary distances between N and substitutions to any amino acids accounting for ancestral amino acids for X.

For each amino acid X and node N, the program also calculates the effect size of distance bias by calculating a number of standard deviations which lie between the real distance and the null-distribution mean (i.e. standard score, or z-score). 

If substitutions to X are biased towards nodes N1 and N2 with different polarity (i.e. X is proximal for one of them and distal for the other), X is considered to have different preferences in phylogenetic vicinity of N1 and N2. 

Applying the method:
Input files: to apply d-test, you need a rooted phylogenetic tree in *.newick format and corresponding amino acid alignment which includes reconstructed ancestral states for internal nodes of the tree in *.fasta format. The names of all internal and external nodes should correspond to each other in alignment and tree.  

As all protein sites are considered independently, the program may treat them in parallel as distinct tasks. For these your machine should be able to proceed arrays of tasks. It will accelerate the calculations substantially, but  the version of script for iterative processing of each site also exists.
Control file FitnessShift.Config.txt:
Contains several lines to be filled by the user:
input=alignment.fasta – path to input fasta file
tree_file=tree.newick – path to input newick file
root_id=node1 – name of the root node of the tree
focal_clades=A,B,C – names of tree regions you are interested in to better orient in outputs
focal_points=node100,node200,node300 – list of focal nodes
significance_threshold=0.01 – significance p-value threshold below which amino acid is considered proximal or distal for the node 
output=Results.txt – final output file


Description of scripts:
1)make_file.pl – makes one-site alignments from the initial *.fasta file. This step is needed only if you are going to make parallel computations.
Output – files *.site for each site of the alignment.
2)FitnessShift.pl – for each amino acid in each protein site, calculates p-value of being “proximal” or “distal” for each focal node as well as corresponding z-scores. The default version is for running as an array of tasks, and the version for full alignment is named “FitnessShift.full_alignment.pl”.

External branches that are longer than 10 medians of external branch length are excluded from the analysis as suspiciously long. You can cut easily find this step in the script and cut it if you do not want it to be done.

Outputs: 
1 - Site.dists_fromCladename_FocalSpecies – for each Site, outputs phylogenetic distances from FocalSpecies to all substitutions in this site. Each line contains:
site ancestral_amino_acid derived_amino_acid distances. Cladename is a name of a part of the phylogenetic tree that corresponds to FocalSpecies (this information is taken from, config file). 
2 – Site.Cladename_FocalSpecies.p_value.10000 – for each amino acid (column1) in Site, outputs p-value of being proximal (column2) or distal (column4) for FocalSpecies. Columns 3 and 5 represent minimal p-values that can be obtained from the null distribution for this amino acid in this site (see paper for more details).  
3 - Site.Cladename_FocalSpecies.z_score.10000 – for each amino acid in Site, this file contains standardized score (z-score) that shows how far mean distance between FocalSpecies and substitutions to this amino acid is from null expectation. It can be shorter (negative values of the score) or further (positive values of the score) than it is expected. File includes 5 columns: amino acid; mean distance between FocalSpecies and substitutions to this amino acid; mean distance between FocalSpecies and substitutions to this amino acid, expected from the null distribution; standard deviation of the null distribution; z-score.
3)ResultingTable.pl – using described intermediate outputs, it makes one resulting table that contains information about amino acids that turned out to be proximal or distal for at least one focal node in a site. The table contains columns:
site – site; AA_number_in_site - number of amino acids in a site; AA – amino acid;
after it, for each focal node, it contains two columns: type of an amino acid for the focal node (“proximal”, “distal” or “none”) and corresponding z-score. 

