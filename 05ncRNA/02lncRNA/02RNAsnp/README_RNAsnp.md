RNAsnp v1.1
===========

CONTENTS:
--------
	1) Introduction	
	2) Installation	
	3) Usage
	4) Input formats
	5) Examples
 	6) References
	7) Contact

1) Introduction
---------------

The program RNAsnp helps for an effective detection of local RNA secondary
structural changes induced by the SNPs. RNAsnp starts with the computation of
base pair probabilities for wild-type and mutant sequence (with SNP) by using
the RNAfold or RNAplfold tools from Vienna RNA package[1,2]. Further, the
structural difference between the wild-type and mutant is assessed by computing
the Pearson correlation coefficient and the Euclidean distance for all sequence
intervals. The sequence interval with maximum distance or minimum correlation is
finally reported with its p-value. This empirical p-value is computed using the
background distribution given the sequence length, G+C content and SNP position. 
	
2) Installation
--------------- 
The package requires the GNU gengetopt to be installed on your machine. If so,
then follow the steps to make the executable, 

./configure
make
make install 

RNAsnp requires an environment variable named RNASNPPATH to run. The location of 
the RNAsnp-1.2 directory needs to be assigned to the RNASNPPATH variable. 
For example in the bash terminal, this can be assigned as,

export RNASNPPATH='<PATH>/RNAsnp-1.2'. 

You can add this line to .bashrc file available in the home directory to avoid
executing the above command every time to start with RNAsnp on new terminal. 

3) Usage 
--------
Usage: RNAsnp -f <seq_file> -s <snp_file> [options]

  -h, --help                    Print help and exit
      --detailed-help           Print help, including all details and hidden
                                  options, and exit
      --full-help               Print help, including hidden options, and exit
  -V, --version                 Print version and exit

Input Options:
  -f, --seq=STRING              File containing the input sequence

  The single input sequence can be provided either in fasta format or linear
  sequence without any gaps

  -s, --snp=STRING              File containing the list of SNP

  The list of SNPs to be tested have to be provided in separate lines, see
  input format section below for more details

  -m, --mode=INT                Select the mode of operation (default=`1')

                                  1 - perform global folding by using RNAfold
                                  and compute the difference in base pair
                                  probabilities for all sequence intervals

                                  2 - perform local folding by using RNAplfold
                                  and compute the difference in base pair
                                  probabilities for all sequence intervals of
                                  fixed length

                                  3 - screen putative structure-disruptive SNPs
                                  in an RNA sequence

  
  Mode 1 is designed to predict the effect of SNPs on short RNA sequences
  (i.e., -w parameter is less than or equal to 500), where the base pair
  probabilities of the wild-type and mutant RNA sequences are calculated using
  the global folding method RNAfold. The structural difference between
  wild-type and mutant is computed using Euclidean distance and Pearson
  correlation measures for all sequence intervals (with minimum size of 50,
  -l). Finally, the interval with maximum base pair distance or minimum
  correlation coefficient and the corresponding p-value is reported.

   Mode 2 is designed to predict the effect of SNPs on large RNA sequence.
  Here, the base pair probabilities are calculated using the local folding
  method RNAplfold (with -W 200 and -L 120 options). As a first step, the
  structural difference is calculated using the Euclidean distance measure for
  all sequence intervals of fixed window length (default: 20, -X) and allowing
  the bases within the window can pair up to a distance of 120 (i.e. the
  maximal span of a base pair, -Y). In the second step, the sequence interval
  [u, v] with maximum base pair distance is selected to re-compute the
  difference for all internal local intervals that starting at u. Finally, the
  interval with maximum base pair distance and the corresponding p-value is
  reported.

   Mode 3, the combination of modes 1 and 2, is designed to screen all possible
  structure-disruptive SNPs in an input sequence using a brute-force approach.
  First, Mode 2 is applied to evaluate the SNP effect for all possible
  substitutions at every nucleotide position. Second, the SNPs with p-value
  less than 0.4 (--pvalue1) are subjected to Mode 1 to re-compute the structure
  effect using a global folding approach. The SNPs that have significant local
  structural effect (p-value less than 0.1, --pvalue2) are finally reported.

  -w, --winsizeFold=INT         length of flanking sequence on either side of
                                  SNP to fold  (default=`200')
  
  By default the program uses +/- 200nts around the SNP position to compute the
  base pair probabilities in all the three modes. This default value can be
  changed between 100 and 800 (inclusive) in multiples of 50 for Mode 1, and
  between 200 and 800 (inclusive) in multiples of 50 for Mode 2 and 3. In order
  to achieve this, however, please make sure that the input sequence is at
  least twice the size of chosen flanking. This restriction is necessary to
  keep the size of parameter tables for the p-value calculations manageable.

   In case the input sequence is less than twice the size of chosen flanking,
  the RNAsnp takes the nts up to the start and end position of the given
  sequence from the SNP position and perform the analysis. However, in this
  case the reporting p-value is not accurate since the input sequence length
  does not match the sequence length available in the pre-computed parameter
  tables.

Additional parameters:
  Please note that the precomputed background scores, which RNAsnp uses to
  estimate p-value, are based on the default value assigned to the following
  parameters. Thus, if the default value is changed for any of the following
  parameters  (except --pvalue1 and --pvalue2), then the reporting p-value is
  not accurate.

  -c, --cutoff=FLOAT            minimum cut-off for the base pair
                                  probabilities. This parameter is applicable
                                  to both Mode 1 and 2  (default=`0.01')

  Base pair probabilities that are above this cut-off are only considered to
  compute the Euclidean distance or correlation coefficient between wild-type
  and mutant.


  Parameters associated with mode -M 1:
  -l, --minLen=INT              minimum length of the sequence interval
                                  (default=`50')

  The structural difference between wild-type and mutant is computed for all
  sequence intervals with the selected minimum length

  Parameters associated with mode -M 2:
  -W, --winsize=INT             Average the pair probabilities over windows of
                                  given size  (default=`200')
  -L, --span=INT                Set the maximum allowed separation of a base
                                  pair to span. i.e. no pairs (i,j) with j-i >
                                  L will be allowed.  (default=`120')
  -X, --regionX=INT             Length of the local structural element that we
                                  expect to have an effect  (default=`20')
  -Y, --regionY=INT             Length of the interval over which the local
                                  structural changes are evaluated, i.e., the
                                  maximal span of a base pair  (default=`120')
  
  The functions of each of these parameters are mentioned in the description of
  mode 2 shown above

  Parameters associated with mode -M 3:
      --pvalue1=FLOAT           p-value threshold to filter SNPs that are
                                  predicted using Mode 2  (default=`0.4')
      --pvalue2=FLOAT           p-value threshold to filter SNPs that are
                                  predicted using Mode 1  (default=`0.1')
  -e, --winsizeExt=INT          size of the flanking region on either side of
                                  SNP that includes the local window returned
                                  by Mode 2. This subsequence is then passed to
                                  Mode 1 for re-computation  (default=`200')

Additional option to compute edist:
  -E, --edist=INT               compute ensemble Euclidean distance between the
                                  distribution of structures between two
                                  sequences  (default=`0')
  -C, --boltzmannPreFactor=DOUBLE
                                Multiply the Boltzmann factor with a prefactor
                                  alpha  (default=`1')

4) Input formats 
----------------
Sequence file must contain one sequence (preferably in FASTA format). 

SNP file must contain the list of SNPs that are given in separate lines. The
SNPs are described as, wild-type nucleotide followed by nucleotide position
followed by mutant nucleotide. In case of multiple SNPs, the SNPs are delimited
by the special character "-". 

Example formats:
Single SNP: A201G
where, A is the wild-type nucleotide in the given sequence, 201 is the sequence
position of wild-type nucleotide and G is the mutant (or SNP). 

Multiple SNPs: A201G-U257A-C260G
The multiple SNPs (which occurs together) are defined next to each other with
the delimiter "-" between them. 

5) Examples 
-----------
-------------------------------------------------------------------
RNAsnp mode 1:
-------------------------------------------------------------------
1) Test for the effect of single SNP with RNAsnp default mode -m 1

$ RNAsnp -f examples/seq1.txt -s examples/snp1.txt
SNP	W   Slen   GC	interval d_max	p-value	interval r_min	p-value
U1013C	200 3344 0.5411	975-1025 0.2432	0.0724	998-1052 0.0615	0.0932

2) Test for the effect of multiple SNPs with RNAsnp default mode -m 1

$ RNAsnp -f examples/seq2.txt -s examples/snp2.txt
SNP  		W   Slen   GC	interval  d_max	 p-value interval  r_min  p-value
C9294A-U9296G	200 9605 0.4814	9261-9310 0.1951 0.0749	 9268-9317 0.2345 0.1213

Output details: 
Column 1: details of the given SNP
Column 2: length of the flanking region considered 
	  on either side of the SNPs to fold
Column 3: length of the given input sequence
Column 4: GC percent of the sequence interval considered
	  for folding
Column 5-7: results of best hit local region detected 
	    by distance measure (d_max)
Column 8-10: results of best hit local region detected 
	     by correlation measure (r_min)

-------------------------------------------------------------------
RNAsnp mode 2:
-------------------------------------------------------------------
1) Test for the effect of single SNP with RNAsnp default mode -m 1

$ RNAsnp -f examples/seq1.txt -s examples/snp1.txt -m 2 
SNP	w   Slen   GC	max_k d_max  p-value interval	d	p-value
U1013C	200 3344 0.5411	994   4.3961 0.2176  994-1019	0.1265	0.1232

2) Test for the effect of single SNP with RNAsnp mode 2 

$ RNAsnp -f examples/seq2.txt -s examples/snp2.txt -m 2
SNP		w   Slen   GC	max_k d_max  p-value interval	d	p-value
C9294A-U9296G	200 9605 0.4814	9270  7.0487 0.0624  9270-9298	0.2463	0.0099

Output details: 
Column 1: details of the given SNP
Column 2: length of the flanking region considered 
	  on either side of the SNPs to fold
Column 3: length of the given input sequence
Column 4: GC percent of the sequence interval considered
	  for folding
Column 5-7: results for the scanning of best local region with
            fixed length. The max_k is the start position where
	    the maximum distance(d_max) is detected.  
Column 8-10: results of best hit local region detected 
	     around the max_k region. 

-------------------------------------------------------------------
RNAsnp mode 3:
-------------------------------------------------------------------
1) Screen possible structure-disruptive SNPs in a sequence with 
   default p-value thresholds (pvalue1<0.4 and pvalue2<0.1)

$ RNAsnp -f examples/seq1.txt -m 3
SNP  w	Slen   GC   interval d  pvalue1 ewin interval d_max	 pvalue2
G1A 200	3344 0.5522  1-39    0.0185 0.2024  200	  1-50	  0.0961 0.0467	
G1C 200	3344 0.5522  1-46    0.0421 0.0755  200	  1-50	  0.1581 0.0183
....
....
....

2) Screen possible structure-disruptive SNPs in a sequence with 
   different p-value thresholds (pvalue1<0.1 and pvalue2<0.1)

$ RNAsnp -f examples/seq1.txt -m 3 --pvalue1 0.1 --pvalue2 0.1
SNP  w	Slen   GC  interval d  pvalue1 ewin	interval d_max	pvalue2
G1C 200	3344 0.5522 1-46    0.0421 0.0755  200   1-50	0.1581	0.0183
G7A 200	3344 0.5556 1-43    0.2236 0.0207  200	 1-50	0.1570	0.0996
....
....
....

Output details: 
Column 1: details of the screened SNPs
Column 2: length of the flanking region considered 
	  on either side of the SNP to fold
Column 3: length of the given input sequence
Column 4: GC percent of the sequence interval considered
	  for folding
Column 5-7: results of initial screen returned by mode 2.
Column 8: ewin is the length of the flanking region considered on either of the SNP
	  that also includes the local interval returned by mode 2. This subsequence
	  is then passed to mode 1 to compute the d_max 
Column 9-11: results of final screen returned by mode 1. 

Note: By default, the ewin(-e) uses the window length of 200. A higher computation
speed is achieved if the ewin value is reduced to 100 or 150, however, there 
will be little difference in the result compare to ewin with 200. In addition,
RNAsnp does automatically increase the length of ewin if the defined length
is not enough to cover the local interval returned by mode 2. 

Graphical representation of RNAsnp mode 3 output
------------------------------------------------
A perl script was developed to parse the output of the RNAsnp scanning mode (or mode 3)
and returns region that are sensitive to SNPs. 

Utils/selectSensRegion.pl -a 0.4  -b 0.1 -l 10 -n 3 -o <file_openingEnergy> < <mode3.output>

Example:
Utils/selectSensRegion.pl -a 0.4  -b 0.1 -l 10 -n 3  \\
  -o examples/sample_output/seq3_openen -j seq3 < examples/sample_output/seq3_mode3.out

To get more details about the options please check,
Utils/selectSensRegion.pl -h 
 

6) References
-------------
[1] Hofacker IL, Fontana W, Stadler PF, Bonhoeffer LS, Tacker M, Schuster P. (1994)
    Fast Folding and Comparison of RNA Secondary Structures.
    Monatshefte f. Chemie 125: 167-188
[2] Lorenz R, Bernhart SH, Honer zu Siederdissen C, Tafer H, Flamm C, Stadler PF,
    Hofacker IL (2011) ViennaRNA Package 2.0. Alg. Mol. Biol. 6:26.

If you find this software useful for your research, please cite the following work:
    Sabarinathan R, Tafer H, Seemann SE, Hofacker IL, Stadler PF, Gorodkin J (2013)
    RNAsnp: Efficient detection of local RNA secondary structure changes induced
    by SNPs. Hum Mutat. 34:546-56.

7) Contact
----------
For any comments or bug reports please contact the authors
Email: sabari@rth.dk, htafer@bioinf.uni-leipzig.de


