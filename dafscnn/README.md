DAFS: simultaneous aligning and folding of RNA sequences by dual decomposition
==============================================================================

Requirements
------------

* [Vienna RNA package](http://www.tbi.univie.ac.at/~ivo/RNA/) (>= 1.8)
* (optional)
  [GNU Linear Programming Kit](http://www.gnu.org/software/glpk/) (>=4.41)
  or [Gurobi Optimizer](http://www.gurobi.com/) (>=2.0)
  or [ILOG CPLEX](http://www.ibm.com/software/products/ibmilogcple/) (>=12.0)

Install
-------

	./configure --with-vienna-rna=/path/to/vienna-rna
	make
	make install

To use "--ipknot" option for pseudoknotted common secondary structure
prediction, build DAFS with an IP solver.

For GLPK,

	./configure --with-vienna-rna=/path/to/vienna-rna --with-glpk

For Gurobi, 

	./configure --with-vienna-rna=/path/to/vienna-rna --with-gurobi=/path/to/gurobi

For CPLEX,

	./configure --with-vienna-rna=/path/to/vienna-rna --with-cplex \
	            --with-cplex-include=/path/to/cplex/include \
		        --with-cplex-lib=/path/to/cplex/lib

Usage
-----

DAFS can take FASTA formatted RNA sequences as input, then produce a structural alignment.

    Usage: dafs [OPTIONS]... [FILES]
    
     -h, --help                Print help and exit
         --full-help           Print help, including hidden options, and exit
     -V, --version             Print version and exit
     -w, --weight=FLOAT        Weight of the expected accuracy score for secondary 
                                 structures  (default=`4.0')
     -m, --max-iter=INT        The maximum number of iteration of the subgradient 
                                 optimization  (default=`600')
     -v, --verbose=INT         The level of verbose outputs  (default=`0')
   
    Options for alignments:
     -a, --align-model=STRING  Alignment model for calcualating matching 
                                 probablities  (possible values="CONTRAlign", 
                                 "ProbCons" default=`ProbCons')
     -u, --align-th=FLOAT      Threshold for matching probabilities  
                                 (default=`0.01')
   
    Options for folding:
     -s, --fold-model=STRING   Folding model for calculating base-pairing 
                                 probablities  (possible values="Boltzmann", 
                                 "Vienna", "CONTRAfold" default=`Boltzmann')
     -t, --fold-th=FLOAT       Threshold for base-pairing probabilities  
                                 (default=`0.2')
     -T, --fold-th1=FLOAT      Threshold for base-pairing probabilities of the 
                                 conclusive common secondary structures
         --ipknot              Use IPknot decoding  (default=off)


Example
-------

	% dafs RF00005:0.fa
	[ 0.0985233 [ 0.585795 [ 0.933469 M68929-1/151018-150946 X00360-1/1-73 ] [ 0.826623 X12857-1/421-494 [ 0.935672 J05395-1/2325-2252 M16863-1/21-94 ] ] ] [ 0.349897 [ 0.780743 J04815-1/3159-3231 [ 0.96716 J01390-1/6861-6932 M20972-1/1-72 ] ] [ 0.74278 K00228-1/1-82 AC009395-7/99012-98941 ] ] ]
	>SS_cons
	(((((((...(((..............))).......(((((..........)))))......(.((((.......))))).))))))).
	> J01390-1/6861-6932
	CAGGUUA-GAGCC-AGGU-GGU-UA--GGCG------UCUUGUUU--GG-GUCAAGA-AAUUGU-UAUGUUCGAAUCAUAA-UAACCUGA
	> J05395-1/2325-2252
	GGUUUCG-UGGUC-UAGUCGGUUAU--GGCA------UCUGCUUA--AC-ACGCAGA-ACGUCC-CCAGUUCGAUCCUGGG-CGAAAUCG
	> K00228-1/1-82
	GGUUGUUUG-GCCGA-GC-GGU-CUAAGGCGCCUGAUUCAAGCUCAGGU-AUCGUAA--GAUGCAAGAGUUCGAAUCUCUU-AGCAACCA
	> AC009395-7/99012-98941
	GGCUCAA-U----------GGU-CUAG-GGGUAUGAUUCUCGCUUUGGG-UGCGAGA--GGUCC-CGGGUUCAAAUCCCGG-UUGAGCCC
	> J04815-1/3159-3231
	AGAGCUU-GCUCC-CAAA-GCU-UG--GGUG------UCUAGCUG--AU-AAUUAGA-CUAUCA-AGGGUUAAAUUCCCUUCAAGCUCUA
	> M20972-1/1-72
	AGGGCUA-UAGCU-CAGC-GGU-AG--AGCG------CCUCGUUU--AC-ACCGAGA-AUGUCU-ACGGUUCAAAUCCGUA-UAGCCCUA
	> M68929-1/151018-150946
	CGCGGGA-UAGAG-UAAUUGGU-AA--CUCG------UCAGGCUC--AU-AAUCUGA-AUGUUG-UGGGUUCGAAUCCGAC-UCCCGCCA
	> X00360-1/1-73
	CCGACCU-UAGCU-CAGUUGGU-AG--AGCG------GAGGACUG---UAGAUCCUU-AGGUCA-CUGGUUCGAAUCCGGU-AGGUCGGA
	> X12857-1/421-494
	GCGGAUG-UAGCC-AAGUGGAUCAA--GGCA------GUGGAUUG--UG-AAUCCACCAUG-CG-CGGGUUCAAUUCCCGU-CAUUCGCC
	> M16863-1/21-94
	GGGCUCG-UAGCU-CAGAGGAUUAG--AGCA------CGCGGCUA--CG-AACCACG-GUGUCG-GGGGUUCGAAUCCCUC-CUCGCCCA


References
----------

* Sato, K., Kato, Y., Akutsu, T., Asai, K., Sakakibara, Y.: DAFS: simultaneous aligning and folding RNA sequences via dual decomposition. *Bioinformatics*, 28(24):3218-3224, 2012.
