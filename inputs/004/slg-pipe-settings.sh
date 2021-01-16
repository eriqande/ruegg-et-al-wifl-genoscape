



### Control What Gets Done
DO_GSI_SIM_STUFF=0;         # 1 means do the MissingDataTable and the GeneClass AssMat
DO_PHYLIP_STUFF=0;          # 1 means build CSE-NJ tree and consensus tree and transfer
                            # bootstrap values back to the distance tree
DO_GENEPOP=0                # must be 1 to do any genepop analyses at all
DO_GENEPOP_HWE=0;           # 1 means do Guo and Thompson MCMC Fisher Exact test for HWE and 0 means do not do it
DO_GENEPOP_LD=0;            # 1 means test for LD between all pairs of loci and 0 means do not do it
PREPARE_COLONY_AREA=0;      # 1 means make the ColonyArea with everything needed to run Colony on each
                            # separate population
PREPARE_FSTAT_DIR=0;        # make the fstat directory which gives you an input file for fstat
PREPARE_GENETIX_DIR=0;      # make the genetix directory with an fstat file to import into genetix
PREPARE_STRUCTURE_AREA=1;   # make a directory for Structure stuff.

DO_ISOL_BY_DIST=0;          # make an isol_by_dist directory and do an isolation by distance analysis
                            # note that if this is done you have to supply the path to the distance matrix
                            # to the next variable
ISOL_DISTMAT=/Users/eriq/Documents/work/prj/NorCalSH_PopStruct/IsolByDist/stlhead_sites_dm_ShNames.txt;
IBD_MANTEL_PERMS=10000;

#### Options for different analyses

## For Phylip
PhylipNBOOT=1000;  # number of bootstrap replicates for phylip consensus tree


## For Colony


## For Structure
StructOutDir="StructureArea";      # the desired name for the output directory.
StructThePops="the_pops.txt";      # the population name file in slg_pipe PopFile format.
DistructNamesColumn=1;             # the column number (tab delimited) in StructThePops to use for the names in distruct
StructNameCol=1;                   # the column number in the PopsFile (which is tab delimited,
                                   # recall) from which you want to grab the names of the pops for labels.
KVals="2 3 4 5 6 7 8 9";                   # is a quoted string (like \"2 3 5 7\") of different K values
                                   # at which to run structure.
StructReps="5 5 5 5 5 5 5 5"                     # should be  either:
                                   #     a) a single number which is the number of reps to run each K at.
                                   #     b) a quoted string (like \"5 5 2 2\") of length exactly equal to  the quoted string
                                   #        giving the number of K values.  In this case, each value of K
                                   #        will have number of reps as given by the quoted string here.

                                   # mainparams template file
StructNumBurnIn=50000;             # number of structure Burn In Sweeps
StructNumSweep=150000;             # number of data collection sweeps.
StructUseLocationPrior=0;          # set this to 1 to use the Location Prior option in Structure (Hubisz et al 2009)
StructMainP=$SLG_PATH/inputs/ForStructure/mainparams;    # template for mainparams
StructExtraP=$SLG_PATH/inputs/ForStructure/extraparams   # template for extraparams
StructBinary=$SLG_PATH/inputs/ForStructure/structure2.3.4  # structure executable
ParallelScript=$SLG_PATH/inputs/parallel;
ClumpAndDistruct=$SLG_PATH/inputs/ForStructure/clump_and_distruct;




