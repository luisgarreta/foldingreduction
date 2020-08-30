#!/usr/bin/Rscript
HELP="
AUTHOR: Luis Garreta 
MAIL: lgarreta@gmail.com

Given a trayectory it uses a fast clustering algorith to
reduce the trayectory to only the main representatives.
"
INPUT="  
   <inDir>      An input directory with the trayectory files 
   <outDir>     An output directory with the results (representative files)
   <threshold>  TM-score threshold for similarity between two protein structures
   <Size Bin>   Number of structures for each bin of the partiioned trajectory
   <N Cores>    Number of cores to use for paralelizing the procedure
"
USAGE  ="
Reduces a folding trajectory using a fast clustering strategy

USAGE	: Rscript folding_reduction.R <inDir> <outDir> <SizeBin> <TMSCORE> <K> <nCores>
EXAMPLE : Rscript folding_reduction.R inputDir outputDir 40 0.5 5 4"

# Default values
#threshold = 0.5   # TM-score threshold for comparisons between protein structures"
#SIZEBIN   = 40	# Number of files for each bin
#NCORES	  = 2	  # Number of cores for multiprocessing
library (parallel)
library (cluster)
library (optparse)

tmscorelib = sprintf ("%s/%s", Sys.getenv("FOLDING_REDUCTION_HOME"), "/tmscorelg.so")
dyn.load (tmscorelib)

#------------------------------------------------------------------
#------------------------------------------------------------------
main <- function () {
	args = commandArgs(trailingOnly = TRUE)
	#args = c("--inDir", "sarscov2", "--threshold", "0.9", "--nCores=7")
	if (length (args) < 1) {
		message (USAGE)
		message ("\nWhere:")
		message (INPUT)
		quit ()
	}
	p = readCheckParameters (args)

	# Create initial dirs for input, outputs
	createDir (p$outDir)
	system (sprintf ("mkdir -p %s/tmp", p$outDir))
	system (sprintf ("ln -s %s/%s %s", getwd(), p$inDir, p$pdbsDir))
	writeLog (p$outDir, p$inDir, p$sizeBin, p$threshold, p$k)

	# Split full trajectory in bins (blocks of pdbs)
	createBins (p$pdbsDir, p$outputDirBins, p$sizeBin)

	# Get Representatives for each bin
	reduceLocal_parallel (p$outputDirBins, p$outDir, p$threshold, p$nCores)

	# Uses K-Medoids to select K medoids for each local bin
	reduceGlobal_parallel (p$inputDirGlobal, p$outputDirGlobal, p$k, p$nCores)
}

#------------------------------------------------------------------
#------------------------------------------------------------------
readCheckParameters <- function (args) {
    indir     = make_option (c("-i", "--indir"), type="character", default=NULL, metavar="character", 
				help="Input directory with the PDB conformations of the trajectory (e.g. pbs2JOF)")
	outdir    = make_option (c("-o", "--outdir"), type="character", default="out", metavar="character", 
				help="Output directory to write reduction results (e.g. out2JOF)")
	threshold = make_option (c("-t", "--threshold"), type="double", default=0.5, metavar="number", 
				help="Threshold to take two protein structures as similar or redundant (e.g. 0.5)")
	nCores    = make_option (c("-c", "--ncores"), type="integer", default=1, metavar="number", 
				help="Number of processing cores to be used in parallel")
	nRepr     = make_option (c("-r", "--nrepresentative"), type="integer", default=1, metavar="number", 
	            help="Number of PDB representatives to select for each partition or bin (e.g. 100)")
	binSize   = make_option (c("-b", "--binsize"), type="integer", default=1, metavar="number", 
	            help="Number of PDB files for each partition or bin (e.g. 500)")

	optLst    = list (indir, outdir, threshold, nCores, nRepr, binSize, nCores)

	optParser = OptionParser (usage=USAGE, option_list=optLst)
	opt       = parse_args (optParser)

	if (is.null (opt$indir)) {
		print_help (optParser)
		quit ()
	}
}
	

#------------------------------------------------------------------
#------------------------------------------------------------------
old_readCheckParameters <- function (args) {
	parser = argparse::ArgumentParser (description="Reduced a protein folding trajectory")
	parser$add_argument ("--inDir", required=TRUE, help="Input directory with the trajectory PDBs (e.g. pdbs2JOF/)")
	parser$add_argument ("--outDir", default="out", help="Output directory to write reduction results (e.g. out2JOF/)")
	parser$add_argument ("--sizeBin", default=10, help="Number of PDB files for each partition or bin (e.g. 500)")
	parser$add_argument ("--k", default=30, help="Number of PDB representatives to select for each partition or bin (e.g. 100)")
	parser$add_argument ("--nCores", default="1", help="Number of cores to run in parallel (e.g. 4)")
	parser$add_argument ("--threshold", default=0.5, help="Threshold to take two protein structures as similar or redundant (e.g. 0.5)")

	args = parser$parse_args (args)

	args$outputDirBins	= sprintf ("%s/tmp/bins", args$outDir)
	args$pdbsDir		= paste0 (args$outDir, "/tmp/pdbs")
	args$inputDirGlobal	= paste0 (args$outDir, "/tmp/binsLocal")
	args$outputDirGlobal = args$outDir 

	numberOfPDBs         = length (list.files (args$inDir))
	args$sizeBin         = args$sizeBin/100. * numberOfPDBs
	args$k               = args$k/100. * args$sizeBin

	message ("Parameters:")
	message ("\t Input dir:          ", args$inDir)
	message ("\t Output dir:         ", args$outDir)
	message ("\t TM-score threshold: ", args$threshold)
	message ("\t Size of bins:       ", args$sizeBin)
	message ("\t K:                  ", args$k)
	message ("\t Num of Cores:       ", args$nCores)
	message ("\n")

	return (args)
}
	
#------------------------------------------------------------------
#------------------------------------------------------------------
writeLog <- function (outDir, proteinName, sizeOfBins, threshold, k) {
	inFile = file (paste0 (outDir, "/params.txt", open="a"))
	writeLines (paste0 ("Protein name: ", basename (proteinName)))
	writeLines (paste0 ("Size of bins : " , sizeOfBins))
	writeLines (paste0 ("TM-score trhreshold : " , threshold))
	writeLines (paste0 ("K representatives : " , k))
	close (inFile)
}


#--------------------------------------------------------
# Split the files of an input directory in bins according to the\n
# size of the bin. The bins are put in an output directory\n
#
# outputDir is the destiny dir for bins
# binSize is the number of file by bin
# sizeFill is the prefix for each bin filename
#--------------------------------------------------------
createBins <- function (inputDir, outputDir, binSize) {
	createDir (outputDir)

	inputFiles  = list.files (inputDir, pattern=".pdb")
	n = length (inputFiles)
	sizeFill = nchar (n)
	binList = splitBins (inputFiles, binSize)

	for (k in 1:length(binList)) {
		lst = binList [[k]]
		binNumber = k

		binName = formatC (k, width=sizeFill, format="d", flag=0)
		binDirname = sprintf ("%s/%s%s", outputDir, "bin", binName)
		message (">>> Creating bin ",  binName, "...")

		binFilename = sprintf ("%s.pdbs" , binDirname)
		filenames = unlist (lst)
		writeLines (filenames, binFilename)
	}
}

#--------------------------------------------------------
# Creates a list of sublist where each sublist corresponds to
# the files of each bin
#--------------------------------------------------------
splitBins <- function (inputFiles, binSize) {
	nSeqs = length (inputFiles)
	#nBins = nSeqs / binSize
	nBins = ceiling (nSeqs / binSize)

	binList = list()
	for (k in 0:(nBins-1)) {
		start = k*binSize
		end   = start + binSize 
		if (k < nBins-1) 
			binList = append (binList, list (inputFiles [(start+1):end]))
		else
			binList = append (binList, list (inputFiles [(start+1):nSeqs]))
	}

	return (binList)
}

#----------------------------------------------------------
# Parallel local reduction (see reduceLocal)
#----------------------------------------------------------
# Make a fast clustering of protein conformations from a 
# trayectory by doing first a fast local clustering 
# INPUT:  
#   <inputDir>       An input directory with the trayectory files 
#   <outputDir>      An output directory with the results (representative files)
#   <RMSD threshold> Threshold for local reduction comparisons that uses RMSD
#   <Size Bin>       Number of structures for each bin of the partiioned trajectory
#   <N Cores>        Number of cores to use for paralelizing the procedure
#
# OUTPUT:            An output dir with the clustering results for each bin
#----------------------------------------------------------
reduceLocal_parallel <- function (INPUTDIR, OUTPUTDIR, THRESHOLD, NCORES) {
	dirBins   = paste (OUTPUTDIR,"/tmp/binsLocal",sep="")
	dirPdbs   = OUTPUTDIR

	createDir (dirBins)

	binPathLst         = list.files (INPUTDIR, pattern=".pdbs", full.names=T)
	clusteringResults  = mclapply (binPathLst, reduceLocal, outputDir=dirBins, 
								        threshold=THRESHOLD,dirPdbs=dirPdbs, mc.cores=NCORES)

	writeClusteringResults (clusteringResults, dirPdbs, "pdbsLocal.pdbs")
}

#----------------------------------------------------------
# Reduction function to reduce a single bin
# Fast clustering following hobbohm algorith
# The first protein in the bin is selected as representative, then
# if the others are differente from it, they are preserved.
# Write the links to the representatives in the output dir
# Return a list with the representative as the first pdb in the group
#----------------------------------------------------------
reduceLocal <- function (inputBinPath, outputDir, threshold, dirPdbs) {
	cat ("\n>>> Local Reducing ", inputBinPath,"\n" )
	# Create the output dir for representatives
	outputBinPath = (paste (getwd(), outputDir, basename (inputBinPath), sep="/"))

	pdbsDir = paste (dirname (dirname (inputBinPath)),"/pdbs",sep="")
	pdbsNames = strsplit (readLines (inputBinPath), split=".pdbs")
	listOfPDBPaths <- paste (pdbsDir, pdbsNames, sep="/")

	# Fast clustering for bin, writes representatives to outputBinPath
	n = length (listOfPDBPaths)
	headProteinPath    = listOfPDBPaths [[n]] # Default head for the first group
	tmscoreValue       = -1                   # To create the link for the first group
	listOfSelectedPdbs = c (basename(headProteinPath))
	k = n - 1
	while (k >= 1) {
		targetProtein = listOfPDBPaths[[k]] 
		results = .Fortran ("gettmscore", pdb1=targetProtein, pdb2=headProteinPath, resTMscore=0.4)
		tmscoreValue = results$resTMscore

		if (tmscoreValue < threshold) {
			listOfSelectedPdbs = append (listOfSelectedPdbs, basename (targetProtein))
			headProteinPath = targetProtein
		}
		k = k - 1
	}
	listOfSelectedPdbs = sort (listOfSelectedPdbs)
	write.table (listOfSelectedPdbs, file=outputBinPath, sep="\n",col.names=F, row.names=F, quote=F)
	return (listOfSelectedPdbs)
}

#----------------------------------------------------------
# Make links of the selected PDBs into the output dir
#----------------------------------------------------------
writeClusteringResults <- function (clusteringResults, outputDir, outFile) {
	listOfPDBs = c ()
	for (binResults in clusteringResults) 
		for (pdbPath in binResults) {
			listOfPDBs = append (listOfPDBs, pdbPath)
		}
	filename = sprintf ("%s/tmp/%s", outputDir, outFile)
	print (filename)
	listOfPDBs = sort (listOfPDBs)
	listOfPDBs = paste ("pdbs/", basename(listOfPDBs),sep="")
	write.table (listOfPDBs, file=filename, sep="\n",col.names=F, row.names=F, quote=F)

	return (listOfPDBs)
}

#----------------------------------------------------------
# Parallel reduceGlobal 
#----------------------------------------------------------
# Makes a detailed global clustering of protein conformations 
# from the representatives resulting from the local clustering 
# INPUT:  inputDir filename with the protein conformations
#         outputDir filename to write the results
#
# OUTPUT: Medoids for each bin cluster and their distance matrix
#----------------------------------------------------------

reduceGlobal_parallel <- function (INPUTDIR, OUTPUTDIR, K, NCORES) {

	listOfBinPaths    = list.files (INPUTDIR, pattern=".pdbs", full.names=T)
	clusteringResults = mclapply (listOfBinPaths, reduceGlobal, K, mc.cores=NCORES)

	listOfPDBs = writeClusteringResults (clusteringResults, OUTPUTDIR, "pdbsGlobal.pdbs")

	# Write selected PDBs to reduction dir
	reductionDir = paste0(OUTPUTDIR, "/reduction")
	createDir (reductionDir)
	print (">>>")
	print (getwd())
	print (reductionDir)
	print (">>>")
	print (listOfPDBs)
	for (p in paste0 ("out/tmp/", listOfPDBs))
		file.copy (p, reductionDir)
	#mclapply (listOfPDBs, file.copy,to=reductionDir) 
}

#----------------------------------------------------------
# Reduction function to reduce a single bin
# Clustering around medoids. Return k medoid for the bin
#----------------------------------------------------------
reduceGlobal <- function (inputBinPath, K) {
	cat ("\n>>> Global Reducing ", inputBinPath )

	# Fast clustering for bin, writes representatives to clusDir
	#listOfPDBPaths <<- list.files (inputBinPath, pattern=".pdb", full.names=T)
	pdbsDir = paste (dirname (dirname (inputBinPath)),"/pdbs",sep="")
	listOfPDBPaths <<- paste (pdbsDir, readLines(inputBinPath), sep="/")


  # Clustering around medoids. Return one medoid for all inputs
	nPdbs = length (listOfPDBPaths)
	if (nPdbs < 2)
		medoids = 1
	else if (nPdbs <= K)
		medoids = seq (nPdbs)
	else {
		binDir = inputBinPath
		cat ("\n>>> Calculating distance matrix", inputBinPath,"\n")
		distanceMatrix <- getTMDistanceMatrix (listOfPDBPaths)
		split          <- -1 * nPdbs / K
		initialMedoids <- round (seq (nPdbs, 1, split))
		pamPDBs        <- pam (distanceMatrix, k=K, diss=F, medoids=initialMedoids)
		medoids        <- pamPDBs$id.med
	}

	medoidName <- listOfPDBPaths [medoids]
	return (medoidName)
}

#--------------------------------------------------------------
# Calculate pairwise using TM-score distance
#--------------------------------------------------------------
getTMDistanceMatrix <- function (listOfPDBPaths) {
	n = length (listOfPDBPaths)
	mat = matrix (seq (1,n))

	distMat = proxy::dist (mat, method=calculateTmscore)
	return (distMat)
}

#----------------------------------------------------------
# Calculate the TM-scores using a external tool "TMscore"
#----------------------------------------------------------
calculateTmscore <- function (targetProtein, referenceProtein) {
	targetProtein    = listOfPDBPaths [[targetProtein]]
	referenceProtein = listOfPDBPaths [[referenceProtein]]

	results = .Fortran ("gettmscore", pdb1=targetProtein, pdb2=referenceProtein,  resTMscore=0.4)
	tmscoreValue = results$resTMscore

	return  (tmscoreValue)
}

#----------------------------------------------------------
# Create dir, if it exists the it is renamed old-XXX
#----------------------------------------------------------
createDir <- function (newDir) {
	checkOldDir <- function (newDir) {
		name  = basename (newDir)
		path  = dirname  (newDir)
		if (dir.exists (newDir) == T) {
			oldDir = sprintf ("%s/old-%s", path, name)
			if (dir.exists (oldDir) == T) {
				checkOldDir (oldDir)
			}

			file.rename (newDir, oldDir)
		}
	}

	checkOldDir (newDir)
	system (sprintf ("mkdir %s", newDir))
}

#--------------------------------------------------------------
#--------------------------------------------------------------
main () 
