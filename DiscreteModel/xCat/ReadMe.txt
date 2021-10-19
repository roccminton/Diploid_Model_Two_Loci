In this folder you can find every file which could be used on xcat.

--General--
Remember here to change the path to your xcat account in the SLiM and the Rscripts!

--Bash files--
For every SLiM file I made a file for running a job, these can be found in the folder "bin". 
Here I first check if one of the result with a special Seed (from 1 to 50)(or a special number of genes or a special lambda) is existing with checkFile.R(or checkGenes.R/checkLambda.R) or another variant of this file. These variants are there to differ between different simulations. For some simulations I also checked for a Burn-In phase with a Rscript called burnIn_False.R or burnIn_True.R. After that the SLiM file starts.
After that I start the slim file itself.

--GetGenes--
The files geneErstellen.R and geneErstellen2.R were made to made genes, which can be load to the SLiM files. The standard genes are load with the files start.txt, genomicType.txt, mutationRate.txt and end.txt to SLiM, which where made originally this "geneErstelle.R".

--SLiM-files--
All SLiM files where at the latest updated on 23.3.2021.
All files with dynasty_ are studying different aspects of the consanguineous mating scheme.
All files with random_ are studying different aspects of the random mating scheme.
The file with western_ was made to simulate a western population which prefer to mate with nonrelatives, but this mating scheme should be review before using it.

In the following I will describe the aspects to study, if there isn't said something different a growth phase starts at generation 500:
_xcat.slim is used to count the mutation load and the prevalence, it's the "standard simulation"

_relate.slim is used to analyse the relationships between the individuals, should also be reviewed

_popsize_smallGenome is used analyse different aspects like population size and mutationload on a simulation which has no population size limit, every couple gets in mean lambda (constant) offspring where lambda can be differ; the genome is much smaller than in the "standard" simulation 

_random_popsize_smallGenome_difLambda is for tests with a much smaller Genome(here 100 genes with each 1 bp) to compair different lambdas (expected value for the litter size of women); lambda is regularized with checkLambda.R

_pois.slim is used to count the mutation load and the prevalence in a WF model with poisson distributed population size, where the total number is poisson distributed; most recent version for the WF model (not done with dynasty)

_nWF.slim is used to count the mutation load and the prevalence with a non Wright-Fisher model

_nWF_populationsize.slim is used to analyse different aspects like population size and mutationload in a non Wright-Fisher model which has no population size limit, every couple gets in mean lambda (constant) offspring where lambda can be differ 

_nWF_nocouples.slim is used to count the mutation load and the prevalence with a non Wright-Fisher model, where a the offspring of a woman have different fathers (no effect on prevalence and mutation load comparing to _nWF.slim)

_nWF_difFitness.slim is used like _nWF.slim but with a difference in the fitness values of the individuals. Here every ill individual will be count!; most recent version of the nonWF model (not done with dynasty)

_mutloadIlls.slim is used to count the mutation load in general and mutation load of all healthy individuals and mutation load of all ill individuals

_litterPois.slim is used to count different aspects like population size, the mutation load and the prevalence in a WF model, which has no population size limit made with a poisson distributed population size, where the litter of each woman is poisson distributed with constant lambda

_Kopie.slim is used get an mean value of the mutation load in a not growing population

_hw.slim is used to count the number of mutations that are in a hardy-weinberg equilibrium

_decline.slim is used to count the mutation load and the prevalence with an additional decline of the population size

_cousins.slim is used to count the mean value of individuals which have the same M grandparents, M can be between 0 and 4 and can be changed in the constant "grandparent"






