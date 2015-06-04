#ifndef _DATA_H
#define _DATA_H

#define ___NA___ (0.0/0.0)

#include "utils/utils.h"
#include "region.h"
#include "utils/portR.h"
#include "utils/auxiliary.h"

class data {
public:
	//INCLUDE/EXCLUDE LISTS
	set < string > sample_inclusion;
	set < string > sample_exclusion;
	set < string > genotype_inclusion;
	set < string > genotype_exclusion;
	set < string > phenotype_inclusion;
	set < string > phenotype_exclusion;
	set < string > covariate_inclusion;
	set < string > covariate_exclusion;

	//REGIONS
	region regionPhenotype;
	region regionGenotype;
	//float cis_window;
    
	//SAMPLES
	int sample_count;									//sample number
	vector < string > sample_id;						//sample IDs

	//GENOTYPES
	int genotype_count;									//variant site number
	vector < vector < float > > genotype_val;			//genotype dosages
	vector < string > genotype_chr;						//variant site chromosome
	vector < string > genotype_id;						//variant site IDs
	vector < int > genotype_pos;						//variant site positions
	vector < double > genotype_sd;						//variabt site standard deviation

	//PHENOTYPES
	int phenotype_count;								//phenotype number
	vector < vector < float > > phenotype_val;			//phenotype values
	vector < string > phenotype_id;						//phenotype ids
	vector < string > phenotype_chr;					//phenotype chromosomes
	vector < int > phenotype_start;						//phenotype start positions
	vector < int > phenotype_end;						//phenotype end positions
    vector < double > phenotype_sd;						//phenotype standard deviation
    vector < vector < int > > phenotype_cluster;		//phenotype cluster (parallel jobs)

	//COVARIATES
	int covariate_count;								//covariate number
	vector < vector < float > > covariate_val;			//covariate values
	vector < string > covariate_id;						//covariate ids

	//PARAMETERS
	bool spearman;

	//CONSTRUCTOR / DESTRUCTOR
	data();
	void clear();

	//READ EXCLUSION/INCLUSION LISTS (IO/readInclusionsExcusions.cpp)
	void readSamplesToExclude(string);
	void readSamplesToInclude(string);
	void readGenotypesToExclude(string);
	void readGenotypesToInclude(string);
	void readPhenotypesToExclude(string);
	void readPhenotypesToInclude(string);
	void readCovariatesToExclude(string);
	void readCovariatesToInclude(string);
	bool checkSample(string &, bool checkDuplicates = true);
	bool checkGenotype(string &);
	bool checkPhenotype(string &);
	bool checkCovariate(string &);

	//REGIONS
	void setPhenotypeRegion(string);
	void setPhenotypeRegion(int);
	string getPhenotypeRegion(int);
	void deduceGenotypeRegion(int);

	//READ GENOTYPE DATA (IO/readGenotypes.cpp)
	void readGenotypesVCF(string);
	void readGenotypesImpute2(string, string);

	//READ PHENOTYPE DATA (IO/readPhenotypes.cpp)
	void readPhenotypes(string);
	void scanPhenotypes(string);

	//READ COVAR DATA
	void readCovariates(string);

	//LINEAR REGRESSION FOR COVARIATES
	void correctGenotypesForCovariates();
	void correctPhenotypesForCovariates();

	//GENERAL MANAGMENT
	void clusterizePhenotypes(int);
	void imputeGenotypes();
    void imputePhenotypes();

	//NORMALISE GENOTYPE AND PHENOTYPE DATA (analysis/basics.cpp)
	void normalisePhenotypes();
	void normaliseGenotypes();
	void computeSDGenotypes();
	void computeSDPhenotypes();
	void normalTranformPhenotypes();
	void rankTranformPhenotypes();
	void rankTranformGenotypes();

	//COMPUTATION METHODS [ALL INLINES FOR SPEED]
	double getCorrelation(int, int);
	double getCorrelation(int, vector < float > &);
	double getBeta(double, double, double);
	double getPvalue(double);

	//ANALYSIS
	void runNominal(string, float);
	void runPermutation(string, float, int, int);

	//COMMANDS
	void writeCommands(string, int, int, char **);
};


//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//

inline double data::getCorrelation(int gidx, int pidx) {
	double corr = 0.0;
	for (int s = 0 ; s < sample_count ; s ++) corr += genotype_val[gidx][s] * phenotype_val[pidx][s];
	return corr;
}

inline double data::getCorrelation(int gidx, vector < float > & phenotypes) {
	double corr = 0.0;
	for (int s = 0 ; s < sample_count ; s ++) corr += genotype_val[gidx][s] * phenotypes[s];
	return corr;
}

inline double data::getBeta(double gidx, double pidx, double corr) {
	if (genotype_sd[gidx] < 1e-16 || phenotype_sd[pidx] < 1e-16) return 0;
	else return corr * phenotype_sd[pidx] / genotype_sd[gidx];
}

inline double data::getPvalue(double corr) {
    return portR::pf((sample_count - 2) * corr * corr / (1 - corr * corr), 1, sample_count - 2,0,0);
}

#endif
