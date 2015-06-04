#include "data.h"
#include "utils/ranker.h"

void data::normalisePhenotypes() {
	LOG.println("\nNormalising phenotypes");
	for (int p = 0; p < phenotype_count ; p ++) {
		double mean = 0.0, sum = 0.0;
		for (int s = 0; s < sample_count ; s ++) mean += phenotype_val[p][s];
		mean /= sample_count;
		for (int s = 0; s < sample_count ; s ++) {
			phenotype_val[p][s] -= mean;
			sum += phenotype_val[p][s] * phenotype_val[p][s];
		}
		sum = sqrt(sum);
		if (sum == 0) sum = 1;
		for (int s = 0; s < sample_count ; s ++) phenotype_val[p][s] /= sum;
	}
}

void data::normaliseGenotypes() {
	LOG.println("\nNormalising genotypes");
	for (int g = 0; g < genotype_count ; g ++) {
		double mean = 0.0, sum = 0.0;
		for (int s = 0; s < sample_count ; s ++) mean += genotype_val[g][s];
		mean /= sample_count;
		for (int s = 0; s < sample_count ; s ++) {
			genotype_val[g][s] -= mean;
			sum += genotype_val[g][s] * genotype_val[g][s];
		}
		sum = sqrt(sum);
		if (sum == 0) sum = 1;
		for (int s = 0; s < sample_count ; s ++) genotype_val[g][s] /= sum;
	}
}

void data::computeSDGenotypes() {
	LOG.println("\nComputing SD for genotypes");
	genotype_sd = vector < double > (genotype_count, 0.0);
	for (int g = 0; g < genotype_count ; g ++) {
		double mean = 0.0;
		for (int s = 0; s < sample_count ; s ++) mean += genotype_val[g][s];
		mean /= sample_count;
		genotype_sd[g] = 0.0;
		for (int s = 0; s < sample_count ; s ++) genotype_sd[g] += (genotype_val[g][s] - mean) * (genotype_val[g][s] - mean);
		genotype_sd[g] = sqrt(genotype_sd[g] / sample_count);
	}
}

void data::computeSDPhenotypes() {
	LOG.println("\nComputing SD for phenotypes");
	phenotype_sd = vector < double > (phenotype_count, 0.0);
	for (int p = 0; p < phenotype_count ; p ++) {
		double mean = 0.0;
		for (int s = 0; s < sample_count ; s ++) mean += phenotype_val[p][s];
		mean /= sample_count;
		phenotype_sd[p] = 0;
		for (int s = 0; s < sample_count ; s ++) phenotype_sd[p] += (phenotype_val[p][s] - mean) * (phenotype_val[p][s] - mean);
		phenotype_sd[p] = sqrt(phenotype_sd[p] / sample_count);
	}
}

void data::normalTranformPhenotypes() {
	string method = "average";
	LOG.println("\nQuantile normalise phenotypes to Normal distribution");
	for (int p = 0; p < phenotype_count ; p ++) {
		vector < float > R;
		rank(phenotype_val[p], R, method);
		double max = 0;
		for (int s = 0 ; s < sample_count ; s ++) {
			R[s] = R[s] - 0.5;
			if (R[s] > max) max = R[s];
		}
		max = max + 0.5;
		for (int s = 0 ; s < sample_count ; s ++) {
			R[s] /= max;
			phenotype_val[p][s] = putils::qnorm(R[s], 0.0, 1.0, 1, 0);
		}
	}
}

void data::rankTranformPhenotypes() {
	string method = "average";
	LOG.println("\nRanking phenotypes");
	for (int p = 0; p < phenotype_count ; p ++) {
		vector < float > R;
		rank(phenotype_val[p], R, method);
		phenotype_val[p] = R;
	}
}

void data::rankTranformGenotypes() {
	string method = "average";
	LOG.println("\nRanking genotypes");
	for (int g = 0; g < genotype_count ; g ++) {
		vector < float > R;
		rank(genotype_val[g], R, method);
		genotype_val[g] = R;
	}
}

