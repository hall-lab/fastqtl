#define _GLOBAL

#include "data.h"

data::data() {
	sample_count = 0;
	genotype_count = 0;
	phenotype_count = 0;
	covariate_count = 0;
	spearman = false;
}

void data::clear() {
	sample_count = 0;
	sample_id.clear();
	genotype_count = 0;
	genotype_val.clear();
	genotype_chr.clear();
	genotype_id.clear();
	genotype_pos.clear();
	genotype_sd.clear();
	phenotype_count = 0;
	phenotype_val.clear();
	phenotype_id.clear();
	phenotype_chr.clear();
	phenotype_start.clear();
	phenotype_end.clear();
	phenotype_sd.clear();
	covariate_count = 0;
	covariate_val.clear();
	covariate_id.clear();
}

bool data::checkSample(string & id, bool checkDuplicate) {
	bool included = ((sample_inclusion.size() == 0)?true:sample_inclusion.count(id));
	bool excluded = ((sample_exclusion.size() == 0)?false:sample_exclusion.count(id));
	if (!included || excluded) return false;
	if (checkDuplicate) for (int s = 0 ; s < sample_id.size() ; s++) if (id == sample_id[s]) LOG.error("Duplicate sample id [" + id +"]");
	return true;
}

bool data::checkGenotype(string & id) {
	bool included = ((genotype_inclusion.size() == 0)?true:genotype_inclusion.count(id));
	bool excluded = ((genotype_exclusion.size() == 0)?false:genotype_exclusion.count(id));
	if (!included || excluded) return false;
	for (int g = 0 ; g < genotype_id.size() ; g++) if (id == genotype_id[g]) LOG.error("Duplicate variant site id [" + id + "]");
	return true;
}

bool data::checkPhenotype(string & id) {
	bool included = ((phenotype_inclusion.size() == 0)?true:phenotype_inclusion.count(id));
	bool excluded = ((phenotype_exclusion.size() == 0)?false:phenotype_exclusion.count(id));
	if (!included || excluded) return false;
	for (int p = 0 ; p < phenotype_id.size() ; p++) if (id == phenotype_id[p]) LOG.error("Duplicate phenotype id [" + id + "]");
	return true;
}

bool data::checkCovariate(string & id) {
	bool included = ((covariate_inclusion.size() == 0)?true:covariate_inclusion.count(id));
	bool excluded = ((covariate_exclusion.size() == 0)?false:covariate_exclusion.count(id));
	if (!included || excluded) return false;
	for (int c = 0 ; c < covariate_id.size() ; c++) if (id == covariate_id[c]) LOG.error("Duplicate covariate id [" + id + "]");
	return true;
}

void data::clusterizePhenotypes(int K) {
	//initialise (1 cluster = 1 chromosome)
	phenotype_cluster = vector < vector < int > > (1, vector < int > (1, 0));
	for (int p = 1 ; p < phenotype_count ; p ++) {
		if (phenotype_chr[p] != phenotype_chr[p-1]) phenotype_cluster.push_back(vector < int > (1, p));
		else phenotype_cluster.back().push_back(p);
	}

	//iterate (split cluster in the middle until K clusters are built)
	int max_idx, max_val, max_mid;
	do {
		max_idx = -1;
		max_val = +1;
		max_mid = -1;
		for (int k = 0 ; k < phenotype_cluster.size() ; k ++) {
			if (phenotype_cluster[k].size() > max_val) {
				max_val = phenotype_cluster[k].size();
				max_idx = k;
			}
		}
		if (max_idx >= 0) {
			max_mid = max_val / 2;
			while (max_mid > 1 && phenotype_end[phenotype_cluster[max_idx][max_mid-1]] >= phenotype_start[phenotype_cluster[max_idx][max_mid]]) max_mid --;
			phenotype_cluster.push_back(vector < int > (phenotype_cluster[max_idx].begin() + max_mid, phenotype_cluster[max_idx].end()));
			phenotype_cluster[max_idx].erase(phenotype_cluster[max_idx].begin() + max_mid, phenotype_cluster[max_idx].end());
		}
	} while (phenotype_cluster.size() < K && max_idx >= 0);
}

void data::setPhenotypeRegion(string reg) {
	regionPhenotype.set(reg);
}

void data::setPhenotypeRegion(int k) {
	regionPhenotype.chr = phenotype_chr[phenotype_cluster[k][0]];
	regionPhenotype.start = phenotype_start[phenotype_cluster[k][0]];
	regionPhenotype.end = phenotype_end[phenotype_cluster[k].back()];
}

string data::getPhenotypeRegion(int k) {
	return string (phenotype_chr[phenotype_cluster[k][0]] + ":" + sutils::int2str(phenotype_start[phenotype_cluster[k][0]]) + "-" + sutils::int2str(phenotype_end[phenotype_cluster[k].back()]));
}

void data::deduceGenotypeRegion(int W) {
	regionGenotype.chr = regionPhenotype.chr;
	regionGenotype.start = regionPhenotype.start - W;
	if (regionGenotype.start < 0) regionGenotype.start = 0;
	regionGenotype.end = regionPhenotype.end + W;
}

void data::imputeGenotypes() {
	LOG.println("\nImputing missing genotypes");
	for (int g = 0; g < genotype_count ; g ++) {
		double mean = 0.0;
		int c_mean = 0;
		for (int s = 0; s < sample_count ; s ++) {
			if (genotype_val[g][s] >= 0) {
				mean += genotype_val[g][s];
				c_mean ++;
			}
		}
		mean /= c_mean;
		for (int s = 0; s < sample_count ; s ++) if (genotype_val[g][s] < 0) genotype_val[g][s] = mean;
	}
}

void data::imputePhenotypes() {
	LOG.println("\nImputing missing phenotypes");
	for (int p = 0; p < phenotype_count ; p ++) {
		double mean = 0.0;
		int c_mean= 0;
		for (int s = 0; s < sample_count; s ++) {
			if (!isnan(phenotype_val[p][s])) {
				mean += phenotype_val [p][s];
				c_mean ++;
			}
		}
		mean /= c_mean;
		for (int s = 0; s < sample_count ; s ++) if (isnan(phenotype_val[p][s])) phenotype_val[p][s] = mean;
	}
}

