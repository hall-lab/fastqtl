#include "data.h"

void data::readCovariates(string fcov) {
	string buffer;
	vector < string > str;
	int n_includedS = 0;
	int n_excludedS = 0;
	int n_includedC = 0;
	int n_excludedC = 0;
	int n_missingS = 0;
	vector < int > mappingS;

	LOG.println("\nReading covariates in [" + fcov + "]");
	ifile fd (fcov);

	//Read samples
	getline(fd, buffer);
	sutils::tokenize(buffer, str);
	for (int t = 1 ; t < str.size() ; t ++) {
		if (checkSample(str[t], false)) {
			int idx_sample = -1;
			for (int i = 0 ; i < sample_count && idx_sample < 0 ; i++) if (sample_id[i] == str[t]) idx_sample = i;
			mappingS.push_back(idx_sample);
			if (idx_sample >= 0) n_includedS ++;
			else n_missingS ++;
		}
	}

	//Read covariates
	while(getline(fd, buffer)) {
		sutils::tokenize(buffer, str);
		if (str.size() < 2) LOG.error("Wrong genotype covariate file format");
		if (checkCovariate(str[0])) {
			covariate_val.push_back(vector < float > (sample_count, 0.0));
			for (int t = 1 ; t < str.size() ; t ++) {
				if (mappingS[t-1] >= 0) {
					covariate_val.back()[mappingS[t-1]] = atof(str[t].c_str());
				}
			}
			n_includedC ++;
		} else n_excludedC ++;
	}

	//Finalise
	covariate_count = n_includedC;
	LOG.println("  * " + sutils::int2str(n_includedS) + " samples included");
	if (n_excludedS > 0) LOG.println("  * " + sutils::int2str(n_excludedS) + " samples excluded");
	if (n_missingS > 0) LOG.println("  * " + sutils::int2str(n_excludedS) + " samples excluded without phenotype data");
	LOG.println("  * " + sutils::int2str(n_includedC) + " covariate(s) included");
	if (n_excludedC > 0) LOG.println("  * " + sutils::int2str(n_excludedC) + " covariate(s) excluded");
	fd.close();
}
