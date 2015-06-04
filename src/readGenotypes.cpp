#include "data.h"
#include "utils/tabix.hpp"

void data::readGenotypesVCF(string fvcf) {
	string buffer;
	vector<string> str, field;
	int n_includedG = 0;
	int n_excludedG = 0;
	int n_excludedF = 0;
	int n_includedS = 0;
	int n_excludedS = 0;
	int n_missingS = 0;
	vector < int > mappingS;

	//Initialise
	LOG.println("\nReading genotype data in [" + fvcf + "] in VCF format");
	if (!futils::isFile(fvcf + ".tbi")) LOG.error("index file missing [" + fvcf + ".tbi]");
	Tabix fd (fvcf);

	//Read samples
	fd.getLastHeader(buffer);
	sutils::tokenize(buffer, str);
	if (str.size() < 10) LOG.error("Wrong VCF file format");
	for (int t = 9 ; t < str.size() ; t ++) {
		if (checkSample(str[t], false)) {
			int idx_sample = -1;
			for (int i = 0 ; i < sample_count && idx_sample < 0 ; i++) if (sample_id[i] == str[t]) idx_sample = i;
			mappingS.push_back(idx_sample);
			if (idx_sample >= 0) n_includedS ++;
			else n_missingS ++;
		} else n_excludedS ++;
	}
	if (n_includedS != sample_count) LOG.error("Genotype data does not overlap with phenotype data, check your files!");

	//Read genotypes
	if (!fd.setRegion(regionGenotype.str())) LOG.error("Failed to get region " + regionGenotype.str() + " in [" + fvcf + "]");
	while (fd.getNextLine(buffer)) {
		sutils::tokenize(buffer, str);
		if (checkGenotype(str[2])) {
			//Check VCF format
			bool gt_field = false;
			int idx_field = -1;
			sutils::tokenize(str[8], field, ":");
			for (int f = 0 ; f < field.size() ; f ++) if (field[f] == "DS") idx_field = f;
			if (idx_field < 0) { for (int f = 0 ; f < field.size() ; f ++) if (field[f] == "GT") idx_field = f; gt_field = true; }
			//Read data is format is correct
			if (idx_field >= 0) {
				genotype_id.push_back(str[2]);
				genotype_chr.push_back(str[0]);
				genotype_pos.push_back(atoi(str[1].c_str()));
				genotype_val.push_back(vector < float > (sample_count, 0.0));
				for (int t = 9 ; t < str.size() ; t ++) {
					if (mappingS[t-9] >= 0) {
						sutils::tokenize(str[t], field, ":");
						if (str[t] == "." || str[t] == "NN" || str[t] == "NA") genotype_val.back()[mappingS[t-9]] = -1.0;
						else if (!gt_field) {
							if (field[idx_field][0] == '.') genotype_val.back()[mappingS[t-9]] = -1.0;
							else {
                                float dosage = atof(field[idx_field].c_str());
                                if (dosage < 0 || dosage > 2) LOG.error("DosareadGenotypesImpute2_yesIndex(fgen, fsam);ges must be between 0 and 2, check: " + field[idx_field]);
                                genotype_val.back()[mappingS[t-9]] = dosage;
							}
						} else {
							if (field[idx_field][0] == '.' || field[idx_field][2] == '.') genotype_val.back()[mappingS[t-9]] = -1.0;
							else {
								int a0 = atoi(field[idx_field].substr(0, 1).c_str());
								int a1 = atoi(field[idx_field].substr(2, 1).c_str());
                                int dosage = a0 + a1;
                                if (dosage < 0 || dosage > 2) LOG.error("Genotypes must be 00, 01, or 11, check: " + field[idx_field]);
                                genotype_val.back()[mappingS[t-9]] = dosage;
							}
						}
					}
				}
				n_includedG ++;
			} else n_excludedF ++;
		} else n_excludedG ++;
	}

	//Finalise
	genotype_count = n_includedG;
	LOG.println("  * region = " + regionGenotype.str());
	LOG.println("  * " + sutils::int2str(n_includedS) + " samples included");
	if (n_excludedS > 0) LOG.println("  * " + sutils::int2str(n_excludedS) + " samples excluded");
	if (n_missingS > 0) LOG.println("  * " + sutils::int2str(n_excludedS) + " samples excluded without phenotype data");
	LOG.println("  * " + sutils::int2str(n_includedG) + " sites included");
	if (n_excludedG > 0) LOG.println("  * " + sutils::int2str(n_excludedG) + " sites excluded");
	if (n_excludedF > 0) LOG.println("  * " + sutils::int2str(n_excludedF) + " sites excluded because of missing GT/DS field");
    if (n_includedG <= 0) LOG.error("No genotypes in this region: " + regionPhenotype.str());
}


void data::readGenotypesImpute2(string fgen, string fsam) {
	string buffer;
	vector<string> str;
	int n_includedG = 0;
	int n_excludedG = 0;
	int n_includedS = 0;
	int n_excludedS = 0;
	int n_missingS = 0;
	vector < int > mappingS;

	//Read samples
	LOG.println("\nReading sample list for genotype data in [" + fsam + "] in Imute2 format");
	ifile fds (fsam);
	getline(fds, buffer, '\n');
	getline(fds, buffer, '\n');
	while (getline(fds, buffer, '\n')) {
		sutils::tokenize(buffer, str);
		if (checkSample(str[1])) {
			int idx_sample = -1;
			for (int i = 0 ; i < sample_count && idx_sample < 0 ; i++) if (sample_id[i] == str[1]) idx_sample = i;
			mappingS.push_back(idx_sample);
			if (idx_sample >= 0) n_includedS ++;
			else n_missingS ++;
		}
	}
	fds.close();
	if (n_includedS != sample_count) LOG.error("Genotype data does not overlap with phenotype data, check your files!");

	//Initialise
	if (!futils::isFile(fgen+".tbi")) LOG.error("index file missing [" + fgen + ".tbi]");
	Tabix fd (fgen);
	if (!fd.setRegion(regionGenotype.str())) LOG.error("Failed to get region " + regionGenotype.str() + " in [" + fgen + "]");

	//Read genotypes
	while (fd.getNextLine(buffer)) {
		sutils::tokenize(buffer, str);
		if (checkGenotype(str[1])) {
			genotype_id.push_back(str[1]);
			genotype_chr.push_back(str[0]);
			genotype_pos.push_back(atoi(str[2].c_str()));
			genotype_val.push_back(vector < float > (sample_count, 0.0));
			for (int t = 5 ; t < str.size() ; t += 3) {
				int idx_sample = (t - 5) / 3;
				if (mappingS[idx_sample] >= 0) genotype_val.back()[mappingS[idx_sample]] = 2.0 * atof(str[t + 2].c_str()) + atof(str[t + 1].c_str());
			}
			n_includedG ++;
		} else n_excludedG ++;
	}

	//Finalise
	genotype_count = n_includedG;
	LOG.println("  * region = " + regionGenotype.str());
	LOG.println("  * " + sutils::int2str(n_includedS) + " samples included");
	if (n_excludedS > 0) LOG.println("  * " + sutils::int2str(n_excludedS) + " samples excluded");
	if (n_missingS > 0) LOG.println("  * " + sutils::int2str(n_excludedS) + " samples excluded without phenotype data");
	LOG.println("  * " + sutils::int2str(n_includedG) + " sites included");
	if (n_excludedG > 0) LOG.println("  * " + sutils::int2str(n_excludedG) + " sites excluded");
    if (n_includedG <= 0) LOG.error("No genotypes in this region: " + regionPhenotype.str() );
}