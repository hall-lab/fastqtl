#include "data.h"
#include "utils/tabix.hpp"

void data::readPhenotypes(string fbed) {
	string buffer;
	vector < string > str;
	int n_includedP = 0;
	int n_excludedP = 0;
	int n_includedS = 0;
	int n_excludedS = 0;
	vector < int > mappingS;

	//Initialise
	LOG.println("\nReading phenotype data in [" + fbed + "]");
	if (!futils::isFile(fbed + ".tbi")) LOG.error("index file missing [" + fbed + ".tbi]");
	Tabix fd (fbed);

	//Read samples ids
	fd.getLastHeader(buffer);
	sutils::tokenize(buffer, str);
    if (str.size() < 5) LOG.error("Wrong BED file format 1");
	for (int t = 4 ; t < str.size() ; t ++) {
		if (checkSample(str[t])) {
			sample_id.push_back(str[t]);
			mappingS.push_back(n_includedS);
			n_includedS ++;
		} else {
			mappingS.push_back(-1);
			n_excludedS ++;
		}
	}
	sample_count = sample_id.size();

	//Read phenotypes
	if (!fd.setRegion(regionPhenotype.str())) LOG.error("Failed to get region " + regionPhenotype.str() + " in [" + fbed + "]");
	while (fd.getNextLine(buffer)) {
		sutils::tokenize(buffer, str);
		if (str.size() != sample_count + 4) LOG.error("Wrong BED file format 2 = " +sutils::int2str(str.size()));
		if (checkPhenotype(str[3])) {
			phenotype_id.push_back(str[3]);
			phenotype_chr.push_back(str[0]);
			phenotype_start.push_back(atoi(str[1].c_str()) + 1); //convert to 1-based, tabix works on 1-based coordinates and all other files are 1-based
			phenotype_end.push_back(atoi(str[2].c_str()) + 1); //convert to 1-based, tabix works on 1-based coordinates and all other files are 1-based
			phenotype_val.push_back(vector < float > (sample_count, 0.0));
			for (int t = 4 ; t < str.size() ; t ++) {
				if (mappingS[t-4] >= 0) {
					if (str[t] == "NA") phenotype_val.back()[mappingS[t-5]] = ___NA___;
					else if (!isNumeric(str[t].c_str())) LOG.error("Phenotype encountered is not a number, check: [" + str[t] + "]");
					else phenotype_val.back()[mappingS[t-4]] = atof(str[t].c_str());
				}
			}
			n_includedP++;
		} else n_excludedP ++;
	}

	//Finalise
	phenotype_count = phenotype_id.size();
	LOG.println("  * region = " + regionPhenotype.str());
	LOG.println("  * " + sutils::int2str(n_includedS) + " samples included");
	if (n_excludedS > 0) LOG.println("  * " + sutils::int2str(n_excludedS) + " samples excluded");
	LOG.println("  * " + sutils::int2str(n_includedP) + " phenotypes included");
	if (n_excludedP > 0) LOG.println("  * " + sutils::int2str(n_excludedP) + " phenotypes excluded");
    if (phenotype_count == 0) LOG.error("No phenotypes in this region");
}

void data::scanPhenotypes(string fbed) {
	string buffer;
	vector < string > str;
	int n_includedP = 0;
	int n_excludedP = 0;

	//Initialise
	LOG.println("\nScanning phenotype data in [" + fbed + "]");
	if (!futils::isFile(fbed + ".tbi")) LOG.error("index file missing [" + fbed + ".tbi]");
	Tabix fd (fbed);

	//Read lines
	fd.getLastHeader(buffer);
	while (fd.getNextLine(buffer)) {
		sutils::tokenize(buffer, str, 4);
		if (str.size() != 4) LOG.error("Wrong BED file format");
		if (checkPhenotype(str[3])) {
			phenotype_id.push_back(str[3]);
			phenotype_chr.push_back(str[0]);
			phenotype_start.push_back(atoi(str[1].c_str()) + 1); //convert to 1-based, tabix works on 1-based coordinates and all other files are 1-based
			phenotype_end.push_back(atoi(str[2].c_str()) + 1); //convert to 1-based, tabix works on 1-based coordinates and all other files are 1-based
			phenotype_val.push_back(vector < float > (sample_count, 0.0));
			n_includedP++;
		} else n_excludedP ++;
	}

	//Finalise
	phenotype_count = phenotype_id.size();
	LOG.println("  * " + sutils::int2str(n_includedP) + " phenotypes included");
	if (n_excludedP > 0) LOG.println("  * " + sutils::int2str(n_excludedP) + " phenotypes excluded");
    if (phenotype_count == 0) LOG.error("No phenotypes in this region");
}
