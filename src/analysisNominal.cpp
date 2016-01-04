//FastQTL: Fast and efficient QTL mapper for molecular phenotypes
//Copyright (C) 2015 Olivier DELANEAU, Alfonso BUIL, Emmanouil DERMITZAKIS & Olivier DELANEAU
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "data.h"


void data::runNominal(string fout, double threshold) {

	//0. Prepare genotypes
	vector < double > genotype_sd = vector < double > (genotype_count, 0.0);
	vector < double > phenotype_sd = vector < double > (phenotype_count, 0.0);
	if (covariate_count > 0) {
		LOG.println("\nCorrecting genotypes & phenotypes for covariates");
		covariate_engine->residualize(genotype_orig);
		covariate_engine->residualize(phenotype_orig);
	}
	for (int g = 0 ; g < genotype_count ; g ++) genotype_sd[g] = RunningStat(genotype_orig[g]).StandardDeviation();
	for (int p = 0 ; p < phenotype_count ; p ++) phenotype_sd[p] = RunningStat(phenotype_orig[p]).StandardDeviation();
	normalize(genotype_orig);
	normalize(phenotype_orig);

	//1. Loop over phenotypes
	ofile fdo (fout);
	for (int p = 0 ; p < phenotype_count ; p ++) {

		LOG.println("\nProcessing gene [" + phenotype_id[p] + "]");

		//1.1. Enumerate all genotype-phenotype pairs within cis-window
		vector < int > targetGenotypes, targetDistances;
		for (int g = 0 ; g < genotype_count ; g ++) {
		  int cisdistance;
		  int startdistance = genotype_pos[g] - phenotype_start[p];
		  int enddistance = genotype_end[g] - phenotype_start[p];

		  // for INVs ignore the span and define the cisdistance
		  // as the distance from the breakpoints to the phenotype_start
		  if (genotype_vartype[g].compare("INV") == 0) {
		    if (abs(startdistance) <= abs(enddistance))
		      cisdistance = startdistance;
		    else
		      cisdistance = enddistance;
		  }

		  // for the variants with span (DEL, DUP, MEI), cisdistance is zero
		  // if the phenotype_start falls within the span, and the distance to
		  // the closest edge otherwise
		  // BNDs get processed here as well, but their END coordinate is the
		  // same as the START coordinate.
		  else {
		    if (startdistance < 0 && enddistance > 0) { // if gene is within SV, then cis distance is 0
		      cisdistance = 0;
		    }
		    else if (startdistance >= 0)
		      cisdistance = startdistance;
		    else
		      cisdistance = enddistance;
		  }

		  if (abs(cisdistance) <= cis_window) {
		    targetGenotypes.push_back(g);
		    targetDistances.push_back(cisdistance);
		  }
		}
		LOG.println("  * Number of variants in cis = " + sutils::int2str(targetGenotypes.size()));

		//1.2. Nominal pass: scan cis-window & compute statistics
		for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
			double corr = getCorrelation(genotype_orig[targetGenotypes[g]], phenotype_orig[p]);
			double pval = getPvalue(corr, sample_count - 2 - covariate_count);
			double slope = getSlope(corr, phenotype_sd[p], genotype_sd[targetGenotypes[g]]);
			if (pval <= threshold ) {
				fdo << phenotype_id[p];
				fdo << " " << genotype_id[targetGenotypes[g]];
				fdo << " " << targetDistances[g];
				fdo << " " << pval;
				fdo << " " << slope;
				fdo << endl;
			}
		}
		LOG.println("  * Progress = " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + "%");
	}
	fdo.close();
}
