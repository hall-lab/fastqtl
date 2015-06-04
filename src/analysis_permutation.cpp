#include "data.h"

// minPermutations = -1 means fixed number of permutations!
void data::runPermutation(string fout, float cis_window, int minPermutations, int maxPermutations) {
	long n_pairs = 0;
	long n_tests = 0;

	LOG.println("\nScan cis-window, permute phenotypes and report nominal, empirical and estimated p-values");
	LOG.println("  * W = [" + sutils::int2str(cis_window) +"]");
	if (minPermutations < 0) LOG.println("  * P = [" + sutils::int2str(maxPermutations) +"] (fixed number)");
	else LOG.println("  * P = [" + sutils::int2str(minPermutations) + ", " + sutils::int2str(maxPermutations) +"] (variable number)");
	LOG.println("  * O = [" + fout + "]");
	LOG.printC("\r  * 0.0 %");

	ofile fdo (fout);
	for (int p = 0 ; p < phenotype_count ; p ++) {

		//STEP1: MAP GENOTYPE-PHENOTYPE PAIRS
		vector < int > targetGenotypes, targetDistances;
		for (int g = 0 ; g < genotype_count ; g ++) {
			int cisdistance = genotype_pos[g] - phenotype_start[p];
			if (abs(cisdistance) <= cis_window) {
				targetGenotypes.push_back(g);
				targetDistances.push_back(cisdistance);
				n_pairs ++;
			}
		}

		//STEP2: COMPUTE CORRELATION COEFFICIENT IN CIS WINDOW
		double bestC = 0.0;
		int bestD = 1000000000, bestI = -1;
		for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
			double corr = getCorrelation(targetGenotypes[g], p);
			if (abs(corr) > abs(bestC) || (abs(corr) == abs(bestC) && abs(targetDistances[g]) < bestD)) {
				bestC = corr;
				bestD = targetDistances[g];
				bestI = targetGenotypes[g];
			}
			n_tests ++;
		}

		//STEP3: PERFORM PERMUTATIONS
		vector < float > phenotype_shuffled = phenotype_val[p];
		int nBetterCorrelation = 0, countPermutations = 0;
		double sumBestPvalues = 0.0;
		double sumBestPvalues2 = 0.0;
		while (countPermutations < maxPermutations && (minPermutations < 0 || (minPermutations > 0 && nBetterCorrelation < minPermutations))) {
			double bestCperm = 0;
			random_shuffle(phenotype_shuffled.begin(), phenotype_shuffled.end());
			for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
				double corr = getCorrelation(targetGenotypes[g], phenotype_shuffled);
				if (abs(corr) > abs(bestCperm)) bestCperm = corr;
				n_tests ++;
			}
			if (abs(bestCperm) >= abs(bestC)) nBetterCorrelation++;
			double pval = getPvalue(bestCperm);
			sumBestPvalues += pval;
			sumBestPvalues2 += pval * pval;
			countPermutations++;
		}

		//STEP3.1: DERIVE STATISTICS FOR BETA APPROXIMATION
		double mean = sumBestPvalues / countPermutations;
		double var = sumBestPvalues2 / countPermutations - mean * mean;
		//double shape1 = ((1 - mean) / var - 1 / mean) * mean * mean;
		double shape1 = 1.0;
		double shape2 = shape1 * (1 / mean - 1);

		//STEP4: WRITE RESULTS
		fdo << phenotype_id[p];
		if (bestI >= 0) {
			double pval = getPvalue(bestC);
			fdo << " " << genotype_id[bestI];
			fdo << " " << bestD;
			if (!spearman) fdo << " " << sutils::double2str(getBeta(bestI, p, bestC), 3);
			fdo << " " << pval;
			fdo << " " << (nBetterCorrelation + 1) * 1.0 / (countPermutations + 1.0);
			fdo << " " << portR::pbeta(pval, shape1, shape2,1,0);
			fdo << endl;
		} else fdo << (spearman?" NA NA NA NA NA":" NA NA NA NA NA NA")  << endl;

		//STEP5: VERBOSE PROGRESS
		LOG.printC("\r  * " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + " %");
		LOG.printL("  * " + sutils::int2str(p+1) + "/" + sutils::int2str(phenotype_count));
	}
	LOG.printlnC("\r  * " + sutils::int2str(n_pairs) + " tested pairs");
	LOG.printlnL("  * " + sutils::int2str(n_pairs) + " tested pairs");
	LOG.println("  * " + sutils::int2str(n_tests) + " association tests performed");
	fdo.close();
}
