#include "data.h"


void data::runNominal(string fout, float cis_window) {
	long n_tests = 0;

	LOG.println("\nNominal pass, no permutation and report all p-values");
	LOG.println("  * W = [" + sutils::int2str(cis_window) + "]");
	LOG.println("  * O = [" + fout + "]");
	LOG.printC("\r  * 0.0 %");
	ofile fdo (fout);
	for (int p = 0 ; p < phenotype_count ; p ++) {
		for (int g = 0 ; g < genotype_count ; g ++) {
			int cisdistance = genotype_pos[g] - phenotype_start[p];
			if (abs(cisdistance) <= cis_window) {
				double corr = getCorrelation(g, p);
				fdo << phenotype_id[p];
				fdo << " " << genotype_id[g];
				fdo << " " << cisdistance;
				if (!spearman) fdo << " " << sutils::double2str(getBeta(g, p, corr), 3);
				fdo << " " << getPvalue(corr) << endl;
				n_tests ++;
			}
		}
		LOG.printC("\r  * " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + " %");
	}
	fdo.close();
	LOG.printlnC("\r  * " + sutils::int2str(n_tests) + " tested pairs");
	LOG.printlnL("  * " + sutils::int2str(n_tests) + " tested pairs");
	LOG.println("  * " + sutils::int2str(n_tests) + " association tests performed");
}
