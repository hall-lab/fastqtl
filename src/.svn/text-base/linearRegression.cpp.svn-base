#include "utils/Eigen/Dense"
#include "utils/Eigen/LU"
#include "utils/auxiliary.h"

#include "data.h"

using namespace Eigen;

void data::correctGenotypesForCovariates() {
    LOG.println("\nRunning linear regression on genotypes with the provided covariates");

    //Initialise the genotypes matrix rows = variants columns = samples
    MatrixXf GenotypeMatrix (genotype_count, sample_count);
    for (int g = 0 ; g < genotype_count ; g++) GenotypeMatrix.row(g) = VectorXf::Map(&genotype_val[g][0], sample_count);

    //Create the covariate matrix with rows = samples and columns = covariates
    MatrixXf CovariateMatrix (sample_count, covariate_count + 1);
    CovariateMatrix.col(0) = VectorXf::Ones(sample_count); //Add a column of 1's for the y-intercept (mean) of the linear regression
    for (int s = 0 ; s < sample_count ; ++s) for (int c = 0 ; c < covariate_count ; ++c) CovariateMatrix(s,c+1) = covariate_val[c][s];

    //Transpose covariate matrix
    MatrixXf TransposedCovariateMatrix = CovariateMatrix.transpose();

    //Actual regression starts here
    MatrixXf mul = TransposedCovariateMatrix * CovariateMatrix;
    MatrixXf inverse = mul.fullPivLu().inverse();
	MatrixXf tmp = inverse * TransposedCovariateMatrix;
	for (int g = 0; g < genotype_count ; g ++) {
		MatrixXf beta = tmp * GenotypeMatrix.row(g).transpose();
		MatrixXf e = GenotypeMatrix.row(g) - (beta.transpose() * TransposedCovariateMatrix);
		for (int s = 0; s < sample_count ; ++s) genotype_val[g][s]= e(0,s);
		LOG.printC("\r  * " + sutils::double2str((g+1) * 100.0 / genotype_count, 1) + " %");
	}
    LOG.printlnC("\r  * " + sutils::int2str(genotype_count) + " genotypes normalized");
	LOG.printlnL("  * " + sutils::int2str(genotype_count) + " genotypes normalized");
}

void data::correctPhenotypesForCovariates() {
    LOG.println("\nRunning linear regression on phenotypes with the provided covariates");

    //Initialise the genotypes matrix rows = variants columns = samples
    MatrixXf PhenotypeMatrix (phenotype_count, sample_count);
    for (int p = 0 ; p < phenotype_count ; p++) PhenotypeMatrix.row(p) = VectorXf::Map(&phenotype_val[p][0], sample_count);

    //Create the covariate matrix with rows = samples and columns = covariates
    MatrixXf CovariateMatrix (sample_count, covariate_count + 1);
    CovariateMatrix.col(0) = VectorXf::Ones(sample_count); //Add a column of 1's for the y-intercept (mean) of the linear regression
    for (int s = 0 ; s < sample_count ; ++s) for (int c = 0 ; c < covariate_count ; ++c) CovariateMatrix(s,c+1) = covariate_val[c][s];

    //Transpose covariate matrix
    MatrixXf TransposedCovariateMatrix = CovariateMatrix.transpose();

    //Actual regression starts here
    MatrixXf mul = TransposedCovariateMatrix * CovariateMatrix;
    MatrixXf inverse = mul.fullPivLu().inverse();
	MatrixXf tmp = inverse * TransposedCovariateMatrix;
	for (int p = 0; p < phenotype_count ; p ++) {
		MatrixXf beta = tmp * PhenotypeMatrix.row(p).transpose();
		MatrixXf e = PhenotypeMatrix.row(p) - (beta.transpose() * TransposedCovariateMatrix);
		for (int s = 0; s < sample_count ; ++s) phenotype_val[p][s]= e(0,s);
		LOG.printC("\r  * " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + " %");
	}
    LOG.printlnC("\r  * " + sutils::int2str(phenotype_count) + " phenotypes normalized");
	LOG.printlnL("  * " + sutils::int2str(phenotype_count) + " phenotypes normalized");
}

