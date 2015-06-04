#include "data.h"

int main(int argc, char ** argv) {
	data D;

	//-------------------------
	// 1.0. DECLARE ALL OPTIONS
	//-------------------------
	bpo::options_description opt_files ("\33[33mInput/Output files\33[0m");
	opt_files.add_options()
		("log,L", bpo::value< string >()->default_value("fastQTL_date_time_UUID.log"), "Screen output is copied in this file.")
		("vcf,V", bpo::value< string >(), "Genotypes in VCF format.")
		("gen,G", bpo::value< vector < string > >()->multitoken(), "Genotypes in Impute2 format.")
		("bed,B", bpo::value< string >(), "Phenotypes in BED format.")
		("cov,C", bpo::value< string >(), "Covariates in TXT format.")
		("region,R", bpo::value< string >(), "Region of interest.")
		("out,O", bpo::value< string >(), "Output file.");

	bpo::options_description opt_exclusion ("\33[33mExclusion/Inclusion files\33[0m");
	opt_exclusion.add_options()
		("exclude-samples", bpo::value< string >(), "List of samples to exclude.")
		("include-samples", bpo::value< string >(), "List of samples to include.")
		("exclude-sites", bpo::value< string >(), "List of sites to exclude.")
		("include-sites", bpo::value< string >(), "List of sites to include.")
		("exclude-phenotypes", bpo::value< string >(), "List of phenotypes to exclude.")
		("include-phenotypes", bpo::value< string >(), "List of phenotypes to include.")
		("exclude-covariates", bpo::value< string >(), "List of covariates to exclude.")
		("include-covariates", bpo::value< string >(), "List of covariates to include.");

	bpo::options_description opt_methods ("\33[33mMethods\33[0m");
	opt_methods.add_options()
		("spearman", "Use Spearman correlation for association testing [Default is Pearson].")
		("normal", "Normal transform the phenotypes.")
		("window", bpo::value< float >()->default_value(1e6), "Cis-window size.")
		("nominal", "Perform a nominal pass on the data without permutations and report all QTL-MP pairs.")
		("permute", bpo::value< vector < int > >()->multitoken(), "Perform permutations and report only best QTL-MP pairs.")
		("seed", bpo::value< int >()->default_value(time(NULL)), "Random number seed. Useful to replicate runs.");

	bpo::options_description opt_parallel ("\33[33mParallelization\33[0m");
	opt_parallel.add_options()
		("chunk", bpo::value< vector < int > >()->multitoken(), "Specify which chunk needs to be processed")
		("commands", bpo::value< vector < string > >()->multitoken(), "Generates all commands");

	bpo::options_description descriptions;
	descriptions.add(opt_files).add(opt_exclusion).add(opt_methods).add(opt_parallel);

	//-------------------
	// 2.0. PARSE OPTIONS
	//-------------------
	bpo::variables_map options;
	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing command line :" << string(e.what()) << endl;
		exit(0);
	}

	//-----------------------
	// 3.0. PRINT HEADER/HELP
	//-----------------------
	cout << endl;
	cout << "\33[33mF\33[0mast \33[33mQTL\33[0m" << endl;
	cout << "  * Authors : Olivier DELANEAU, Halit ONGEN, Alfonso BUIL & Manolis DERMITZAKIS" << endl;
	cout << "  * Contact : olivier.delaneau@gmail.com" << endl;
	cout << "  * Webpage : http://funpopgen.unige.ch/" << endl;
	cout << "  * Version : v1.0" << endl;
	if (options.count("help")) { cout << descriptions<< endl; exit(1); }

	//--------------
	// 4.0. LOG FILE
	//--------------
	struct timeval start_time, stop_time;
	gettimeofday(&start_time, 0);
	START_DATE = time(0);
	time(&START_DATE);
	string logfile = "fastQTL_" + sutils::date2str(&START_DATE, "%d%m%Y_%Hh%Mm%Ss") + "_" + putils::getRandomID() + ".log";
	if (!options["log"].defaulted()) logfile = options["log"].as < string > ();
	if (!LOG.open(logfile)) {
		cerr << "Impossible to open log file[" << options["log"].as < string > () << "] check writing permissions!" << endl;
		exit(1);
	}

	//--------------------------------
	// 4.0. CHECK /AND VERBOSE OPTIONS
	//--------------------------------

	// 4.1 CHECK INPUT/OUTPUT FILES
	int nInputGenotypes = options.count("vcf") + options.count("impute");
	if (nInputGenotypes != 1) LOG.error("Genotype data needs to be specified with --vcf [file.vcf] OR --gen [file.gen] [file.sample]");
	if (options.count("vcf") && !futils::isFile(options["vcf"].as < string > ())) LOG.error(options["vcf"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("gen")) {
		vector < string > ifiles = options["gen"].as < vector < string > > ();
		if (ifiles.size() != 2) LOG.error("--gen requires 2 arguments: [file.gen] [file.sample]");
		if (!futils::isFile(ifiles[0])) LOG.error(ifiles[0] + " is impossible to open, check file existence or reading permissions");
		if (!futils::isFile(ifiles[1])) LOG.error(ifiles[1] + " is impossible to open, check file existence or reading permissions");
	}
	if (!options.count("bed")) LOG.error("Phenotype data needs to be specified with --bed [file.bed]");
	if (!futils::isFile(options["bed"].as < string > ())) LOG.error(options["bed"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("cov") && !futils::isFile(options["cov"].as < string > ())) LOG.error(options["cov"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (!options.count("out")) LOG.error("Output needs to be specified with --out [file.out]");
	if (!options.count("commands") && !futils::createFile(options["out"].as < string > ())) LOG.error(options["out"].as < string > () + " is impossible to create, check file existence or writing permissions");

	// 4.2. CHECK INCLUSION/EXCLUSION FILES
	if (options.count("exclude-samples") && !futils::isFile(options["exclude-samples"].as < string > ())) LOG.error(options["exclude-samples"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("include-samples") && !futils::isFile(options["include-samples"].as < string > ())) LOG.error(options["include-samples"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("exclude-sites") && !futils::isFile(options["exclude-sites"].as < string > ())) LOG.error(options["exclude-sites"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("include-sites") && !futils::isFile(options["include-sites"].as < string > ())) LOG.error(options["include-sites"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("exclude-phenotypes") && !futils::isFile(options["exclude-phenotypes"].as < string > ())) LOG.error(options["exclude-phenotypes"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("include-phenotypes") && !futils::isFile(options["include-phenotypes"].as < string > ())) LOG.error(options["include-phenotypes"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("exclude-covariates") && !futils::isFile(options["exclude-covariates"].as < string > ())) LOG.error(options["exclude-covariates"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("include-covariates") && !futils::isFile(options["include-covariates"].as < string > ())) LOG.error(options["include-covariates"].as < string > () + " is impossible to open, check file existence or reading permissions");

	// 4.3. CHECK METHODS/PARAMETERS [I REALLY DON'T LIKE THIS SECTION, TOO MESSY]
	LOG.println("\nParameters:");
	int nCisAnalysis = options.count("nominal") + options.count("permute");
	if (nCisAnalysis != 1) LOG.error("Analysis type is wrong/unspecified. This can be either --nominal or --permute");
	else if (options.count("nominal")) LOG.println("  * All SNP-MP pairs' statistics are reported");
	else LOG.println("  * Only best SNP-MP pairs' statistics are reported");
	if (options.count("permute")) {
		vector < int > nPerm = options["permute"].as < vector < int > > ();
		if (nPerm.size() > 2 || nPerm.size() < 1) LOG.error ("Option --permute takes 1 or 2 arguments");
		if (nPerm.size() == 1) {
			if (nPerm[0] <= 0) LOG.error("Permutation number needs to be positive");
			LOG.println("  * " + sutils::int2str(nPerm[0]) + " permutations");
		} else {
			if (nPerm[0] <= 0 || nPerm[1] <= 0) LOG.error("Permutation number needs to be positive");
			LOG.println("  * Between " + sutils::int2str(nPerm[0]) + " and " + sutils::int2str(nPerm[1]) + " permutations");
		}
	}
	if (options["seed"].as < int > () < 0) LOG.error("Random number generator needs a positive seed value");
	else srand(options["seed"].as < int > ());
	LOG.println("  * Seed = " + sutils::int2str(options["seed"].as < int > ()));
	if (options["window"].as < float > () <= 0) LOG.error ("Incorrect argument for option --cis-windows (null or negative value)");
	LOG.println("  * Cis-window = " + sutils::int2str(options["window"].as < float > ()));
	if (options.count("spearman")) LOG.println("  * Association testing using Spearman correlation");
	else LOG.println("  * Association testing using linear regression");
	int nParallel = options.count("chunk") + options.count("commands") + options.count("region");
	if (nParallel != 1) LOG.error("one option [--region, --chunk, --commands] needs to be specified");
	else {
		if (options.count("chunk")) {
			vector < int > nChunk = options["chunk"].as < vector < int > > ();
			if (nChunk.size() != 2) LOG.error ("--chunk needs 2 integer arguments");
			LOG.println ("  * Chunk processed " + sutils::int2str(nChunk[0]) + " / " + sutils::int2str(nChunk[1]));
		} else if (options.count("commands")) {
			vector < string > nCommands = options["commands"].as < vector < string > > ();
			if (nCommands.size() != 2) LOG.error ("--commands needs 2 arguments");
			LOG.println ("  * " + nCommands[0] + " commands output in [" + nCommands[1] +"]");
		}
	}

	//----------------------------------
	// 5.0. READ EXCLUDE / INCLUDE FILES
	//----------------------------------
	if (options.count("exclude-samples")) D.readSamplesToExclude(options["exclude-samples"].as < string > ());
	if (options.count("include-samples")) D.readSamplesToInclude(options["include-samples"].as < string > ());
	if (options.count("exclude-sites")) D.readGenotypesToExclude(options["exclude-sites"].as < string > ());
	if (options.count("include-sites")) D.readGenotypesToInclude(options["include-sites"].as < string > ());
	if (options.count("exclude-phenotypes")) D.readPhenotypesToExclude(options["exclude-phenotypes"].as < string > ());
	if (options.count("include-phenotypes")) D.readPhenotypesToInclude(options["include-phenotypes"].as < string > ());
	if (options.count("exclude-covariates")) D.readCovariatesToExclude(options["exclude-phenotypes"].as < string > ());
	if (options.count("include-covariates")) D.readCovariatesToInclude(options["include-phenotypes"].as < string > ());

	if (options.count("commands")) {
		//-----------------------
		// 6.0. GENERATE COMMANDS
		//-----------------------
		int nChunks = atoi(options["commands"].as < vector < string > > ()[0].c_str());
		D.scanPhenotypes(options["bed"].as < string > ());
		D.clusterizePhenotypes(nChunks);
		D.writeCommands(options["commands"].as < vector < string > > ()[1], nChunks, argc, argv);
	} else {
		//----------------
		// 7.0. SET REGION
		//----------------
		if (options.count("chunk")) {
			D.scanPhenotypes(options["bed"].as < string > ());
			D.clusterizePhenotypes(options["chunk"].as < vector < int > > ()[1]);
			D.setPhenotypeRegion(options["chunk"].as < vector < int > > ()[0] - 1);
			D.clear();
		} else D.setPhenotypeRegion(options["region"].as < string > ());
		D.deduceGenotypeRegion(options["window"].as < float > ());

		//----------------
		// 8.0. READ FILES
		//----------------
		D.readPhenotypes(options["bed"].as < string > ());
		if (options.count("vcf")) D.readGenotypesVCF(options["vcf"].as < string > ());
		if (options.count("imp")) D.readGenotypesImpute2(options["imp"].as < vector < string > > ()[0], options["imp"].as < vector < string > > ()[1]);
		if (options.count("cov")) D.readCovariates(options["cov"].as < string > ());

		//-------------------------
		// 9.0. INITIALIZE ANALYSIS
		//-------------------------
		D.imputeGenotypes();
		D.imputePhenotypes();
		if (options.count("cov")) {
			D.correctGenotypesForCovariates();
			D.correctPhenotypesForCovariates();
		}
		if (options.count("spearman")) {
			D.spearman = true;
			D.rankTranformPhenotypes();
			D.rankTranformGenotypes();
		}
		if (options.count("normal")) D.normalTranformPhenotypes();
		D.computeSDGenotypes();
		D.computeSDPhenotypes();
		D.normalisePhenotypes();
		D.normaliseGenotypes();

		//-------------------
		// 10.0. RUN ANALYSIS
		//-------------------
		if (options.count("nominal")) D.runNominal(options["out"].as < string > (), options["window"].as < float > ());
		if (options.count("permute")) {
			vector < int > nPerm = options["permute"].as < vector < int > > ();
			if (nPerm.size() == 1) D.runPermutation(options["out"].as < string > (), options["window"].as < float > (), -1, nPerm[0]);
			else D.runPermutation(options["out"].as < string > (), options["window"].as < float > (), nPerm[0], nPerm[1]);
		}
	}

	//------------------
	// 11.0. TERMINATION
	//------------------
	D.clear();
	gettimeofday(&stop_time, 0);
	int n_seconds = (int)floor(stop_time.tv_sec - start_time.tv_sec);
	LOG.println("\nRunning time: " + sutils::int2str(n_seconds) + " seconds");
	LOG.close();
}
