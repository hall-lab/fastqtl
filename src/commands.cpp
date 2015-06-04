#include "data.h"

void data::writeCommands(string fout, int nChunks, int argc, char ** argv) {
	LOG.println("\nGenerate " + sutils::int2str(nChunks) + " commands in [" + fout + "]");

	vector < string > args;
	for (int a = 0 ; a < argc ; a ++) args.push_back(argv[a]);

	//map commands argument
	int idx_commands_arg = -1;
	for (int a = 0 ; a < argc ; a ++) if (args[a] == "--commands") idx_commands_arg = a;
	args.erase (args.begin() + idx_commands_arg, args.begin() + idx_commands_arg + 3);

	//map output argument
	int idx_output_arg = -1;
	string str_output_arg = "";
	for (int a = 0 ; a < argc ; a ++) {
		if (args[a] == "--out") {
			idx_output_arg = a;
			str_output_arg = string(args[a+1]);
		}
	}

	//write commands [loop over regions]
	ofile fd(fout);
	for (int c = 0 ; c < nChunks ; c ++) {
		setPhenotypeRegion(c);
		string region = getPhenotypeRegion(c);

		for (int a = 0 ; a < args.size() ; a ++ )
			if (a == idx_output_arg + 1) fd << " " << str_output_arg + "." + region;
			else fd << " " << args[a];

		fd << " --region " << region << endl;
	}
}
