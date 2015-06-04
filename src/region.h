#ifndef _REGION_H
#define _REGION_H

#define POS_MIN 0
#define POS_MAX 1000000000

class region {
public:
	string chr;
	int start;
	int end;

	region() {
		chr = "chrNA";
		start = POS_MIN;
		end = POS_MAX;
	}

	string str() {
		ostringstream s2( stringstream::out );
		s2 << chr << ":" << start << "-" << end;
		return s2.str();
	}

	bool check (string _chr, int _pos) {
		if (_chr == chr && _pos >= start && _pos < end) return true;
		else return false;
	}

	bool set(string r) {
		vector < string > chr_split, pos_split;
		sutils::tokenize(r, chr_split, ":");
		if (chr_split.size() == 2) {
			chr = chr_split[0];
			sutils::tokenize(chr_split[1], pos_split, "-");
			if (pos_split.size() == 2) {
				start = atoi(pos_split[0].c_str());
				end = atoi(pos_split[1].c_str());
				return true;
			} else return false;
		} else if (chr_split.size() == 1) {
			chr = chr_split[0];
			start = POS_MIN;
			end = POS_MAX;
			return true;
		} else return false;
	}
};

#endif
