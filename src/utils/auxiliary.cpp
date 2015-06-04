/*Copyright (C) 2012  Halit Ongen, Emmanouil T. Dermitzakis

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include "auxiliary.h"

string getTime(){
	time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
	stringstream time;
	time << setw(2) << setfill('0') << now->tm_hour << ":" << setw(2) << setfill('0') << now->tm_min << ":" << setw(2) << setfill('0') << now->tm_sec;
	return time.str();
}

string extractDirectory( const string& path ){
	return path.substr( 0, path.find_last_of( '/' ) +1 );
}

string extractFilename( const string& path ){
	return path.substr( path.find_last_of( '/' ) +1 );
}

string changeExtension( const string& path, const string& ext ){
	string filename = extractFilename( path );
	return extractDirectory( path ) +filename.substr( 0, filename.find_last_of( '.' ) ) +ext;
}

string removeExtension( const string& path ){
	string filename = extractFilename( path );
	return extractDirectory( path ) +filename.substr( 0, filename.find_last_of( '.' ) );
}

void movefile(const string& from, const string& to) {
	copyfile(from,to);
	remove(from.c_str());
}

void copyfile(const string& from, const string& to) {
	ifstream ifs(from.c_str(), ios::in | ios::binary);
	if(!ifs.is_open()){
        cerr<<"ERROR: Cannot open done file: " << from << endl;
        exit(1);
	}
	ofstream ofs(to.c_str(), ios::out | ios::binary);
	if(!ofs.is_open()){
        cerr<<"ERROR: Cannot open done file: " << to << endl;
        exit(1);
	}
	ofs << ifs.rdbuf();
	ifs.close();
	ofs.close();
}

bool dirExists(const string& path){
	struct stat sb;
	return (stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) ? true : false;
}

bool fileExists(const string& filename){
	ifstream ifile(filename.c_str());
	return ifile;
}

bool isPositiveInt(char * input){
	int i = 0;
	do{
		if(!isdigit(input[i]))
			return false;
		i++;
	}while(input[i] != '\0');
	return true;
}

bool isPositiveInt(const char * input){
	int i = 0;
	do{
		if(!isdigit(input[i]))
			return false;
		i++;
	}while(input[i] != '\0');
	return true;
}

bool isNumber(char * input){
	int i = 0;
	int ecount = 0;
	int dashcount=0;
	int dotcount=0;
	do{
		if(!isdigit(input[i]) && input[i]!= '.' && input[i]!= '-' && input[i]!= '+' && input[i]!= 'E' && input[i]!= 'e')
			return false;
		if((input[i]== 'E' || input[i]== 'e') && (i ==0 || input[i+1] == '\0' || !isdigit(input[i-1]) || ( !isdigit(input[i+1]) && input[i+1]!='-' &&  input[i+1]!='+')))
			return false;
		if((input[i]== '+' || input[i]== '-' || input[i]=='.' ) && (input[i+1] == '\0' || !isdigit(input[i+1])))
			return false;
		if(input[i]== 'E' || input[i]== 'e') ecount++;
		if(input[i]== '+' || input[i]== '-') dashcount++;
		if(input[i]=='.') dotcount++;
		i++;
	}while(input[i] != '\0');
	if (ecount > 1 || dotcount > 1 || dashcount>2)
		return false;
	return true;
}

bool isNumber(const char * input){
	int i = 0;
	int ecount = 0;
	int dashcount=0;
	int dotcount=0;
	do{
		if(!isdigit(input[i]) && input[i]!= '.' && input[i]!= '-' && input[i]!= '+' && input[i]!= 'E' && input[i]!= 'e')
			return false;
		if((input[i]== 'E' || input[i]== 'e') && (i ==0 || input[i+1] == '\0' || !isdigit(input[i-1]) || ( !isdigit(input[i+1]) && input[i+1]!='-' &&  input[i+1]!='+')))
			return false;
		if((input[i]== '+' || input[i]== '-' || input[i]=='.' ) && (input[i+1] == '\0' || !isdigit(input[i+1])))
			return false;
		if(input[i]== 'E' || input[i]== 'e') ecount++;
		if(input[i]== '+' || input[i]== '-') dashcount++;
		if(input[i]=='.') dotcount++;
		i++;
	}while(input[i] != '\0');
	if (ecount > 1 || dotcount > 1 || dashcount>2)
		return false;
	return true;
}

int FileRead( istream & is, vector <char> & buff ) {
    is.read( &buff[0], buff.size() );
    return is.gcount();
}

uint CountLines( const vector <char> & buff, int sz ) {
    uint newlines = 0;
    const char * p = &buff[0];
    for ( int i = 0; i < sz; i++ ) {
    	if ( p[i] == '\n' ) {
    		newlines++;
    	}
    }
    return newlines;
}

string concurrentCopyTo (string & file, string & dir){
	string tempName = dir + "/" + extractFilename(file);
	string tempDone = tempName + ".cpdone";
	//Check if file exists
	if(fileExists(tempName)){
		uint seconds_waited =0;
		while(!fileExists(tempDone)){
			seconds_waited += 60;
			sleep(60);
			if (seconds_waited > 86400){
				cerr << "ERROR: We waited over 24 hrs for copying of " << file << " to " << dir << " exiting" << endl;
				exit(200);
			}
		}
	}else{
		copyfile(file,tempName);
        // read, write, and execute permissions for the owner, and no permissions for group and others
        my_chmod(tempName.c_str(),S_IRWXU);
		ofstream doneFile(tempDone.c_str());
		if (doneFile.is_open()){
			doneFile << "DO NOT DELETE THIS BEFORE DELETING " << tempName << endl;
			doneFile.close();
		}else{
			cerr<<"ERROR: Cannot open done file: " << doneFile << endl;
			exit(1);
		}
	}
	return tempName;
}

bool isNumeric( const char* pszInput, int nNumberBase )
{
	istringstream iss( pszInput );
    
	if ( nNumberBase == 10 )
	{
		double dTestSink;
		iss >> dTestSink;
	}
	else if ( nNumberBase == 8 || nNumberBase == 16 )
	{
		int nTestSink;
		iss >> ( ( nNumberBase == 8 ) ? oct : hex ) >> nTestSink;
	}
	else
		return false;
    
	// was any input successfully consumed/converted?
	if ( ! iss )
		return false;
    
	// was all the input successfully consumed/converted?
	return ( iss.rdbuf()->in_avail() == 0 );
}

