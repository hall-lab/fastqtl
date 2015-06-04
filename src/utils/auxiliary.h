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

#ifndef __AUXILIARY_H_INCLUDED__
#define __AUXILIARY_H_INCLUDED__


#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <vector>
#include <cstddef>
#include <ctype.h>
#include <utility>
#include <algorithm>
#include <map>
#include <unistd.h>
#include <iostream>
#include <locale>

using namespace std;
typedef unsigned int uint;

#ifndef S_ISDIR
#define S_ISDIR(mode)  (((mode) & S_IFMT) == S_IFDIR)
#endif

//From http://stackoverflow.com/questions/592448/c-how-to-set-file-permissions-cross-platform
#ifdef _WIN32
#   include <io.h>

typedef int mode_t;

/// @Note If STRICT_UGO_PERMISSIONS is not defined, then setting Read for any
///       of User, Group, or Other will set Read for User and setting Write
///       will set Write for User.  Otherwise, Read and Write for Group and
///       Other are ignored.
///
/// @Note For the POSIX modes that do not have a Windows equivalent, the modes
///       defined here use the POSIX values left shifted 16 bits.

static const mode_t S_ISUID      = 0x08000000;           ///< does nothing
static const mode_t S_ISGID      = 0x04000000;           ///< does nothing
static const mode_t S_ISVTX      = 0x02000000;           ///< does nothing
static const mode_t S_IRUSR      = mode_t(_S_IREAD);     ///< read by user
static const mode_t S_IWUSR      = mode_t(_S_IWRITE);    ///< write by user
static const mode_t S_IXUSR      = 0x00400000;           ///< does nothing
#   ifndef STRICT_UGO_PERMISSIONS
static const mode_t S_IRGRP      = mode_t(_S_IREAD);     ///< read by *USER*
static const mode_t S_IWGRP      = mode_t(_S_IWRITE);    ///< write by *USER*
static const mode_t S_IXGRP      = 0x00080000;           ///< does nothing
static const mode_t S_IROTH      = mode_t(_S_IREAD);     ///< read by *USER*
static const mode_t S_IWOTH      = mode_t(_S_IWRITE);    ///< write by *USER*
static const mode_t S_IXOTH      = 0x00010000;           ///< does nothing
#   else
static const mode_t S_IRGRP      = 0x00200000;           ///< does nothing
static const mode_t S_IWGRP      = 0x00100000;           ///< does nothing
static const mode_t S_IXGRP      = 0x00080000;           ///< does nothing
static const mode_t S_IROTH      = 0x00040000;           ///< does nothing
static const mode_t S_IWOTH      = 0x00020000;           ///< does nothing
static const mode_t S_IXOTH      = 0x00010000;           ///< does nothing
#   endif
static const mode_t MS_MODE_MASK = 0x0000ffff;           ///< low word

static inline int my_chmod(const char * path, mode_t mode)
{
    int result = _chmod(path, (mode & MS_MODE_MASK));

    if (result != 0)
    {
        result = errno;
    }

    return (result);
}
#else
static inline int my_chmod(const char * path, mode_t mode)
{
    int result = chmod(path, mode);

    if (result != 0)
    {
        result = errno;
    }

    return (result);
}
#endif


string extractDirectory( const string& );
string extractFilename( const string& );
string changeExtension( const string& , const string& );
string removeExtension( const string& );
bool dirExists(const string& );
bool fileExists(const string& );
void movefile(const string& , const string&);
void copyfile(const string& , const string&);
string concurrentCopyTo (string &, string &);
string getTime();
bool isPositiveInt(char *);
bool isPositiveInt(const char *);
bool isNumber(char * );
bool isNumber(const char * input);
bool isNumeric( const char* pszInput, int nNumberBase=10 );
int FileRead( istream &, vector <char> &  );
uint CountLines( const vector <char> & , int );
template <typename M> void FreeClear( M & avec ) {
	for ( typename M::iterator it = avec.begin(); it != avec.end(); ++it ) {
		delete (*it);
	}
	avec.clear();
};

template <typename M> void FreeClearMap( M & amap ) {
	for ( typename M::iterator it = amap.begin(); it != amap.end(); ++it ) {
		delete it->second;
	}
	amap.clear();
};

template <typename Container> Container& split( Container &result, const typename Container::value_type& s, const typename Container::value_type& delimiters, bool empties = false ){
	result.clear();
	size_t current;
	size_t next = -1;
	do
	{
		if (!empties)
		{
			next = s.find_first_not_of( delimiters, next + 1 );
			if (next == Container::value_type::npos) break;
			next -= 1;
		}
		current = next + 1;
		next = s.find_first_of( delimiters, current );
		result.push_back( s.substr( current, next - current ) );
	}
	while (next != Container::value_type::npos);
	return result;
};

template <typename Map>
bool map_compare (Map const &lhs, Map const &rhs) {
    // No predicate needed because there is operator== for pairs already.
    return lhs.size() == rhs.size()
        && std::equal(lhs.begin(), lhs.end(),
                      rhs.begin());
}

template<class T>
string FormatWithCommas(T value)
{
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << value;
    return ss.str();
}

#endif
