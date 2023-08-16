
#ifndef __PARAMS_INCLUDED__
#define __PARAMS_INCLUDED__

#include "ray_tracer.h"

// The ParamReader class makes the ReadString method of each Plugin very
// easy to write.  Given the parameters as a string, each "[]" operator
// attempts to strip off the given item, and returns "true" if successful.

class ParamReader {
	public:
		ParamReader( string s ) { params = s; }
		bool operator[]( string field );
		bool operator[]( Interval &I );
		bool operator[]( Color    & );
		bool operator[]( Vec3     & );
		bool operator[]( Vec2     & );
		bool operator[]( double   & );
		bool operator[]( unsigned &x );
		bool getString ( string &field );
		bool returnCoordOfFace ( Object & object, int & numOfTriangle, vector<vector<Vec3>> & );
		bool getPoint ( Vec3 & v );
		bool getTexture ( Vec2 & v );
		string returnString () { return params ; } 
	private:
		void SkipBlanks();
		string params;
};

inline void addUniqueNormals ( vector<vector<Vec3>> & adjMatrix, int & index, Vec3 & _normal ){

	for ( auto &  normal : adjMatrix[index] ) {
		if ( Length(normal - _normal ) < Epsilon ) return ;
	}

	adjMatrix[index].push_back ( _normal );
}

#endif

