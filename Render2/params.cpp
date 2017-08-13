
#include "params.h"
#include <algorithm>

using namespace std;

void ParamReader::SkipBlanks()
	{
	while( params[0] == ' ' || params[0] == '\t' ) params.erase(0,1);
	}

bool ParamReader::operator[]( string field ){
	SkipBlanks();
	if( params.find( field, 0 ) == 0 )
		{
		params.erase( 0, field.length() );
		return true;
		}
	return false;
	}

bool ParamReader::getString (string & field )
{
	SkipBlanks();
	char c[200];
	if ( sscanf ( params.c_str(), "%s", c ) == 1 )
	{
		field = c;
		params.erase( 0, field.length() ) ;
		return true;
	}
	return false;
}

bool ParamReader::getPoint ( Vec3 & v )
{
	SkipBlanks();
	if( sscanf( params.c_str(), "%lf %lf %lf", &v.x, &v.y, &v.z ) == 3 )
	{
		//int len = 1 + params.find( ")", 0 );
		//params.erase( 0, len );
		return true;
	}
	return false;
}

bool ParamReader::getTexture ( Vec2 & v )
{
	SkipBlanks();
	if( sscanf( params.c_str(), "%lf %lf", &v.x, &v.y ) == 2 )
	{
		//int len = 1 + params.find( ")", 0 );
		//params.erase( 0, len );
		return true;
	}
	return false;
}


bool ParamReader::operator[]( Interval &I )
{
	if( sscanf( params.c_str(), " ( %lf , %lf )", &I.min, &I.max ) == 2 )
		{
		int len = 1 + params.find( ")", 0 );
		params.erase( 0, len );
		return true;
		}
	return false;
}

bool ParamReader::operator[]( Color &c )
	{
	if( sscanf( params.c_str(), " [ %lf , %lf , %lf ]", &c.red, &c.green, &c.blue ) == 3 )
		{
		int len = 1 + params.find( "]", 0 );
		params.erase( 0, len );
		return true;
		}
	return false;
	}

bool ParamReader::operator[]( Vec3 &v )
	{
	SkipBlanks();
	if( sscanf( params.c_str(), "(%lf , %lf , %lf)", &v.x, &v.y, &v.z ) == 3 )
		{
			int len = 1 + params.find( ")", 0 );
			params.erase( 0, len );
			
		return true;
		}
	return false;
	}

bool ParamReader::operator[]( Vec2 &v )
	{
	if( sscanf( params.c_str(), " ( %lf , %lf )", &v.x, &v.y ) == 2 )
		{
		int len = 1 + params.find( ")", 0 );
		params.erase( 0, len );
		return true;
		}
	return false;
	}

bool ParamReader::operator[]( double &x )
{
	SkipBlanks();
	if( sscanf( params.c_str(), "%lf", &x ) == 1 )
		{
		int len = 1 + params.find( " ", 0 );
		if( len > 1 ) params.erase( 0, len );
		return true;
		}
	return false;
}

bool ParamReader::operator[]( unsigned &x )
{
	SkipBlanks();
	if( sscanf( params.c_str(), "%u", &x ) == 1 )
		{
		int len = 1 + params.find( " ", 0 );
		if( len > 1 ) params.erase( 0, len );
		return true;
		}
	return false;
}

// Record the triangle in the object. 
// Also if the face has an order > 3 then triangulate.
bool ParamReader::returnCoordOfFace ( Object & object, int & numOfTriangle, vector<vector<Vec3>> & adjMatrix ) 
{
	SkipBlanks();
	
	char * pch;

	char * writable = new char[params.size() + 1]; // Params is from the Parameter class constructor, stores exactly the line read in the file.
	std::copy(params.begin(), params.end(), writable);
	writable[params.size()] = '\0'; // don't forget the terminating 0

	vector<char * > vecOut ;
	vector < vector < int > > faceData ; //2d vector storing information per face
	
	pch = strtok (writable," ");  // Split the face string
	while (pch != NULL){
		vecOut.push_back( pch );
		pch = strtok (NULL, " ,.-");
	}

	int numOfFieldsPerFace = 0;
	int x = 0;
	bool hasTextureVert = true;
	bool hasVertexNormal = false ;
	//TODO : each face can have n order number of vertices, not a fixed amount
	while ( vecOut[0][x] != '\0') {
		if ( vecOut[0][x] == '/') {
			if ( vecOut[0][x+1] == '/' ) { // Either v or v//vn or v/vt/vn, find out if there is a '/' fight after a '/'
				numOfFieldsPerFace --;
				hasTextureVert = false;
				hasVertexNormal = true;
			}
			x++;
			numOfFieldsPerFace++;
		}
		else {
			numOfFieldsPerFace++;
		}
		x++;
	}

	//store each face data and it's elements in a 2d array
	for ( int i = 0 ; i < vecOut.size() ; ++ i ) {
		faceData.push_back(vector<int>());
		for ( int j = 0 ; j <= numOfFieldsPerFace ; ++ j ) {	
			if ( j == 0 ) pch = strtok ( vecOut [i] , "/" );
			else pch = strtok ( NULL , "/" );
			if ( pch != NULL ) faceData.at(i).push_back (atoi(pch) - 1); // Store the index of the vertex of the face
		}
	}

	object.triangles.push_back ( BSP_Node(numOfTriangle) );
	object.triangles.back().triVert1 = object.points[faceData[0][0]]; //first point
	object.triangles.back().triVert2 = object.points[faceData[1][0]]; //second point
	object.triangles.back().triVert3 = object.points[faceData[2][0]]; //third point
	object.triangles.back().triIndx1 = faceData[0][0];
	object.triangles.back().triIndx2 = faceData[1][0];
	object.triangles.back().triIndx3 = faceData[2][0];

	//CALC NORMAL OF TRIANGLE, MAKE SURE IT IS IN THE SAME DIR AS THE NORMAL OF THE MIDDLE POINT, (B - A) x (C - A)
	Vec3 calcNormal = Unit((object.triangles.back().triVert2 - object.triangles.back().triVert1) ^
					  (object.triangles.back().triVert3 - object.triangles.back().triVert2));

	if ( hasVertexNormal ) {
		object.triangles.back().triNormal = object.vertexNormals[faceData[2][numOfFieldsPerFace - 1] ]; //get vertex normal of second point
		if ( (calcNormal * (object.triangles.back().triNormal)) < 0.0 )
			calcNormal = -1.0 * calcNormal ;
	}

	object.triangles.back().triNormal = calcNormal;

	//adjMatrix[object.triangles.back().triIndx1].push_back(calcNormal);
	//adjMatrix[object.triangles.back().triIndx2].push_back(calcNormal);
	//adjMatrix[object.triangles.back().triIndx3].push_back(calcNormal);

	addUniqueNormals ( adjMatrix, object.triangles.back().triIndx1, calcNormal );
	addUniqueNormals ( adjMatrix, object.triangles.back().triIndx2, calcNormal );
	addUniqueNormals ( adjMatrix, object.triangles.back().triIndx3, calcNormal );

	// TODO Vertex textures and interpolation
	// TODO Trapezoidal decomp
	if ( vecOut.size() == 4 ) //the face is a quad, or of order n. Triangulate
	{
		object.triangles.push_back ( BSP_Node(++numOfTriangle) ); //make a new triangle for the fourth point
		object.triangles.back().triVert1 = object.points[faceData[2][0]] ; //vert 3 of prior triangle
		object.triangles.back().triVert2 = object.points[faceData[3][0]];  //vert 4 
		object.triangles.back().triVert3 = object.points[faceData[0][0]] ; //vert 1 of prior triangle
		object.triangles.back().triNormal = object.vertexNormals[faceData[3][numOfFieldsPerFace-1]]; //get vertex normal of fourth point

		//CALC NORMAL OF TRIANGLE, MAKE SURE IT IS IN THE SAME DIR AS THE NORMAL OF THE MIDDLE POINT
		Vec3 calcNormal = Unit((object.triangles.back().triVert1 - object.triangles.back().triVert2) ^
						  (object.triangles.back().triVert3 - object.triangles.back().triVert2));

		if ( (calcNormal * (object.triangles.back().triNormal)) < 0.0 )
			calcNormal = -1.0 * calcNormal ;
		
		object.triangles.back().triNormal = calcNormal;
	}

	delete[] writable;

	return true; 
}
