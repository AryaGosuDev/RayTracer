
// misc. utilities such as predicates on materials & objects.      *
#include <random>
#include "ray_tracer.h"
#include "util.h"


//static std::list< PrimitiveObject * > * all_light_types = NULL;

double min( double x, double y, double z )
	{
	if( x <= y && x <= z ) return x;
	if( y <= z ) return y;
	return z;
	}

double max( double x, double y, double z )
	{
	if( x >= y && x >= z ) return x;
	if( y >= z ) return y;
	return z;
	}

AABB GetBox( const PrimitiveObject &obj ){
	AABB box;
	box.X = obj.GetSlab( Vec3( 1, 0, 0 ) );
	box.Y = obj.GetSlab( Vec3( 0, 1, 0 ) );
	box.Z = obj.GetSlab( Vec3( 0, 0, 1 ) );
	return box;
	}

// return an axis aligned bounding box for a polygon
AABB GetBoxPolygon ( const Object * _object ) {

	vector<Vec3>::const_iterator pointsIter = _object->points.cbegin();
	double HiX , LowX = 0.0;
	double HiY , LowY = 0.0;
	double HiZ , LowZ = 0.0;
	bool isFirst = true  ;

	for ( ; pointsIter != _object->points.cend(); ++pointsIter ) {

		if ( isFirst ) { 
			HiX = LowX = pointsIter->x ;
			HiY = LowY = pointsIter->y ;
			HiZ = LowZ = pointsIter->z ;
			isFirst = !isFirst ;

		}
		else {
			if ( pointsIter->x > HiX )  HiX = pointsIter->x ;
			else if ( pointsIter->x < LowX )  LowX = pointsIter->x ;

			if ( pointsIter->y > HiY )  HiY = pointsIter->y ;
			else if ( pointsIter->y < LowY )  LowY = pointsIter->y ;

			if ( pointsIter->z > HiZ )  HiZ = pointsIter->z ;
			else if ( pointsIter->z < LowZ )  LowZ = pointsIter->z ;
		}
	}

	AABB tempBox ( Interval ( LowX, HiX ), Interval ( LowY, HiY ), Interval ( LowZ, HiZ ));

	return tempBox ;
}

bool AABBIntersect ( const AABB & _a, const AABB & _b ){

	double diff1 = ((_a.X.max - _a.X.min) / 2.0) + ((_b.X.max - _b.X.min ) / 2.0 ) ;
	if ( diff1 < abs( ( (_a.X.max + _a.X.min) / 2.0 ) - ( (_b.X.max + _b.X.min ) / 2.0 ) ) ) return false;
	double diff2 = ((_a.Y.max - _a.Y.min) / 2.0) + ((_b.Y.max - _b.Y.min ) / 2.0 ) ;
	if ( diff2 < abs( ( (_a.Y.max + _a.Y.min) / 2.0 ) - ( (_b.Y.max + _b.Y.min ) / 2.0 ) ) ) return false;
	double diff3 = ((_a.Z.max - _a.Z.min) / 2.0) + ((_b.Z.max - _b.Z.min ) / 2.0 ) ;
	if ( diff3 < abs( ( (_a.Z.max + _a.Z.min) / 2.0 ) - ( (_b.Z.max + _b.Z.min ) / 2.0 ) ) ) return false;
	return true;
}

double rand( double a, double b ){
	double x = float(rand()) / RAND_MAX;
	if( x < 0.0 ) x = -x;
	return a + x * ( b - a );
	}

bool operator==( const Material &a, const Material &b )
	{
	return
		a.diffuse      == b.diffuse      &&
		a.specular     == b.specular     &&
		a.emission     == b.emission     &&
		a.ambient      == b.ambient      &&
		a.reflectivity == b.reflectivity &&
		a.translucency == b.translucency &&
		a.Phong_exp    == b.Phong_exp    &&
		a.ref_index    == b.ref_index    &&
		a.type         == b.type; 
	}

int returnHighestValueIndx ( double * _a, int & n ) {
	double val = 0.0 ;
	int indx = 0 ;
	//int sizeOfArray = sizeof ( _a ) / sizeof ( double ) ;

	for ( int i = 0 ; i < n ; ++ i ) {
		if ( val < _a[i] ) {
			val = _a[i] ;
			indx = i ;
		}
	}

	return indx ;
}

inline Vec3 returnRandomPointOnElement( QuadTreeNode * i ) {
	
	try {
		std::random_device rd;
		std::mt19937 gen(rd());

		double u = 0.0, v = 0.0 ;

		do {
			std::uniform_real_distribution<> dis (0.0, 1.0 ) ;

			u = dis ( gen ) ;
			v = dis ( gen ) ;

		} while ( u + v > 1.0 ) ;

		Vec3 P ;
		P.x = (( 1.0 - u - v ) * i->triVert1.x) + (u * i->triVert2.x) + (v * i->triVert3.x) ;
		P.y = (( 1.0 - u - v ) * i->triVert1.y) + (u * i->triVert2.y) + (v * i->triVert3.y) ;
		P.z = (( 1.0 - u - v ) * i->triVert1.z) + (u * i->triVert2.z) + (v * i->triVert3.z) ;

		return P ;
	}
	catch ( std::exception ex ) {
		cout <<"Error in inline Vec3 returnRandomPointOnElement : " << ex.what()  << endl ;
	}

	return Vec3 ( 0.0, 0.0, 0.0 );
}

inline std::pair<double, double> returnUVofTriangle( const Vec3 & _p , QuadTreeNode * _Qnode ) {

	Mat3x2 A ;
	A(0,0) = ( _Qnode->triVert2 - _Qnode->triVert1 ).x ;
	A(0,1) = ( _Qnode->triVert3 - _Qnode->triVert1 ).x ;
	A(1,0) = ( _Qnode->triVert2 - _Qnode->triVert1 ).y ;
	A(1,1) = ( _Qnode->triVert3 - _Qnode->triVert1 ).y ;
	A(2,0) = ( _Qnode->triVert2 - _Qnode->triVert1 ).z ;
	A(2,1) = ( _Qnode->triVert3 - _Qnode->triVert1 ).z ;

	Vec3 _A, _B  ;

	if ( calculatePseudoInverse ( A, _A, _B ) ) {

		double u = _A.x * (_p - _Qnode->triVert1 ).x + _A.y * (_p - _Qnode->triVert1 ).y + _A.z * (_p - _Qnode->triVert1 ).z  ;
		double v = _B.x * (_p - _Qnode->triVert1 ).x + _B.y * (_p - _Qnode->triVert1 ).y + _B.z * (_p - _Qnode->triVert1 ).z  ;

		return std::pair<double, double> ( u, v ) ;
	}

	return std::pair<double, double> ( 0.0, 0.0 ) ;
}