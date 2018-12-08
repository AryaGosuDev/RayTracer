
// misc. utilities such as predicates on materials & objects.      *

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

