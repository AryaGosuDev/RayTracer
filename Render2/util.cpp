/***************************************************************************
* util.cpp                                                                 *
*                                                                          *
* Miscellaneous utilities, such as predicates on materials & objects.      *
*                                                                          *
* History:                                                                 *
*   10/16/2005  Added ToString function for plugin_type.                   *
*   12/11/2004  Initial coding.                                            *
*                                                                          *
***************************************************************************/
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

