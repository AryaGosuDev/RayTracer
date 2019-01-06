
#ifndef __UTIL_INCLUDED__
#define __UTIL_INCLUDED__

#include <list>
#include "ray_tracer.h"
#include "aabb.h"


extern bool RegisterPlugin(
	PrimitiveObject *pObj
	);

#define REGISTER_PLUGIN( class_name ) static bool dummy_variable_##class_name = RegisterPlugin( new class_name );

inline double sqr( double x )
	{
	return x * x;
	}

inline double min( double x, double y )
	{
	return x <= y ? x : y;
	}

inline double max( double x, double y )
	{
	return x >= y ? x : y;
	}

inline bool Emitter( const Object *obj )
	{ 
	//return ( obj != NULL ) & ( obj->material->emission != 0.0 );
		return obj->material->emission != 0.0 ;
	}

inline bool Emitter( const Material &mat )
	{ 
	return mat.emission != 0.0;
	}

inline bool Reflective( const Object *obj )
	{ 
   // return ( obj != NULL ) & ( obj->material->reflectivity != 0.0 );
	}

inline bool Translucent( const Object *obj )
	{ 
	//return ( obj != NULL ) & ( obj->material->translucency != 0.0 );
	}

extern AABB GetBox(  // Construct a box from three slabs.
	const PrimitiveObject &obj
	);

extern AABB GetBoxPolygon ( // Construct a bounding box of a polygon
        const Object * 
 );

extern int returnHighestValueIndx ( double *, int & ) ;

extern inline Vec3 returnRandomPointOnElement( QuadTreeNode *  ) ;

extern bool operator==(
	const Material &a,
	const Material &b
	);

extern double min(
	double x,
	double y,
	double z
	);

extern double max(
	double x,
	double y,
	double z
	);

extern double rand(    // Return a random number uniformly distributed in [a,b].
	double a,
	double b
	);

extern inline std::pair<double, double> returnUVofTriangle(const Vec3 &, QuadTreeNode * _Qnode ) ;


#endif

