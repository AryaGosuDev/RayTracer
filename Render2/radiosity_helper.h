#ifndef __RADIOSITY_HELPER_INCLUDED__  // Include this file only once.
#define __RADIOSITY_HELPER_INCLUDED__

#include "base.h"         // Includes system headers and defines constants.
#include "vec3.h"         // Defines the Vec3 class, which are points in R3.
#include "vec2.h"         // Defines the Vec2 class, which are points in R2.
#include "mat3x3.h"       // Defines 3x3 matrices.
#include "mat3x4.h"       // Defines 3x4 matrices; i.e. affine transforms.
#include "color.h"        // Defines the Color class; real RGB values.
#include "ray.h"          // Defines rays in 3-space: origin, direction, etc.
#include "interval.h"     // Defines a (min,max) interval of the real line.
#include "params.h"

struct Radiosity_Helper : public Radiosity { 

	Radiosity_Helper() ;
	Radiosity_Helper(Scene *) ;
	virtual ~Radiosity_Helper() ;

	void detectTriangleIntersections ( Object * , Object *) ;

	Scene * scene ;

};

#endif