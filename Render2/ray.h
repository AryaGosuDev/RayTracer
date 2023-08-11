

#ifndef __RAY_INCLUDED__
#define __RAY_INCLUDED__

#define ERROR_SIDE_NUMBER -90     // Indicates an error in determining the side that one object lie with another object

#include "base.h"
#include "vec3.h"

// Possible ray "types" that may be set and used by a shader and/or an
// acceleration method.
enum ray_type {         // A flag that may affect the processing of a ray.
	undefined_ray = 0,  // Nothing known about this ray.
	primary_ray   = 1,  // A ray cast from the eye by the rasterizer.
	generic_ray   = 2,  // Any ray; No special meaning.
	shadow_ray    = 3,  // A ray cast for shadow detection only.
	indirect_ray  = 4,  // A ray cast from a surface to sample illumination.
	light_ray     = 5,  // A ray cast from a light source (e.g. for photon mapping).
	special_ray   = 6   // A ray used for some other special purpose.
};

struct Ray {            // A ray in R3.
	inline Ray();
	inline Ray( const Ray &r );
	virtual ~Ray() {}
	Vec3 origin;         // The ray originates from this point.
	Vec3 direction;      // Unit vector indicating direction of ray.
	ray_type  type;      // Different rays may be processed differently.
	unsigned generation; // How deep in the ray tree.  1 == generated from eye.
	const Object *from;  // The object from which the ray was cast.
};

inline Ray::Ray()
{
	generation = 1;
	type = generic_ray;
	from = NULL;
}

inline Ray::Ray( const Ray &r )
{
	origin     = r.origin;
	direction  = r.direction;
	generation = r.generation;
	type       = r.type;
	from       = r.from;
}

inline int sideTest3d (  Vec3 a,  Vec3 b ,  Vec3 c,  Vec3 x )
{
		double result = 
		((b-a).x * ((c-a).y * (x-a).z - (c-a).z * (x-a).y)) - 
		((c-a).x * ((b-a).y * (x-a).z - (b-a).z * (x-a).y)) +
		((x-a).x * ((b-a).y * (c-a).z - (b-a).z * (c-a).y));
		
		if ( result < -Epsilon3)
			return -1;
		else if (result >= -Epsilon3 && result <= Epsilon3)
			return 0;
		else return 1;
		
		return ERROR_SIDE_NUMBER;
}

inline int sideOfLine3D (Vec3 & _Point, Vec3 & _LineNormal ) {

	double result = _LineNormal * _Point ;

	if ( result == 0.0 ) return 0;
	else if ( result > 0.0 ) return 1;
	else return -1;
}

/*
 //TODO
inline int sideOfLine3D ( Vec3 & _QLine, Vec3 & _VLine, Vec3 & _Point ) {
	return 0;
}
*/

// Compute the reflected ray given the incident ray (i.e. directed
// toward the surface), and the normal to the surface.  The normal
// may be directed away from or into the surface.  Both the surface
// normal and the ray direction vector are assumed to be normalized.
inline Vec3 Reflect( const Ray &r, const Vec3 &N )
{
	Vec3 U( r.direction );
	return U - ( 2.0 * ( U * N ) ) * N;
}

// Compute the refracted ray given the incident ray (i.e. directed
// toward the surface), the normal to the surface, and the ratio of the
// refractive indices of the material containing the ray and the
// material into which the ray is refracted.  The normal
// may be directed away from or into the surface.  Both the surface
// normal and the ray direction vector are assumed to be normalized.
inline Vec3 Refract( const Ray &r, const Vec3 &N, double eta1_over_eta2 ) {
	Vec3 U( r.direction );
	return U - ( 2.0 * ( U * N ) ) * N;
	}

// An output method, useful for debugging.
inline ostream &operator<<( ostream &out, const Ray &r )
	{
	out << "[ray:"
		<< " org "  << r.origin
		<< " dir "  << r.direction
		<< " type " << r.type
		<< " gen "  << r.generation
		<< "]";
	return out;
	}

#endif

