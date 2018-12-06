
/*
  A sphere of given radius and center.                                     
         
//From Book
To intersect a ray with a sphere with the given center and radius, we    
solve the following equation for s: || Q + sR - C || = radius, where Q   
is the ray origin, R is the ray direction, and C is the center of the    
sphere.  This is equivalent to  ( A + sR )^T ( A + sR ) = radius^2,      
where A = Q - C.  Expanding, A.A + 2 s A.R + s^2 R.R = radius^2,         
where "." denotes the dot product.  Since R is a unit vercor, R.R = 1.   
Rearranging, we have s^2 + (2 A.R) s + (A.A - radius^2) = 0, which is    
a quadratic equation in s, the distance along the ray to the point of    
intersection.  If this equation has complex roots, then the ray misses   
the sphere.  Otherwise, we must determine whether either of the roots    
falls on the positive part of the ray, and if so, which is closer.       
*/

#include "ray_tracer.h"
#include "util.h"
#include "params.h"

// Register the new object with the toytracer.  When this module is linked in, 
// the toytracer will automatically recognize the sphere object.

REGISTER_PLUGIN( SphereLight );

SphereLight::SphereLight( const Vec3 &cent, double rad )
	{
	center  = cent;
	radius  = rad;
	radius2 = rad * rad;
	}

PrimitiveObject *SphereLight::ReadString( const string &params ) // Reads params from a string.
	{
	Vec3 cent;
	double r;
	ParamReader get( params );
	if( get[MyName()] && get[cent] && get[r] ){
		return new SphereLight( cent, r );
	}
	return NULL;
}

Interval SphereLight::GetSlab( const Vec3 &v ) const
	{
	const double len = Length(v);
	const double dot = ( v * center ) / len;
	return Interval( dot - radius, dot + radius ) / len;
	}

bool SphereLight::Intersect( const Ray &ray, HitInfo &hitinfo ) const
	{
	const Vec3 A( ray.origin - center ); // (5, -2, 2.79 ) - ( 0, 0, 0 )
	const Vec3 R( ray.direction ); // based upon screen resolution, etc.... ray.direction = Unit( O + (j + 0.5) * dR - (i + 0.5) * dU  ); //du and dr are gradations based on the resolution, O = origin
	const double b = 2.0 * ( A * R ); 
	const double discr = b * b - 4.0 * ( A * A - radius2 );  // The discriminant.

	// If the discriminant is negative, the quadratic equation has negative
	// roots, and the ray misses the sphere.

	if( discr < 0.0 ) return false;

	const double radical = sqrt( discr );

	// First try the smaller of the two roots.  If this is positive, it is the
	// closest intersection.

	double s = 0.5 * ( -b - radical );
	if( s > 0.0 )
	{
		// If the closest intersection is too far away, report a miss.
		if( s > hitinfo.distance ) return false;
	}
	else
	{
		// Now try the other root, since the smallest root puts the
		// point of intersection "behind"us.  If the larger root is now
		// positive, it means we are inside the sphere.
		s = 0.5 * ( -b + radical );
		if( s <= 0 ) return false;
		if( s > hitinfo.distance ) return false;
	}

	// We have an actual hit.  Fill in all the geometric information so
	// that the shader can shade this point.

	hitinfo.distance = s;
	hitinfo.point    = ray.origin + s * R;
	hitinfo.normal   = Unit( hitinfo.point - center );
	//hitinfo.object   = this;
	return true;
	}

int SphereLight::GetSamples( const Vec3 &P, const Vec3 &N, Sample *samples, int n ) const
	{
	int count = 0;
	Vec3   R( center - P );
	Vec3   Z = Unit( R );
	Vec3   U = Unit( OrthogonalTo( Z ) );
	Vec3   V = Unit( R ^ U );
	double d = Length( R );
	double a = sqrt( d * d - radius * radius );
	double sa = TwoPi * ( d - a ) / d;
	double w  = sa / ( n * n );
	for( int i = 0; i < n; i++ )
	for( int j = 0; j < n; j++ )
		{
		double x1 = ( i + rand(0,1) ) / n;
		double x2 = ( j + rand(0,1) ) / n;
		double z  = x1 + ( 1.0 - x1 ) * a / d;
		double r  = sqrt( 1.0 - z * z );
		double s  = r * sin( x2 * TwoPi );
		double c  = r * cos( x2 * TwoPi );
		samples[ count ].P = P + d * ( z * Z + s * U + c * V );
		samples[ count ].w = w;
		count++;
		}
	return n * n;
	}

