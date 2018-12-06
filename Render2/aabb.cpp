                                                                      
// Encodes a three-dimensional axis-aligned box as a collection of three intervals, one for each axis.  
// These boxes are ideal as bounding boxes as they are equipped with a very efficient ray-box intersector.  
// Methods and functions that expand the box, transform, surface area, etc...             

#include "aabb.h"
#include "util.h"

// This operator allows a 3x3 matrix to transform a bounding box.
// The result is another axis-aligned box guaranteed to contain all
// of the the original vertices transformed by the given matrix.

AABB operator*( const Mat3x3 &M, const AABB &box )
    {
    double new_min[] = { 0, 0, 0 };
    double new_max[] = { 0, 0, 0 };
    double old_min[] = { box.X.min, box.Y.min, box.Z.min };
    double old_max[] = { box.X.max, box.Y.max, box.Z.max }; 

    // Find the extreme points by considering the product
    // of the min and max with each component of M.
                     
    for( int i = 0; i < 3; i++ )
    for( int j = 0; j < 3; j++ )
        {
        double a = M(i,j) * old_min[j];
        double b = M(i,j) * old_max[j];
        new_min[i] += min( a, b ); 
        new_max[i] += max( a, b );
        }

    // Return the smallest bounding box enclosing the extremal points.
    AABB b;
    b.X = Interval( new_min[0], new_max[0] );
    b.Y = Interval( new_min[1], new_max[1] );
    b.Z = Interval( new_min[2], new_max[2] );
    return b;
    }

// Determine whether the given ray intersects the box.  If there is an intersection,
// but it's farther than max_dist, then regard it as a miss.

bool Hit( const Ray &ray, const AABB &A, double max_dist )
    {
    Vec3 R = ray.direction;
    Vec3 Q = ray.origin;
    double r, s, t;
    double min = 0.0;
    double max = max_dist;

    if( R.x > 0 ) // Looking down positive X-axis.
        {
        if( Q.x > A.X.max ) return false;
        r = 1.0 / R.x;
        s = ( A.X.min - Q.x ) * r; if( s > min ) min = s;
        t = ( A.X.max - Q.x ) * r; if( t < max ) max = t;
        }
    else if( R.x < 0.0 ) // Looking down negative X-axis.
        {
        if( Q.x < A.X.min ) return false;
        r = 1.0 / R.x;
        s = ( A.X.max - Q.x ) * r; if( s > min ) min = s;
        t = ( A.X.min - Q.x ) * r; if( t < max ) max = t;
        }
    else if( Q.x < A.X.min || Q.x > A.X.max ) return false;

    if( min > max ) return false; // Degenerate interval.

    if( R.y > 0 ) // Looking down positive Y-axis.
        {
        if( Q.y > A.Y.max ) return false;
        r = 1.0 / R.y;
        s = ( A.Y.min - Q.y ) * r; if( s > min ) min = s;
        t = ( A.Y.max - Q.y ) * r; if( t < max ) max = t;
        }
    else if( R.y < 0 ) // Looking down negative Y-axis.
        {
        if( Q.y < A.Y.min ) return false;
        r = 1.0 / R.y;
        s = ( A.Y.max - Q.y ) * r; if( s > min ) min = s;
        t = ( A.Y.min - Q.y ) * r; if( t < max ) max = t;
        }
    else if( Q.y < A.Y.min || Q.y > A.Y.max ) return false;

    if( min > max ) return false; // Degenerate interval.

    if( R.z > 0 ) // Looking down positive Z-axis.
        {
        if( Q.z > A.Z.max ) return false;
        r = 1.0 / R.z;
        s = ( A.Z.min - Q.z ) * r; if( s > min ) min = s;
        t = ( A.Z.max - Q.z ) * r; if( t < max ) max = t;
        }
    else if( R.z < 0 ) // Looking down negative Z-axis.
        {
        if( Q.z < A.Z.min ) return false;
        r = 1.0 / R.z;
        s = ( A.Z.max - Q.z ) * r; if( s > min ) min = s;
        t = ( A.Z.min - Q.z ) * r; if( t < max ) max = t;
        }
    else if( Q.z < A.Z.min || Q.z > A.Z.max ) return false;

    // There is a hit if and only if the intersection interval [min,max]
    // is not degenerate.
    return min <= max;
    }

bool Inside( const Vec3 &p, const AABB &A )
    {
    return Inside( p.x, A.X ) && Inside( p.y, A.Y ) && Inside( p.z, A.Z );
    }

double SurfaceArea( const AABB &A )
    {
    if( IsNull( A ) ) return 0.0;
    return 2.0 * ( Len(A.X) * ( Len(A.Y) + Len(A.Z) ) + Len(A.Y) * Len(A.Z) );
    }

double Volume( const AABB &A )
    {
    if( IsNull( A ) ) return 0.0;
    return Len(A.X) * Len(A.Y) * Len(A.Z);
    }

bool AABB::Expand( const AABB &b )
    {
    // Expand this box to enclose box "b" if it does not already do so.
    // Return true iff the box actually had to be expanded.
    bool expand = false;
    if( b.X.min < X.min ) { expand = true; X.min = b.X.min; } 
    if( b.X.max > X.max ) { expand = true; X.max = b.X.max; } 
    if( b.Y.min < Y.min ) { expand = true; Y.min = b.Y.min; } 
    if( b.Y.max > Y.max ) { expand = true; Y.max = b.Y.max; } 
    if( b.Z.min < Z.min ) { expand = true; Z.min = b.Z.min; } 
    if( b.Z.max > Z.max ) { expand = true; Z.max = b.Z.max; } 
    return expand;
    }

AABB operator*( const Mat3x4 &M, const AABB &box ) 
    {
    // Applies an affine transformation (encoded as a 3x4 matrix) to a
    // 3D box and returns the tightest-fitting axis-aligned enclosing
    // the result.
    AABB new_box = M.mat * box;
    return new_box += M.vec;
    }