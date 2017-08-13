
#ifndef __AABB_INCLUDED__
#define __AABB_INCLUDED__

#include "base.h"
#include "ray.h"
#include "interval.h"
#include "mat3x4.h"

struct AABB {  // A box in R3.  Useful for bounding boxes.
    AABB() {}
    AABB( const Interval &A, const Interval &B, const Interval &C ) { X = A; Y = B; Z = C; }
    inline Vec3 MinCorner() const { return Vec3( X.min, Y.min, Z.min ); }
    inline Vec3 MaxCorner() const { return Vec3( X.max, Y.max, Z.max ); }
    bool   Expand( const AABB & );
    static inline AABB Null();
    Interval X;
    Interval Y;
    Interval Z;
    };

inline AABB AABB::Null()
    {
    return AABB( Interval::Null(), Interval::Null(), Interval::Null() );
    }

inline bool IsNull( const AABB &A )
    { 
    return IsNull( A.X ) || IsNull( A.Y ) || IsNull( A.Z );
    }

inline Vec3 Center( const AABB &A )
    {
    return 0.5 * ( A.MinCorner() + A.MaxCorner() );
    }

inline AABB operator<<( AABB &A, const AABB &B ) 
    {
    // Expand box A to enclose box B.
    A.X << B.X;
    A.Y << B.Y;
    A.Z << B.Z;
    return A;
    };

inline AABB operator<<( AABB &A, const Vec3 &p ) 
    {
    // Expand box A to enclose the point p.
    A.X << p.x;
    A.Y << p.y;
    A.Z << p.z;
    return A;
    }

inline AABB operator+=( AABB &A, const Vec3 &p ) 
    {
    // Translate the box A by p.
    A.X += p.x;
    A.Y += p.y;
    A.Z += p.z;
    return A;
    }

inline AABB operator-=( AABB &A, const Vec3 &p ) 
    {
    // Translate the box A by p.
    A.X -= p.x;
    A.Y -= p.y;
    A.Z -= p.z;
    return A;
    }

inline ostream &operator<<( ostream &out, const AABB &box )
    {
    out << "(" << box.X << "," << box.Y << "," << box.X << ")";
    return out;
    }

extern bool Hit(   // Does the ray r intersect the box A? 
    const Ray &r,
    const AABB &A,
    double max_dist = Infinity
    );

extern bool Inside(  // Is the point P inside (or on) the box A?
    const Vec3 &P,
    const AABB &A
    );

extern double Volume(
    const AABB &
    );

extern double SurfaceArea(
    const AABB &
    );

extern AABB operator*(  // Transforms a bounding box by a 3x3 matrix.
    const Mat3x3 &M,
    const AABB &box
    ); 

extern AABB operator*(  // Transforms a bounding box by a 3x4 matrix.
    const Mat3x4 &M,
    const AABB &box
    );

#endif


