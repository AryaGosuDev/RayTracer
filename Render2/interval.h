

#ifndef __INTERVAL_INCLUDED__
#define __INTERVAL_INCLUDED__

#include "base.h"

struct Interval { // An interval of real numbers.
    inline Interval() { min = Infinity; max = -Infinity; }
    inline Interval( double a ) { min = a; max = a; }
    inline Interval( double a, double b ) { min = a; max = b; }
    static inline Interval Null();
    double min;
    double max;
    };

inline double Len( const Interval &I )
    {
    // Return the length of interval I.  If I is degenerate, the length is zero.
    double len = I.max - I.min;
    return len > 0.0 ? len : 0.0;
    }

inline bool IsNull( const Interval &I )
    { 
    // Is the interval I degenerate?
    return I.max < I.min;
    }

inline Interval Interval::Null()
    {
    // Return an interval that is maximally degenerate.
    return Interval( Infinity, -Infinity );
    }

inline Interval &operator<<( Interval &A, const Interval &B ) 
    {
    // Expand interval A to include interval B.
    if( B.min < A.min ) A.min = B.min;
    if( B.max > A.max ) A.max = B.max;
    return A;
    };

inline Interval &operator<<( Interval &I, const double &x ) 
    {
    // Expand interval I to include the real value x.
    if( x < I.min ) I.min = x;
    if( x > I.max ) I.max = x;
    return I;
    };

inline Interval &operator+=( Interval &I, const double &x ) 
    {
    // Translate the interval I by x.
    I.min += x;
    I.max += x;
    return I;
    }

inline Interval &operator-=( Interval &I, const double &x )
    {
    // Translate the interval I by -x.
    I.min -= x;
    I.max -= x;
    return I;
    }

inline Interval operator+( double x, const Interval &I )
    {
    // Return a new interval that is I translated by x.
    return Interval( I.min + x, I.max + x );
    }

inline Interval operator+( const Interval &I, double x ) 
    {
    // Return a new interval that is I translated by x.
    return Interval( I.min + x, I.max + x );
    }

inline Interval operator-( const Interval &I, double x )
    {
    // Return a new interval that is I translated by -x.
    return Interval( I.min - x, I.max - x );
    }

inline Interval operator/( const Interval &I, double c )
    {
    // Scale all elements of the interval I by 1/c.
    double a = I.min / c;
    double b = I.max / c;
    return c >= 0.0 ? Interval( a, b ) : Interval( b, a );
    }

inline bool Inside( double x, const Interval &I )
    {
    // Is the point x within the interval I?
    return ( I.min <= x ) && ( x <= I.max );
    }

inline ostream &operator<<( ostream &out, const Interval &I )
    {
    // An output operator, useful for debugging.
    out << "[" << I.min << "," << I.max << "]";
    return out;
    }

#endif


