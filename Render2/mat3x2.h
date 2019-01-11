#ifndef __MAT3X2_INCLUDED__
#define __MAT3X2_INCLUDED__

#include "vec3.h"
#include "vec2.h"

struct Mat3x2 {
    inline Mat3x2();
    inline Mat3x2( const Mat3x2 &M ) { *this = M; }
    inline ~Mat3x2() {}
    inline       double &operator()( int i, int j )       { return m[i][j]; }
    inline const double &operator()( int i, int j ) const { return m[i][j]; }
    inline static Mat3x2 Identity();
    inline Vec2 row( int i ) const { return Vec2( m[i][0], m[i][1] ); }
    inline Vec2 col( int j ) const { return Vec2( m[0][j], m[1][j] ); }
    double m[3][2];

    };

inline Mat3x2::Mat3x2()
    {
    for( int i = 0; i < 3; i++ )
    for( int j = 0; j < 2; j++ ) m[i][j] = 0.0;
    }

inline Mat3x2 operator+( const Mat3x2 &A, const Mat3x2 &B )
    {
    Mat3x2 C;
    for( int i = 0; i < 3; i++ )
    for( int j = 0; j < 2; j++ ) C(i,j) = A(i,j) + B(i,j);
    return C;
    }

inline Mat3x2 operator-( const Mat3x2 &A, const Mat3x2 &B )
    {
    Mat3x2 C;
    for( int i = 0; i < 3; i++ )
    for( int j = 0; j < 2; j++ ) C(i,j) = A(i,j) - B(i,j);
    return C;
    }

inline Mat3x2 operator*( double c, const Mat3x2 &M )
    {
    Mat3x2 A;
    for( int i = 0; i < 3; i++ )
    for( int j = 0; j < 2; j++ ) A(i,j) = c * M(i,j);
    return A;    
    }

inline Mat3x2 operator*( const Mat3x2 &M, double c )
    {
    return c * M;
    }

/*
inline Vec2 operator*( const Mat3x2 &M, const Vec2 &A )
    {
    // Transform the column vector A, multiplying on the left by matrix M.
    return Vec2(
        M(0,0) * A.x + M(0,1) * A.y + M(0,2) * A.z,
        M(1,0) * A.x + M(1,1) * A.y + M(1,2) * A.z,
        M(2,0) * A.x + M(2,1) * A.y + M(2,2) * A.z
        );
    }

inline Vec2 operator^( const Mat3x2 &M, const Vec2 &A )
    {
    // Multiply the vector A on the left by the TRANSPOSE of the matrix.
    return Vec2(
        M(0,0) * A.x + M(1,0) * A.y + M(2,0) * A.z,
        M(0,1) * A.x + M(1,1) * A.y + M(2,1) * A.z,
        M(0,2) * A.x + M(1,2) * A.y + M(2,2) * A.z
        );
    }
	*/
inline Mat3x2 operator*( const Mat3x2 &A, const Mat3x2 &B )
    {
    Mat3x2 C;
    for( int i = 0; i < 3; i++ )
    for( int j = 0; j < 2; j++ )
        C(i,j) = A(i,0) * B(0,j) + A(i,1) * B(1,j) + A(i,2) * B(2,j);
    return C;
    }

inline Mat3x2 operator/( const Mat3x2 &M, double c )
    {
    return (1/c) * M;
    }

inline Mat3x2 Adjoint( const Mat3x2 &M )  
    {
    // The Adjoint matrix is the inverse transpose times the determinant.
    Mat3x2 A;
    A(0,0) = M(1,1) * M(2,2) - M(1,2) * M(2,1);
    A(0,1) = M(1,2) * M(2,0) - M(1,0) * M(2,2);
    A(0,2) = M(1,0) * M(2,1) - M(1,1) * M(2,0);
 
    A(1,0) = M(0,2) * M(2,1) - M(0,1) * M(2,2);
    A(1,1) = M(0,0) * M(2,2) - M(0,2) * M(2,0);
    A(1,2) = M(0,1) * M(2,0) - M(0,0) * M(2,1);

    A(2,0) = M(0,1) * M(1,2) - M(0,2) * M(1,1);
    A(2,1) = M(0,2) * M(1,0) - M(0,0) * M(1,2);
    A(2,2) = M(0,0) * M(1,1) - M(0,1) * M(1,0);
    return A;
    }

inline double det( const Mat3x2 &M )  // Determinant.
    {
    return
        M(0,0) * ( M(1,1) * M(2,2) - M(1,2) * M(2,1) )
      - M(0,1) * ( M(1,0) * M(2,2) - M(1,2) * M(2,0) )
      + M(0,2) * ( M(1,0) * M(2,1) - M(1,1) * M(2,0) );
    }

inline Mat3x2 Transpose( const Mat3x2 &M )
    {
    Mat3x2 W;
    for( int i = 0; i < 3; i++ )
    for( int j = 0; j < 3; j++ ) W(i,j) = M(j,i);
    return W;
    }

inline Mat3x2 Inverse( const Mat3x2 &M )
    {
    double d = det( M );
    return Transpose( Adjoint( M ) ) / d;
    }

inline Mat3x2 Mat3x2::Identity()
    {
    Mat3x2 I;
    I(0,0) = 1.0;
    I(1,1) = 1.0;
    I(2,2) = 1.0;
    return I;
    }
/*
inline Mat3x2 Rotate_X( double angle )
    {
    // Rotation about the X-axis.
    Mat3x2 M = Mat3x2::Identity();
    M(1,1) = cos( angle );  M(1,2) = -sin( angle );
    M(2,1) = sin( angle );  M(2,2) =  cos( angle );
    return M;
    }

inline Mat3x2 Rotate_Y( double angle )
    {
    // Rotation about the Y-axis.
    Mat3x2 M = Mat3x2::Identity();
    M(0,0) = cos( angle );  M(0,2) = -sin( angle );
    M(2,0) = sin( angle );  M(2,2) =  cos( angle );
    return M;
    }

inline Mat3x2 Rotate_Z( double angle )
    {
    // Rotation about the Z-axis.
    Mat3x2 M = Mat3x2::Identity();
    M(0,0) = cos( angle );  M(0,1) = -sin( angle );
    M(1,0) = sin( angle );  M(1,1) =  cos( angle );
    return M;
    }

inline Mat3x2 Scale( double x, double y, double z )
    {
    // Allows non-uniform scaling along all three coordinate axes.
    Mat3x2 M;
    M(0,0) = x;
    M(1,1) = y;
    M(2,2) = z;
    return M;
    }
*/

inline ostream &operator<<( ostream &out, const Mat3x2 &M )
    {
    out << "\n";
    out << "| " << M(0,0) << " " << M(0,1) << " " << M(0,2) << " |\n";
    out << "| " << M(1,0) << " " << M(1,1) << " " << M(1,2) << " |\n";
    out << "| " << M(2,0) << " " << M(2,1) << " " << M(2,2) << " |\n";
    out << endl;
    return out;
    }

inline bool calculatePseudoInverse ( Mat3x2 & _M, Vec3 & _row1, Vec3 & _row2  ) {

	// A^T * A
	Vec2 AtACol_1( _M(0,0) * _M(0,0) + _M(1,0) * _M(1,0) + _M(2,0) * _M(2,0) ,
		           _M(0,1) * _M(0,0) + _M(1,1) * _M(1,0) + _M(2,1) * _M(2,0) );
	Vec2 AtACol_2( _M(0,0) * _M(0,1) + _M(1,0) * _M(1,1) + _M(2,0) * _M(2,1) ,
		           _M(0,1) * _M(0,1) + _M(1,1) * _M(1,1) + _M(2,1) * _M(2,1) );

	// coeff of A^T * A
	double x = AtACol_1.x ;
	double z = AtACol_1.y ;
	double y = AtACol_2.x ;
	double w = AtACol_2.y ;

	// calc roots of A-(lambda)(I)
	double eigenV_1 = ( x + w ) + sqrt ( pow(( x + w ), 2) - ( 4.0 * ( ( x * w ) - ( y * z ) ) ) ) ;
	eigenV_1 /= 2.0;
	double eigenV_2 = ( x + w ) - sqrt ( pow(( x + w ), 2) - ( 4.0 * ( ( x * w ) - ( y * z ) ) ) ) ;
	eigenV_2 /= 2.0;

	// not an invertable matrix
	if ( eigenV_1 == 0.0 || eigenV_2 == 0.0 ) return false;

	// if both roots are the same then these are already the SVD
	if ( eigenV_1 == eigenV_2  ) {

		_row1.x =  _M(0,0) / eigenV_1 ;
		_row1.y =  _M(1,0) / eigenV_1 ;
		_row1.z =  _M(2,0) / eigenV_1 ;

		_row2.x =  _M(0,1) / eigenV_1;
		_row2.y =  _M(1,1) / eigenV_1;
		_row2.z =  _M(2,1) / eigenV_1;

		return true ;
	}

	if (  z == 0.0 && y == 0.0 ) {

		_row1.x =  _M(0,0) / eigenV_1 ;
		_row1.y =  _M(1,0) / eigenV_1 ;
		_row1.z =  _M(2,0) / eigenV_1 ;

		if ( abs(( _row1.x * _M(0,0) + _row1.y * _M(1,0) + _row1.z * _M(2,0) ) - 1.0 ) > Epsilon ) {
			_row1.x =  _M(0,0) / eigenV_2 ;
			_row1.y =  _M(1,0) / eigenV_2 ;
			_row1.z =  _M(2,0) / eigenV_2 ;

			_row2.x =  _M(0,1) / eigenV_1;
			_row2.y =  _M(1,1) / eigenV_1;
			_row2.z =  _M(2,1) / eigenV_1;

			return true;
		}

		_row2.x =  _M(0,1) / eigenV_2;
		_row2.y =  _M(1,1) / eigenV_2;
		_row2.z =  _M(2,1) / eigenV_2;

		return true ;
	}

	// calc V
	double lambda_1_x2 = abs ( AtACol_1.x - eigenV_1) ;
	double lambda_1_x1 = ((-1.0 * AtACol_2.x) / ( AtACol_1.x - eigenV_1)) * lambda_1_x2 ;

	double lambda_2_x2 = abs ( AtACol_1.x - eigenV_2) ;
	double lambda_2_x1 = ((-1.0 * AtACol_2.x) / ( AtACol_1.x - eigenV_2)) * lambda_2_x2 ; 

	Vec2 V1 = Unit(Vec2 ( lambda_1_x1, lambda_1_x2 )) ;
	Vec2 V2 = Unit(Vec2 ( lambda_2_x1, lambda_2_x2 )) ;

	Vec3 U1 = Vec3 (  _M(0,0) * V1.x + _M(0,1) * V1.y , _M(1,0) * V1.x + _M(1,1) * V1.y , _M(2,0) * V1.x + _M(2,1) * V1.y ) / sqrt ( eigenV_1 ) ;  
	Vec3 U2 = Vec3 (  _M(0,0) * V2.x + _M(0,1) * V2.y , _M(1,0) * V2.x + _M(1,1) * V2.y , _M(2,0) * V2.x + _M(2,1) * V2.y ) / sqrt ( eigenV_2 ) ;

	U1 = Unit ( U1 );
	U2 = Unit ( U2 );

	Vec2 interMedPseudo1 ( V1.x * (1.0 / sqrt ( eigenV_1 )), V2.x * ( 1.0 / sqrt (eigenV_2 ) ) ) ;  
	Vec2 interMedPseudo2 ( V1.y * (1.0 / sqrt ( eigenV_1 )), V2.y * ( 1.0 / sqrt (eigenV_2 ) ) ) ; 

	_row1.x =  interMedPseudo1.x * U1.x + interMedPseudo1.y * U2.x;
	_row1.y =  interMedPseudo1.x * U1.y + interMedPseudo1.y * U2.y;
	_row1.z =  interMedPseudo1.x * U1.z + interMedPseudo1.y * U2.z;

	_row2.x =  interMedPseudo2.x * U1.x + interMedPseudo2.y * U2.x;
	_row2.y =  interMedPseudo2.x * U1.y + interMedPseudo2.y * U2.y;
	_row2.z =  interMedPseudo2.x * U1.z + interMedPseudo2.y * U2.z;

	return true;
}

#endif


