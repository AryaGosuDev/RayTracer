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
#include "util.h"

struct Radiosity_Helper : public Radiosity { 

	Radiosity_Helper() ;
	Radiosity_Helper(Scene *) ;
	virtual ~Radiosity_Helper() ;

	void detectTriangleIntersections ( QuadTreeNode * , QuadTreeNode *) ;
	void returnFilledElementsOfObject ( QuadTreeNode * _v , vector<QuadTreeNode * > & _a ) ;
	Color trace_ray ( Ray & , vector<QuadTreeNode *> & );

	Scene * scene ;

};

inline bool isInView ( QuadTreeNode * _A, QuadTreeNode * _B, Radiosity * _this ) {
	int hits = 0 ;
	int formFactorLoopIterations = 30 ;

	//shoot formFactorLoopIterations # of rays from element i to element j
	for ( int k = 0 ; k < formFactorLoopIterations ; ++k ) {
		Vec3 xi = returnRandomPointOnElement ( _A ) ;
		Vec3 xj = returnRandomPointOnElement ( _B ) ;

		if ( _this->isVisibleXiToXj ( _A, _B, xi, xj ) ) hits++; 
	}
	if ( hits >= formFactorLoopIterations * 0.95 ) return true ;
	return false;
}


inline bool doesPointLieOnLine (const Vec3 & _P, Vec3 & _Q, Vec3 & _V ) {
	Vec3 P = _Q - _P ;
	double tx = ( _V.x == 0.0 ? 0.0 : P.x / _V.x   ) ;
	double ty = ( _V.y == 0.0 ? 0.0 : P.y / _V.y   ) ;
	double tz = ( _V.z == 0.0 ? 0.0 : P.z / _V.z   ) ;

	if ( tx != 0.0 ) {
		if ( tx * _V.y == P.y && tx * _V.z == P.z ) return true;
		return false;
	}
	if ( ty != 0.0 ) {
		if ( ty * _V.x == P.x && ty * _V.z == P.z ) return true;
		return false;
	}
	if ( tz != 0.0 ) {
		if ( tz * _V.x == P.x && tz * _V.y == P.y ) return true;
		return false;
	}
	return true;
}


inline Vec3 centerOfTriangle (QuadTreeNode * _Qnode ) {

	Vec3 QLine1 = _Qnode->triVert1 ;
	Vec3 VLine1 = Unit(((_Qnode->triVert3 - _Qnode->triVert2) / 2.0) - _Qnode->triVert1) ; 

	Vec3 QLine2 = _Qnode->triVert2 ;
	Vec3 VLine2 = Unit(((_Qnode->triVert1 - _Qnode->triVert3) / 2.0) - _Qnode->triVert2) ; 

	Mat3x2 A;
	A(0,0) = VLine1.x ;
	A(0,1) = -VLine2.x ;
	A(1,0) = VLine1.y ;
	A(1,1) = -VLine2.y ;
	A(2,0) = VLine1.z ;
	A(2,1) = VLine2.z ;

	Vec3 _A, _B  ;

	calculatePseudoInverse ( A, _A, _B ) ;

	double t = _A.x * (QLine2 - QLine1).x + 
			   _A.y * (QLine2 - QLine1).y +
			   _A.z * (QLine2 - QLine1).z ;

	return (t * VLine1) + QLine1 ;
}


inline Color returnQuadTreeColor ( Vec3 & _p , QuadTreeNode * _Qnode ) {

	Vec3 centerPointOfTriangle = centerOfTriangle ( _Qnode ) ;

	int numberOfPolygonPoints = 0;

	vector<std::pair<Vec3, double>> pointsInsideTriangle ;
	vector<std::pair<Vec3, double>> pointsOnTriangle ;
	std::map<Vec3, double> interpolationPoints ;
	std::map<Vec3, int> numberOfSamples;

	
	for ( const auto & i : _Qnode->nextAdj ) {
		for ( const auto &j : i.second ) {

			if ( doesPointLieOnLine ( j, _Qnode->triVert1, _Qnode->triVert2 - _Qnode->triVert1 ) || 
				 doesPointLieOnLine ( j, _Qnode->triVert2, _Qnode->triVert3 - _Qnode->triVert2 ) ||
				 doesPointLieOnLine ( j, _Qnode->triVert3, _Qnode->triVert1 - _Qnode->triVert3 ) ) {
					 interpolationPoints.insert ( std::pair<Vec3, double> (j, 0.0) );
			}
			else pointsInsideTriangle.emplace_back ( std::pair<Vec3, double> ( j , 0.0 ) ) ;
		}
	}
	pointsInsideTriangle.emplace_back ( std::pair<Vec3, double> ( centerPointOfTriangle, _Qnode->radiosityValue )) ;

	for ( const auto &v : pointsInsideTriangle ) {

		for ( const auto & adj : _Qnode->nextAdj ) {

			for ( const auto &adjPoints : adj.second ) {

				if ( doesPointLieOnLine ( adjPoints, _Qnode->triVert1, _Qnode->triVert2 - _Qnode->triVert1 ) || 
				     doesPointLieOnLine ( adjPoints, _Qnode->triVert2, _Qnode->triVert3 - _Qnode->triVert2 ) ||
				     doesPointLieOnLine ( adjPoints, _Qnode->triVert3, _Qnode->triVert1 - _Qnode->triVert3 ) ) {
						 double fullLength = Length ( centerPointOfTriangle - adjPoints ) + Length ( adjPoints - centerOfTriangle(adj.first ));

						 interpolationPoints[adjPoints] += ((1.0 - Length ( centerPointOfTriangle - adjPoints )) / fullLength ) * _Qnode->radiosityValue +
							                              ((1.0 - Length ( adjPoints - centerOfTriangle(adj.first ))) / fullLength ) * _Qnode->radiosityValue ;

						 numberOfSamples[adjPoints] ++;
				}
			}

		}
	}
	
	double w1 = (coTangent( _p, _Qnode->triVert1, _Qnode->triVert2 ) + coTangent( _p, _Qnode->triVert1, _Qnode->triVert3 )) / (pow(Length( _p - _Qnode->triVert1 ), 2 ));
	w1 /= 3.0 ;
	double w2 = (coTangent( _p, _Qnode->triVert2, _Qnode->triVert1 ) + coTangent( _p, _Qnode->triVert2, _Qnode->triVert3 )) / (pow(Length( _p - _Qnode->triVert2 ), 2 ));
	w2 /= 3.0 ;
	double w3 = (coTangent( _p, _Qnode->triVert3, _Qnode->triVert1 ) + coTangent( _p, _Qnode->triVert3, _Qnode->triVert2 )) / (pow(Length( _p - _Qnode->triVert3 ), 2 ));
	w3 /= 3.0 ;

	double finalRadiosityValue = 0.0 ; 

	for ( const auto &v : interpolationPoints ) {

		if ( Length(v.first - _Qnode->triVert1) < Epsilon ) finalRadiosityValue += (w1 * v.second) ;
		else if ( Length(v.first - _Qnode->triVert2) < Epsilon ) finalRadiosityValue += (w2 * v.second) ;
		else if ( Length(v.first - _Qnode->triVert3) < Epsilon ) finalRadiosityValue += (w3 * v.second) ;

		else cout << "Error in inline Color returnQuadTreeColor " << endl ;
	}

	return Color ( finalRadiosityValue, finalRadiosityValue, finalRadiosityValue );
}

#endif