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
	Color trace_ray ( const Ray & , vector<QuadTreeNode *> & );

	//Scene * scene ;

};

struct EdgeList {

	bool isEdge ;
	bool isVertex ;

	Vec3 vTx ;
	Vec3 LineQ;
	Vec3 LineV;

	EdgeList * next ;
	EdgeList * prev ;

	double pointRadiosity ;
	int nodeNumber;
};


inline EdgeList * createEdgeList ( const BSP_Node & _bsp_N ) {

	EdgeList * newEdgeListHead = new EdgeList () ;
	EdgeList * edgeListCurrent = newEdgeListHead ;
	EdgeList * edgeListTail = newEdgeListHead ;

	vector<Vec3> triPoints;
	triPoints.emplace_back( _bsp_N.triVert1 );
	triPoints.emplace_back( _bsp_N.triVert2 );
	triPoints.emplace_back( _bsp_N.triVert3 );
	triPoints.emplace_back( _bsp_N.triVert1 );

	for ( int i = 0 ; i < 3 ; ++ i ) {
		edgeListCurrent->isVertex = true ;
		edgeListCurrent->vTx = triPoints[i] ;
		edgeListCurrent->nodeNumber = i;
		edgeListCurrent->next = new EdgeList() ;

		edgeListCurrent->next->isEdge = true ;
		edgeListCurrent->next->LineQ = edgeListCurrent->vTx ;
		edgeListCurrent->next->LineV = triPoints[i+1] - triPoints[i] ;
		edgeListCurrent->next->prev = edgeListCurrent ;
		edgeListCurrent = edgeListCurrent->next ;

		if ( i == 2 ) {
			edgeListCurrent->next = newEdgeListHead ;
			newEdgeListHead->prev = edgeListCurrent;
		}
		else {
			edgeListCurrent->next = new EdgeList();
			edgeListCurrent->next->prev = edgeListCurrent ;
			edgeListCurrent = edgeListCurrent->next ;
		}
	}
	return newEdgeListHead ;
}

inline void deleteEdgeList ( EdgeList * e ) {
	EdgeList * current = e ;
	EdgeList * currentNext = e->next ;

	while ( current != NULL ) {
		if ( current->prev != NULL ) current->prev->next = NULL;
		if ( current->next != NULL ) current->next->prev = NULL;
		delete current ;
		current = NULL ;
		current = currentNext ;
		if ( currentNext != NULL ) currentNext = currentNext->next ;
	}
}

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
	Vec3 P = _P - _Q ;
	double tx = ( _V.x == 0.0 ? 0.0 : P.x / _V.x   ) ;
	double ty = ( _V.y == 0.0 ? 0.0 : P.y / _V.y   ) ;
	double tz = ( _V.z == 0.0 ? 0.0 : P.z / _V.z   ) ;

	return (tx * _V.x == P.x && ty * _V.y == P.y && tz * _V.z == P.z) ;
}

inline bool isPointNearLine (const Vec3 & _P, Vec3 & _Q, Vec3 & _V ) {
	Vec3 P = _P - _Q ;
	double tx = ( _V.x == 0.0 ? 0.0 : P.x / _V.x   ) ;
	double ty = ( _V.y == 0.0 ? 0.0 : P.y / _V.y   ) ;
	double tz = ( _V.z == 0.0 ? 0.0 : P.z / _V.z   ) ;

	double threshold = Epsilon * 100.0 ;

	return (abs (tx * _V.x - P.x) <= threshold && abs(ty * _V.y - P.y) <= threshold && abs(tz * _V.z - P.z) <= threshold) ;
}

inline Vec3 centerOfTriangle (QuadTreeNode * _Qnode ) {

	Vec3 QLine1 = _Qnode->triVert1 ;
	Vec3 tempAdjLinePoint = (0.5 * (_Qnode->triVert3 - _Qnode->triVert2)) + _Qnode->triVert2 ;
	//Vec3 VLine1 = Unit(((_Qnode->triVert3 - _Qnode->triVert2) / 2.0) - _Qnode->triVert1) ; 
	Vec3 VLine1 = tempAdjLinePoint -  _Qnode->triVert1 ;

	Vec3 QLine2 = _Qnode->triVert2 ;
	tempAdjLinePoint = ( 0.5 * ( _Qnode->triVert1 - _Qnode->triVert3 )) + _Qnode->triVert3 ;
	Vec3 VLine2 = tempAdjLinePoint - _Qnode->triVert2 ;

	Mat3x2 A;
	A(0,0) = VLine1.x ;
	A(0,1) = -VLine2.x ;
	A(1,0) = VLine1.y ;
	A(1,1) = -VLine2.y ;
	A(2,0) = VLine1.z ;
	A(2,1) = -VLine2.z ;

	Vec3 _A, _B  ;

	calculatePseudoInverse ( A, _A, _B ) ;

	double t = _A.x * (QLine2 - QLine1).x + 
			   _A.y * (QLine2 - QLine1).y +
			   _A.z * (QLine2 - QLine1).z ;

	return (t * VLine1) + QLine1 ;
}

// check to make sure that the vec3 we are adding to a quadnode's adj list doesn't already exist.
inline bool isAlreadyInAdjList ( Vec3 & _P, QuadTreeNode * _A, QuadTreeNode * _B  )  {

	for ( auto & unordered_set_vector : _A->nextAdj[_B] ) {
		if ( relativeClose ( unordered_set_vector, _P, 0.5 )) 
			return true;
	}
	return false;

}

inline Color returnQuadTreeColor ( Vec3 & _p , QuadTreeNode * _Qnode ) {

	EdgeList * edge_list = createEdgeList ( *dynamic_cast<BSP_Node *>(_Qnode) );
	EdgeList * currentEdgeListNode = edge_list ;

	Vec3 centerPointOfTriangle = centerOfTriangle ( _Qnode ) ;

	int numberOfPolygonPoints = 0;

	vector<std::pair<Vec3, double>> pointsInsideTriangle ;
	vector<std::pair<Vec3, double>> pointsOnTriangle ;
	//std::map<Vec3, double> interpolationPoints ;
	std::unordered_map<Vec3, double> interpolationPoints ;
	// UNUSED FOR NOW
	std::unordered_map<Vec3, int> numberOfSamples;

	for ( const auto & i : _Qnode->nextAdj ) {
		for ( const auto &j : i.second ) {

			if ( doesPointLieOnLine ( j, _Qnode->triVert1, _Qnode->triVert2 - _Qnode->triVert1 ) || 
				 doesPointLieOnLine ( j, _Qnode->triVert2, _Qnode->triVert3 - _Qnode->triVert2 ) ||
				 doesPointLieOnLine ( j, _Qnode->triVert3, _Qnode->triVert1 - _Qnode->triVert3 ) ) {
					 interpolationPoints.insert ( std::pair<Vec3, double> (j, 0.0) );
					 //cout << "Added vec3 to interpolationPoints : " << j << endl;
			}
			else{
				//pointsInsideTriangle.emplace_back ( std::pair<Vec3, double> ( j , 0.0 ) ) ;
				//cout << "Didnt add vec3 to interpolationPoints : " << j << endl; 
				}

		}
	}
	pointsInsideTriangle.emplace_back ( std::pair<Vec3, double> ( centerPointOfTriangle, _Qnode->radiosityValue )) ;

	// for all points not on the lines of the triangle
	for ( const auto &v : pointsInsideTriangle ) {
		//for all triangles adj to center triangle
		for ( const auto & adj : _Qnode->nextAdj ) {
			// for all the adj points that are shared between the center triangle and the adj triangle
			for ( const auto &adjPoints : adj.second ) {

				if ( doesPointLieOnLine ( adjPoints, _Qnode->triVert1, _Qnode->triVert2 - _Qnode->triVert1 ) || 
				     doesPointLieOnLine ( adjPoints, _Qnode->triVert2, _Qnode->triVert3 - _Qnode->triVert2 ) ||
				     doesPointLieOnLine ( adjPoints, _Qnode->triVert3, _Qnode->triVert1 - _Qnode->triVert3 ) ) {
						 double fullLength = Length ( centerPointOfTriangle - adjPoints ) + Length ( adjPoints - centerOfTriangle(adj.first ));

						 interpolationPoints[adjPoints] += ((1.0 - (Length ( centerPointOfTriangle - adjPoints ) / fullLength) ) * _Qnode->radiosityValue) +
							                              ((1.0 - (Length ( adjPoints - centerOfTriangle(adj.first )) / fullLength) ) * (adj.first)->radiosityValue) ;

						 numberOfSamples[adjPoints] ++;
				}
			}
		}
	}

	for ( auto & v : interpolationPoints ) {
		v.second /= numberOfSamples[v.first];
	}

	double w1, w2, w3;

	//boundary case : is the point near the side of the triangle
	if ( isPointNearLine ( _p, _Qnode->triVert1, _Qnode->triVert2 - _Qnode->triVert1 ) || 
		 isPointNearLine ( _p, _Qnode->triVert2, _Qnode->triVert3 - _Qnode->triVert2 ) ||
	     isPointNearLine ( _p, _Qnode->triVert3, _Qnode->triVert1 - _Qnode->triVert3 ) ) {
			 double areaOfQuadNode = Area ( _Qnode->triVert1, _Qnode->triVert2, _Qnode->triVert3 );
			 w1 = Area ( _p, _Qnode->triVert2, _Qnode->triVert3 ) / areaOfQuadNode ;
			 w2 = Area ( _p, _Qnode->triVert1, _Qnode->triVert3 ) / areaOfQuadNode ;
			 w3 = Area ( _p, _Qnode->triVert1, _Qnode->triVert2 ) / areaOfQuadNode ;
	}
	else {
			w1 = (coTangent( _p, _Qnode->triVert1, _Qnode->triVert2 ) + coTangent( _p, _Qnode->triVert1, _Qnode->triVert3 )) / (pow(Length( _p - _Qnode->triVert1 ), 2 ));
			w1 /= 3.0 ;
			w2 = (coTangent( _p, _Qnode->triVert2, _Qnode->triVert1 ) + coTangent( _p, _Qnode->triVert2, _Qnode->triVert3 )) / (pow(Length( _p - _Qnode->triVert2 ), 2 ));
			w2 /= 3.0 ;
			w3 = (coTangent( _p, _Qnode->triVert3, _Qnode->triVert1 ) + coTangent( _p, _Qnode->triVert3, _Qnode->triVert2 )) / (pow(Length( _p - _Qnode->triVert3 ), 2 ));
			w3 /= 3.0 ;
	}

	double areaOfQuadNode = Area ( _Qnode->triVert1, _Qnode->triVert2, _Qnode->triVert3 );
	w1 = Area ( _p, _Qnode->triVert2, _Qnode->triVert3 ) / areaOfQuadNode ;
	w2 = Area ( _p, _Qnode->triVert1, _Qnode->triVert3 ) / areaOfQuadNode ;
	w3 = Area ( _p, _Qnode->triVert1, _Qnode->triVert2 ) / areaOfQuadNode ;

	double triWeights[3] = { w1, w2, w3 };

	double finalRadiosityValue = 0.0 ; 

	do {

		bool found = false;

		for ( const auto &v : interpolationPoints ) {

			if ( Length(v.first - currentEdgeListNode->vTx) <= Epsilon ) { finalRadiosityValue += (triWeights[currentEdgeListNode->nodeNumber] * v.second) ; found = true ; break; }
			
			//else cout << "error in Color returnQuadTreeColor " << endl;
		}

		if ( !found ) {
			finalRadiosityValue +=  triWeights[currentEdgeListNode->nodeNumber] * _Qnode->radiosityValue ;	
		}

		currentEdgeListNode = currentEdgeListNode->next->next ;
	
	} while ( currentEdgeListNode != edge_list );

	deleteEdgeList ( edge_list );

	//finalRadiosityValue /= 3.0;

	return Color ( finalRadiosityValue, finalRadiosityValue, finalRadiosityValue );
}



#endif