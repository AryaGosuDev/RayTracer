// Class and implementation file used to assist in radioisty calculations
// especially with triangle - trianlge intersection and re-triangulation.
// Also rasterization of 

#include "ray_tracer.h"
#include "radiosity_helper.h"

inline 

inline int returnNumOfIntersectionsWithTriangle ( Vec3 & _Q, Vec3 & _V, const BSP_Node & _bsp_N ) {
	
	/*//  for testing
	Mat3x2 A ;
	A(0,0) = -1.0;
	A(0,1) = -1.0;
	A(1,0) = 0.0;
	A(1,1) = 1.0;
	A(2,0) = 0.0;
	A(2,1) = 0.0;
	Vec3 _A ;
	Vec3 _B ;
	calculatePseudoInverse ( A, _A, _B );
	*/
	/*
	//  for testing
	Mat3x2 A ;
	A(0,0) = -1.0;
	A(0,1) = -2.0;
	A(1,0) = 0.0;
	A(1,1) = 0.0;
	A(2,0) = 0.0;
	A(2,1) = 0.0;
	Vec3 _A ;
	Vec3 _B ;
	calculatePseudoInverse ( A, _A, _B );
	*/
	
	int intersections = 0;

	Vec3 QLine = _bsp_N.triVert1 ;
	Vec3 VLine = _bsp_N.triVert2 - _bsp_N.triVert1 ;

	Mat3x2 A ;
	A(0,0) = _V.x ;
	A(0,1) = -1.0 * VLine.x ;
	A(1,0) = _V.y ;
	A(1,1) = -1.0 * VLine.y ;
	A(2,0) = _V.z ;
	A(2,1) = -1.0 * VLine.z ;

	Vec3 delQ ( QLine.x - _Q.x,  QLine.y - _Q.y,  QLine.z - _Q.z ) ;

	Vec3 _A, _B;

	// line is the same line as a side of the triangle
	if ( !calculatePseudoInverse ( A, _A, _B ) ) return -1;

	double tP = _A.x * delQ.x + _A.y * delQ.y + _A.z * delQ.z  ;
	double tTriangle = _B.x * delQ.x + _B.y * delQ.y + _B.z * delQ.z  ;
	
	Vec3 pointOfIntersection ( tTriangle * VLine.x + QLine.x ,  tTriangle * VLine.y + QLine.y,  tTriangle * VLine.z + QLine.z ); 



	return 0 ;
}

Radiosity_Helper::Radiosity_Helper() {

}

Radiosity_Helper::Radiosity_Helper( Scene * _scene) : scene ( _scene ) {

	
	returnNumOfIntersectionsWithTriangle ( Vec3 ( 0.0, 0.0, 0.0 ) ,  Vec3 ( 0.0, 0.0, 0.0 ), NULL ) ;

}

Radiosity_Helper::~Radiosity_Helper() {

}

// Discovered that these two objects' AABB interesect. Now find which triangles intersect and re-triangulate the triangles.
// Avoid mesh discontinuities. Make sure all appropriate adjacent triangles are detected that span different objects.
// Find the intersecting parametric line between two planes P(t) = Q + t(N(1) X N(2)). 
void Radiosity_Helper::detectTriangleIntersections ( Object * _1 , Object * _2) {
	for ( int i = 0 ; i < _1->triangles.size() ; ++ i) {
		for ( int j = 0 ; j < _2->triangles.size() ; j++ ) {
			Vec3 V = _1->triangles[i].triNormal ^ _2->triangles[j].triNormal  ;
			// the triangles are co-planar
			if (V == Vec3 ( 0.0, 0.0, 0.0 ) ) {
			}
			else { // the triangles are not co-planar
				// Find parametric line between both planes. The two planes are formed from the normals of the two triangles.
				Mat3x3 Q;
				Q(0,0) = _1->triangles[i].triNormal.x ;
				Q(0,1) = _1->triangles[i].triNormal.y ;
				Q(0,2) = _1->triangles[i].triNormal.z ;
				Q(1,0) = _2->triangles[j].triNormal.x ;
				Q(1,1) = _2->triangles[j].triNormal.y ;
				Q(1,2) = _2->triangles[j].triNormal.z ;
				Q(2,0) = V.x;
				Q(2,1) = V.y;
				Q(2,2) = V.z;
				Q = Inverse(Q);
				Vec3 D ( (-1.0) * ( (-1.0 * _1->triangles[i].triNormal ) * _1->triangles[i].triVert1 ),
					     (-1.0) * ( (-1.0 * _2->triangles[j].triNormal ) * _2->triangles[j].triVert1 ),
						 0.0 );

				Vec3 QVec = Q * D ; 

				int numOfIntersections_1 = returnNumOfIntersectionsWithTriangle ( QVec, V, _1->triangles[i] );
				int numOfIntersections_2 = returnNumOfIntersectionsWithTriangle ( QVec, V, _2->triangles[j] );


			}
		}
	}
}


