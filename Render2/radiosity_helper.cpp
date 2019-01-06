// Class and implementation file used to assist in radioisty calculations
// especially with triangle - trianlge intersection and re-triangulation.
// Also rasterization of 

#include "ray_tracer.h"
#include "radiosity_helper.h"
#include <stack>


struct EdgeList {

	bool isEdge ;
	bool isVertex ;

	Vec3 vTx ;
	Vec3 LineQ;
	Vec3 LineV;

	EdgeList * next ;
	EdgeList * prev ;
};

struct FormFactorStackNode {

	FormFactorStackNode ( QuadTreeNode * _quadTreeNode  ) {
		quadTreeNode = _quadTreeNode;
		if ( quadTreeNode->children.size() > 0 ) internal = true ; 
		else internal = false;

		visited = false;
	}

	QuadTreeNode * quadTreeNode;
	bool visited ;
	bool internal; 
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
	EdgeList * head = e ;
	EdgeList * current = e ;
	EdgeList * currentNext = e->next ;

	while ( current != NULL ) {
		current->next = current->prev = NULL;
		delete current ;
		current = currentNext ;
		currentNext = currentNext->next ;
	}
}

// find the intersections of 2 adj triangles. Adds them to the triangle adj list along with intersections.
inline void findTriangleAdjPoints ( QuadTreeNode * _1 , QuadTreeNode * _2, EdgeList * _EL1, EdgeList * _EL2  ) {

	EdgeList * EL1Head = _EL1 ;
	EdgeList * EL2Head = _EL2 ;
	EdgeList * EL1curr = _EL1 ;
	EdgeList * EL2curr = _EL2 ;

	do  {

		EL1curr = EL1curr->next;

		Vec3 QLine = EL1curr->LineQ;
		Vec3 VLine = EL1curr->LineV ;

		do  {
			EL2curr = EL2curr->next;

			Mat3x2 A ;
			A(0,0) = EL1curr->LineV.x ;
			A(0,1) = -1.0 * EL2curr->LineV.x ;
			A(1,0) = EL1curr->LineV.y ;
			A(1,1) = -1.0 * EL2curr->LineV.y ;
			A(2,0) = EL1curr->LineV.z ;
			A(2,1) = -1.0 * EL2curr->LineV.z ;

			//Vec3 delQ ( QLine.x - EL2curr->LineQ.x,  QLine.y - EL2curr->LineQ.y,  QLine.z - EL2curr->LineQ.z ) ;
			Vec3 delQ ( EL2curr->LineQ.x - QLine.x , EL2curr->LineQ.y - QLine.y, EL2curr->LineQ.z - QLine.z ) ;

			Vec3 _A, _B;

			// line is the same line as a side of the triangle
			if ( calculatePseudoInverse ( A, _A, _B ) ) {

				double tP = _A.x * delQ.x + _A.y * delQ.y + _A.z * delQ.z  ;
				double tTriangle = _B.x * delQ.x + _B.y * delQ.y + _B.z * delQ.z  ;

				if ( tTriangle >= 0.0 && tTriangle <= 1.0 ) {
					//_1->nextAdj.insert( _2,  ) ;
					_1->nextAdj[_2].insert ( Vec3 ( tTriangle * EL2curr->LineV.x + EL2curr->LineQ.x ,  tTriangle * EL2curr->LineV.y + EL2curr->LineQ.y,  tTriangle * EL2curr->LineV.z + EL2curr->LineQ.z ) ) ;
					//_2->nextAdj.emplace( _1 );
					_2->nextAdj[_1].insert ( Vec3 ( tTriangle * EL2curr->LineV.x + EL2curr->LineQ.x ,  tTriangle * EL2curr->LineV.y + EL2curr->LineQ.y,  tTriangle * EL2curr->LineV.z + EL2curr->LineQ.z ) ) ;
				}
			}
			EL2curr = EL2curr->next;
		} while ( EL2curr != EL2Head ) ;
		EL1curr = EL1curr->next ;
	} while ( EL1curr != EL1Head ) ;
}

// find the intersections as one triangle goes into another triangle
inline void linePlaneIntersectionAdj  ( QuadTreeNode * _1 , QuadTreeNode * _2, EdgeList * _EL1, EdgeList * _EL2  ) {

	EdgeList * EL1Head = _EL1 ;
	EdgeList * EL2Head = _EL2 ;
	EdgeList * EL1curr = _EL1 ;
	EdgeList * EL2curr = _EL2 ;

	do  {

		EL1curr = EL1curr->next;

		Vec3 QLine = EL1curr->LineQ;
		Vec3 VLine = EL1curr->LineV ;

			EL2curr = EL2curr->next;

			double D = (_2->triNormal * -1.0) * _2->triVert1 ;
			double Denom = _2->triNormal * VLine ;

			// line and plane are not parallel
			if ( Denom > Epsilon && Denom < -Epsilon ) {
				double t = (-1.0 * ( _2->triNormal * QLine + D )) / Denom ;
				Vec3 P =  (VLine * t) + QLine ;

				int result = sideOfLine3D ( P, EL2curr->LineV ^ _2->triNormal );
				EL2curr = EL2curr->next->next;
				int result2 = sideOfLine3D ( P, EL2curr->LineV ^ _2->triNormal );
				EL2curr = EL2curr->next->next;
				int result3 = sideOfLine3D ( P, EL2curr->LineV ^ _2->triNormal );

				if ( result >= 0 && result2 >= 0 && result3 >= 0 ) {
					//_2->nextAdj.emplace(_1);
					_2->nextAdj[_1].insert(P );
				}
				if ( result <= 0 && result2 <= 0 && result3 <= 0 ) {
					//_2->nextAdj.emplace(_1);
					_2->nextAdj[_1].insert(P );
				}
			}
			EL1curr = EL1curr->next ;
		} while ( EL1curr != EL1Head );

	EL1curr = _EL1 ;
	EL2curr = _EL2 ;

	do  {

		EL2curr = EL2curr->next;

		Vec3 QLine = EL2curr->LineQ;
		Vec3 VLine = EL2curr->LineV ;

			EL1curr = EL1curr->next;

			double D = (_1->triNormal * -1.0) * _1->triVert1 ;
			double Denom = _1->triNormal * VLine ;

			// line and plane are not parallel
			if ( Denom > Epsilon && Denom < -Epsilon ) {
				double t = (-1.0 * ( _1->triNormal * QLine + D )) / Denom ;
				Vec3 P =  (VLine * t) + QLine ;

				int result = sideOfLine3D ( P, EL1curr->LineV ^ _1->triNormal );
				EL1curr = EL1curr->next->next;
				int result2 = sideOfLine3D ( P, EL1curr->LineV ^ _1->triNormal );
				EL1curr = EL1curr->next->next;
				int result3 = sideOfLine3D ( P, EL1curr->LineV ^ _1->triNormal );

				if ( result >= 0 && result2 >= 0 && result3 >= 0 ) {
					//_1->nextAdj.emplace(_2);
					_1->nextAdj[_2].insert(P );
				}
				if ( result <= 0 && result2 <= 0 && result3 <= 0 ) {
					//_1->nextAdj.emplace(_2);
					_1->nextAdj[_2].insert(P );
				}
			}
			EL2curr = EL2curr->next ;
		} while ( EL2curr != EL2Head );
}


inline int returnNumOfIntersectionsWithTriangle ( Vec3 & _Q, Vec3 & _V, const BSP_Node & _bsp_N, vector<Vec3> & newIntersections ) {
	
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
	int numOfIntersections = 0;

	vector<Vec3> triPoints;
	triPoints.emplace_back( _bsp_N.triVert1 );
	triPoints.emplace_back( _bsp_N.triVert2 );
	triPoints.emplace_back( _bsp_N.triVert3 );
	triPoints.emplace_back( _bsp_N.triVert1 );

	for ( int i = 0 ; i < 3 ; i ++ ) {

		//Vec3 QLine = _bsp_N.triVert1 ;
		//Vec3 VLine = _bsp_N.triVert2 - _bsp_N.triVert1 ;
		Vec3 QLine = triPoints[i];
		Vec3 VLine = triPoints[i+1] - triPoints[i] ;

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

		if ( tTriangle > 1.0 || tTriangle < 0.0 ) return 0;

		//if ( tTriangle >= 0.0 && tTriangle <= 1000.0 * Epsilon ) return 1;
		//if ( tTriangle <= 1.0 && tTriangle >= 1000.0 * Epsilon ) return 1;

		if ( tTriangle >= 0.0 && tTriangle <= 1.0 ) numOfIntersections++;
		
		Vec3 pointOfIntersection ( tTriangle * VLine.x + QLine.x ,  tTriangle * VLine.y + QLine.y,  tTriangle * VLine.z + QLine.z ); 

		newIntersections.emplace_back( pointOfIntersection );
	}
	//return 2 ;
	return numOfIntersections ;
}

Radiosity_Helper::Radiosity_Helper() {
}

Radiosity_Helper::Radiosity_Helper( Scene * _scene) : scene ( _scene ) {
	//returnNumOfIntersectionsWithTriangle ( Vec3 ( 0.0, 0.0, 0.0 ) ,  Vec3 ( 0.0, 0.0, 0.0 ), NULL ) ;
}

Radiosity_Helper::~Radiosity_Helper() {
}



// Discovered that these two objects' AABB interesect. Now find which triangles intersect and re-triangulate the triangles.
// Avoid mesh discontinuities. Make sure all appropriate adjacent triangles are detected that span different objects.
// Find the intersecting parametric line between two planes P(t) = Q + t(N(1) X N(2)). 
void Radiosity_Helper::detectTriangleIntersections ( QuadTreeNode * _1 , QuadTreeNode * _2) {
	
	EdgeList * edgeL_1 = createEdgeList( *(dynamic_cast<BSP_Node *>(_1)));
	EdgeList * edgeL_2 = createEdgeList( *(dynamic_cast<BSP_Node *>(_2)));
	vector<Vec3> newIntersections1 ;
	vector<Vec3> newIntersections2 ;
	Vec3 V = _1->triNormal ^ _2->triNormal  ;
	// the triangles are co-planar
	if (V == Vec3 ( 0.0, 0.0, 0.0 ) ) {
		findTriangleAdjPoints ( _1, _2 , edgeL_1, edgeL_2 );
	} 
	else { // the triangles are not co-planar
		// Find parametric line between both planes. The two planes are formed from the normals of the two triangles.
		Mat3x3 Q;
		Q(0,0) = _1->triNormal.x ;
		Q(0,1) = _1->triNormal.y ;
		Q(0,2) = _1->triNormal.z ;
		Q(1,0) = _2->triNormal.x ;
		Q(1,1) = _2->triNormal.y ;
		Q(1,2) = _2->triNormal.z ;
		Q(2,0) = V.x;
		Q(2,1) = V.y;
		Q(2,2) = V.z;
		Q = Inverse(Q);
		Vec3 D ( (-1.0) * ( (-1.0 * _1->triNormal ) * _1->triVert1 ),
					(-1.0) * ( (-1.0 * _2->triNormal ) * _2->triVert1 ),
					0.0 );

		Vec3 QVec = Q * D ; 

		int numOfIntersections_1 = returnNumOfIntersectionsWithTriangle ( QVec, V, *(dynamic_cast<BSP_Node *>(_1)), newIntersections1 );
		int numOfIntersections_2 = returnNumOfIntersectionsWithTriangle ( QVec, V, *(dynamic_cast<BSP_Node *>(_2)), newIntersections2 );

		//triangles do not intersect
		if ( numOfIntersections_1 == 0 || numOfIntersections_2 == 0 ) return ;

		if ( numOfIntersections_1 == -1 && numOfIntersections_2 == -1 && isInView ( _1, _2, dynamic_cast<Radiosity *> (this) ))
			findTriangleAdjPoints ( _1, _2, edgeL_1, edgeL_2 ); 

		else if ( ((numOfIntersections_1 == -1 && numOfIntersections_2 == 1) || (numOfIntersections_1 == 1 && numOfIntersections_2 == -1)) && 
			        isInView ( _1, _2, dynamic_cast<Radiosity *> (this) ) )
			findTriangleAdjPoints ( _1, _2 , edgeL_1 , edgeL_2 );

		else if ( numOfIntersections_1 == 1 && numOfIntersections_2 == 1 && isInView ( _1, _2, dynamic_cast<Radiosity *> (this) ) )
			findTriangleAdjPoints ( _1, _2, edgeL_1, edgeL_2 );

		else if ( numOfIntersections_1 == -1 && numOfIntersections_2 == 2 ) {
			linePlaneIntersectionAdj  ( _1 , _2,edgeL_1, edgeL_2  ) ;
			for ( auto &v : _2->nextAdj[_1] ) {
				//_1->nextAdj.emplace( _2) ;
				_1->nextAdj[_2].insert(v );
			}
		}

		else if ( numOfIntersections_1 == 2 && numOfIntersections_2 == -1 ) {
			linePlaneIntersectionAdj  ( _1 , _2,edgeL_1, edgeL_2  ) ;
			for ( auto &v : _1->nextAdj[_2] ) {
				//_2->nextAdj.emplace( _1) ;
				_2->nextAdj[_1].insert(v );
			}
		} 

		else if ( numOfIntersections_1 == 1 && numOfIntersections_2 == 2 ) {
			linePlaneIntersectionAdj  ( _1 , _2,edgeL_1, edgeL_2  ) ;
			for ( auto &v : _2->nextAdj[_1] ) {
				//_1->nextAdj.emplace( _2) ;
				_1->nextAdj[_2].insert(v );
			}
		}

		else if ( numOfIntersections_1 == 2 && numOfIntersections_2 == 1 ) {
			linePlaneIntersectionAdj  ( _1 , _2,edgeL_1, edgeL_2  ) ;
			for ( auto &v : _1->nextAdj[_2] ) {
				//_2->nextAdj.emplace( _1) ;
				_2->nextAdj[_1].insert(v );
			}
		}  

		else if ( numOfIntersections_1 == 2 && numOfIntersections_2 == 2 ) {
			linePlaneIntersectionAdj  ( _1 , _2,edgeL_1, edgeL_2  ) ;
		}  

		/*
	for ( int i = 0 ; i < _1->object->triangles.size() ; ++ i) {
		for ( int j = 0 ; j < _2->object->triangles.size() ; j++ ) {

			EdgeList * edgeL_1 = createEdgeList( _1->object->triangles[i] );
			EdgeList * edgeL_2 = createEdgeList( _2->object->triangles[j] );
			vector<Vec3> newIntersections1 ;
			vector<Vec3> newIntersections2 ;
			Vec3 V = _1->object->triangles[i].triNormal ^ _2->object->triangles[j].triNormal  ;
			// the triangles are co-planar
			if (V == Vec3 ( 0.0, 0.0, 0.0 ) ) {
				if ( findTriangleAdjPoints ( _1, _2 ) ) { 
					_1->nextAdj.emplace_back( _2 ) ;
					_2->nextAdj.emplace_back( _1 ) ;
				}
			} 
			else { // the triangles are not co-planar
				// Find parametric line between both planes. The two planes are formed from the normals of the two triangles.
				Mat3x3 Q;
				Q(0,0) = _1->object->triangles[i].triNormal.x ;
				Q(0,1) = _1->object->triangles[i].triNormal.y ;
				Q(0,2) = _1->object->triangles[i].triNormal.z ;
				Q(1,0) = _2->object->triangles[j].triNormal.x ;
				Q(1,1) = _2->object->triangles[j].triNormal.y ;
				Q(1,2) = _2->object->triangles[j].triNormal.z ;
				Q(2,0) = V.x;
				Q(2,1) = V.y;
				Q(2,2) = V.z;
				Q = Inverse(Q);
				Vec3 D ( (-1.0) * ( (-1.0 * _1->object->triangles[i].triNormal ) * _1->object->triangles[i].triVert1 ),
					     (-1.0) * ( (-1.0 * _2->object->triangles[j].triNormal ) * _2->object->triangles[j].triVert1 ),
						 0.0 );

				Vec3 QVec = Q * D ; 

				int numOfIntersections_1 = returnNumOfIntersectionsWithTriangle ( QVec, V, _1->object->triangles[i], newIntersections1 );
				int numOfIntersections_2 = returnNumOfIntersectionsWithTriangle ( QVec, V, _2->object->triangles[j], newIntersections2 );

				//triangles do not intersect
				if ( numOfIntersections_1 == 0 || numOfIntersections_2 == 0 ) return ; 

				// the line between the triangles is the same line as the side of the triangle, both triangles share a side
				if ( numOfIntersections_1 == -1 && numOfIntersections_2 == -1  ) {
					    //do they share points?
						if ( findTriangleAdjPoints ( _1, _2 ) ) { 
							_1->nextAdj.emplace_back( _2 ) ;
							_2->nextAdj.emplace_back( _1 ) ;
						}
						else cout << "Case in detectTriangleIntersections not detected" << endl; 
				}
				if ( numOfIntersections_1 == -1 && numOfIntersections_2 == 1 ) {
					if ( relativeClose ( newIntersections2[0] , _1->triVert1 ) || 
						 relativeClose ( newIntersections2[0] , _1->triVert2 ) ||
						 relativeClose ( newIntersections2[0] , _1->triVert3 ) ) {
							_1->nextAdj.emplace_back( _2 ) ;
							_2->nextAdj.emplace_back( _1 ) ;
					}
					else {
						EdgeList * currentEdgeListNode = edgeL_1->next ;
						while ( !doesPointLieOnLine ( newIntersections2[0], currentEdgeListNode->LineQ, currentEdgeListNode->LineV ) ) {
							currentEdgeListNode = currentEdgeListNode->next->next ;
						}
						EdgeList * incidentEdge = currentEdgeListNode ;
						_1->children.emplace_back( new QuadTreeNode() ) ;
						_1->children.back()->triVert1 = newIntersections2[0] ;
						_1->children.back()->triVert2 = incidentEdge->next->vTx ; 
						_1->children.back()->triVert3 = incidentEdge->next->next->next->vTx ;
						_1->children.back()->triNormal = _1->triNormal;
						_1->children.emplace_back( new QuadTreeNode() ) ;
						_1->children.back()->triVert1 = newIntersections2[0] ;
						_1->children.back()->triVert2 = incidentEdge->prev->vTx ; 
						_1->children.back()->triVert3 = incidentEdge->prev->prev->prev->vTx ;
						_1->children.back()->triNormal = _1->triNormal;
						_1->children.back()->nextAdj.emplace_back ( _1->children.front
					}

				}
				if ( numOfIntersections_1 == 1 && numOfIntersections_2 == -1 ) {
					if ( relativeClose ( newIntersections1[0] , _2->triVert1 ) || 
						 relativeClose ( newIntersections1[0] , _2->triVert2 ) ||
						 relativeClose ( newIntersections1[0] , _2->triVert3 ) ) {
							_1->nextAdj.emplace_back( _2 ) ;
							_2->nextAdj.emplace_back( _1 ) ;
					}
					else {
						EdgeList * currentEdgeListNode = edgeL_2->next ;
						while ( !doesPointLieOnLine ( newIntersections1[0], currentEdgeListNode->LineQ, currentEdgeListNode->LineV ) ) {
							currentEdgeListNode = currentEdgeListNode->next->next ;
						}
						EdgeList * incidentEdge = currentEdgeListNode ;
						_2->children.emplace_back( new QuadTreeNode() ) ;
						_2->children.back()->triVert1 = newIntersections2[0] ;
						_2->children.back()->triVert2 = incidentEdge->next->vTx ; 
						_2->children.back()->triVert3 = incidentEdge->next->next->next->vTx ;
						_2->children.back()->triNormal = _2->triNormal;
						_2->children.emplace_back( new QuadTreeNode() ) ;
						_2->children.back()->triVert1 = newIntersections2[0] ;
						_2->children.back()->triVert2 = incidentEdge->prev->vTx ; 
						_2->children.back()->triVert3 = incidentEdge->prev->prev->prev->vTx ;
						_2->children.back()->triNormal = _2->triNormal;
					}

				}
				if ( numOfIntersections_1 == -1 && numOfIntersections_2 == 2 ) {

				}
			}
			deleteEdgeList ( edgeL_1 );
			deleteEdgeList ( edgeL_2 );
		}
	}
	*/		
     }
}

void Radiosity_Helper::returnFilledElementsOfObject ( QuadTreeNode * _v , vector<QuadTreeNode * > & _a ) {

	std::stack<FormFactorStackNode> stackFormFactor ;

	//for all the initial elements of the object node
	//for ( auto &elements : _v->children ) {
		
		for ( vector<QuadTreeNode *>::const_reverse_iterator iterQuads = _v->children.crbegin() ; iterQuads != _v->children.crend() ; ++iterQuads )
			stackFormFactor.push ( FormFactorStackNode(*iterQuads));
	//}

	while ( !stackFormFactor.empty() ) {
		FormFactorStackNode temp = stackFormFactor.top();
		stackFormFactor.pop();
		if ( temp.quadTreeNode->children.size() > 0 ) 
			for ( vector<QuadTreeNode *>::const_reverse_iterator iterQuads = temp.quadTreeNode->children.crbegin() ;
				iterQuads != temp.quadTreeNode->children.crend() ; 
				++iterQuads ) 

					stackFormFactor.push ( FormFactorStackNode(*iterQuads));
		else {
			_a.emplace_back ( temp.quadTreeNode ) ;
			//_a.emplace_back ( stackFormFactor.top().quadTreeNode ) ;
			//stackFormFactor.pop() ;
		}
	}
}



Color Radiosity_Helper::trace_ray ( Ray & _ray , vector<QuadTreeNode *> & _a  ) {

	Color c ;
	double t ;
	 
	double dist = 0.0 ;
	QuadTreeNode * hitNode = NULL ;
	Vec3 pointOfIntersection;

	// check for the intersection of rays with planes extended by the respective triangle
	for ( int i = 0 ; i < _a.size(); ++ i ) {

		double denom = _a[i]->triNormal * _ray.direction ;

		if ( denom != 0 ) {
			t = -1.0 * ( (_a[i]->triNormal * _ray.origin ) + ( (-1.0 * _a[i]->triNormal) * _a[i]->triVert1 ) ) / denom ;

			pointOfIntersection = t * _ray.direction + _ray.origin ;

			std::pair<double, double> thisUV = returnUVofTriangle( pointOfIntersection , _a[i] ) ;
			/*
			Mat3x2 A ;
			A(0,0) = (_a[i]->triVert2 - _a[i]->triVert1 ).x ;
			A(0,1) = (_a[i]->triVert3 - _a[i]->triVert1 ).x ;
			A(1,0) = (_a[i]->triVert2 - _a[i]->triVert1 ).y ;
			A(1,1) = (_a[i]->triVert3 - _a[i]->triVert1 ).y ;
			A(2,0) = (_a[i]->triVert2 - _a[i]->triVert1 ).z ;
			A(2,1) = (_a[i]->triVert3 - _a[i]->triVert1 ).z ;

			Vec3 _A, _B ;

			if ( calculatePseudoInverse ( A, _A, _B ) ) {
			
				double u = _A.x * pointOfIntersection.x + _A.y * pointOfIntersection.y + _A.z * pointOfIntersection.z  ;
				double v = _B.x * pointOfIntersection.x + _B.y * pointOfIntersection.y + _B.z * pointOfIntersection.z  ;
				*/
			double u = thisUV.first ;
			double v = thisUV.second ;
				// this ray hit the triangle
				if ( u + v <= 1.0 ) {
					hitNode = _a[i] ;
					if ( dist == 0.0 )
						dist = Length ( pointOfIntersection - _ray.origin ) ;
					else if ( dist > Length ( pointOfIntersection - _ray.origin ) ) 
						dist = Length ( pointOfIntersection - _ray.origin ) ;
				}
			//}
		}
	}

	// nothing was hit
	if ( hitNode == NULL ) return Green ;

	return returnQuadTreeColor ( pointOfIntersection, hitNode )  ;
}
