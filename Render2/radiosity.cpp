#include "ray_tracer.h"
#include "util.h"
#include <random>
#include <stack>
#include "radiosity_helper.h"

struct ObjectNode {

	struct LinkListAdjNode ;

	ObjectNode() {
		object = NULL ;
		aabbDoesIntersect = false ;
		finishedWithIntersections = false;
		AdjLL = NULL  ;
	}

	void addAdj ( ObjectNode * _adj ) {
		if ( AdjLL == NULL ) {
			AdjLL = new LinkListAdjNode () ;
			AdjLL->oN = _adj ;
			AdjLL->next = NULL ;
			tailAdjLL = AdjLL;
		}
		else {
			tailAdjLL->next = new LinkListAdjNode() ;
			tailAdjLL = tailAdjLL->next ;
			tailAdjLL->oN = _adj ;
			tailAdjLL->next = NULL ;
		}
	}

	Object * object;
	AABB aabbInfo ;
	bool aabbDoesIntersect ;
	bool finishedWithIntersections;

	LinkListAdjNode * AdjLL ;
	LinkListAdjNode * tailAdjLL ;

	QuadTreeNode * qNode ;

	struct LinkListAdjNode {

		LinkListAdjNode () {
		}

		void add (ObjectNode * _oN) {
			LinkListAdjNode * temp = this;
		}

		ObjectNode * oN ;
		LinkListAdjNode * next ;
	};
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

// objects in the scene excluding light
vector<ObjectNode *> radObjectNodes ;
// form factor matrix
double ** FF; 
// K matrix K = M - PFF where M = Identity, P = Reflectiviy and FF = form factor
double ** K ; 
// radiosities
double * B ;

// find all the AABB intersections of every object with each other
inline void findIntersections () {

	for ( auto &v : radObjectNodes ) {
		for ( auto &o : radObjectNodes ) {
			if (v != o ) {
				if  ( AABBIntersect ( v->aabbInfo, o->aabbInfo ) ) {
					v->aabbDoesIntersect = true;
					v->addAdj ( o ) ;
				}
			}
		}
	}
}

inline bool findTriangleAdj ( BSP_Node * _a, BSP_Node * _b ) {
	if ( relativeClose(_a->triVert1, _b->triVert1 ) || relativeClose(_a->triVert1, _b->triVert2 ) || relativeClose(_a->triVert1, _b->triVert3 ) ) return true ; 
	if ( relativeClose(_a->triVert2, _b->triVert1 ) || relativeClose(_a->triVert2, _b->triVert2 ) || relativeClose(_a->triVert2, _b->triVert3 ) ) return true ;
	if ( relativeClose(_a->triVert3, _b->triVert1 ) || relativeClose(_a->triVert3, _b->triVert2 ) || relativeClose(_a->triVert3, _b->triVert3 ) ) return true ;
	return false;
}

inline void returnFilledElementsOfObject ( QuadTreeNode * _v , vector<QuadTreeNode * > & _a ) {

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

inline void setRadiosityForAllElements ( QuadTreeNode * _qR ) {

	vector<QuadTreeNode * > tempQuadVector ;

	//for every object
	for ( auto &v : _qR->children ) 
		returnFilledElementsOfObject ( v, tempQuadVector ) ;

	int i = 0 ;

	for ( auto &v : tempQuadVector ) 
		v->radiosityValue = B[i++];
}

Radiosity::Radiosity( Scene * _scene, Camera * _cam ) :
	scene ( _scene ), cam ( _cam ),  numOfElements ( 0 ), quadTreeRoot ( NULL ) {

	try {

		radiosityHelper = new Radiosity_Helper( scene ) ;
		vector<Object *>::const_iterator sceneObjIter = scene->sceneObjects.cbegin ();
		tempQuadVector = new vector<QuadTreeNode * >(0);

		//radiosityHelper
		// find AABB for all objects in the scene to check for poly-poly intersection.
		// in which case, we must allow meshing of objects that are not part of the same polygon
		// aka correcting mesh disontinuity, for example if there is a column that is attached to the floor
		// and that they are different polygons
		for ( ; sceneObjIter != scene->sceneObjects.cend() ; ++sceneObjIter ) {

			ObjectNode * tempObjNode = new ObjectNode() ;

			radObjectNodes.push_back ( tempObjNode );
			tempObjNode->object = *sceneObjIter;
			tempObjNode->aabbInfo = GetBoxPolygon ( *sceneObjIter );
		}
		findIntersections () ;

		//Now we shall build the initial quadtree
		buildInitialQuadTree() ;

		// find form factor between all elements in the scene
		findFormFactor() ;

		// find the solution to the linear equation and remesh if needed
		progressiveRefinement () ;

		//Camera camera;
		setRadiosityForAllElements ( this->quadTreeRoot ) ;

		scene->rasterize->Radiosity_Raster ( "Radiosity" ,*cam, this, radiosityHelper ) ; 
	}
	catch ( std::exception ex ) {
		cout << "Error in Radiosity::Radiosity( Scene * _scene ) : "  << ex.what() << endl;
	}
}

Radiosity::~Radiosity () {

	int n = numOfElements + scene->NumLights() ;

	for ( auto &v : radObjectNodes ) delete v ;

	for ( int i = 0 ; i < n ; ++ i )
		free ( FF[i] ) ;
	free ( FF ) ;

	for ( int i = 0 ; i < n ; ++ i )
		free ( K[i] ) ;
	free ( K ) ;

	free ( B );

	// TODO : post order deletion of quadtree recursively
}

// build the initial structure of the quadtree
void Radiosity::buildInitialQuadTree () {

	if ( quadTreeRoot == NULL ) quadTreeRoot = new QuadTreeNode () ;

	// for all the objects in the scene
	for ( auto &v : radObjectNodes ) {
		quadTreeRoot->children.push_back(new QuadTreeNode() );
		quadTreeRoot->children.back()->object = v->object ;
		v->qNode = quadTreeRoot->children.back();

		// add all the triangles in the BSP of an object to a vector to check for adjacent triangles
		vector<QuadTreeNode * > tempQuadVector ;

		BSP_Node * tempBSPNode = v->object->root ;
			
		//tempQuadVector.emplace_back( dynamic_cast<QuadTreeNode *>(tempBSPNode) ) ;
		tempQuadVector.emplace_back( new QuadTreeNode ( tempBSPNode)  ) ;

		numOfElements += v->object->triangles.size() ;

		// add all the triangles into a vector 
		for ( int i = 0 ; i < v->object->triangles.size() ; ++ i ) {
			if ( tempQuadVector[i]->left != NULL )tempQuadVector.emplace_back ( new QuadTreeNode ( tempQuadVector[i]->left )  ) ;
			if ( tempQuadVector[i]->right != NULL ) tempQuadVector.emplace_back ( new QuadTreeNode ( tempQuadVector[i]->right ) ) ; 
		}

		QuadTreeNode * currentQuadObject = quadTreeRoot->children.back();

		// add triangles as new children of object
		for ( int i = 0 ; i < tempQuadVector.size() ; ++ i ) 
			currentQuadObject->children.emplace_back (tempQuadVector[i] ) ;

		currentQuadObject = quadTreeRoot->children.back();
		// make the quad tree nodes and check for adjacency with other triangles, checking the similarity of vertices within the same polygon
		for ( int i = 0 ; i < currentQuadObject->children.size() ; ++ i ) {
			for ( int j = 0 ; j < currentQuadObject->children.size() ; ++ j ) {
				if ( i != j ) {
					radiosityHelper->detectTriangleIntersections ( currentQuadObject->children[i], currentQuadObject->children[j] );
				}
			}
		}
	}

	// intersection detection and re-triangulation of overlapping/colliding triangles
	for ( int i = 0 ; i < quadTreeRoot->children.size() ; i++ ) {
		if ( radObjectNodes[i]->aabbDoesIntersect ) {
			ObjectNode::LinkListAdjNode * curr = radObjectNodes[i]->AdjLL ;
			while ( curr != NULL ) {
				QuadTreeNode * otherNode = curr->oN->qNode ;
				if ( curr->oN->finishedWithIntersections == false ) {
					for ( int q1 = 0 ; q1 < radObjectNodes[i]->qNode->children.size() ; ++ q1 ) {
						for ( int q2 = 0 ; q2 < otherNode->children.size() ; ++ q2 ) {
							radiosityHelper->detectTriangleIntersections ( radObjectNodes[i]->qNode->children[q1] , otherNode->children[q2] ) ;
						}
					}
				}
				curr = curr->next ;
			}
			radObjectNodes[i]->finishedWithIntersections = true ;
		}
	} 
}

bool Radiosity::castFromElementToElement ( const Ray &ray, HitInfo &hitinfo, QuadTreeNode * i, QuadTreeNode * j, const Vec3 &xi, const Vec3 &yi  ) {

	double cmprDistance = 0.0 ;
	HitInfo tempInfo ;

	for ( std::vector<Object * >::const_iterator objsBegin = scene->sceneObjects.begin() ; 
		  objsBegin != scene->sceneObjects.end() && ( *objsBegin != NULL || *objsBegin != hitinfo.ignore ) ; 
		  objsBegin++ ){

				if( (*objsBegin)->Intersect( ray, hitinfo ) ){
					if ( cmprDistance == 0 ){
						cmprDistance = hitinfo.distance ;
						hitinfo.object = (*objsBegin) ;
						tempInfo = hitinfo ;
					}
					else if ( cmprDistance > hitinfo.distance ){
				
						cmprDistance = hitinfo.distance ;
						hitinfo.object = (*objsBegin);
						tempInfo = hitinfo ;
					}
					hitinfo.ray = ray; // Save the ray in world coordinates.
				}
	}
	hitinfo = tempInfo ;

	return ((hitinfo.object != 0) ?  true :  false);
}

bool Radiosity::isVisibleXiToXj ( QuadTreeNode * i, QuadTreeNode * j, const Vec3 &xi, const Vec3 &yi ) {

	try {
		Ray ray;
		ray.direction = Unit (yi - xi ) ;
		ray.origin     = xi + ( Epsilon * ray.direction ) ; // cast from above the element surface
		ray.type       = generic_ray; // These rays are given no special meaning.
		ray.generation = 1;           // Rays cast from the element are first-generation.

		HitInfo hitinfo;             
		hitinfo.ignore = NULL;       // Don't ignore any objects.
		hitinfo.distance = Infinity; // Follow the full ray.

		if ( castFromElementToElement ( ray, hitinfo, i, j, xi, yi ) ){
			if( hitinfo.object == NULL ){ /* cout << "Error : hitinfo in isVisibleXitoXj  did not return an collision with object" << endl ;*/ }

			else {

				return ((hitinfo.object != 0 &&
				 //hitinfo.object == j->object &&
				 Length(hitinfo.point - yi) < Epsilon &&
				 ray.direction * j->triNormal < 0  ) ?  true :  false);
			}
		}
		//else cout << "isVisibleXitoXj detected no collisions between elements" << endl ; 
	}
	catch ( std::exception ex ) {
		cout << "Error in Radiosity::isVisibleXiToXj " << ex.what() << endl ;
		
	}
	return false;
}

double Radiosity::findFormFactorTermWithLight ( QuadTreeNode * i ) {
	
	double formFactorij ;
	Vec3 xi;
	int formFactorLoopIterations = 500 ;

	const PrimitiveObject *light = scene->GetLight(0);
	AABB box = GetBox( *light );
	Vec3 centerEx = Center(box);

	for ( int j = 0 ; j < formFactorLoopIterations ; ++ j ) {

		formFactorij = 0.0 ;
		xi = returnRandomPointOnElement ( i ) ;

		Ray ray;
		ray.direction  = Unit (centerEx - xi);
		ray.origin     = xi + ( Epsilon * i->triNormal ) ;
		ray.type       = generic_ray; // These rays are given no special meaning.
		ray.generation = 1;           // Rays cast from the element are first-generation.

		HitInfo hitinfo;             
		hitinfo.ignore = NULL;       // Don't ignore any objects.
		hitinfo.distance = Infinity; // Follow the full ray.

		// if some element is hit on the way from element i to a light, then clearly element i can't see the light
		if ( !castFromElementToElement ( ray, hitinfo, i, NULL, xi, centerEx) ) {
			//double tempFormFactorij = (( 1.0 / ( 4.0 * Pi ) ) * (( ray.direction * i->triNormal ) / (pow(Length( centerEx - xi ), 2 )))) / 
				//Area ( i->triVert1, i->triVert2, i->triVert3 ) ;
			double tempFormFactorij = (( 1.0 / ( 4.0 * Pi ) ) * (( ray.direction * i->triNormal ) / (pow(Length( centerEx - xi ), 2 ))));
				
			if (tempFormFactorij > 0.0 )
				formFactorij += tempFormFactorij ;

			//return ((formFactorij > 0.0 ) ? formFactorij : 0.0) ;
		}

	}
	return (formFactorij / ( (double ) formFactorLoopIterations ))   ;
}

double Radiosity::findFormFactorTerm ( QuadTreeNode * i , QuadTreeNode * j ) {
	
	double formFactorij = 0 ;
	int formFactorLoopIterations = 500 ;

	//shoot formFactorLoopIterations # of rays from element i to element j
	for ( int k = 0 ; k < formFactorLoopIterations ; ++k ) {
		Vec3 xi = returnRandomPointOnElement ( i ) ;
		Vec3 xj = returnRandomPointOnElement ( j ) ;

		if ( isVisibleXiToXj ( i, j, xi, xj ) ) {
			double distanceSquared = ( xi - xj ) * ( xi - xj ) ;
			Vec3 cos1 = Unit ( xj - xi ) ;
			double cosF1 = cos1 * i->triNormal ;
			Vec3 cos2 = Unit ( xi - xj ) ;
			double cosF2 = cos2 * j->triNormal ;
			double cosFactor =  (Unit( xj - xi ) * i->triNormal) * ( Unit ( xi - xj ) * j->triNormal ) ; 
			double delF = (((Unit( xj - xi ) * i->triNormal) < 0.0 || ( Unit ( xi - xj ) * j->triNormal ) < 0.0 ) ?  
				0.0 : 
				cosFactor / (( Pi * distanceSquared ) + (Area(j->triVert1, j->triVert2, j->triVert3 ) / formFactorLoopIterations))) ; 
			if ( delF > 0.0 ) {
				formFactorij += delF ;
				//cout << xi << " : " << xj << endl ;
			}
		}
	}
	formFactorij = formFactorij / ((double)formFactorLoopIterations) ;
	return formFactorij * Area(j->triVert1, j->triVert2, j->triVert3 ) ;
}

void Radiosity::findFormFactor() {

	// Size of K matrix.
	int n = numOfElements + scene->NumLights() ;

	// K matrix = M - P * FF
	K = ( double ** ) malloc ( n * sizeof ( double * )) ;

	for ( int i = 0 ; i < n ; ++ i )
		K[i] = ( double * ) malloc ( n * sizeof ( double ));

	// form factor matrix
	FF = ( double ** ) malloc ( n * sizeof ( double * )) ;

	for ( int i = 0 ; i < n ; ++ i )
		FF[i] = ( double * ) malloc ( n * sizeof ( double ));

	// identity matrix
	double ** I = ( double ** ) malloc ( n * sizeof ( double * )) ;

	for ( int i = 0 ; i < n ; ++ i )
		I[i] = ( double * ) malloc ( n * sizeof ( double )); 

	// Find the form factor matrix.
	// Iterate through all the triangles and find the FF with all other triangles in the scene.
	// Use the structure of the QuadTree to iterate through all the triangles.
	vector<QuadTreeNode * > tempQuadVector ;

	//for every object
	for ( auto &v : quadTreeRoot->children ) {
		returnFilledElementsOfObject ( v, tempQuadVector ) ;
	}

	for ( int i = 0 ; i < numOfElements ; ++i ) {
		for (  int j = 0 ; j < numOfElements  ; ++j ) {
			//cout << "j : " << j << endl ;
			if ( i == j ) FF [i][j] = 0.0 ; 
			else FF[i][j] = findFormFactorTerm ( tempQuadVector[i], tempQuadVector[j] ) ;
		}
	}

	// Now find the form factor for the light source.
	// Point light source.

	// Fill row
	for ( int i = 0; i < n ; i ++ ) {
		FF[n-1][i] = 0.0 ;
	}

	// Fill ith col. Find the form factor between the light and all elements.
	// Does the element hit any other elements along it's path?
	for ( int i = 0 ; i < numOfElements ; ++i ) {
		//for ( int j = 0 ; j < numOfElements ; ++ j ) {
			//FF[i][n-1] = findFormFactorTermWithLight ( tempQuadVector[i], tempQuadVector[j] ) ;
			FF[i][n-1] = findFormFactorTermWithLight ( tempQuadVector[i] ) ;
			cout << "FF light value : " << FF[i][n-1] << endl;
		//}	
	}

	// Fill I matrix
	for ( int i = 0 ; i < n ; ++ i ) {
		for ( int j = 0 ; j < n ; ++ j ) {
			if ( i == j ) I[i][j] = 1.0 ;
			else I[i][j] = 0.0 ;
		}
	}

	// calculate K matrix = M - P ( FF ) // p is a scalar
	for ( int i = 0 ; i < numOfElements; ++i ){
		for ( int j = 0 ; j < numOfElements ; ++j ) {
			K[i][j] = I[i][j] - ( REFLECTIVITY_INDEX *  FF[i][j] ) ;
		}
	}

	// Calcuate last col ( light col ) minus the last row
	for ( int i = 0 ; i < numOfElements ; ++ i ) {
		K [i][n-1] = I[i][n-1] - ( REFLECTIVITY_INDEX * FF[i][n-1] ) ; 
	}

	// Calculate last row ( light row )
	for ( int i = 0 ; i < n ; ++ i ) {
		K[n-1][i] = I[n-1][i] - FF[n-1][i] ;
	}

	string outFileString = "FF.txt" ;
	std::ofstream myfile;
	myfile.open (outFileString.c_str());
	string _1 ;
	//logging FF
	for ( int i = 0 ; i < n ; ++i ) {
		_1 = "";
		for ( int j = 0 ; j < n ; ++j ) {
			_1 += "[" + std::to_string(i) + "][" + std::to_string(j) + "]" + ":" + std::to_string(FF[i][j]) + "  " ;
		}
		_1 += "\n";
		myfile << _1 ;
	}
	myfile.close();

	outFileString = "K.txt" ;
	myfile.open (outFileString.c_str());
	_1 ;
	//logging FF
	for ( int i = 0 ; i < n ; ++i ) {
		_1 = "";
		for ( int j = 0 ; j < n ; ++j ) {
			_1 += "[" + std::to_string(i) + "][" + std::to_string(j) + "]" + ":" + std::to_string(K[i][j]) + "  " ;
		}
		_1 += "\n";
		myfile << _1 ;
	}
	myfile.close();

    // delete I
	for (int i = 0 ; i < n ; ++ i ) 
		//delete[] I[i] ;
			free ( I[i] );

	//delete [] I ;
	free ( I );

}

bool Radiosity::checkIfConverged ( double * _E ) {

	int n = numOfElements + scene->NumLights() ;

	for ( int i = 0 ; i < n ; i ++ ) {
		double sum = 0.0 ;
		for ( int j = 0 ; j < n ; j ++ ) {
			sum += K[i][j] * B[j] ;
		}
		if ( abs (sum - _E[i] ) > Epsilon ) return false ; 
	}
	return true ;
}

void Radiosity::progressiveRefinement() {

	int n = numOfElements + scene->NumLights() ;
	vector<QuadTreeNode * > tempQuadVector ;

	//for every object
	for ( auto &v : quadTreeRoot->children ) {
		returnFilledElementsOfObject ( v, tempQuadVector ) ;
	}

	// emissitivity
	double * E = ( double * ) malloc ( n * sizeof ( double ) ) ;

	for ( int i = 0 ; i < n - 1 ; i ++ )
		//E[i] = SOLAR_RADIANT_FLUX ;
			E[i] = 0.0 ;

	E[n-1] = SOLAR_RADIOSITY_POWER  ;

	// radiosities 
	B = ( double * ) malloc ( n * sizeof ( double ) ) ;

	double * dB = ( double * ) malloc ( n * sizeof ( double )) ;

	for ( int i = 0 ; i < n ; ++ i ) {
		dB[i] = E[i] ;
		B[i]  = E[i] ;
	}

	double * tempdBArea = ( double * ) malloc ( n * sizeof ( double ) ) ;

	double dRad = 0.0 ;

	bool firstShot = true ;

	while ( !checkIfConverged( E ) ){

		int indx ;
		if ( firstShot ) {
			firstShot = !firstShot ;
			for ( int i = 0 ; i < n - 1 ; ++ i ) {
				dRad = dB[n-1] * FF[i][n-1] * REFLECTIVITY_INDEX ;
				dB[i] += dRad;
				B[i] += dRad ;
			}
			indx = n - 1 ;
		}
		else {
			for ( int i = 0 ; i < n - 1 ; i ++ )
				tempdBArea[i] = dB[i] * Area ( tempQuadVector[i]->triVert1, tempQuadVector[i]->triVert2, tempQuadVector[i]->triVert3 ) ;

			indx = returnHighestValueIndx (tempdBArea, n ) ;

			for ( int i = 0 ; i < n ; ++ i ) {
				dRad = dB[indx] * FF[i][indx] * REFLECTIVITY_INDEX ;
				dB[i] += dRad ;
				B[i] += dRad ;
			}
		}
		dB[indx] = 0.0 ;
	}

	// make an output B solutions file
	string outFileString = "B.txt" ;
	std::ofstream myfile;
	myfile.open (outFileString.c_str());
	string _1 ;
	//logging FF
	for ( int i = 0 ; i < numOfElements + scene->NumLights() ; ++i ) {
			
		_1 = "[" + std::to_string(i) + "] :" + std::to_string(B[i]) + "  " ;
		_1 += "\n";
			
		myfile << _1 ;
	}
		
	myfile.close();

	//delete [] E ;
	//delete [] dB ;
	//delete [] tempdBArea ;
	free ( E );
	free ( dB);
	free (tempdBArea ) ;
}

/*
bool Radiosity::adapdtiveMeshSubDivision() {
}
*/