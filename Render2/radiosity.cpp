#include "ray_tracer.h"
#include "util.h"
#include <random>
#include <stack>

struct ObjectNode {

	ObjectNode() {
		object = NULL ;
		intersectsWith = NULL;
		aabbDoesIntersect = false ;
	}

	Object * object;
	Object * intersectsWith ;
	AABB aabbInfo ;
	bool aabbDoesIntersect ;
};

// create initial quadtree
// must sustain the adjacency list of all triangles with other triangles
struct QuadTreeNode : public BSP_Node {
	
	QuadTreeNode () {
		object = NULL;
		parentObject = NULL ;
		visited = false;
		delVisited = false;
	}

	QuadTreeNode ( BSP_Node * _BSP_Node ) {

		object = NULL;
		visited = false;
		delVisited = false;

		triVert1 = _BSP_Node->triVert1;
		triVert2 = _BSP_Node->triVert2;
		triVert3 = _BSP_Node->triVert3;
		triIndx1 = _BSP_Node->triIndx1;
		triIndx2 = _BSP_Node->triIndx2;
		triIndx3 = _BSP_Node->triIndx3;
		vertNormal1 = _BSP_Node->vertNormal1;
		vertNormal2 = _BSP_Node->vertNormal2;
		vertNormal3 = _BSP_Node->vertNormal3;
		triNormal = _BSP_Node->triNormal;
		textVert = _BSP_Node->textVert;
		parentObject = _BSP_Node->parentObject;
	}

	vector<QuadTreeNode *> children ;
	// maintain QuadTreeNode vector of all adjacent nodes
	vector<QuadTreeNode *> nextAdj;

	Object * object;
	Object * parentObject;

	double radiosityValue;
	double reflectivityIndex;
	double emissitivity;
	double formFactor ;

	bool visited;
	bool delVisited;
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
QuadTreeNode * quadTreeRoot = NULL ;
// form factor matrix
double ** FF; 
// K matrix K = M - PFF where M = Identity, P = Reflectiviy and FF = form factor
double ** K ; 

// find all the AABB intersections of every object with each other
inline void findIntersections () {

	for ( auto &v : radObjectNodes ) {
		for ( auto &o : radObjectNodes ) {
			if (v != o && v->intersectsWith != o->object  ) {
				if  ( AABBIntersect ( v->aabbInfo, o->aabbInfo ) ) {
					v->aabbDoesIntersect = true;
					v->intersectsWith = o->object ;
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


Radiosity::Radiosity( Scene * _scene ) :
	scene ( _scene ), numOfElements ( 0 ) {

	try {

		vector<Object *>::const_iterator sceneObjIter = scene->sceneObjects.cbegin ();

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
	}
	catch ( std::exception ex ) {
		cout << "Error in Radiosity::Radiosity( Scene * _scene ) : "  << ex.what() << endl;
	}
}

Radiosity::~Radiosity () {

	for ( auto &v : radObjectNodes ) delete v ;

	// TODO : post order deletion of quadtree recursively


}

// build the initial structure of the quadtree
void Radiosity::buildInitialQuadTree () {

	if ( quadTreeRoot == NULL ) quadTreeRoot = new QuadTreeNode () ;

	// for all the objects in the scene
	for ( auto &v : radObjectNodes ) {
		quadTreeRoot->children.push_back(new QuadTreeNode() );
		quadTreeRoot->children.back()->object = v->object ;

		// add all the triangles in the BSP of an object to a vector to check for adjacent triangles
		vector<QuadTreeNode * > tempQuadVector ;

		BSP_Node * tempBSPNode = v->object->root ;
			
		tempQuadVector.emplace_back( tempBSPNode ) ;

		numOfElements += v->object->triangles.size() ;

		// add all the triangles into a vector 
		for ( int i = 0 ; i < v->object->triangles.size() ; ++ i ) {
			if ( tempQuadVector[i]->left != NULL ) tempQuadVector.emplace_back ( tempQuadVector[i]->left ) ;
			if ( tempQuadVector[i]->right != NULL ) tempQuadVector.emplace_back ( tempQuadVector[i]->right ) ; 
		}

		QuadTreeNode * currentQuadObject = quadTreeRoot->children.back();

		// add triangles as new children of object
		for ( int i = 0 ; i < tempQuadVector.size() ; ++ i ) 
			currentQuadObject->children.emplace_back (new QuadTreeNode(tempQuadVector[i]) ) ;

		currentQuadObject = quadTreeRoot->children.back();
		// make the quad tree nodes and check for adjacency with other triangles
		for ( int i = 0 ; i < currentQuadObject->children.size() ; ++ i ) {
			for ( int j = 0 ; j < currentQuadObject->children.size() ; ++ j ) {
				if ( i != j ) {
					if ( findTriangleAdj ( currentQuadObject->children[i], currentQuadObject->children[j] ) ) 
						currentQuadObject[i].nextAdj.emplace_back ( currentQuadObject->children[j] ) ;
				}
			}
		}
	}

	// TODO : finish later
	// problem : the two intersecting objects might not share the same triangle endpoints, find other way of intersection detection
	// now if any of the objects' AABB intersect each other, need to join their triangles
	for ( int i = 0 ; i < quadTreeRoot->children.size() ; i++ ) {
		if ( radObjectNodes[i]->aabbDoesIntersect ) {
			for ( int j = 0 ; j < quadTreeRoot->children.size(); j++ ) {
				if ( i != j && radObjectNodes[i]->intersectsWith == quadTreeRoot->children[j]->object ) {
				}
			}
		}
	}
}


inline void returnFilledElementsOfObject ( QuadTreeNode * _v , vector<QuadTreeNode * > & _a ) {

	std::stack<FormFactorStackNode> stackFormFactor ;

	//for all the initial elements of the object node
	for ( auto &elements : _v->children ) {

		for ( vector<QuadTreeNode *>::const_iterator iterQuads = v->children.cend() - 1 ; iterQuads >= v->children.cbegin() ; --iterQuads )
			stackFormFactor.push ( FormFactorStackNode(*iterQuads));
	}

	while ( !stackFormFactor.empty() ) {
		FormFactorStackNode temp = stackFormFactor.top();
		stackFormFactor.pop();
		if ( temp.quadTreeNode->children.size() > 0 ) 
			for ( vector<QuadTreeNode *>::iterator iterQuads = temp.quadTreeNode->children.end - 1 ;
				iterQuads >= temp.quadTreeNode->children.begin() ; 
				--iterQuads ) 

					stackFormFactor.push ( FormFactorStackNode(*iterQuads));
		else {
			_a.emplace_back ( stackFormFactor.top().quadTreeNode ) ;
			stackFormFactor.pop() ;
		}
	}
}

bool Radiosity::castFromElementToElement ( const Ray &ray, HitInfo &hitinfo, QuadTreeNode * i, QuadTreeNode * j, const Vec3 &xi, const Vec3 &yi  ) {

	double cmprDistance = 0.0 ;
	HitInfo tempInfo ;

	for ( std::vector<Object * >::const_iterator objsBegin = scene->sceneObjects.begin ; 
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

	return ((hitinfo.object != 0 &&
		     hitinfo.object == j->object &&
		     Length(hitinfo.point - yi) < Epsilon &&
			 ray.direction * j->triNormal < 0  ) ?  true :  false);
}

bool Radiosity::isVisibleXiToXj ( QuadTreeNode * i, QuadTreeNode * j, const Vec3 &xi, const Vec3 &yi ) {

	Ray ray;
	ray.direction = Unit (yi - xi ) ;
	ray.origin     = xi + ( Epsilon * ray.direction ) ; // cast from above the element surface
	ray.type       = generic_ray; // These rays are given no special meaning.
	ray.generation = 1;           // Rays cast from the element are first-generation.

	HitInfo hitinfo;             
	hitinfo.ignore = NULL;       // Don't ignore any objects.
	hitinfo.distance = Infinity; // Follow the full ray.

	if ( castFromElementToElement ( ray, hitinfo ) ){
		if( hitinfo.object == NULL ) cout << "Error : hitinfo in isVisibleXitoXj  did not return an collision with object" << endl ; 

		else {

		}
	}
	else cout << "isVisibleXitoXj detected no collisions between elements" << endl ; 
}

inline Vec3 returnRandomPointOnElement( QuadTreeNode * i ) {
	
	try {
		std::random_device rd;
		std::mt19937 gen(rd());

		double u = 0.0, v = 0.0 ;

		do {
			std::uniform_real_distribution<> dis (0.0, 1.0 ) ;

			u = dis ( gen ) ;
			v = dis ( gen ) ;

		} while ( u + v <= 1.0 ) ;

		Vec3 P ;
		P.x = (( 1.0 - u - v ) * i->triVert1.x) + (u * i->triVert2.x) + (v * i->triVert3.x) ;
		P.y = (( 1.0 - u - v ) * i->triVert1.y) + (u * i->triVert2.y) + (v * i->triVert3.y) ;
		P.z = (( 1.0 - u - v ) * i->triVert1.z) + (u * i->triVert2.z) + (v * i->triVert3.z) ;

		return P ;
	}
	catch ( std::exception ex ) {
		cout <<"Error in inline Vec3 returnRandomPointOnElement : " << ex.what()  << endl ;
	}
}


double Radiosity::findFormFactorTerm ( QuadTreeNode * i , QuadTreeNode * j ) {
	
	double formFactorij = 0 ;

	//shoot 5 rays from element i to element j
	for ( int k = 0 ; k < 5 ; ++k ) {
		Vec3 xi = returnRandomPointOnElement ( i ) ;
		Vec3 xj = returnRandomPointOnElement ( j ) ;

		bool isVisibleXitoXj ( i, j, xi, xj ) ;
	}
}

void Radiosity::findFormFactor() {

	// Size of K matrix.
	int n = numOfElements + scene->NumLights() ;

	K = ( double ** ) malloc ( n * sizeof ( double * )) ;

	for ( int i = 0 ; i < n ; ++ i )
		K[i] = ( double * ) malloc ( n * sizeof ( double ));

	FF = ( double ** ) malloc ( n * sizeof ( double * )) ;

	for ( int i = 0 ; i < n ; ++ i )
		FF[i] = ( double * ) malloc ( n * sizeof ( double ));

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
			if ( i == j ) FF [i][j] = 1.0 ; 
			else FF[i][j] = findFormFactorTerm ( tempQuadVector[i], tempQuadVector[j] ) ;
		}
	}
}
