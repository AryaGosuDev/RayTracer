#include "ray_tracer.h"
#include "util.h"

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

void Radiosity::findFormFactor() {

	K = ( double ** ) malloc ( numOfElements * sizeof ( double * )) ;

	for ( int i = 0 ; i < numOfElements ; ++ i )
		K[i] = ( double * ) malloc ( numOfElements * sizeof ( double ));

	FF = ( double ** ) malloc ( numOfElements * sizeof ( double * )) ;

	for ( int i = 0 ; i < numOfElements ; ++ i )
		FF[i] = ( double * ) malloc ( numOfElements * sizeof ( double ));

	// Find the form factor matrix.
	// Iterate through all the triangles and find the FF with all other triangles in the scene.
	// Use the structure of the QuadTree to iterate through all the triangles.









}
