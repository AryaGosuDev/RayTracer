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
struct QuadTreeNode {
	
	QuadTreeNode () {
		object = NULL;
		parentObject = NULL ;
	}

	vector<QuadTreeNode *> children ;

	Object * object;
	Object * parentObject;

	// maintain QuadTreeNode of all adjacent nodes
	QuadTreeNode * nextAdj;

	Vec3 triVert1;
	Vec3 triVert2;
	Vec3 triVert3;
	int triIndx1;
	int triIndx2;
	int triIndx3;
	Vec3 vertNormal1;
	Vec3 vertNormal2;
	Vec3 vertNormal3;
	Vec3 triNormal;

	double radiosityValue;
	double reflectivityIndex;
	double emissitivity;
	double formFactor ;
};

vector<ObjectNode *> radObjectNodes ;
QuadTreeNode * root = NULL ;

// find all the AABB intersections of every object with each other
void findIntersections () {

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





Radiosity::Radiosity( Scene * _scene ) :
	scene ( _scene ) {

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
		buildInitialQuadTree();



	}
	catch ( std::exception ex ) {
		cout << "Error in Radiosity::Radiosity( Scene * _scene ) : "  << ex.what() << endl;
	}
}

Radiosity::~Radiosity () {

	for ( auto &v : radObjectNodes ) delete v ;

	// post order deletion of quadtree


}

bool Radiosity::buildInitialQuadTree () {

	if ( root == NULL ) root = new QuadTreeNode () ;

	//for ( auto &v : radObjectNodes ) {
		//root->children.insert(
	//}

	return 1;


}
