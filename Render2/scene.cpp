#include <list>
#include "ray_tracer.h"
#include "util.h"
#include <iostream>
#include <fstream> 
#include <random>

//Object * object ;

using std::string ;

static std::list< PrimitiveObject * > * all_light_types = NULL;  // All registered plug-ins.

bool RegisterPlugin( PrimitiveObject *pObj ) {
	if( all_light_types == NULL ) all_light_types = new std::list< PrimitiveObject* >(); 
	all_light_types->push_back( pObj );
	return true;
}

static const Color
	DEFAULT_BACKGROUND_COLOR  = Color( 0.15, 0.25, 0.35 ),
	DEFAULT_TERMINATION_COLOR = Color( 0.0, 0.0, 0.0 );

static Material * Copy( Material * mat, const Material &material ){
	if( mat == NULL || !(*mat == material) ) mat = new Material( material );
	return mat;
}

static bool Skip( char *line )
{
	for( char c = *line; c != '\0' && c != '\n' && c != '#'; c = *++line )
	{
		if( c != ' ' && c != '\t' ) return false;
	}
	return true;
}

// Make the DoT file to visualize the tree
void makeDotFile ( BSP_Node * v, string &nodeString, string &connectionString ) {

	// Only called once at the beginning of the function. Output the normal of the face.
	if ( nodeString == "" ) {
		//std::ofstream myfile;
		//myfile.open ("outGraph.dot");
		nodeString = "digraph g { graph []; node [ fontsize=\"16\" shape = \"ellipse\"]; edge [];" ;
		nodeString += "\"node" + std::to_string(v->node) + "\" [ label = \"" + std::to_string(v->node) + " | " + 
								 std::to_string(v->triNormal.x) + " " + 
								 std::to_string(v->triNormal.y) + " " +
								 std::to_string(v->triNormal.z) + " " + "\" shape = \"record\" ];";
		if ( v->left != NULL ){
			connectionString += "\"node" + std::to_string(v->node) + "\"-> \"node" + std::to_string(v->left->node) + "\" [ label = \"left\"];";
			makeDotFile ( v->left, nodeString, connectionString );
		}
		if ( v->right != NULL ){
			connectionString += "\"node" + std::to_string(v->node) + "\"-> \"node" + std::to_string(v->right->node) + "\" [ label = \"right\"];";
			makeDotFile ( v->right, nodeString, connectionString);
		}

		connectionString += "}";
	
		//myfile << nodeString;
		//myfile << connectionString;
		//myfile.close();
	}
	else{

		nodeString += "\"node" + std::to_string(v->node) + "\" [ label = \"" + std::to_string(v->node) + " | " +
								 std::to_string(v->triNormal.x) + " " + 
								 std::to_string(v->triNormal.y) + " " +
								 std::to_string(v->triNormal.z) + " " + "\" shape = \"record\" ];";

		if ( v->left != NULL ) {
			connectionString += "\"node" + std::to_string(v->node) + "\"-> \"node" + std::to_string(v->left->node) + "\" [ label = \"left\"];";
			makeDotFile ( v->left, nodeString, connectionString );
		}
		if ( v->right != NULL ) {
			connectionString += "\"node" + std::to_string(v->node) + "\"-> \"node" + std::to_string(v->right->node) + "\" [ label = \"right\"];";
			makeDotFile ( v->right, nodeString, connectionString);
		}
	}
}

Scene::Scene() {
	rasterize = NULL;
	radiosity = NULL ;
	max_tree_depth = default_max_tree_depth;

	shader = new Shader();
}

/*
//Will return true if the point that is passed is already inside some object
bool Scene::IsAlreadyInside ( const Vec3 &p, HitInfo &hitinfo ) const
{
	//TODO : Figure this out
	if( object == NULL || object == hitinfo.ignore ) return false;

	if( object->Inside ( p))
	{
	   
		return true;
	}
	return false;

}
*/

bool Scene::calculateVertexNormals ( Object & obj, vector<vector<Vec3>> &  adjMatrix ) {

	if ( adjMatrix.size() == 0 ) return false;

	for ( std::vector<vector<Vec3>>::iterator adjHashIterator = adjMatrix.begin() ; adjHashIterator != adjMatrix.end() ; ++ adjHashIterator ) {
		Vec3 tempVec ( 0.0, 0.0, 0.0);
		for  ( int Vec3InHash = 0 ; Vec3InHash < (*adjHashIterator).size() ; ++Vec3InHash ) 
			tempVec += (*adjHashIterator)[Vec3InHash];
	
		(*adjHashIterator)[0] = Unit(tempVec);
	}

	for ( std::vector<BSP_Node>::iterator triangleIterator = obj.triangles.begin() ; triangleIterator != obj.triangles.end() ; ++triangleIterator ) {
		(*triangleIterator).vertNormal1 = adjMatrix[(*triangleIterator).triIndx1][0];
		(*triangleIterator).vertNormal2 = adjMatrix[(*triangleIterator).triIndx2][0];
		(*triangleIterator).vertNormal3 = adjMatrix[(*triangleIterator).triIndx3][0];
	}
	return true;
}

// Assigns a color to a ray. Cast a ray and attempt to hit any of the objects. If a collision is made, shade the scene.
// Must be recursive to generate generations of rays. 
// End the generation with a generation limit.
// This is called from rasterizer.
Color Scene::Trace( const Ray &ray ) const {

	Color   color;               
	HitInfo hitinfo;             
	hitinfo.ignore = NULL;       // Don't ignore any objects.
	hitinfo.distance = Infinity; // Follow the full ray.

	if( ray.generation > max_tree_depth ){
		// Reached the end of the generations
		color = DEFAULT_TERMINATION_COLOR;
	}
	else if( Cast( ray, hitinfo ) ) {
		if( hitinfo.object == NULL ) return Green;
		
		if( this->shader != NULL )
			 color = this->shader->Shade( *this, hitinfo ); // Shade the scene with the scene's shader and parameters. TODO : add shading params
		else 
			cout << "Scene::Trace : there is no shader defined" << endl ; 
	}
	else {
		// The ray has failed to collide with an object.  Use the environment map
		// associated with the scene to determine the color.
		// If no envmap was specified, use the default background color.
		 color = DEFAULT_BACKGROUND_COLOR;
	}
	
	return color;
}

// Finds the closest point of intersection among all the objects.
// The information about the collision is stored in hitinfo.
bool Scene::Cast( const Ray &ray, HitInfo &hitinfo ) const {
	double cmprDistance = 0 ;
	HitInfo tempInfo ;

	// Iterate through all the objects, lights are not considered objects at the moment.
	for ( std::vector<Object * >::const_iterator objsBegin = sceneObjects.begin() ; 
				objsBegin != sceneObjects.end() && ( *objsBegin != NULL || *objsBegin != hitinfo.ignore ) ; 
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

bool Scene::BuildScene ( string fileName, string fileObj, Camera &camera ) {
	int line_num = 0;
	char input_line[512];
	//Plugin    *dat = NULL;  // A container for arbitrary data.
	PrimitiveObject    *obj = NULL;  // The current object, which was the last object created.
	//Shader    *shd = NULL;  // The current shader. 
	//Aggregate *agg = NULL;  // Current aggregate object, to which all objects are now added.
	//Envmap    *env = NULL;  // Current environment map.
	Material  *mat = NULL;  // Current material pointer.
	Material   material;    // Current material.
	Transformation tempTrnsf;
	string	   tempObjectString;
	Object * object;

	//scene.object    = NULL;
	
	// Attempt to open the input file.
	fileName += ".sdf";
	std::ifstream fin;

	fin.open( fileName.c_str() );
	if( fin.fail() ){
		cerr << "Error: Could not open file " << fileName << endl;
		return false;  // Report failure.
	}
	cout << "Reading " << fileName << "... " << endl;
	cout.flush();

	// Set some defaults.
	material.diffuse      = White;
	material.emission     = Black;
	material.specular     = White;
	material.ambient      = Black;
	material.reflectivity = Black;
	material.translucency = Black;
	material.ref_index    = 0.0;
	material.Phong_exp    = 0.0;
	material.microFacetsRoughness = 0.0;

	bool objectFlipYZFlip = false;
	bool objectScale	  = false;
	bool objectRotate	  = false;
	bool objectTranslate  = false;

	camera.x_res = default_image_width;
	camera.y_res = default_image_height;
	camera.x_win = Interval( -1.0, 1.0 );
	camera.y_win = Interval( -1.0, 1.0 );

	// Process lines until the end of file is reached.
	// Print a warning for all lines that are unrecognizable.

	while( fin.getline( input_line, 512 ) ) {
		line_num++;
		if( Skip( input_line ) ) continue;

		if( all_light_types == NULL ) return NULL;
		PrimitiveObject *pObj = NULL;
		std::list< PrimitiveObject* >::iterator iter;
		for( iter = all_light_types->begin(); iter != all_light_types->end(); iter++ ) { // check to see if the input object is a light
			pObj = (*iter)->ReadString( input_line );
			if( pObj != NULL ) break;
		}

		if ( pObj != NULL )  {  // A light has been entered, add it to the list
			obj = (PrimitiveObject*) pObj ;
			obj->material = Copy ( mat, material );
			if ( Emitter ( material ) ) this->lights.push_back  ( obj );
			continue;
		}
		
		// Now look for all the other stuff...  materials, camera, lights, etc.

		ParamReader get( input_line );
 
		if( get["eye"]          && get[camera.eye]            ) continue;
		if( get["lookat"]       && get[camera.lookat]         ) continue;            
		if( get["up"]           && get[camera.up]             ) continue;            
		if( get["vpdist"]       && get[camera.vpdist]         ) continue;                   
		if( get["x_res"]        && get[camera.x_res]          ) continue;
		if( get["y_res"]        && get[camera.y_res]          ) continue;
		if( get["x_win"]        && get[camera.x_win]          ) continue;
		if( get["y_win"]        && get[camera.y_win]          ) continue;
		
		if( get["ambient"]      && get[material.ambient]      ) continue;            
		if( get["diffuse"]      && get[material.diffuse]      ) continue;
		if( get["specular"]     && get[material.specular]     ) continue;
		if( get["emission"]     && get[material.emission]     ) continue;
		if( get["reflectivity"] && get[material.reflectivity] ) continue;
		if( get["translucency"] && get[material.translucency] ) continue;
		if( get["Phong_exp"]    && get[material.Phong_exp]    ) continue;
		if( get["ref_index"]    && get[material.ref_index]    ) continue;
		if( get["roughness"]    && get[material.microFacetsRoughness] ) continue;

		if( get["flipYZ"]		&& get.getString ( tempObjectString ) ) {
			if ( tempObjectString == "0" ||  tempObjectString == "1" )
				objectFlipYZFlip = std::stoi( tempObjectString ) ;
			else cerr << "Error in assigning flipYZ String" << endl;
			continue ;	
		}

		if ( get["translate"]    && get[tempTrnsf.translation] ) {
			objectTranslate = true ;
			continue;
		}
		if ( get["rotation"]    && get[tempTrnsf.rotation] ) {
			objectRotate = true ;
			continue;
		}
		if ( get["scale"]    && get[tempTrnsf.scale] ) {
			objectScale = true ;
			continue;
		}

		if( get["object"]		&& get.getString ( tempObjectString )) { // encountered an object in the scene file, instantiate the object and assign the material

			cout << "Encountered an object named : " << tempObjectString << endl ;

			object = new Object ();
			this->sceneObjects.push_back ( object );
			object->nameOfObject = tempObjectString;
			object->material = Copy ( mat, material );
			object->trnsf = tempTrnsf ;
			object->flipYZ = objectFlipYZFlip ;
			object->trnsf.Istranslation = objectTranslate;
			object->trnsf.Isrotation = objectRotate;
			object->trnsf.Isscale = objectScale;

			// Set defaults.
			material.diffuse      = White;
			material.emission     = Black;
			material.specular     = White;
			material.ambient      = Black;
			material.reflectivity = Black;
			material.translucency = Black;
			material.ref_index    = 0.0;
			material.Phong_exp    = 0.0;
			material.microFacetsRoughness = 0.0 ;
			objectFlipYZFlip	  = false; 
			objectTranslate		  = false;
			objectRotate		  = false;
			objectScale			  = false;
			tempTrnsf.reset() ;

			continue;
		}
		
		// If no object is defined at this point, it's an error.
		if( !get["end"] ){ 
				cerr << "Warning: Unrecognized or ill-formed command at line "
			 << line_num << ": "
			 << input_line << endl;
		}
	}
	cout << "Done with reading Scene File(s)." << endl;

	/************ Object file input ************ */
	for ( std::vector<Object *>::iterator currentIteratorObject = this->sceneObjects.begin() ; currentIteratorObject != this->sceneObjects.end() ; ++currentIteratorObject ) {

		cout << "Reading Object file(s). " << endl ;

		fileObj = DefaultScene + (*currentIteratorObject)->nameOfObject + ".obj" ; //read the name of the object which should be the name of the object file
		std::ifstream finObj;
		finObj.clear();
		finObj.open ( fileObj.c_str() );

		if ( finObj.fail () ){
			cerr << "Error: Could not open file " << fileObj << endl;
			continue;
		}

		cout << "Reading " << fileObj << "... " << endl ;
		cout.flush();

		line_num = 0 ;
		int numOfObj = 0 ;
		int numOfTriangle = 0;

		Vec3 tempVec ;
		Vec2 tempVec2;
		string tempObjectName;
		vector<vector<Vec3>> adjMatrix;

		while ( finObj.getline ( input_line , 512 ) ){
			
			if( Skip( input_line ) ) continue;

			ParamReader get( input_line );

			// Make sure the first line of the object file has the same name as the object name in the scene file.
			if ( line_num == 0 ) {
				if ( get["o"] && get.getString(tempObjectName) && tempObjectName == (*currentIteratorObject)->nameOfObject ) {
					numOfObj ++ ;
					line_num++;
				}
				else break;
			}
			else  {
				// materials file
				if ( get["mtllib"] ){

				}
				else if ( get ["vt"] ){
					if ( !get.getTexture (tempVec2 ))
						cout << "Error retrieving texture" << endl;
					else
						(*currentIteratorObject)->vertexTexture.push_back ( Vec2( tempVec2.x, tempVec2.y ) );
				}
				else if ( get ["vn"] ){
					if ( !get.getPoint (tempVec ))
						cout << "Error retrieving normal" << endl;
					else
						(*currentIteratorObject)->vertexNormals.push_back (Vec3(tempVec.x, tempVec.y, tempVec.z ));
				}	
				else if ( get["v"] ){
					if ( !get.getPoint (tempVec ))
						cout << "Error retrieving point at line " << line_num << endl;
					else	
						((*currentIteratorObject)->flipYZ ) ? (*currentIteratorObject)->points.push_back  ( Vec3(tempVec.x, tempVec.z, tempVec.y )) :
															  (*currentIteratorObject)->points.push_back  ( Vec3(tempVec.x, tempVec.y, tempVec.z ));
					if ( (*currentIteratorObject)->trnsf.Istranslation ){
						((*currentIteratorObject)->points[(*currentIteratorObject)->points.size() - 1]).x += (*currentIteratorObject)->trnsf.translation.x ;
						((*currentIteratorObject)->points[(*currentIteratorObject)->points.size() - 1]).y += (*currentIteratorObject)->trnsf.translation.y ;
						((*currentIteratorObject)->points[(*currentIteratorObject)->points.size() - 1]).z += (*currentIteratorObject)->trnsf.translation.z ;
					}
				}
				else if ( get["f"] ){
					adjMatrix.resize((*currentIteratorObject)->points.size());
					if ( !get.returnCoordOfFace ( *(*currentIteratorObject), numOfTriangle,adjMatrix  )) //TODO assign each node to it's parent object
						cout << "Error retrieving face at line " << line_num << endl;
					numOfTriangle++;
				}
				else
					cout << "Unidentified ID in obj file @ line : " << line_num << " : " << get.returnString() <<  endl;	
				line_num++;
			}
		}
		if ( !calculateVertexNormals( **currentIteratorObject, adjMatrix ) ) {
			cerr << "Error in calculating vertex normals." << endl ;
		}
	}
	if ( line_num > 0 ) {
		cout << "Successfully read Object file(s)." << endl;
		return true;
	}
	else {
		cerr << "Unsuccesfully read Object File." << endl << "Either the object file did not start with the object name or it was empty." << endl;
		return false;
	}
	return true;
}

// Build each objects BSP that occupies the scene
bool Scene::BuildBSP () {

	for ( std::vector< Object * >::iterator sceneObjectIterator = sceneObjects.begin() ; sceneObjectIterator != sceneObjects.end() ; ++sceneObjectIterator ) {
		
		//TODO : Turned off random triangles for Radiosity debugging
		//randomizeTriangles ( (*sceneObjectIterator)->triangles);
		vector < BSP_Node * > tris ;

		for ( unsigned int i = 0 ; i < (*sceneObjectIterator)->triangles.size(); i++ )
			tris.push_back( &(*sceneObjectIterator)->triangles[i] ) ;

		(*sceneObjectIterator)->root = buildBSPTree ( tris );
		
		string _1 = "", _2 = "" ;

		string outFileString = (*sceneObjectIterator)->nameOfObject + "outGraph.dot" ;
		std::ofstream myfile;
		myfile.open (outFileString.c_str());

		//makeDotFile ( (*sceneObjectIterator)->root , _1, _2 );

		myfile << _1;
		myfile << _2;
		myfile.close();
	}	
	return true;
}

//Fisher–Yates shuffle
void  Scene::randomizeTriangles ( vector <BSP_Node> & triangles) {

	//srand ((unsigned int ) time(NULL));
	try {
	std::random_device rd;
	std::mt19937 gen(rd());
	//std::uniform_real_distribution<int> dis(0.0, 2.0);

	for ( int i = (int) triangles.size()  ; i >=2 ; i-- ){
		//int k = rand () % i + 1 ;
		std::uniform_real_distribution<> dis(1, i);
		int k = dis ( gen ) ;

		if ( i != k ){
			BSP_Node temp; 

			temp = triangles[ i - 1 ];
			triangles[ i - 1 ] = triangles[ k - 1 ];
			triangles[ k - 1 ] = temp;
		} 
		} }
	catch ( std::exception ex ) {
		cout << ex.what()  << endl ;
		}
}

// Build the Binary space partition tree recursivly using auto partitioning.
// If a triangle's vertice are between 0 and 3 then it is on the left side of the Hyperplane.
// If a triangle's vertices are between 0 and -3 then it is on the right side of the Hyperplane.
// Otherwise the triangle is coplaner with the Hyperplane. In which case, randomly add to the right or left side.
// Base case : the left or right items in the list is either 0 or 1
BSP_Node * Scene::buildBSPTree ( vector<BSP_Node *>  & v ) {
	if ( v.size() == 0 )		
		return NULL;
	if ( v.size() <= 1){ //leaf node
		return v[0] ; //return the node that is by itself.
	}
	else{
		vector<BSP_Node *> vRight ;
		vector<BSP_Node *> vLeft;
		for ( unsigned int i = 1  ; i < v.size() ; i ++ ){
			
			int result = sideTest3d(v[0]->triVert1, v[0]->triVert2, v[0]->triVert3, v[i]->triVert1 );
			result += sideTest3d(v[0]->triVert1, v[0]->triVert2, v[0]->triVert3, v[i]->triVert2 );
			result += sideTest3d(v[0]->triVert1, v[0]->triVert2, v[0]->triVert3, v[i]->triVert3 );

			if ( result == 0 ){
				if ( rand() % 2 == 0 ) vLeft.push_back( v[i] ); //the point is on the plane, randomly choose a side to place the triangle. It will eventually be visited in the painters algorithm.
				else vRight.push_back( v[i] );
			}
			else if ( result <= 3 && result > 0 ){ //TODO : CHECK FOR TRIANGLE SPLIT
				vLeft.push_back( v[i] );
			}
			else if ( result >= -3 && result < 0 ) {  //TODO : CHECK FOR TRIANGLE SPLIT
				vRight.push_back( v[i] );
			}
			else {
				cout << "Error buildBSPTree" << endl ;
				return NULL;
			}
		}

		v[0]->left = buildBSPTree ( vLeft );
		v[0]->right = buildBSPTree ( vRight );

		return v[0];
	}
}