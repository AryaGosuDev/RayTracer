#include "ray_tracer.h"
#include <chrono>

int BSP_Node::numOfTriangle = 1;

int main( int argc, char *argv[] ) {
	Scene  scene;
	Camera camera;
	scene.rasterize = new Rasterizer();
	
	cout << "Tracer Begin" << endl ;

	string fname = DefaultScene;
	string fnameObject = DefaultScene;
	string fnameOutput = DefaultScene;

	if( argc > 1 ){
		fname += argv[1];
		fnameObject += argv[2];
		fnameOutput += argv[3];
	}
	
	//TODO : Each object's BSP should only be attribted to their respective object, not the entire scene
	if ( !scene.BuildScene ( fname , fnameObject , camera ))
		cout << "Build did not work." << endl ;
	
	if ( !scene.BuildBSP ())
		cout << "Could not build BSP." << endl;
	
	//auto start = std::chrono::steady_clock::now();

	scene.radiosity = new Radiosity(&scene, &camera);

	if( !scene.rasterize->Rasterize( fnameOutput, camera, scene)){
		cerr << "Error encountered while rasterizing." << endl;
		return error_rasterizing_image;
	}

	//auto end = std::chrono::steady_clock::now();

	//auto diff = end - start;

	//cout <<  "Shading time performance : " << std::chrono::duration_cast<std::chrono::seconds> (diff).count() << " seconds" << endl; 

	delete scene.rasterize;
	//delete scene.radiosity;

	return 0;
}
