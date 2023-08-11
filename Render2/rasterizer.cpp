#include <thread>
#include "ray_tracer.h"
#include "ppm_image.h"
#include "params.h"
#include "radiosity_helper.h"

int Sample::debug_line = 0;

using namespace std;

const struct RasterDetails {

	double xmin;
	double ymax; 
	double width;
	double height;
	
	Vec3 G ;          // Gaze direction.
	Vec3 U ;          // Up vector.
	Vec3 R ;          // Right vector.
	Vec3 O ;          // "Origin" of the 3D raster.
	Vec3 dR;          // Right increments.
	Vec3 dU;          // Up increments.

	static const double Lwidth ;
	static const double Lheight ;

	int DistLensDim ;
	Vec3 XGradation ;
	Vec3 YGradation ;
	//const Vec3 XGradation( Lwidth * R /cam.x_res);
	//const Vec3 YGradation ( Lheight * U / cam.y_res) ;
	int DistLensGradations  ;
	string img_file_name ;
	Scene * scene ;
	Camera * cam ;
	PPM_Image * I ;
	
	RasterDetails() ;
	RasterDetails(Scene & _scene, Camera & _cam, string & _in_file_name, PPM_Image & _I  ) {

			xmin = _cam.x_win.min ;
			ymax = _cam.y_win.max ;
			width  = Len( _cam.x_win );
			height = Len( _cam.y_win );

			// compute increments based on camera geometry. Will be used to define the ray direction at each pixel.

			G =  Unit( _cam.lookat - _cam.eye ) ;           // Gaze direction.
			U =  Unit( _cam.up / G ) ;                      // Up vector.
			R =  Unit( G ^ U ) ;                            // Right vector.
			O =  _cam.vpdist * G + 
						 xmin * R + 
						 ymax * U ;                 // "Origin" of the 3D raster.
			dR = width * R / _cam.x_res ;           // Right increments.
			dU =  height * U / _cam.y_res ;         // Up increments.

			DistLensDim = 7 ;
			XGradation = Lwidth * R / DistLensDim;
			YGradation = Lheight * U / DistLensDim ;
			DistLensGradations = DistLensDim / 2 ;

			scene = &_scene;
			cam = &_cam;
			I = &_I ;

			img_file_name = _in_file_name + ".ppm";
	}
};

const double RasterDetails::Lwidth = 0.375;
const double RasterDetails::Lheight = 0.3;
const int threadDivisionsInX = 8 ;
const int threadDivisionsInY = 6 ;

// Maps values that are in [0,1] and maps them to integers between 0 and 255.  Values above 1 are truncated.
// Doesn't implement HDR
// TODO HDR
static Pixel ToneMap( const Color &color )
{
	int red   = (int)floor( 256 * color.red   );
	int green = (int)floor( 256 * color.green );
	int blue  = (int)floor( 256 * color.blue  );
	channel r = (channel)( red   >= 255 ? 255 : red   ); 
	channel g = (channel)( green >= 255 ? 255 : green ); 
	channel b = (channel)( blue  >= 255 ? 255 : blue  );
	return Pixel( r, g, b );
}

static Pixel ToneMapRadiosity( const Color &color ) {
    /*
    int red   = (int)floor( 256 * color.red   );
	int green = (int)floor( 256 * color.green  );
	int blue  = (int)floor( 256 * color.blue    );
	channel r = (channel)( red   >= 255 ? 255 : red   ); 
	channel g = (channel)( green >= 255 ? 255 : green ); 
	channel b = (channel)( blue  >= 255 ? 255 : blue  );
	return Pixel( r, g, b );
	*/

    if ( color == Green ) return Pixel ( 0, 255, 0 );

	// or we can evade the bottom code and do a more simpler tone mapping sampling
	// with this active code, the bottom code is unreachable
	Color mappedColor;
	channel finalRGBAlt = 255;
	//if (finalRGBAlt == 255) return Pixel(finalRGBAlt, finalRGBAlt, finalRGBAlt);
	finalRGBAlt = (channel)((double)finalRGBAlt * (color.red / (color.red + 1.0)));
	return Pixel(finalRGBAlt, finalRGBAlt, finalRGBAlt);

	//constants
	static const double LUMINANCE_DISPLAY = 200.0 ;
	static const double CONTRAST_RATIO = 1000000.0;
	static const double GAMMA_FACTOR = 2.1 ;

	double a_rw = .41 * log10 (color.red) + 2.92;
	double a_disp = .41 * log10 (LUMINANCE_DISPLAY ) + 2.92 ;
	double b_rw = -.41 * pow ( log10 ( color.red ) , 2 ) + ( -2.584 * log10 ( color.red ) ) ;
	double b_disp = -.41 * pow ( log10 ( LUMINANCE_DISPLAY ) , 2 ) + ( -2.584 * log10 ( LUMINANCE_DISPLAY ) ) ;

	double inverseContrast = 1.0 / CONTRAST_RATIO ;

	double enumerator = pow ( color.red , a_rw / a_disp ) ;
	enumerator *= pow ( 10.0 , ( b_rw - b_disp) / a_disp ) ;
	enumerator /= LUMINANCE_DISPLAY ;
	enumerator -= inverseContrast ;
	enumerator = pow ( enumerator , 1.0 / GAMMA_FACTOR ) ;

	int col   = (int)floor( 256 * enumerator );
	channel finalRGB = (channel)( col  >= 255 ? 255 : col  );

	return Pixel( finalRGB, finalRGB, finalRGB );
}

void _thread_function_to_call_ ( const Rasterizer * _this ,
								 void (Rasterizer:: * fptr)( RasterDetails &, const int , const int ),
								 RasterDetails & rasterD,
								 const int idx,
								 const int idy
								 ) {

			  ((const_cast<Rasterizer*>(_this))->*fptr) ( rasterD, idx, idy ) ;
			  // ^ pointer to member function
}
								 
// Rasterize casts all the initial rays starting from the eye.
// casts one ray per pixel, in raster order, then writes the pixels out to a file.
bool Rasterizer::Rasterize( string file_name_out, const Camera &cam, const Scene &scene ) const
{
    try {
		// Create an image of the given resolution.
		PPM_Image I( cam.x_res, cam.y_res );
			
		//const Vec3 OLens ( cam.vpdist * G + xmin * LLens + ymax * ULens  

		Scene & tempScene = const_cast<Scene&>(scene);
		Camera & tempCamera = const_cast<Camera&>(cam);
		// Initiate struct of raster ray increments and axes
		RasterDetails rasterD(tempScene, tempCamera, file_name_out, I );
		
		/*
		std::thread * tt = new std::thread[threadDivisionsInX * threadDivisionsInY];
		int totalThreads = 0;

		for (int i = 0; i < threadDivisionsInX ; ++i) {
			for ( int j = 0 ; j < threadDivisionsInY ; ++j ) {
					 
				tt[totalThreads++] = std::thread(_thread_function_to_call_, this, &Rasterizer::Normal_Raster, rasterD, i, j);
				//tt[totalThreads++] = std::thread(_thread_function_to_call_, this, &Rasterizer::Anti_Aliasing, rasterD, i, j);
				//tt[totalThreads++] = std::thread(_thread_function_to_call_, this, &Rasterizer::Depth_Of_Field_Effect, rasterD, i, j);
				//tt[totalThreads++] = std::thread(_thread_function_to_call_, this, &Rasterizer::Radiosity_Raster, rasterD, i, j);
			}
		}

		for (int i = 0; i < threadDivisionsInX * threadDivisionsInY ; ++i) {
			tt[i].join();
		}
		delete [] tt ;
		cout << "Total threads started : " << totalThreads << endl ;
		*/
		
		Serial_Normal_Raster( rasterD ) ;
	}
	catch ( std::exception ex ) {
		cout << ex.what() << endl;
	}

	return true;
}

void Rasterizer::Serial_Normal_Raster ( RasterDetails & rasterD ) const {

	try {
		    PPM_Image & I = *rasterD.I ;
			// Find the form factor matrix.
			// Iterate through all the triangles and find the FF with all other triangles in the scene.
			// Use the structure of the QuadTree to iterate through all the triangles.
			//vector<QuadTreeNode * > tempQuadVector ;

			//for every object
			//for ( auto &v : _rad->quadTreeRoot->children ) {
				//_rad_helper->returnFilledElementsOfObject ( v, *_rad->tempQuadVector ) ;
			//}
		    // Create an image of the given resolution.
		    //PPM_Image I( rasterD.cam->x_res, rasterD.cam->y_res ); 
			
			//Camera & tempCamera = const_cast<Camera&>(rasterD.cam);
			// Initiate struct of raster ray increments and axes
			//RasterDetails rasterD(*_rad->scene, tempCamera, _fname, I );

			//PPM_Image & I = *rasterD.I ;

			// Initialize all the fields of the first-generation ray except for "direction".
			Ray ray;
			ray.origin     = rasterD.cam->eye;     // All initial rays originate from the eye.
			ray.type       = generic_ray; // These rays are given no special meaning.
			ray.generation = 1;           // Rays cast from the eye are first-generation.

			for ( unsigned int i = 0 ; i < rasterD.cam->y_res ; ++ i ) {
				    // Overwrite the line number written to the console.
					cout << rubout( i ) << (i+1);
					cout.flush();
				for ( unsigned int j = 0 ; j < rasterD.cam->x_res ; ++ j ) {
					ray.direction = Unit( rasterD.O + (j + 0.5) * rasterD.dR - (i + 0.5) * rasterD.dU  );
					I(i,j) = ToneMap( rasterD.scene->Trace ( ray ) );
				}
			}

			I.Write( rasterD.img_file_name );
			cout << "done." << endl;	
	}
	catch ( std::exception ex ) {
		cout << ex.what() << endl;
	}
}

void Rasterizer::Normal_Raster (RasterDetails & rasterD , const int idx, const int idy   ) {

	PPM_Image & I = *rasterD.I ;

	// Initialize all the fields of the first-generation ray except for "direction".
	Ray ray;
	ray.origin     = rasterD.cam->eye;     // All initial rays originate from the eye.
	ray.type       = generic_ray; // These rays are given no special meaning.
	ray.generation = 1;           // Rays cast from the eye are first-generation.

	// Loop over the entire image, casting a single ray per pixel.

	double workItemsPerThreadX = (double)rasterD.cam->x_res / (double)threadDivisionsInX ;
	double workItemsPerThreadY = (double)rasterD.cam->y_res / (double)threadDivisionsInY ;

	int endingWorkItemsX = idx * workItemsPerThreadX + workItemsPerThreadX ;
	int endingWorkItemsY = idy * workItemsPerThreadY + workItemsPerThreadY ;
	
	/********* NORMAL CODE **********/
								 
	for( unsigned int i = idy * workItemsPerThreadY ; i < endingWorkItemsY; i++ ){

		// Overwrite the line number written to the console.
		//cout << rubout( i ) << (i+1);
		//cout << spaceout( idy * threadDivisionsInX  + idx  ) << rubout( i ) << (i+1);
		//cout.flush();

		Sample::debug_line = i;

		for( unsigned int j = idx * workItemsPerThreadX ; j < endingWorkItemsX ; j++ ){

			// on irfanview, j = X and i = Y
			if ( j == 250 && i == 300 ) {
				//cout << " now -> " << endl ;
				ray.direction = Unit( rasterD.O + (j + 0.5) * rasterD.dR - (i + 0.5) * rasterD.dU  );
			}
			else ray.direction = Unit( rasterD.O + (j + 0.5) * rasterD.dR - (i + 0.5) * rasterD.dU  );

			//ray.direction = Unit (( O + (j + 0.5) * dR - (i + 0.5) * dU  ) - cam.eye);

			I(i,j) = ToneMap( rasterD.scene->Trace( ray ) );
		}
	}
	
	// Thus far the image exists only in memory.  Now write it out to a file.

	//cout << "\nWriting image file " << rasterD.img_file_name << "... ";
	//cout.flush();
	I.Write( rasterD.img_file_name );
	cout << "done." << endl;
	//return true;

}

void Rasterizer::Anti_Aliasing (RasterDetails & rasterD , const int idx, const int idy)
{
    PPM_Image & I = *rasterD.I ;

	// Initialize all the fields of the first-generation ray except for "direction".
	Ray ray;
	ray.origin     = rasterD.cam->eye;     // All initial rays originate from the eye.
	ray.type       = generic_ray; // These rays are given no special meaning.
	ray.generation = 1;           // Rays cast from the eye are first-generation.

	// Loop over the entire image, casting a single ray per pixel.

	double workItemsPerThreadX = (double)rasterD.cam->x_res / (double)threadDivisionsInX ;
	double workItemsPerThreadY = (double)rasterD.cam->y_res / (double)threadDivisionsInY ;

	int endingWorkItemsX = idx * workItemsPerThreadX + workItemsPerThreadX ;
	int endingWorkItemsY = idy * workItemsPerThreadY + workItemsPerThreadY ;

	/***********************    ANTI ALIASING CODE ***************/
	//cout << "Rendering line 0";
	for( unsigned int i = idy * workItemsPerThreadY ; i < endingWorkItemsY; i++ ){

		// Overwrite the line number written to the console.
		//cout << rubout( i ) << (i+1);
		//cout.flush();

		for( unsigned int j = idx * workItemsPerThreadX ; j < endingWorkItemsX ; j++ ){

			Color contributedColor ( 0.0, 0.0, 0.0 );

			for(float fragmenty = 1; fragmenty < 4; fragmenty += 2)
			{
				for(float fragmentx = 1; fragmentx < 4; fragmentx += 2)
				{
					float coef = 0.25f; 
					
					//ray.direction = Unit( O + (j + 0.5) * dR - (i + 0.5) * dU  );

					ray.direction = Unit( rasterD.O + ((j + ( fragmentx *  0.25)) * rasterD.dR) - ((i + (fragmenty *0.25)) * rasterD.dU));

					contributedColor += rasterD.scene->Trace( ray );
					
					//I(i,j) = ToneMap( scene.Trace( ray ) );
				}
			}

			contributedColor = contributedColor / 4.0;

			I(i,j) = ToneMap( contributedColor ); 
		}
	}
	I.Write( rasterD.img_file_name );
}

void Rasterizer::Depth_Of_Field_Effect(RasterDetails & rasterD , const int idx, const int idy)
{
	Color contributedDOF (0, 0, 0);

	PPM_Image & I = *rasterD.I ;

	// Initialize all the fields of the first-generation ray except for "direction".
	Ray ray;
	ray.origin     = rasterD.cam->eye;     // All initial rays originate from the eye.
	ray.type       = generic_ray; // These rays are given no special meaning.
	ray.generation = 1;           // Rays cast from the eye are first-generation.

	// Loop over the entire image, casting a single ray per pixel.

	double workItemsPerThreadX = (double)rasterD.cam->x_res / (double)threadDivisionsInX ;
	double workItemsPerThreadY = (double)rasterD.cam->y_res / (double)threadDivisionsInY ;

	int endingWorkItemsX = idx * workItemsPerThreadX + workItemsPerThreadX ;
	int endingWorkItemsY = idy * workItemsPerThreadY + workItemsPerThreadY ;

	/**************  Depth of Field code ***************/
	for( unsigned int i = idy * workItemsPerThreadY ; i < endingWorkItemsY; i++ ){
		//cout << rubout( i ) << (i+1);
		//cout.flush();

		 for( unsigned int j = idx * workItemsPerThreadX ; j < endingWorkItemsX ; j++ ){

			 for ( int yLens = rasterD.DistLensGradations ; yLens >= -rasterD.DistLensGradations ; yLens--){
				 
				 for ( int xLens = -rasterD.DistLensGradations ; xLens <= rasterD.DistLensGradations  ; xLens ++ ){

					 
					 Vec3 LensOrigin (  rasterD.cam->eye + ( xLens  * rasterD.XGradation) + ( yLens * rasterD.YGradation));
					 Vec3 GLens ( Unit( rasterD.cam->lookat - LensOrigin ) );
					 Vec3 O ( rasterD.cam->vpdist * GLens + rasterD.xmin * rasterD.R + rasterD.ymax * rasterD.U );
					 Vec3 focalPoint (  O + (j + 0.5) * rasterD.dR - (i + 0.5) * rasterD.dU);

					 ray.direction = Unit ( focalPoint);
					 ray.origin = LensOrigin;
					 //ray.origin = cam.eye;
					 ray.generation = 1;
					 ray.type = generic_ray ;
					 //contributedDOF += scene.Trace( ray );
					 contributedDOF += rasterD.scene->Trace( ray );
				 }
			 }
			 //cout << "here" << endl;
			 contributedDOF /= ( rasterD.DistLensDim * rasterD.DistLensDim );
			 I(i,j) = ToneMap( contributedDOF );
			 contributedDOF.red = contributedDOF.green = contributedDOF.blue = 0.0;
		 }
	}
	cout << "done." << endl;
	I.Write( rasterD.img_file_name );
	
}

bool Rasterizer::Radiosity_Raster ( string _fname , const Camera & _camera, Radiosity * _rad, Radiosity_Helper * _rad_helper ) {

	try {

			// Find the form factor matrix.
			// Iterate through all the triangles and find the FF with all other triangles in the scene.
			// Use the structure of the QuadTree to iterate through all the triangles.
			//vector<QuadTreeNode * > tempQuadVector ;

			//for every object
			//for ( auto &v : _rad->quadTreeRoot->children ) {
				//_rad_helper->returnFilledElementsOfObject ( v, *_rad->tempQuadVector ) ;
			//}
		    // Create an image of the given resolution.
		    PPM_Image I( _camera.x_res, _camera.y_res ); 

			Camera & tempCamera = const_cast<Camera&>(_camera);
			// Initiate struct of raster ray increments and axes
			RasterDetails rasterD(*_rad->scene, tempCamera, _fname, I );

			//PPM_Image & I = *rasterD.I ;

			// Initialize all the fields of the first-generation ray except for "direction".
			Ray ray;
			ray.origin     = rasterD.cam->eye;     // All initial rays originate from the eye.
			ray.type       = generic_ray; // These rays are given no special meaning.
			ray.generation = 1;           // Rays cast from the eye are first-generation.

			for ( unsigned int i = 0 ; i < rasterD.cam->y_res ; ++ i ) {
				    // Overwrite the line number written to the console.
					cout << rubout( i ) << (i+1);
					cout.flush();

				for ( unsigned int j = 0 ; j < rasterD.cam->x_res ; ++ j ) {
					
					// on irfanview, j = X and i = Y
					if ( i == 232 && j == 357 ) {
						int fdgfdg = 4 ;
						//ray.direction = Unit( rasterD.O + (j + 0.5) * rasterD.dR - (i + 0.5) * rasterD.dU  );
						//I(i,j) = ToneMapRadiosity( _rad_helper->trace_ray( ray, tempQuadVector ) );
					}

					//cout << i << " : " << j << endl ; 

					ray.direction = Unit( rasterD.O + (j + 0.5) * rasterD.dR - (i + 0.5) * rasterD.dU  );
					I(i,j) = ToneMapRadiosity( _rad_helper->trace_ray( ray, *_rad->tempQuadVector ) );
				}
			}

			I.Write( rasterD.img_file_name );
			cout << "done." << endl;	
	}
	catch ( std::exception ex ) {
		cout << ex.what() << endl;
	}

	return true;
}
