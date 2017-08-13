
#include "ray_tracer.h"
#include "ppm_image.h"
#include "params.h"

int Sample::debug_line = 0;

using namespace std;


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


// Rasterize casts all the initial rays starting from the eye.
// casts one ray per pixel, in raster
// order, then writes the pixels out to a file.
bool Rasterizer::Rasterize( string file_name_out, const Camera &cam, const Scene &scene ) const
{
	file_name_out += ".ppm";

	
	// Create an image of the given resolution.

	PPM_Image I( cam.x_res, cam.y_res );

	// Initialize all the fields of the first-generation ray except for "direction".

	Ray ray;
	ray.origin     = cam.eye;     // All initial rays originate from the eye.
	ray.type       = generic_ray; // These rays are given no special meaning.
	ray.generation = 1;           // Rays cast from the eye are first-generation.

	const double xmin   = cam.x_win.min;
	const double ymax   = cam.y_win.max;
	const double width  = Len( cam.x_win );
	const double height = Len( cam.y_win );

	// Compute increments etc. based on the camera geometry.  These will be used
	// to define the ray direction at each pixel.

	const Vec3 G ( Unit( cam.lookat - cam.eye ) );          // Gaze direction.
	const Vec3 U ( Unit( cam.up / G ) );                    // Up vector.
	const Vec3 R ( Unit( G ^ U ) );                         // Right vector.
	const Vec3 O ( cam.vpdist * G + xmin * R + ymax * U );  // "Origin" of the 3D raster.
	const Vec3 dR( width  * R / cam.x_res );                // Right increments.
	const Vec3 dU( height * U / cam.y_res );                // Up increments.

	const double Lwidth = 0.375;
	const double Lheight = 0.3;
	
	const int DistLensDim = 7 ;
	const Vec3 XGradation = Lwidth * R /DistLensDim;
	const Vec3 YGradation = Lheight * U / DistLensDim ;
	//const Vec3 XGradation( Lwidth * R /cam.x_res);
	//const Vec3 YGradation ( Lheight * U / cam.y_res) ;
	const int DistLensGradations = DistLensDim /2 ;
	Color contributedDOF (0, 0, 0);
	
	//const Vec3 OLens ( cam.vpdist * G + xmin * LLens + ymax * ULens  

	// Loop over the entire image, casting a single ray per pixel.
	
	/********* NORMAL CODE **********/
	
	cout << "Rendering line 0";
	for( unsigned int i = 0; i < cam.y_res; i++ )
	{
		// Overwrite the line number written to the console.
		cout << rubout( i ) << (i+1);
		cout.flush();

		Sample::debug_line = i;

		if ( i == 252 ) {
			int fdgfdg = 4 ;
		}

		for( unsigned int j = 0; j < cam.x_res; j++ )
		{
			// i == y
			// j == x
			if ( j == 177 && i == 152 ) {
				int fdgfd = 43;
			}

			if ( j == 206 && i == 147 ) {
				int fdgfd = 43;
			}

			if ( j == 236 && i == 142 ) {
				int fdgfdg = 5 ;
			}

			ray.direction = Unit( O + (j + 0.5) * dR - (i + 0.5) * dU  );

			//ray.direction = Unit (( O + (j + 0.5) * dR - (i + 0.5) * dU  ) - cam.eye);

			I(i,j) = ToneMap( scene.Trace( ray ) );
		}
	}
	
	// Thus far the image exists only in memory.  Now write it out to a file.

	cout << "\nWriting image file " << file_name_out << "... ";
	cout.flush();
	I.Write( file_name_out );
	cout << "done." << endl;
	return true;
}

bool Rasterizer::Anti_Aliasing  ()
{
		/***********************    ANTI ALIASING CODE ***************/
	/*
	cout << "Rendering line 0";
	for( unsigned int i = 0; i < cam.y_res; i++ )
	{
		// Overwrite the line number written to the console.
		cout << rubout( i ) << (i+1);
		cout.flush();

		for( unsigned int j = 0; j < cam.x_res; j++ )
		{
			Color contributedColor ( 0, 0, 0 );

			for(float fragmenty = 1; fragmenty < 4; fragmenty += 2)
			{
				for(float fragmentx = 1; fragmentx < 4; fragmentx += 2)
				{
					float coef = 0.25f; 
					
					//ray.direction = Unit( O + (j + 0.5) * dR - (i + 0.5) * dU  );

					ray.direction = Unit( O + ((j + ( fragmentx *  0.25)) * dR) - ((i + (fragmenty *0.25)) * dU)  );

					contributedColor += scene.Trace( ray );
					
					//I(i,j) = ToneMap( scene.Trace( ray ) );
				}

			}

			contributedColor = contributedColor / 4;

			I(i,j) = ToneMap( contributedColor ); 
		}
	}
	*/
	return true;
}

bool Rasterizer::Depth_Of_Field_Effect()
{
	/**************  Depth of Field code ***************/
	/*
	for ( unsigned int i = 0; i < cam.y_res ; i++ )
	{
		cout << rubout( i ) << (i+1);
		cout.flush();

		 for( unsigned int j = 0; j < cam.x_res; j++ )
		 {
			 for ( int yLens = DistLensGradations ; yLens >= -DistLensGradations ; yLens--)
			 {
				 for ( int xLens = -DistLensGradations ; xLens <= DistLensGradations  ; xLens ++ )
				 {
					 //cout << "here" << endl;
					 Vec3 LensOrigin (  cam.eye + ( xLens  * XGradation) + ( yLens * YGradation));
					 Vec3 GLens ( Unit( cam.lookat - LensOrigin ) );
					 Vec3 O ( cam.vpdist * GLens + xmin * R + ymax * U );
					 Vec3 focalPoint (  O + (j + 0.5) * dR - (i + 0.5) * dU);

					// Vec3 focalL ( focalPoint
					 

					 //ray.direction = Unit ( focalPoint - cam.eye );
					 //ray.direction = Unit (focalPoint - LensOrigin) ;

					 ray.direction = Unit ( focalPoint);
					 ray.origin = LensOrigin;
					 //ray.origin = cam.eye;
					 ray.generation = 1;
					 ray.type = generic_ray ;
					 //contributedDOF += scene.Trace( ray );
					 contributedDOF += scene.Trace( ray );


				 }
			 }
			 contributedDOF = contributedDOF / ( DistLensDim * DistLensDim );
			 I(i,j) = ToneMap( contributedDOF );
			 contributedDOF.red = 0;
			 contributedDOF.green = 0;
			 contributedDOF.blue = 0;
		 }
	}
	*/
	return true ;
}

