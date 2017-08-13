
#include <fstream>
#include <sstream>
#include "ppm_image.h"

using std::string;
using std::ifstream;
using std::ofstream;
using std::stringstream;

PPM_Image::PPM_Image( int x_res, int y_res )
{
    width  = x_res;
    height = y_res;
    pixels = new Pixel[ width * height ];
    Pixel *p = pixels;
    for( int i = 0; i < width * height; i++ ) *p++ = Pixel(0,0,0);
}

bool PPM_Image::Write( string file_name )
{
    ofstream fout( file_name.c_str(), ofstream::binary );
    if( !fout.is_open() ) return false;
    fout << "P6\n" 
         << width << " " << height << "\n"
         << 255 << "\n";
    fout.write( (const char *)pixels, 3 * height * width );
    fout.close();
    return true;
}

bool PPM_Image::Read( string file_name )
{
    string type;
    int depth;
    char buff[512];
    ifstream fin( file_name.c_str(), ifstream::binary );
    if( !fin.is_open() ) return false;

    fin.getline( buff, 512 );  // P6
    stringstream line1( buff );

    fin.getline( buff, 512 );  // width and height
    stringstream line2( buff );

    fin.getline( buff, 512 ); // Image depth.
    stringstream line3( buff );

    line1 >> type;
    line2 >> width >> height;
    line3 >> depth;

    if( type  != "P6" ) return false;
    if( depth !=  255 ) return false;

    delete pixels;
    pixels = new Pixel[ width * height ];
    Pixel *p = pixels;

    for( int i = 0; i < height; i++ )
    for( int j = 0; j < width ; j++ )
        {
        channel r = fin.get();
        channel g = fin.get();
        channel b = fin.get();
        *p++ = Pixel( r, g, b );
        }

    fin.close();
    return true;
}


