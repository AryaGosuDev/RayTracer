
#ifndef __PPM_IMAGE_INCLUDED__
#define __PPM_IMAGE_INCLUDED__

#include <string>

typedef unsigned char channel;  

struct Pixel {
    Pixel() { r = 0; g = 0; b = 0; }
    Pixel( channel R, channel G, channel B ) { r = R; g = G; b = B; }
    channel r;
    channel g;
    channel b;
};

struct PPM_Image {
    PPM_Image( int x_res, int y_res );
   ~PPM_Image() { delete[] pixels; }
    bool Read ( std::string file_name );
    bool Write( std::string file_name );
    inline  Pixel &operator()( int i, int j ) { return *( pixels + ( i * width + j ) ); }  
    Pixel  *pixels;
    int     width;
    int     height;
};

inline Pixel operator*( const float x, Pixel a ){
    return Pixel( static_cast<channel> (x * a.r), static_cast<channel> (x * a.g), static_cast<channel> (x * a.b) );
}

inline Pixel operator+( Pixel a, Pixel b ){
	return Pixel( a.r + b.r, a.g + b.g, a.b + b.b );
}

#endif


