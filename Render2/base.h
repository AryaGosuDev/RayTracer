


#ifndef __BASE_INCLUDED__    // Include this file only once.
#define __BASE_INCLUDED__

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <stdlib.h>
#include <time.h>     

using std::cout;
using std::cerr;
using std::endl; 
using std::ostream;
using std::string;
using std::vector;

// We declare all the main structures here so that we can include
// pointers to these objects before their full definitions are provided.

struct Object;     // Anything that can be ray traced.
struct Shader;     // Shaders that can be associated with surfaces.
struct Material;   // Surface parameters available to shaders.
struct Camera;     // The parameters of the camera, image resolution, etc.
struct Scene;      // The camera, lights, object(s), etc.
struct Rasterizer; // The function that casts primary rays & creates an image.
struct Builder;    // Builds the scene, usually by reading a file (e.g. sdf).
struct BSP_Node;   // Node of the tree

// Path constants.
static const string DefaultScene = "scenes/";

// Miscellaneous numerical constants.

static const double 
    Pi          = 3.14159265358979,
    TwoPi       = 2.0 * Pi,
    FourPi      = 4.0 * Pi,
    DegToRad    = Pi / 180.0,  // Convert degrees to radians.
    Infinity    = 1.0E20,      // Should suffice for "infinity"
    Epsilon     = 1.0E-5,      // A reasonable value for the "machine epsilon".
    OnePlusEps  = 1.0 + Epsilon,
    OneMinusEps = 1.0 - Epsilon;

// Miscellaneous default values.

static const int
    default_image_width  = 400,  // Default image width (x resolution).
    default_image_height = 400,  // Default image height (y resolution).
    default_max_tree_depth = 2;  // Default cap on ray tree depth.

enum raytracer_error {
    no_errors = 0,
    error_opening_input_file,
    error_reading_input_file,
    error_opening_image_file,
    error_no_builder,
    error_building_scene,
    error_no_rasterizer,
    error_rasterizing_image
    };


static const char *rubout( int i )
{
    if( i < 0 ) i = -10 * i; // Add a backspace for the negative sign.
    if( i < 10   ) return "\b";
    if( i < 100  ) return "\b\b";
    if( i < 1000 ) return "\b\b\b";
    return "\b\b\b\b";
}

#endif

