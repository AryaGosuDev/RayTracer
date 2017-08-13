#ifndef __TOYTRACER_INCLUDED__  // Include this file only once.
#define __TOYTRACER_INCLUDED__

#include "base.h"         // Includes system headers and defines constants.
#include "vec3.h"         // Defines the Vec3 class, which are points in R3.
#include "vec2.h"         // Defines the Vec2 class, which are points in R2.
#include "mat3x3.h"       // Defines 3x3 matrices.
#include "mat3x4.h"       // Defines 3x4 matrices; i.e. affine transforms.
#include "color.h"        // Defines the Color class; real RGB values.
#include "ray.h"          // Defines rays in 3-space: origin, direction, etc.
#include "interval.h"     // Defines a (min,max) interval of the real line.
#include "params.h"



struct Sample {           // A point and weight returned from a sampling algorithm.
	Vec3   P;
	double w;
	static int debug_line ;
	};

struct Material {         // Surface material for shading.
	Color  diffuse;       // Diffuse color.
	Color  specular;      // Color of highlights.
	Color  emission;      // Emitted light.
	Color  ambient;       // Ambient light (from all directions).
	Color  reflectivity;  // Weights for refleced light, between 0 and 1.
	Color  translucency;  // Weights for refracted light, between 0 and 1.
	double Phong_exp;     // Phong exponent for specular highlights.
	double ref_index;     // Refractive index.
	long   type;          // Reserved for future use.
	};

struct HitInfo {          // Records all info at ray-object intersection.
	const Object *ignore; // One object in scene to ignore (used by Cast).
	const Object *object; // The object that was hit (set by Intersect).
	BSP_Node * tri; // The triangle that has been hit by the ray
	double  distance;     // Distance to hit (used & reset by Intersect).
	Vec3    point;        // ray-object intersection point (set by Intersect).
	Vec3    normal;       // Surface normal (set by Intersect).
	Vec2    uv;           // Texture coordinates (set by intersect).
	Ray     ray;          // The ray that hit the surface (set by Cast).
	float   BayCentricU;  // Baycentric parameters of the triangle
	float	BayCentricV;  
};

struct Camera {           // Defines the position of the eye/camera.
	Vec3     eye;         // Position of eye.
	Vec3     lookat;      // The point we are looking toward.
	Vec3     up;          // A vector not parallel to the direction of gaze.
	double   vpdist;      // Distance to the view plane.
	Interval x_win;       // Horizontal extent of view window (typically [-1,1]).
	Interval y_win;       // Vertical extent of view window (typically [-1,1]).
	unsigned x_res;       // Horizontal image resolution in pixels.
	unsigned y_res;       // Vertical image resolution in pixels.
	};

struct PrimitiveObject ;

struct Scene {
	Scene();
	virtual ~Scene() { lights.clear(); }
	Color Trace( const Ray &ray ) const;
	bool  Cast ( const Ray &ray, HitInfo &hitinfo ) const;
	bool  IsAlreadyInside ( const Vec3 &p, HitInfo &hitinfo ) const ;
	bool  BuildScene ( string file, string fileObj, Camera &camera );
	bool  BuildBSP () ;
	void  randomizeTriangles ( vector <BSP_Node> & );
	bool calculateVertexNormals ( Object & obj, vector<vector<Vec3>> &  adjMatrix ) ;
	BSP_Node * buildBSPTree (vector<BSP_Node *>  &v);
	virtual const PrimitiveObject *GetLight( unsigned i ) const { return lights[i]; } 
	virtual unsigned NumLights() const { return (unsigned ) lights.size(); }
	Rasterizer *rasterize;   // This casts all primary rays & makes the image.
	vector<PrimitiveObject*> lights;  // All objects that are emitters.
	vector<Object*> sceneObjects ;

	Shader * shader ;

	unsigned max_tree_depth; // Limit on depth of the ray tree.
};

struct Shader {  
	Shader() {}
	virtual ~Shader() {}
	virtual Color Shade( const Scene &scene, const HitInfo &hitinfo ) const ;
	virtual Color generateNormalShadows (const Scene &, const HitInfo &, Color &color ) const ; 
	virtual Color generateSoftShadows (const Scene &, const HitInfo &, Color &color ) const ; 
	virtual double findOcclusion (const Scene &, const HitInfo &, Color &color ) const ; 

	
	virtual string MyName() const { return "shader"; }
	virtual bool Default() const { return true; }
	};

struct Envmap { // Each object can have its own environment map.
	Envmap() {}
	virtual ~Envmap() {}
	virtual Color Shade( const Ray &ray ) const = 0;
	};

struct Rasterizer  {  // The rasterizer creates all the primary rays.
	Rasterizer() {}
	virtual ~Rasterizer() {}
	bool Rasterize(
		string fname,         // File to write to (must include the extension).
		const Camera &camera, // Defines the view.
		const Scene &scene    // Global scene description: object, envmap, etc.
		) const;

	virtual bool Anti_Aliasing();
	virtual bool Depth_Of_Field_Effect();
	};

struct Object{
	Object () {}
	virtual ~Object () {}
	virtual bool Intersect (const Ray & ray, HitInfo & hitinfo ) const ;
	virtual double Cost() const { return 1.0; }
	string nameOfObject;
	vector <BSP_Node> triangles ;
	vector <Vec3> points;
	vector <Vec3> vertexNormals;
	vector <Vec2> vertexTexture;
	BSP_Node * root ;
	Material * material;
};

struct BSP_Node{
	BSP_Node () {
		right = 0;
		left = 0;
	};
	BSP_Node (int _node) {
		right = 0;
		left = 0;
		node = _node;
	}
	virtual ~BSP_Node () {}
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
	Vec2 textVert;
	int node;
	BSP_Node * right ;
	BSP_Node * left  ;

	Object * parentObject ;
};



struct PrimitiveObject{
	PrimitiveObject() { material = 0; shader = 0; envmap = 0;  }
	virtual ~PrimitiveObject() {}
	virtual PrimitiveObject *ReadString( const string & ) = NULL; // Instance from string.
	virtual bool Intersect( const Ray &ray, HitInfo & ) const = 0;
	virtual unsigned GetSamples( const Vec3 &P, const Vec3 &N, Sample *samples, unsigned n ) const { return 0; }
	virtual bool Inside( const Vec3 & ) const = 0;
	virtual Interval GetSlab( const Vec3 & ) const = 0;
	virtual double Cost() const { return 1.0; }
	Material  *material;
	Shader    *shader;
	Envmap    *envmap;
};

struct Light : public PrimitiveObject {
	Light() {}
	virtual ~Light() {}
};

struct SphereLight : public Light {
	SphereLight() {}
	SphereLight( const Vec3 &center, double radius );
	virtual bool Intersect( const Ray &ray, HitInfo & ) const;
	virtual bool Inside( const Vec3 &P ) const { return dist( P, center ) <= radius; } 
	virtual Interval GetSlab( const Vec3 & ) const;
	virtual int GetSamples( const Vec3 &P, const Vec3 &N, Sample *samples, int n ) const;
	virtual PrimitiveObject *ReadString( const string &params );
	virtual string MyName() const { return "spherelight"; }
	Vec3   center;
	double radius;
	double radius2;
};

struct QuadLight : public Light {
	QuadLight() {}
	QuadLight( const Vec3 &A, const Vec3 &B, const Vec3 &C, const Vec3 &D );
	virtual ~QuadLight() {}
	virtual bool Intersect( const Ray &ray, HitInfo & ) const;
	virtual bool Inside( const Vec3 & ) const { return false; }
	virtual Interval GetSlab( const Vec3 & ) const;
	virtual int GetSamples( const Vec3 &P, const Vec3 &N, Sample *samples, int n ) const;
	virtual PrimitiveObject *ReadString( const string &params );
	virtual string MyName() const { return "quadlight"; }
	Vec3   N;    // Normal to plane of quad;
	double d;    // Distance from origin to plane of quad.
	double area; // Area of the quad.
	Vec3   Eab, Ebc, Ecd, Eda;
	Vec3   A, B, C, D;  
};

#endif