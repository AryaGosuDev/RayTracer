#ifndef __TOYTRACER_INCLUDED__  // Include this file only once.
#define __TOYTRACER_INCLUDED__

#include "base.h"         // Includes system headers and defines constants.
#include "vec3.h"         // Defines the Vec3 class, which are points in R3.
#include "vec2.h"         // Defines the Vec2 class, which are points in R2.
#include "mat3x2.h"       // Defines 3x2 matrices. for line-line intersections.
#include "mat3x3.h"       // Defines 3x3 matrices.
#include "mat3x4.h"       // Defines 3x4 matrices; i.e. affine transforms.
#include "color.h"        // Defines the Color class; real RGB values.
#include "ray.h"          // Defines rays in 3-space: origin, direction, etc.
#include "interval.h"     // Defines a (min,max) interval of the real line.
#include "params.h"       // Parsing an object file.

const struct RasterDetails;

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

struct Transformation {
	Vec3 translation;
	Vec3 rotation;
	Vec3 scale ;
	bool Istranslation;
	bool Isrotation;
	bool Isscale;

	Transformation() : translation ( 0.0, 0.0, 0.0 ),
					   rotation ( 0.0, 0.0, 0.0) ,
					   scale ( 0.0, 0.0, 0.0 ),
					   Istranslation ( false ),
					   Isrotation ( false ),
					   Isscale ( false ){
	}

	void reset() {
		translation = Vec3 ( 0.0, 0.0, 0.0 );
		rotation = Vec3 ( 0.0, 0.0, 0.0 );
		scale = Vec3 ( 0.0, 0.0, 0.0 );
		Istranslation = false ;
		Isrotation = false ;
		Isscale = false ;
	}
};

struct HitInfo {          // Records all info at ray-object intersection.
	const Object *ignore; // One object in scene to ignore (used by Cast).
	const Object *object; // The object that was hit (set by Intersect).
	BSP_Node * tri;		  // The triangle that has been hit by the ray
	double  distance;     // Distance to hit (used & reset by Intersect).
	Vec3    point;        // ray-object intersection point (set by Intersect).
	Vec3    normal;       // Surface normal (set by Intersect).
	Vec2    uv;           // Texture coordinates (set by intersect).
	Ray     ray;          // The ray that hit the surface (set by Cast).
	float   BayCentricU;  // Baycentric parameters of the triangle
	float	BayCentricV;  

	HitInfo() : ignore(0) , object(0), tri(0) {}

	HitInfo & operator= ( const HitInfo & rhs ) {
		ignore = rhs.ignore;
		object = rhs.object;
		tri = rhs.tri;
		distance = rhs.distance;
		point = rhs.point;
		normal = rhs.normal;
		uv = rhs.uv;
		ray = rhs.ray;
		BayCentricU = rhs.BayCentricU;
		BayCentricV = rhs.BayCentricV;

		return *this;
	}
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
	Radiosity *radiosity;    // Pointer to the radiosity object. Will initiate and calculate the radiosity of the scene.
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
	virtual double findOcclusionNew (const Scene &, const HitInfo &, Color &color, double radiusOfHemisphere, int  ) const ; 

	virtual string MyName() const { return "shader"; }
	virtual bool Default() const { return true; }

	int ambientOcclusionSamples ;
	int ambientOcclusionHemisphereRadius ;

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

	void Normal_Raster(RasterDetails & , const int, const int);
	void Anti_Aliasing(RasterDetails & , const int, const int);
	void Depth_Of_Field_Effect(RasterDetails & , const int, const int);
	bool Radiosity_Raster ( string , const Camera &camera, Radiosity *, Radiosity_Helper *  ) ;
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
	bool flipYZ ;
	BSP_Node * root ;
	Material * material;
	Transformation trnsf ;
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


	

// create initial quadtree
// must sustain the adjacency list of all triangles with other triangles
struct QuadTreeNode : public BSP_Node {
	
	QuadTreeNode () {
		object = NULL;
		parentObject = NULL ;
		visited = false;
		delVisited = false;
	}

	~QuadTreeNode() {}

	QuadTreeNode ( BSP_Node * _BSP_Node ) {

		object = _BSP_Node->parentObject;
		parentObject = _BSP_Node->parentObject;
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
		left = _BSP_Node->left;
		right = _BSP_Node->right ;
	}

	

	vector<QuadTreeNode *> children ;
	// maintain QuadTreeNode vector of all adjacent nodes
	//std::map<QuadTreeNode *, std::set<Vec3>> nextAdj;
	std::map<QuadTreeNode *, std::unordered_set<Vec3>> nextAdj;
	

	Object * object;
	Object * parentObject;

	double radiosityValue;
	double reflectivityIndex;
	double emissitivity;
	double formFactor ;

	bool visited;
	bool delVisited;
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

struct Radiosity {
	Radiosity() {}
	Radiosity(  Scene * _scene, Camera * ) ;
	virtual ~Radiosity() ;

	virtual void buildInitialQuadTree();
	void findFormFactor();
	double findFormFactorTerm ( QuadTreeNode * , QuadTreeNode * );
	double findFormFactorTermWithLight ( QuadTreeNode * );
	bool isVisibleXiToXj ( QuadTreeNode * , QuadTreeNode * , const Vec3 &, const Vec3 & ) ;
	bool castFromElementToElement ( const Ray &, HitInfo &, QuadTreeNode *, QuadTreeNode *, const Vec3 &, const Vec3 & ) ;
	void progressiveRefinement();
	bool checkIfConverged ( double *  ) ;
	//bool adapdtiveMeshSubDivision() ;

	Scene * scene ;
	Camera * cam;
	int numOfElements;
	QuadTreeNode * quadTreeRoot ;
	Radiosity_Helper * radiosityHelper ;

};

// include template specialization for std::set<Vec3> in the Radiosity class
namespace std {
		template<>
		struct hash<Vec3> {
			typedef size_t result_type;
			typedef Vec3 argument_type ;
			size_t operator() ( const Vec3 & v ) const ;
		};
		
		inline size_t hash<Vec3>::operator()(const Vec3 & v ) const {
			return hash<double>()( v.x ) ^ hash<double>() ( v.y ) ^ hash<double>() ( v.z ) ;
		}
	}

#endif