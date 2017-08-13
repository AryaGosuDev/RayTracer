#include "ray_tracer.h"

struct InterceptInfo { // Intermediate data structure for the intersect function

	bool found ;
	const Ray * ray ;
	HitInfo * hitinfo ;
	BSP_Node * root ;
	Scene tempScene;
	float u,v;
};

inline void updateHitInfo ( double & distanceToIntersection, InterceptInfo & info, BSP_Node * current ) {
	info.found = true;
	info.hitinfo->ray = *(info.ray) ;
	info.hitinfo->normal = current->triNormal;
	info.hitinfo->distance = distanceToIntersection ;
	info.hitinfo->uv = current->textVert;
	info.hitinfo->BayCentricU = info.u;
	info.hitinfo->BayCentricV = info.v;
	info.hitinfo->tri = current;
	//info.hitinfo->object = current ;
	info.hitinfo->point = info.ray->origin + (info.ray->direction * ((current->triVert1 - info.hitinfo->ray.origin ) * current->triNormal )/(current->triNormal * info.hitinfo->ray.direction)) ;
}

// orig and dir defines the ray. v0, v1, v2 defines the triangle.
// returns the distance from the ray origin to the intersection or 0.
// Möller-Trumbore algorithm
double triangle_intersection(const Vec3& orig, const Vec3& dir, const Vec3 & planeNormal, const Vec3& v0, const Vec3& v1, const Vec3& v2, float &u, float & v ) {
	Vec3 e1 = v1 - v0;
	Vec3 e2 = v2 - v0;
	// Calculate planes normal vector
	Vec3 pvec = dir ^ e2 ;
	double det = e1 * pvec ;

	// Ray is parallel to plane
	if (det < 1e-8 && det > -1e-8) {
		return 0;
	}

	// Check if the triange is behind the ray
	double inv_det = 1 / det;
	Vec3 tvec = orig - v0;
	u = (tvec * pvec ) * inv_det;
	if (u < 0 || u > 1) 
		return 0;
	
	Vec3 qvec = tvec ^ e1;
	v = (dir * qvec ) * inv_det;
	if (v < 0 || u + v > 1)
		return 0;

	//determing U and V alt method

	/*

	float denom = planeNormal * planeNormal;

	 Vec3 P = orig + (dir * ((v0 - orig ) * planeNormal )/(planeNormal * dir)) ;

	 // edge 0
	Vec3 edge0 = v1 - v0; 
	Vec3 vp0 = P - v0; 
	Vec3 C = edge0^ vp0; 

	 // edge 1
	Vec3 edge1 = v2 - v1; 
	Vec3 vp1 = P - v1; 
	C = edge1 ^ vp1; 
	if ((u = planeNormal*C) < 0)  return 0; // P is on the right side 
 
	// edge 2
	Vec3 edge2 = v0 - v2; 
	Vec3 vp2 = P - v2; 
	C = edge2^vp2; 
	if ((v = planeNormal*C) < 0) return 0; // P is on the right side; 
 
	u /= denom; 
	v /= denom; 

	*/
	
	return (e2 * qvec) * inv_det;
}

double rayTriangleIntersect( 
	const Vec3 &orig, const Vec3 &dir, const Vec3 & planeNormal,
	const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, 
	float &u, float &v) 
{ 
	double t;
	// compute plane's normal
	Vec3 v0v1 = v1 - v0; 
	Vec3 v0v2 = v2 - v0; 
	// no need to normalize
	Vec3 N = v0v1 ^ v0v2; // N 
	float denom = N * N; 
 
	// Step 1: finding P
 
	// check if ray and plane are parallel ?
	float NdotRayDirection = N * dir; 
	if (fabs(NdotRayDirection) < 1e-8) // almost 0 
		return false; // they are parallel so they don't intersect ! 
 
	// compute d parameter using equation 2
	double d = N * v0; 
 
	// compute t (equation 3)
	t = ((N*orig) + d) / NdotRayDirection; 
	// check if the triangle is in behind the ray
	if (t < 0) return false; // the triangle is behind 
 
	// compute the intersection point using equation 1
	Vec3 P = orig + t * dir; 
 
	// Step 2: inside-outside test
	Vec3 C; // vector perpendicular to triangle's plane 
 
	// edge 0
	Vec3 edge0 = v1 - v0; 
	Vec3 vp0 = P - v0; 
	C = edge0^ vp0; 
	if (N*C < 0) return false; // P is on the right side 
 
	// edge 1
	Vec3 edge1 = v2 - v1; 
	Vec3 vp1 = P - v1; 
	C = edge1 ^ vp1; 
	if ((u = N*C) < 0)  return false; // P is on the right side 
 
	// edge 2
	Vec3 edge2 = v0 - v2; 
	Vec3 vp2 = P - v2; 
	C = edge2^vp2; 
	if ((v = N*C) < 0) return false; // P is on the right side; 
 
	u /= denom; 
	v /= denom; 
 
	return d; // this ray hits the triangle 
} 

// Recursive depth first search intersect test ray with triangle
bool intersectTreeTraversal (InterceptInfo & info, BSP_Node * current ){ 

	double distanceToIntersection = 0;
	
	//leaf node, recursion base case
	if ( !info.found && current != NULL && current->left == NULL && current->right == NULL ) { 
		distanceToIntersection = triangle_intersection( info.ray->origin, info.ray->direction, current->triNormal, current->triVert1, current->triVert2, current->triVert3, info.u, info.v );
		//distanceToIntersection = rayTriangleIntersect( info.ray->origin, info.ray->direction, current->triNormal, current->triVert1, current->triVert2, current->triVert3, info.u, info.v );
		if ( distanceToIntersection > 0.0 ){
			updateHitInfo ( distanceToIntersection, info, current );
			return true;
		}
	}
	else if ( !info.found && current != NULL ) {
		
		//inner node, find the side of the hyperplane that the ray origin resides
		int side = sideTest3d( current->triVert1, current->triVert2, current->triVert3, info.ray->origin );

		if ( side == 0 ) {
			side = sideTest3d( current->triVert1, current->triVert2, current->triVert3, info.ray->origin + (info.ray->direction * 1.1 ) );
			if ( side == -1 ){
				intersectTreeTraversal ( info , current->right );
				intersectTreeTraversal ( info , current->left );
			}
			else if ( side == 1 ){
				intersectTreeTraversal ( info , current->left );
				intersectTreeTraversal ( info , current->right );
			}
			else {
				cout << "error in intersectTreeTraversal, side test result error" << endl ;
				exit(0);
			}
		}
		else if ( side == -1 ){
			
			intersectTreeTraversal ( info , current->right );
			intersectTreeTraversal ( info , current->left );

			distanceToIntersection = triangle_intersection( info.ray->origin, info.ray->direction, current->triNormal, current->triVert1, current->triVert2, current->triVert3, info.u, info.v );
			//distanceToIntersection = rayTriangleIntersect( info.ray->origin, info.ray->direction, current->triNormal, current->triVert1, current->triVert2, current->triVert3, info.u, info.v );
		
			if ( distanceToIntersection > 0.0 && info.hitinfo->distance > distanceToIntersection  ) {
				updateHitInfo ( distanceToIntersection, info, current );
				return true;
			}
		}
		else if ( side == 1 ){
			
			intersectTreeTraversal ( info , current->left );
			intersectTreeTraversal ( info , current->right );

			distanceToIntersection = triangle_intersection( info.ray->origin, info.ray->direction, current->triNormal, current->triVert1, current->triVert2, current->triVert3, info.u, info.v );
			//distanceToIntersection = rayTriangleIntersect( info.ray->origin, info.ray->direction, current->triNormal, current->triVert1, current->triVert2, current->triVert3, info.u, info.v );
		
			if ( distanceToIntersection > 0.0 && info.hitinfo->distance > distanceToIntersection ) {
				updateHitInfo ( distanceToIntersection, info, current );
				return true;
			}
		}
		else {
			cout << "error in intersectTreeTraversal, side test result error" << endl ;
			exit(0); 
		}
	}
	return false;
}

bool Object::Intersect (const Ray & ray, HitInfo & hitinfo ) const {

	InterceptInfo _info ;
	_info.found = false;
	Scene _thisScene;
	_info.root = this->root;
	_info.tempScene = _thisScene;
	//_info.ray = &( const_cast < Ray & > (ray) ) ;
	_info.ray = &ray ;
	_info.hitinfo = &hitinfo ;

	intersectTreeTraversal ( _info , _info.root ) ;

	if ( _info.found == true ) {
		//cout << Sample::debug_line << endl;
	}

	return _info.found;

}

