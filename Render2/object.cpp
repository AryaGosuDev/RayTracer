#include "ray_tracer.h"

struct InterceptInfo { // Intermediate data structure for the intersect function

	bool found ;
	const Ray * ray = NULL ;
	HitInfo * hitinfo = NULL ;
	BSP_Node * root = NULL ;
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
double triangle_intersection(const Vec3& orig, 
							 const Vec3& dir, 
							 const Vec3 & planeNormal, 
							 const Vec3& v0, 
							 const Vec3& v1, 
							 const Vec3& v2, 
							 float &u, float & v ) {
	Vec3 e1 = v1 - v0;
	Vec3 e2 = v2 - v0;
	// Calculate planes normal vector
	Vec3 pvec = dir ^ e2 ;
	double det = e1 * pvec ;

	// Ray is parallel to plane
	if (det < 1e-8 && det > -1e-8) {
		return 0.0;
	}

	// Check if the triange is behind the ray
	double inv_det = 1.0 / det;
	Vec3 tvec = orig - v0;
	u = (tvec * pvec ) * inv_det;
	if (u < 0.0f || u > 1.0f) 
		return 0.0;
	
	Vec3 qvec = tvec ^ e1;
	v = (dir * qvec ) * inv_det;
	if (v < 0.0f || u + v > 1.0f)
		return 0.0;
	
	return (e2 * qvec) * inv_det;
}

// Recursive depth first search intersect test ray with triangle
bool intersectTreeTraversal (InterceptInfo & info, BSP_Node * current ){ 

	double distanceToIntersection = 0.0;
	BSP_Node* hyperPlane = current;
	//leaf node, recursion base case
	if ( !info.found && hyperPlane != NULL && hyperPlane->left == NULL && hyperPlane->right == NULL ) {
		distanceToIntersection = triangle_intersection( info.ray->origin, info.ray->direction, hyperPlane->triNormal, hyperPlane->triVert1, hyperPlane->triVert2, hyperPlane->triVert3, info.u, info.v );
		if ( distanceToIntersection > 0.0 ){
			updateHitInfo ( distanceToIntersection, info, hyperPlane);
			return true;
		}
	}
	else if ( !info.found && hyperPlane != NULL ) {
		
		
		//inner node, find the side of the hyperplane that the ray origin resides
		int side = sideTest3d(hyperPlane->triVert1, hyperPlane->triVert2, hyperPlane->triVert3, info.ray->origin );

		if ( side == 0 ) {
			side = sideTest3d(hyperPlane->triVert1, hyperPlane->triVert2, hyperPlane->triVert3, info.ray->origin + (info.ray->direction * 1.1));
			if ( side == -1 ){
				intersectTreeTraversal ( info , hyperPlane->right );
				intersectTreeTraversal ( info , hyperPlane->left );
			}
			else {
				intersectTreeTraversal(info, hyperPlane->left);
				intersectTreeTraversal(info, hyperPlane->right);
			}
			/*
			else if ( side == 1 ){
				intersectTreeTraversal ( info , hyperPlane->left );
				intersectTreeTraversal ( info , hyperPlane->right );
			}
			else {
				//cout << "error in intersectTreeTraversal, side test result error : the ray being shot is co-planar with the plane shooting the ray " << endl << endl ;
				return false;
			}
			*/
		}
		else if ( side == -1 ){
			if (intersectTreeTraversal(info, hyperPlane->right)) {
				distanceToIntersection = triangle_intersection(info.ray->origin, info.ray->direction, hyperPlane->triNormal, hyperPlane->triVert1, hyperPlane->triVert2, hyperPlane->triVert3, info.u, info.v);
				if (distanceToIntersection > 0.0 && info.hitinfo->distance > distanceToIntersection) {
					updateHitInfo(distanceToIntersection, info, hyperPlane);
					return true;
				}
			}
			if (intersectTreeTraversal(info, hyperPlane->left)) {
				distanceToIntersection = triangle_intersection(info.ray->origin, info.ray->direction, hyperPlane->triNormal, hyperPlane->triVert1, hyperPlane->triVert2, hyperPlane->triVert3, info.u, info.v);
				if (distanceToIntersection > 0.0 && info.hitinfo->distance > distanceToIntersection) {
					updateHitInfo(distanceToIntersection, info, hyperPlane);
					return true;
				}
			}
		}
		else if ( side == 1 ){
			
			if (intersectTreeTraversal(info, hyperPlane->left)) {
				distanceToIntersection = triangle_intersection(info.ray->origin, info.ray->direction, hyperPlane->triNormal, hyperPlane->triVert1, hyperPlane->triVert2, hyperPlane->triVert3, info.u, info.v);
				if (distanceToIntersection > 0.0 && info.hitinfo->distance > distanceToIntersection) {
					updateHitInfo(distanceToIntersection, info, hyperPlane);
					return true;
				}
			}
			if (intersectTreeTraversal(info, hyperPlane->right)) {
				distanceToIntersection = triangle_intersection(info.ray->origin, info.ray->direction, hyperPlane->triNormal, hyperPlane->triVert1, hyperPlane->triVert2, hyperPlane->triVert3, info.u, info.v);

				if (distanceToIntersection > 0.0 && info.hitinfo->distance > distanceToIntersection) {
					updateHitInfo(distanceToIntersection, info, hyperPlane);
					return true;
				}
			}
		}
		else {
			throw std::logic_error ( "error in intersectTreeTraversal, side test result error");
		}
	}
	return false;
}

// Recursive depth first search intersect test ray with triangle
bool intersectTreeTraversalNoAccel(InterceptInfo& info, BSP_Node* current) {
	double distanceToIntersection = 0;

	//leaf node, recursion base case
	if ( current != NULL && current->left == NULL && current->right == NULL) {
		distanceToIntersection = triangle_intersection(info.ray->origin, info.ray->direction, current->triNormal, current->triVert1, current->triVert2, current->triVert3, info.u, info.v);
		if (distanceToIntersection > 0.0 && info.hitinfo->distance > distanceToIntersection) {
			updateHitInfo(distanceToIntersection, info, current);
			return true;
		}
	}
	else if ( current != NULL) {

		intersectTreeTraversalNoAccel(info, current->left);
		intersectTreeTraversalNoAccel(info, current->right);
		/*
		distanceToIntersection = triangle_intersection(info.ray->origin, info.ray->direction, current->triNormal, current->triVert1, current->triVert2, current->triVert3, info.u, info.v);

		if (distanceToIntersection > 0.0 && info.hitinfo->distance > distanceToIntersection) {
			updateHitInfo(distanceToIntersection, info, current);
			return true;
		}
		*/
	}
	return false;
}

bool Object::Intersect (const Ray & ray, HitInfo & hitinfo ) const {

	InterceptInfo info ;
	Scene thisScene;
	info.found = false;
	info.root = this->root;
	info.tempScene = thisScene;
	info.ray = &ray ;
	info.hitinfo = &hitinfo ;

	//intersectTreeTraversal (info, info.root ) ;
	intersectTreeTraversalNoAccel(info, info.root);

	return info.found;
}

