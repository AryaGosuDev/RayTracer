
#include <vector>
#include "ray_tracer.h"
#include "radiosity_helper.h"
#include "util.h"
#include "params.h"
#include <random>

inline Vec3 getInterpolatedNormal ( const HitInfo & hitinfo ) {
	return hitinfo.BayCentricU * hitinfo.tri->vertNormal2 + hitinfo.BayCentricV * hitinfo.tri->vertNormal3 + ( 1 - hitinfo.BayCentricU - hitinfo.BayCentricV ) * hitinfo.tri->vertNormal1 ;
}

double ToneMapRadiosity( const Color &color )
{
    if ( color == Green ) return 1.0;

	//constants
	static const double LUMINANCE_DISPLAY = 400.0 ;
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
	//channel finalRGB = (channel)( col  >= 255 ? 255 : col  );

	return enumerator;
}

Color Shader::generateSoftShadows( const Scene &scene, const HitInfo &hit, Color &color ) const
{
	const PrimitiveObject *obj = scene.GetLight(0);
	const SphereLight *Sobj    = dynamic_cast<SphereLight*>(const_cast<PrimitiveObject *>(obj)) ;

	AABB lightBox = GetBox( *obj ); //3d box

	Vec3 temp;
	Color colorSoftShadow ( 0 , 0 , 0 );

	std::vector<Vec3> spherePoints;

	Vec3 LightCenter = Center ( lightBox ) ;

	const int spherePointsSize = 100 ;
	int divisionsInX = spherePointsSize;
	int divisionsInY = 1;
	//int radiusOfSphere = Sobj->radius ;

	while ( divisionsInX % 2 == 0 ) {
		divisionsInX = divisionsInX / 2 ;
		divisionsInY *= 2 ;
	}

	//create distribution of lights on a rectangular plane cutting through the great circle of the sphere light

	for ( int i = 0 ; i < divisionsInX ; ++ i ) {
		for ( int j = 0 ; j < divisionsInY  ; ++ j ) {
				spherePoints.push_back(Vec3(LightCenter.x - 1 + (( j * 2 + 1 ) /  (double)divisionsInY)  , LightCenter.y - 1 + ((i * 2 + 1) / (double)divisionsInX), LightCenter.z)) ;
		}
	}

	int lightContributors = 0;

	for ( int i = 0 ; i < spherePoints.size() ; ++ i ){
		Vec3 P = hit.point; //point where ray hit object

		HitInfo hitinfoshadow;             // Holds info to pass to shader.
		//hitinfoshadow.ignore = obj;       // Don't ignore any objects.
		hitinfoshadow.ignore = NULL;       // Don't ignore any objects.
		hitinfoshadow.distance = Infinity; // Follow the full ray.e point on the object back to the light source

		//if we hit an object along the way to the light, this should be a shadow
		Ray shadowray;
		shadowray.origin     = P;     // All shadow rays originate from the point on the shape.
		shadowray.type       = shadow_ray; // These rays are given no special meaning.
		shadowray.from       = hit.object;
		shadowray.generation = 1;           // Rays cast from the eye are first-generation.

		shadowray.direction = Unit ( spherePoints [ i ] - P ) ;
		shadowray.origin = P + (shadowray.direction * 0.05) ;

		if ( scene.Cast ( shadowray, hitinfoshadow ) && hit.object != hitinfoshadow.object ) {
			lightContributors++;
		}

	}

	if ( lightContributors > 0 )
		color = color * ( 1.0 - ((double)lightContributors)/((double)spherePoints.size()) );

	return color;
}

Color Shader::generateNormalShadows( const Scene &scene, const HitInfo &hit, Color & color ) const{
	
	Vec3   P = hit.point; // Point where ray hit object in the scene
	double colorDivisor = 1.3;
	
	//make a ray from the object            
	HitInfo hitinfoshadow;             // Holds info to pass to shader.
	hitinfoshadow.ignore = NULL;       // Don't ignore any objects.
	hitinfoshadow.distance = Infinity; // Follow the full ray.e point on the object back to the light source

	// If we hit an object along the way to the light, this should be a shadow
	Ray shadowray;
	shadowray.origin     = P;     // All shadow rays originate from the point on the shape.
	shadowray.type       = shadow_ray; // Shadow ray
	shadowray.from = hit.object;
	const PrimitiveObject *light = scene.GetLight(0);
	AABB box = GetBox( *light );
	shadowray.direction = Unit (Center(box) - P);

	Vec3 centerEx = Center(box);
	
	//shadowray.origin = P + (shadowray.direction * 0.001) ;
	shadowray.origin = P + (hit.normal * 0.001) ;

	// Determine whether the ray cast to the light has intersected with a face facing the ray.
	 //if ( scene.Cast ( shadowray, hitinfoshadow ) && shadowray.direction * hitinfoshadow.normal <= 0 ) // && hit.object != hitinfoshadow.object ) //&& !hit.object->Inside(shadowray.origin))// && hit.object != hitinfoshadow.object ) 
	 if ( scene.Cast ( shadowray, hitinfoshadow ) ) //&& hit.object != hitinfoshadow.object )
     {
		 //color = color / colorDivisor ;
	    //color *= 0.1;
	     color = Black ;

	 }
	 return color;
}

double Shader::findOcclusionNew ( const Scene & scene, const HitInfo & hit , Color & color, double radiusOfHemisphere, int AORaysToGenerate ) const {

	if ( AORaysToGenerate == 0 ) return 1;

	int occlusions = 0;
	Vec3 newAxesZ, newAxesX, newAxesY, t ;

	//rng
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.0, 2.0);

	//define new axes based on the tangent of the hit surface normal
	newAxesZ = Unit (hit.normal ) ;
	t =  newAxesZ;

	double smallestComponentZ = abs ( newAxesZ.z ),
		   smallestComponentX = abs ( newAxesZ.x ),
		   smallestComponentY = abs ( newAxesZ.y );

	if ( smallestComponentZ <= smallestComponentX ) {
		if ( smallestComponentZ <= smallestComponentY )
			t.z = 1;
	}
	else if ( smallestComponentX <= smallestComponentZ ) {
		if ( smallestComponentX <= smallestComponentY ) 
			t.x = 1 ;
	}
	else t.y = 1 ;

	newAxesX = Unit ( t ^ newAxesZ ) ;
	newAxesY = newAxesZ ^ newAxesX ; 

	//Unit Square to unit disc code

	for ( int AORayIteration = 0 ; AORayIteration < AORaysToGenerate ; ++AORayIteration ) {

		Vec2 squarePoint ( dis(gen) - 1, dis(gen) - 1 );
		double phi, r, u, v ;
		double a = squarePoint.x ;
		double b = squarePoint.y;

		if ( a > -b ){
			if ( a > b ) {
				r = a ;
				phi = ( Pi / 4.0 ) * ( b / a ) ;
			}
			else {
				r = b ;
				phi = ( Pi / 4.0 ) * ( 2.0 - (a/b));
			}
		}
		else {
			if ( a < b ) {
				r = -a ;
				phi = (Pi / 4.0 ) * ( 4.0 + ( b / a ) );
			}
			else {
				r = -b ;
				if ( b != 0 )
					phi = (Pi / 4.0 ) * ( 6.0 - (a/b));
				else
					phi = 0 ;
			}
		}

		u = r * cos ( phi ) ;
		v = r * sin ( phi ) ;

		// finding z intercept of the ambient occlusion hemisphere
		Vec3 VectorComponentX = u * newAxesX, VectorComponentY = v * newAxesY, VectorComponentZ  ;
		Vec3 VectorBaseZ = VectorComponentX + VectorComponentY ;    // the base vector of the ray that intercepts the hemisphere
		double rayLegLength = (( radiusOfHemisphere * radiusOfHemisphere ) - ( Length( VectorBaseZ ) * Length( VectorBaseZ ) ));

		Vec3 RayFromCenterOfHemisphere = ((rayLegLength / radiusOfHemisphere) * newAxesZ ) + VectorBaseZ ;

		//Shoot rays from the hemisphere in all directions and find a hit
		Ray AOray;
		AOray.direction = Unit(RayFromCenterOfHemisphere );
		AOray.generation = 1;
		AOray.origin = hit.point + ( AOray.direction * .0001 ) ;

		HitInfo hitinfoAO;             // Holds info to pass to shader.
		hitinfoAO.ignore = NULL;       // Don't ignore any objects.
		hitinfoAO.distance = Infinity; // Follow the full ray.e point on the object back to the light source

		if ( scene.Cast ( AOray, hitinfoAO) && hitinfoAO.distance < radiusOfHemisphere ) {
			occlusions ++;
		}
	}

	return  (1.0 - (occlusions / (double)AORaysToGenerate));
}

double Shader::findOcclusion( const Scene &scene, const HitInfo &hit, Color &color ) const
{
	//color.red = 1;
	//color.green = 1;
	//color.blue = 1;

	double occlusions = 0;

	double x = 0, y = 0 , z = 0;
	double numberOfDivisions = 50;
	double distanceAO = 0.5;
	double totalRayCount = 0;

	Vec3 normal1 ;
	
	normal1.x = 0;
	normal1.y = 0;
	normal1.z = -1;

	double startingTheta = acos ( hit.normal.z /  sqrt ( hit.normal.x * hit.normal.x + hit.normal.y * hit.normal.y + hit.normal.z * hit.normal.z ));
	double startingPhi = atan2 ( hit.normal.y, hit.normal.x );
	//double startingTheta = acos ( normal1.z /  sqrt ( normal1.x * normal1.x + normal1.y * normal1.y + normal1.z * normal1.z ));
	//double startingPhi = atan2 ( normal1.y, normal1.x );

	double angleRange = 0;
	

	if ( abs ( hit.normal.z ) > abs ( hit.normal.x ) ){

		angleRange =  Pi / numberOfDivisions  ; // the entire range
		double theta = startingTheta + (Pi / 2.0) - angleRange ;
		double phi = angleRange;

		for ( int outerTheta = 0; outerTheta < (numberOfDivisions - 1) ; ++ outerTheta)
		{ 
			phi = angleRange;

			for ( int innerPhi = 0 ; innerPhi < (numberOfDivisions - 1) ; ++innerPhi )
			{
				totalRayCount++;
				x = hit.point.x + cos ( phi ) * sin ( theta );
				y = hit.point.y + sin ( phi ) * sin ( theta );
				z = hit.point.z + cos ( theta );

				Ray AOray;
				AOray.direction = Unit( Vec3 ( x, y, z ) - hit.point );
				AOray.generation = 1;
				AOray.origin = hit.point + ( AOray.direction * .0001 ) ;

				HitInfo hitinfoAO;             // Holds info to pass to shader.
				//hitinfoshadow.ignore = obj;       // Don't ignore any objects.
				hitinfoAO.ignore = NULL;       // Don't ignore any objects.
				hitinfoAO.distance = Infinity; // Follow the full ray.e point on the object back to the light source

				phi +=  angleRange;


				if ( scene.Cast ( AOray, hitinfoAO) == true && hitinfoAO.distance < distanceAO )
				{
					occlusions ++;
				}
				//else if ( scene.IsAlreadyInside ( AOray.origin ,  hitinfoAO ) )
				//{
					//occlusions ++;

				//}
				
			}

			phi = -angleRange ;

			for ( int innerPhi = 0 ; innerPhi < (numberOfDivisions - 1) ; ++innerPhi )
			{
				totalRayCount++;
				x = hit.point.x + cos ( phi ) * sin ( theta );
				y = hit.point.y + sin ( phi ) * sin ( theta );
				z = hit.point.z + cos ( theta );

				Ray AOray;
				AOray.direction = Unit( Vec3 ( x, y, z ) - hit.point );
				AOray.generation = 1;
				AOray.origin = hit.point + ( AOray.direction * .0001 ) ;

				HitInfo hitinfoAO;             // Holds info to pass to shader.
				//hitinfoshadow.ignore = obj;       // Don't ignore any objects.
				hitinfoAO.ignore = NULL;       // Don't ignore any objects.
				hitinfoAO.distance = Infinity; // Follow the full ray.e point on the object back to the light source

				phi -=  angleRange;
				if ( scene.Cast ( AOray, hitinfoAO) == true && hitinfoAO.distance < distanceAO )
				{
					occlusions ++;
				}
				//else if ( scene.IsAlreadyInside ( AOray.origin ,  hitinfoAO ) )
				//{
					//occlusions ++;

				//}
			}
			theta -= angleRange;
		}
		//color = color *  (1 - (occlusions / totalRayCount))


	}
	else
	{
		 
		angleRange =  Pi / numberOfDivisions  ; // the entire range
		double theta = angleRange;
		double phi = (startingPhi - (Pi/2.0)) + angleRange;
		
		//for ( double theta = angleRange ; theta < Pi ; theta += angleRange )
		for ( int outerTheta = 0; outerTheta < (numberOfDivisions - 1) ; ++ outerTheta)
		{
			phi = (startingPhi - (Pi/2.0)) + angleRange;
			
			//for ( double phi = (startingPhi - (Pi/2.0)) + angleRange ; phi < (startingPhi + (Pi/2.0 )) ; phi +=  angleRange )
			for ( int innerPhi = 0 ; innerPhi < (numberOfDivisions - 1) ; ++innerPhi )
			{
				totalRayCount++;
				x = hit.point.x + cos ( phi ) * sin ( theta );
				y = hit.point.y + sin ( phi ) * sin ( theta );
				z = hit.point.z + cos ( theta );

				Ray AOray;
				AOray.direction = Unit( Vec3 ( x, y, z ) - hit.point );
				AOray.generation = 1;
				AOray.origin = hit.point + ( AOray.direction * .0001 ) ;

				HitInfo hitinfoAO;             // Holds info to pass to shader.
				//hitinfoshadow.ignore = obj;       // Don't ignore any objects.
				hitinfoAO.ignore = NULL;       // Don't ignore any objects.
				hitinfoAO.distance = Infinity; // Follow the full ray.e point on the object back to the light source


				phi +=  angleRange;
				if ( scene.Cast ( AOray, hitinfoAO) == true && hitinfoAO.distance < distanceAO )
				{
					occlusions ++;
				}
				//else if ( scene.IsAlreadyInside ( AOray.origin ,  hitinfoAO ) )
				//{
					//occlusions ++;

				//}
				
			}
			theta += angleRange;
			

		}
		//color = color *  (1 - (occlusions / totalRayCount));
	}
	if ( totalRayCount == 0 ) return 1;
	return  (1.0 - (occlusions / totalRayCount));
	//return color;

}

// Hit an object, now we are going to paint that point on the object called from scene
Color Shader::Shade( const Scene &scene, const HitInfo &hit ) const {
		
	Vec3 L ; // light direction
	Vec3 RR;
	HitInfo otherhit;
	
	//TODO : WE HAVE TO MAKE SURE IF WE CAN HIT A LIGHT SOURCE, DRAW THE LIGHT SOURCE
	//if( Emitter( hit.object ) ) return hit.object->material->emission;
	
	Material *mat   = hit.object->material;
	Color  diffuse  = mat->diffuse;
	Color  specular = mat->specular;
	Color  ambient  = mat->ambient;
	//Color  color    = mat->ambient * diffuse; 
	Color  color    = diffuse;
	Color  reflective = mat->reflectivity;
	Color  transl   = mat->translucency ;
	Vec3   O = hit.ray.origin; //ray origin
	Vec3   P = hit.point; //point where ray hit object
	Vec3   SurfaceNonInterpolatedN = hit.normal; //normal of the surface
	Vec3   N = Unit (getInterpolatedNormal ( hit )) ; // Interpolated normal of the surface
	Vec3   E = Unit( O - P );  // direction of the ray
	Vec3   R = Unit( ( 2.0 * ( E * N ) ) * N - E );
	
	double e = mat->Phong_exp;
	double k = mat->ref_index;

	double diffuse_intensity = .2;

	//if( E * N < 0.0 ) N = -N;  // Flip the normal if necessary.

	color.blue = 0;
	color.red = 0;
	color.green = 0;

	if ( scene.radiosity != NULL ) {

		double radiosityAmbient = ToneMapRadiosity (scene.radiosity->radiosityHelper->trace_ray( hit.ray, *scene.radiosity->tempQuadVector ));
		double occlusion = findOcclusionNew ( scene, hit, color, 2.0, 100 ) ;

		HitInfo hitinfoOcclusion;             // Holds info to pass to shader.
		hitinfoOcclusion.ignore = NULL;       // Don't ignore any objects.
		hitinfoOcclusion.distance = Infinity; // Follow the full ray.e point on the object back to the light source

		// If we hit an object along the way to the light, this should be a shadow
		Ray Occlusionray;
		Occlusionray.origin     = P;     // All shadow rays originate from the point on the shape.
		Occlusionray.type       = shadow_ray; // Shadow ray
		Occlusionray.from = hit.object;
		
		for( unsigned i = 0; i < scene.NumLights(); i++ ){
			const PrimitiveObject *light = scene.GetLight(i);
			Color emission = light->material->emission;
			AABB box = GetBox( *light ); //3d box
			Vec3 LightPos( Center( box ) ); 

			// AMBIET + DIFFUSE
			L = Unit ( LightPos - P) ;

			Occlusionray.direction = L;
			Occlusionray.origin = P + (N * 0.001) ;
			if ( scene.Cast ( Occlusionray, hitinfoOcclusion ) ) {
				color += emission * ( occlusion * ( diffuse * .3 * radiosityAmbient ) ) ;
				//color += emission * ( ( diffuse * .3 * radiosityAmbient ) ) ;
			}
			else { 
				color += emission * ( occlusion * ( diffuse * radiosityAmbient ) ) ;
			}

			//SPECULATIVE
			RR = (-1.0 * L) + (2.0 * ( L * N ) * N);

			double temp1 = max ( 0, E * RR ) ;
			double temp2 = pow ( max ( 0, E * RR ) , e);

			color += emission * ( specular * ( pow ( max ( 0.0, E * RR ) , e) ) ) ;
		}

		//generateSoftShadows ( scene, hit, color );
		//generateNormalShadows ( scene, hit, color );
		return color ;
	}
	else {

		double occlusion = findOcclusionNew ( scene, hit, color, 2.0, 10 ) ;
	
		//return color;

		for( unsigned i = 0; i < scene.NumLights(); i++ ){
			const PrimitiveObject *light = scene.GetLight(i);
			Color emission = light->material->emission;
			AABB box = GetBox( *light ); //3d box
			Vec3 LightPos( Center( box ) ); 

			// AMBIET + DIFFUSE
			L = Unit ( LightPos - P) ;

			color += emission * ( (  (.85 *occlusion) * ambient) + ( diffuse * max ( 0, N * L) ) ) ;
			//color += emission * ( ( ambient) + ( diffuse * max ( 0, N * L) ) ) ;

			//SPECULATIVE
			RR = (-1.0 * L) + (2.0 * ( L * N ) * N);

			double temp1 = max ( 0, E * RR ) ;
			double temp2 = pow ( max ( 0, E * RR ) , e);

			color += emission * ( specular * ( pow ( max ( 0.0, E * RR ) , e) ) ) ;
		}

		//generateSoftShadows ( scene, hit, color );
		generateNormalShadows ( scene, hit, color );
		//return color;
	}
	
	/*****************  REFLECTION CODE ***********/
	Ray reflectiveRay;
	
	reflectiveRay.generation = hit.ray.generation + 1 ;

	//R =  (- 1.0 * E) + (2.0 * ( E * N ) * N);
	R =  (- 1.0 * E) + (2.0 * ( E * SurfaceNonInterpolatedN ) * SurfaceNonInterpolatedN);
	
	reflectiveRay.direction = R;
	reflectiveRay.origin = P + ( reflectiveRay.direction * Epsilon)   ;
	reflectiveRay.type = indirect_ray;
	reflectiveRay.from = hit.object;

	//return color;

	/*****************  REFRACTION CODE ***********/
	if ( transl.blue > 0.0 && transl.green > 0.0 && transl.red > 0.0 ) {
		Ray transLRay ;
		double refractionCoef = 0;
		Color beerRefract;
		double invK = 1.0 / k ;
		double dnRefract;
		Vec3 cosRefract;
		double cosSqRefract;
		Vec3 tRefr;

		double beerRefractr, beerRefractg, beerRefractb ;
		/*
		if( E * N < 0.0 )
			R =  (- 1.0 * E) + (2.0 * ( E * ( -1 * N )) * ( -1 * N));
		else 
			R =  (- 1.0 * E) + (2.0 * ( E * N ) * N);

		if ( (( -1 * E ) * N ) < 0.0 ) {
			dnRefract = (-1 * E) * N ; 
			cosRefract = ( ( -1 * E ) - ( N * dnRefract )   ) / k ;
			cosSqRefract = 1 - ((1 -  ( dnRefract * dnRefract )) / ( k * k));

			tRefr = cosRefract - ( N * sqrt ( cosSqRefract ) ) ;

			refractionCoef = E * N ;
			*/

		// **** USE NON INTERPOLATED NORMAL ***

		if( E * SurfaceNonInterpolatedN < 0.0 )
			R =  (- 1.0 * E) + (2.0 * ( E * ( -1 * SurfaceNonInterpolatedN )) * ( -1 * SurfaceNonInterpolatedN));
		else 
			R =  (- 1.0 * E) + (2.0 * ( E * SurfaceNonInterpolatedN ) * SurfaceNonInterpolatedN);

		if ( (( -1 * E ) * SurfaceNonInterpolatedN ) < 0.0 ) {
			dnRefract = (-1 * E) * SurfaceNonInterpolatedN ; 
			cosRefract = ( ( -1 * E ) - ( SurfaceNonInterpolatedN * dnRefract )   ) / k ;
			cosSqRefract = 1 - ((1 -  ( dnRefract * dnRefract )) / ( k * k));

			tRefr = cosRefract - ( SurfaceNonInterpolatedN * sqrt ( cosSqRefract ) ) ;

			refractionCoef = E * SurfaceNonInterpolatedN ;

			transl.blue = 1;
			transl.green = 1;
			transl.red = 1;
		}
		else {
			beerRefractr = exp ( -1.0 * transl.red * hit.distance ); 
			beerRefractg = exp ( -1.0 * transl.green * hit.distance );
			beerRefractb = exp ( -1.0 * transl.blue * hit.distance );

			transl.red = beerRefractr;
			transl.green = beerRefractg;
			transl.blue = beerRefractb;

			//dnRefract = (-1 * E) * ( -1 * N) ; 
			//cosRefract = (  ( -1 * E ) - ( ( -1 * N) * dnRefract )   ) / k ;
			//cosSqRefract = 1 - ((1 -  ( dnRefract * dnRefract  )) / ( k * k));

			// **** USE NON INTERPOLATED NORMAL ***

			dnRefract = (-1 * E) * ( -1 * SurfaceNonInterpolatedN) ; 
			cosRefract = (  ( -1 * E ) - ( ( -1 * SurfaceNonInterpolatedN) * dnRefract )   ) / k ;
			cosSqRefract = 1 - ((1 -  ( dnRefract * dnRefract  )) / ( k * k));

			if ( cosSqRefract < 0.0 ) {
				transLRay.origin = P + ( R * Epsilon );
				transLRay.direction = R ;
				transLRay.generation = hit.ray.generation + 1;
				transLRay.type = indirect_ray ;

				return ( transl * scene.Trace (transLRay)); 
			}
			else {
				//tRefr = cosRefract - ( ( -1 * N) * sqrt ( cosSqRefract ) ) ;
				//refractionCoef = tRefr * N ;

				// **** USE NON INTERPOLATED NORMAL ***

				tRefr = cosRefract - ( ( -1 * SurfaceNonInterpolatedN) * sqrt ( cosSqRefract ) ) ;
				refractionCoef = tRefr * SurfaceNonInterpolatedN ;
			}
		}

		double beerR0 = (( 1 - k ) * ( 1 - k )) / (( 1 + k ) * ( 1 + k ));
		double beerR = beerR0 + ( 1 - beerR0 ) * pow ( 1 - refractionCoef , 5 ) ;

		reflectiveRay.direction = R;
		reflectiveRay.origin = P + ( R * Epsilon)   ;
		reflectiveRay.type = indirect_ray;
		reflectiveRay.generation = hit.ray.generation + 1;
	
		transLRay.origin     = P  + (tRefr * Epsilon );
		transLRay.direction  = tRefr ;
		transLRay.generation = hit.ray.generation + 1;
		transLRay.type = special_ray ;

		return ( transl * ( ( beerR * scene.Trace (reflectiveRay)) + (( 1 - beerR ) * scene.Trace ( transLRay  )))) ;
	}
	else {
		return ( color + ( reflective * scene.Trace (reflectiveRay)));
	}
}
