
#include <vector>
#include "ray_tracer.h"
#include "util.h"
#include "params.h"

inline Vec3 getInterpolatedNormal ( const HitInfo & hitinfo ) {
	return hitinfo.BayCentricU * hitinfo.tri->vertNormal2 + hitinfo.BayCentricV * hitinfo.tri->vertNormal3 + ( 1 - hitinfo.BayCentricU - hitinfo.BayCentricV ) * hitinfo.tri->vertNormal1 ;
	//return hitinfo.BayCentricU * hitinfo.tri->vertNormal1 + hitinfo.BayCentricV * hitinfo.tri->vertNormal3 + ( 1 - hitinfo.BayCentricU - hitinfo.BayCentricV ) * hitinfo.tri->vertNormal2 ;
}

Color Shader::generateSoftShadows( const Scene &scene, const HitInfo &hit, Color &color ) const
{
	const PrimitiveObject *obj = scene.GetLight(0);
	
	AABB lightBox = GetBox( *obj ); //3d box

	Vec3 temp;
	Color colorSoftShadow ( 0 , 0 , 0 );

	//COMMENTED OUT
	/*
	for ( int i = 0 ; spherePoints[99].z == 0 ; ++ i )
	{
		
		spherePoints[i] =  Vec3 ( rand ( lightBox.X.min, lightBox.X.max ) , 
									   rand ( lightBox.Y.min, lightBox.Y.max ) ,
									   rand ( lightBox.Z.min, lightBox.Z.max )) ;
	}
	*/

	
	double lightContributors = 0.0;

	//COMMENTED OUT
	/*

	for ( int i = 0 ; i < spherePoints.size() ; ++ i )
	{
		Vec3   P = hit.point; //point where ray hit object

		HitInfo hitinfoshadow;             // Holds info to pass to shader.
		//hitinfoshadow.ignore = obj;       // Don't ignore any objects.
		hitinfoshadow.ignore = obj;       // Don't ignore any objects.
		hitinfoshadow.distance = Infinity; // Follow the full ray.e point on the object back to the light source

		//if we hit an object along the way to the light, this should be a shadow
		Ray shadowray;
		shadowray.origin     = P;     // All shadow rays originate from the point on the shape.
		shadowray.type       = shadow_ray; // These rays are given no special meaning.
		//shadowray.from = hit.object;
		shadowray.generation = 1;           // Rays cast from the eye are first-generation.

		shadowray.direction = Unit ( spherePoints [ i ] - P ) ;
		shadowray.origin = P + (shadowray.direction * 0.05) ;

		if ( scene.Cast ( shadowray, hitinfoshadow ) == true && hit.object != hitinfoshadow.object ) 
		{
			//color = color / 6 ;
			lightContributors++;
		}

	}
	*/

	lightContributors = lightContributors / 100.0; 

	if ( lightContributors > 0.0 )
	{
		//if ( lightContributors <= ( 1.0/21.0) ) lightContributors = 1.0/19.0;
		//color = color *  ( 1.0 / ( 20.0 *  lightContributors ));
		color = color * ( 1 - lightContributors );
	}

	return color;
}

Color Shader::generateNormalShadows( const Scene &scene, const HitInfo &hit, Color &color ) const{
	
	Vec3   P = hit.point; // Point where ray hit object in the scene
	double colorDivisor = 2.0;
	
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
	
	shadowray.origin = P + (shadowray.direction * 0.001) ;

	// Determine whether the ray cast to the light has intersected with a face facing the ray.
	 if ( scene.Cast ( shadowray, hitinfoshadow ) && shadowray.direction * hitinfoshadow.normal <= 0 ) // && hit.object != hitinfoshadow.object ) //&& !hit.object->Inside(shadowray.origin))// && hit.object != hitinfoshadow.object ) 
	 {
		//color.blue = 0;
		//color.red = 0 ;
		//color.green = .9;
		 color = color / colorDivisor ;
	 }
	 return color;
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
	

	if ( abs ( hit.normal.z ) > abs ( hit.normal.x ) )
	{

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
	Color  color    = mat->ambient * diffuse; 
	Color  reflective = mat->reflectivity;
	Color  transl   = mat->translucency ;
	Vec3   O = hit.ray.origin; //ray origin
	Vec3   P = hit.point; //point where ray hit object
	//Vec3   N = hit.normal; //normal of the surface
	Vec3   N = getInterpolatedNormal ( hit ) ; // Interpolated normal of the surface
	Vec3   E = Unit( O - P );  // direction of the ray
	Vec3   R = Unit( ( 2.0 * ( E * N ) ) * N - E );
	
	double e = mat->Phong_exp;
	double k = mat->ref_index;

	double diffuse_intensity = .2;

	double occlusion = 0;

	//if( E * N < 0.0 ) N = -N;  // Flip the normal if necessary.

	color.blue = 0;
	color.red = 0;
	color.green = 0;

	//occlusion =  findOcclusion ( scene, hit, color ) ;

	//return color;

	for( unsigned i = 0; i < scene.NumLights(); i++ ){
		const PrimitiveObject *light = scene.GetLight(i);
		Color emission = light->material->emission;
		AABB box = GetBox( *light ); //3d box
		Vec3 LightPos( Center( box ) ); 

		// AMBIET + DIFFUSE
		L = Unit ( LightPos - P) ;

		//color += emission * ( (occlusion*   ambient) + ( diffuse * max ( 0, N * L) ) ) ;
		color += emission * ( ( ambient) + ( diffuse * max ( 0, N * L) ) ) ;

		//SPECULATIVE
		RR = (-1.0 * L) + (2.0 * ( L * N ) * N);

		double temp1 = max ( 0, E * RR ) ;
		double temp2 = pow ( max ( 0, E * RR ) , e);

		color += emission * ( specular * ( pow ( max ( 0.0, E * RR ) , e) ) ) ;
	}

	//generateSoftShadows ( scene, hit, color );
	generateNormalShadows ( scene, hit, color );
	return color;
	
	Ray reflectiveRay;
	
	reflectiveRay.generation = hit.ray.generation + 1 ;

	R =  (- 1.0 * E) + (2.0 * ( E * N ) * N);
	
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

			dnRefract = (-1 * E) * ( -1 * N) ; 
			cosRefract = (  ( -1 * E ) - ( ( -1 * N) * dnRefract )   ) / k ;
			cosSqRefract = 1 - ((1 -  ( dnRefract * dnRefract  )) / ( k * k));

			if ( cosSqRefract < 0.0 ) {
				transLRay.origin = P + ( R * Epsilon );
				transLRay.direction = R ;
				transLRay.generation = hit.ray.generation + 1;
				transLRay.type = indirect_ray ;

				return ( transl * scene.Trace (transLRay)); 
			}
			else {
				tRefr = cosRefract - ( ( -1 * N) * sqrt ( cosSqRefract ) ) ;
				refractionCoef = tRefr * N ;
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
		//color = color * occlusion ;
		return ( color + ( reflective * scene.Trace (reflectiveRay ) )  );
	}
}
