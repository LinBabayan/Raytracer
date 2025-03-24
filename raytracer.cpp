#include "image.h"
#include "kdtree.h"
#include "ray.h"
#include "raytracer.h"
#include "scene_types.h"
#include <stdio.h>

#include <glm/gtc/epsilon.hpp>

/// acne_eps is a small constant used to prevent acne when computing
/// intersection
//  or boucing (add this amount to the position before casting a new ray !
const float acne_eps = 1e-4;

bool intersectPlane(Ray *ray, Intersection *intersection, Object *obj) {

  //! \todo : compute intersection of the ray and the plane object
  
  // Get the plane's normal and distance from origin
    vec3 normal = obj->geom.plane.normal;
    float distance = obj->geom.plane.dist;

    // Compute the denominator 
    float denom = dot(ray->dir, normal);

    // If denominator is 0, the ray is parallel to the plane (no intersection)
    if (denom == 0) {
        return false;
    }

    // Calculate t(parameter of the ray at the intersection point)
    float t = -(dot(ray->orig, normal) + distance) / denom;

    // If t is outside the valid range [tmin, tmax], no intersection
    if (t < ray->tmin || t > ray->tmax) {
        return false;
    }

    // Calculate the intersection point
    intersection->position = ray->orig + t * ray->dir;

    intersection->normal = normal;

    intersection->mat = &obj->mat;

    return true;

  
}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj) {

//! \todo : compute intersection of the ray and the sphere object

// Get the sphere's center and radius
    vec3 center = obj->geom.sphere.center;
    float radius = obj->geom.sphere.radius;

    // Vector from ray origin to sphere center
    vec3 oc = ray->orig - center;

    // Calculate the coefficients of the quadratic equation
    float a = dot(ray->dir, ray->dir);  // dot product of ray direction with itself
    float b = 2.0f * dot(oc, ray->dir);  // 2 * dot product of (ray origin - sphere center) with ray direction
    float c = dot(oc, oc) - radius * radius;  // dot product of (ray origin - sphere center) with itself - radius^2

    // Calculate the discriminant
    float discriminant = b * b - 4.0f * a * c;

    // If the discriminant is negative, no intersection
    if (discriminant < 0) {
        return false;
    }

    // Calculate the two possible values of t (solutions to the quadratic equation)
    float sqrtDiscriminant = sqrt(discriminant);
    float t1 = (-b - sqrtDiscriminant) / (2.0f * a);
    float t2 = (-b + sqrtDiscriminant) / (2.0f * a);

    // Pick the smallest positive t (the closest intersection)
    float t = -1.0f;  // Default to an invalid t value

    if (t1 > 0 && t2 > 0) {
        t = (t1 < t2) ? t1 : t2;  // Choose the smaller positive t
    } else if (t1 > 0) {
        t = t1;  // Choose t1 if t2 is negative
    } else if (t2 > 0) {
        t = t2;  // Choose t2 if t1 is negative
    }

    // If no valid t is found, return false
    if (t < ray->tmin || t > ray->tmax) {
        return false;
    }

    // Calculate the intersection position
    intersection->position = ray->orig + t * ray->dir;

    // Calculate the normal at the intersection point
    intersection->normal = normalize(intersection->position - center);  // Normal is the direction from center to intersection

    // Assign the material of the object to the intersection
    intersection->mat = &obj->mat;

    return true;


}

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {
  
//!\todo loop on each object of the scene to compute intersection

  bool hasIntersection = false;
    size_t objectCount = scene->objects.size();
    float closestT = ray->tmax;  // Initialize closest intersection distance

    // Loop through all objects in the scene
    for (size_t i = 0; i < objectCount; i++) {
        Object *obj = scene->objects[i];  // Get current object
        Intersection tempIntersection;  // Temporary intersection storage

        bool hit = false;

        // Check object type and call the appropriate intersection function
        if (obj->geom.type == SPHERE) {
            hit = intersectSphere(ray, &tempIntersection, obj);
        } else if (obj->geom.type == PLANE) {
            hit = intersectPlane(ray, &tempIntersection, obj);
        }

        // If an intersection was found and is the closest so far, update the main intersection
        if (hit && tempIntersection.position != ray->orig) {
            float t = length(tempIntersection.position - ray->orig);
            if (t < closestT) {
                closestT = t;
                *intersection = tempIntersection;
                hasIntersection = true;
            }
        }
    }

    return hasIntersection;
}

/*
bool intersect_kd_tree(const Ray& ray, KdTreeNode* node) {
    if (!ray_intersects_box(ray, node->box)) {
        return false;  
    }

    bool hit = false;
    for (Object* obj : node->objects) {
        hit |= intersect(ray, *obj);  // Test intersection with each object in the node
    }

    if (node->left) hit |= intersect_kd_tree(ray, node->left);  // Check left child
    if (node->right) hit |= intersect_kd_tree(ray, node->right);  // Check right child
    return hit;
}

*/
/* ---------------------------------------------------------------------------
 */
/*
 *	The following functions are coded from Cook-Torrance bsdf model
 *description and are suitable only
 *  for rough dielectrics material (RDM. Code has been validated with Mitsuba
 *renderer)
 */

// Shadowing and masking function. Linked with the NDF. Here, Smith function,
// suitable for Beckmann NDF

float RDM_chiplus(float c) { return (c > 0.f) ? 1.f : 0.f; }

/** Normal Distribution Function : Beckmann
 * NdotH : Norm . Half
 */
 


float RDM_Beckmann(float NdotH, float alpha) {


  //! \todo compute Beckmann normal distribution
	float res = 0;

	if(NdotH > 0){
		float tan = (1-NdotH*NdotH)/(NdotH*NdotH);
		tan = -tan/(alpha*alpha);
		float dem = M_PI*alpha*alpha*NdotH*NdotH*NdotH*NdotH;
		res = pow(M_E,tan);
		res = res/dem;
	}
	
	return res;

}

// Fresnel term computation. Implantation of the exact computation. we can use
// the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR) {

	//! \todo compute Fresnel term
	float s,p,f;
	float ior = extIOR/intIOR;
	float sin2 = ior*ior*(1-LdotH*LdotH);
	
	if(sin2 > 1.0f){
		return 1.0f;
	}
	
	float cos = sqrt(1-sin2);
	s = powf(extIOR*LdotH-intIOR*cos,2) / (powf(extIOR*LdotH + intIOR*cos,2));
	p = powf(extIOR*cos-intIOR*LdotH,2) / (powf(extIOR*cos + intIOR*LdotH,2));
	f = (s+p)/2;
	
	return f;

}

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha) {

  //! \todo compute G1 term of the Smith fonction
	float k = DdotH/DdotN;
	
	if(k <= 0){
		return 0.0f;
	}
	
	float tan = sqrt(1-DdotN*DdotN)/DdotN;
	float b = 1/(alpha*tan);
	float res = 1;
	
	if(b < 1.6f){
		float b2 = b*b;
		res = (3.535*b + 2.181*b2)/(1+2.276*b+2.577*b2);
	}
	
	return res;

}

// LdotH : Light . Half
// LdotN : Light . Norm
// VdotH : View . Half
// VdotN : View . Norm

float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN,
                float alpha) {

  //! \todo the Smith fonction
	float G1l = RDM_G1(LdotH, LdotN,alpha);
	float G1V = RDM_G1(VdotH, VdotN,alpha);
	return G1l*G1V;

}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm


color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, Material *m) {

	//! \todo specular term of the bsdf, using D = RDB_Beckmann, F = RDM_Fresnel, G
	//! = RDM_Smith
	
	float alpha = m->roughness;
	float G = RDM_Smith(LdotH,LdotN, VdotH, VdotN,alpha);
	float D = RDM_Beckmann(NdotH, alpha);
	float F = RDM_Fresnel(LdotH,1,m->IOR);
	
	color3 res = (D*F*G*(m->specularColor))/(4*LdotN*VdotN);
	
	return res;

}


// diffuse term of the cook torrance bsdf


color3 RDM_bsdf_d(Material *m) {

  //! \todo compute diffuse component of the bsdf
	  return m->diffuseColor/3.14f;

}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                Material *m) {

	//! \todo compute bsdf diffuse and specular term
	
	color3 d = RDM_bsdf_d(m);
	color3 s = RDM_bsdf_s(LdotH,NdotH, VdotH, LdotN,VdotN,m);
	return d+s;
	return color3(0.f);

}




color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat) {
	color3 ret = color3(0.f);
	vec3 h = normalize(v+l);
	float LdotH = dot(l,h);
	float NdotH = dot(n,h);
	float VdotH = dot(v,h);
	float LdotN = dot(l,n);
	float VdotN = dot(v,n);
	
	if(LdotN < 0) {
		return ret;
	}
	//! \todo compute bsdf, return the shaded color taking into account the
	//! lightcolor

	ret = lc*RDM_bsdf(LdotH,NdotH,VdotH,LdotN,VdotN,mat) * LdotN;
	return ret;
}




//! if tree is not null, use intersectKdTree to compute the intersection instead
//! of intersect scene

color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree) {
  color3 ret = color3(0, 0, 0);
  color3 zero = color3(0, 0, 0);
  vec3 lum;
  Intersection intersection;
  if(intersectScene(scene,ray,&intersection)){
  
	vec3 v = -(ray->dir);
	size_t lightsCount = scene->lights.size();
	
	for(size_t i = 0;i < lightsCount;i++){
	
		lum = normalize(scene->lights[i]->position - intersection.position);
		Ray r;
		rayInit(&r,intersection.position+acne_eps*lum ,lum,acne_eps,length(scene->lights[i]->position - intersection.position),ray->depth);
		Intersection inter;

		if(!(intersectScene(scene,&r,&inter))){
		
			ret += shade(normalize(intersection.normal),v,lum,scene->lights[i]->color,intersection.mat);
			if(ray->depth < 10){
			
				vec3 r = reflect(ray->dir,normalize(intersection.normal));
				Ray ref;
				rayInit(&ref,intersection.position+acne_eps*normalize(intersection.normal) ,r,acne_eps,FLT_MAX,ray->depth+1);
				color3 refColor = trace_ray(scene,&ref,tree);
				float f = RDM_Fresnel(dot(lum,normalize(lum + v)),1.0f,intersection.mat->IOR);
				color3 sp = intersection.mat->specularColor;
				ret += f*refColor*sp;
			}
		}
	}
  }
  else{
  ret = scene->skyColor;
  }

  return ret;
}
/*
color3 trace_with_anti_aliasing(int x, int y, int num_samples) {
    color3 color = Color(0, 0, 0);
    for (int i = 0; i < num_samples; ++i) {
        float offsetX = random() - 0.5;  // Random offset for sub-pixel
        float offsetY = random() - 0.5;
        Ray ray = generate_ray(x + offsetX, y + offsetY);
        color += trace_ray(ray);
    }
    return color / num_samples;  // Average the results
}
*/


void renderImage(Image *img, Scene *scene) {

	float aspect = 1.f / scene->cam.aspect;

	KdTree *tree = NULL;

	float delta_y = 1.f / (img->height * 0.5f);   //! one pixel size
	vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step
	vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) *
		     aspect * scene->cam.ydir;

	float delta_x = 1.f / (img->width * 0.5f);
	vec3 dx = delta_x * scene->cam.xdir;
	vec3 ray_delta_x = (0.5f - img->width * 0.5f) / (img->width * 0.5f) * scene->cam.xdir;


	for (size_t j = 0; j < img->height; j++) {

		if (j != 0){
			printf("\033[A\r");
		}

		float progress = (float)j / img->height * 100.f;
		printf("progress\t[");
		int cpt = 0;
		for (cpt = 0; cpt < progress; cpt += 5)
			printf(".");
		for (; cpt < 100; cpt += 5)
			printf(" ");
		printf("]\n");

		#pragma omp parallel for
		for (size_t i = 0; i < img->width; i++) {
			color3 *ptr = getPixelPtr(img, i, j);
			color3 pixel_color = {0, 0, 0};

			vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
				     float(i) * dx + float(j) * dy;

			Ray rx;
			rayInit(&rx, scene->cam.position, normalize(ray_dir));

			*ptr = trace_ray(scene, &rx, tree);

		}
	}
}
