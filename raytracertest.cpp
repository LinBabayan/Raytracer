
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

  float sc = dot(ray->dir,obj->geom.plane.normal);
if(sc == 0){
return false;
}
vec3 o = vec3(ray->orig.x,ray->orig.y,ray->orig.z);
float t = -(obj->geom.plane.dist + dot(o,obj->geom.plane.normal))/sc;
if(t < ray->tmin || t > ray->tmax){
return false;
}
ray->tmax = t;
intersection->position.x = ray->orig.x + t*ray->dir.x;
intersection->position.y = ray->orig.y + t*ray->dir.y;
intersection->position.z = ray->orig.z + t*ray->dir.z;
intersection->normal = obj->geom.plane.normal;
intersection->mat = &(obj->mat);
return true;
}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj) {
vec3 oc = vec3(ray->orig.x - obj->geom.sphere.center.x ,ray->orig.y-obj->geom.sphere.center.y ,ray->orig.z -obj->geom.sphere.center.z );
float c = dot(oc,oc) - (obj->geom.sphere.radius*obj->geom.sphere.radius);
float b = 2*dot(oc,ray->dir);
float dist = b*b - 4*c;

if(dist < 0){
return false;
}
float t;
float tmin = ray->tmin;
float tmax = ray->tmax;
if (dist == 0)
{
t = -b/2;
if(t < tmin || t > tmax){
return false;
}
}
else{
float t1 = (-b + sqrt(dist))/2;
float t2 = (-b - sqrt(dist))/2;
if(t1 < tmin || t1 > tmax){
t = t2;
if(t < tmin || t > tmax){
return false;
}
}
if(t2 < tmin || t2 > tmax){
t = t1;
}
else{
if(t2 > t1){
t = t1;
}
else{
t = t2;
}
}
}

//! \todo : compute intersection of the ray and the sphere object
ray->tmax = t;
intersection->position.x = ray->orig.x + t*ray->dir.x;
intersection->position.y = ray->orig.y + t*ray->dir.y;
intersection->position.z = ray->orig.z + t*ray->dir.z;
//vec3 v = vec3(intersection->position.x-)
vec3 v = intersection->position - obj->geom.sphere.center;
intersection->normal = normalize(v);
intersection->mat = &(obj->mat);
return true;
}

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {
  bool hasIntersection = false;
  size_t objectCount = scene->objects.size();

for(size_t i=0; i < objectCount;i++){
if(scene->objects[i]->geom.type == PLANE){
if(intersectPlane(ray,intersection,(scene->objects[i]))){
hasIntersection = true;
};
}
else{
if(intersectSphere(ray,intersection,(scene->objects[i]))){
hasIntersection = true;
};
}
}
/*if(hasIntersection){
scene->skyColor.x = 0.5*intersection->normal.x + 0.5;
scene->skyColor.y = 0.5*intersection->normal.y + 0.5;
scene->skyColor.z = 0.5*intersection->normal.z + 0.5;
}*/

//!\todo loop on each object of the scene to compute intersection

return hasIntersection;
}

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
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN,float VdotN, Material *m) {
  
  float alpha = m->roughness;
  float G = RDM_Smith(LdotH,LdotN, VdotH, VdotN,alpha);
  float D = RDM_Beckmann(NdotH, alpha);
  float F = RDM_Fresnel(LdotH,1,m->IOR);
  color3 res = (D*F*G*(m->specularColor))/(4*LdotN*VdotN);
  return res;

}
// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m) {

  
  return m->diffuseColor/3.14f;

}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                Material *m) {

  color3 d = RDM_bsdf_d(m);
  color3 s = RDM_bsdf_s(LdotH,NdotH, VdotH, LdotN,VdotN,m);
  return d+s;

}



color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat) {
  color3 ret = color3(0.f);
  vec3 h = normalize(v+l);
  float LdotN = dot(l,n);
  float LdotH = dot(l,h);
  float NdotH = dot(n,h);
  float VdotH = dot(v,h);
  float VdotN = dot(v,n);
  if(LdotN < 0) {
    return ret;
  }
//! \todo compute bsdf, return the shaded color taking into account the
//! lightcolor
    //ret = mat->diffuseColor * LdotN * lc/3.14f;
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
  //ret = color3(0.5*intersection.normal.x + 0.5,0.5*intersection.normal.y + 0.5,0.5*intersection.normal.z + 0.5);
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

void renderImage(Image *img, Scene *scene) {

  //! This function is already operational, you might modify it for antialiasing
  //! and kdtree initializaion
  float aspect = 1.f / scene->cam.aspect;

  KdTree *tree = NULL;


//! \todo initialize KdTree

  float delta_y = 1.f / (img->height * 0.5f);   //! one pixel size
  vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step
  vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) *
                     aspect * scene->cam.ydir;

  float delta_x = 1.f / (img->width * 0.5f);
  vec3 dx = delta_x * scene->cam.xdir;
  vec3 ray_delta_x =
      (0.5f - img->width * 0.5f) / (img->width * 0.5f) * scene->cam.xdir;


  for (size_t j = 0; j < img->height; j++) {
    if (j != 0)
      printf("\033[A\r");
    float progress = (float)j / img->height * 100.f;
    printf("progress\t[");
    int cpt = 0;
    for (cpt = 0; cpt < progress; cpt += 5)
      printf(".");
    for (; cpt < 100; cpt += 5)
      printf(" ");
    printf("]\n");
    //int samples_per_pixel = 4;
#pragma omp parallel for
    for (size_t i = 0; i < img->width; i++) {
      color3 *ptr = getPixelPtr(img, i, j);
      color3 pixel_color = {0, 0, 0};
      /*for (int s = 0; s < samples_per_pixel; s++) {
        float u = (float(i) + (rand() / (float)RAND_MAX)) / img->width;
        float v = (float(j) + (rand() / (float)RAND_MAX)) / img->height;
        vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                       (u - 0.5f) * 2.0f * aspect * scene->cam.xdir +
                       (v - 0.5f) * 2.0f * scene->cam.ydir;*/
      vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                     float(i) * dx + float(j) * dy;

      Ray rx;
      rayInit(&rx, scene->cam.position, normalize(ray_dir));
      //pixel_color = pixel_color + trace_ray(scene, &rx, tree);
      //}
      *ptr = trace_ray(scene, &rx, tree);
      //pixel_color = pixel_color * (1.0f / samples_per_pixel);
      //*ptr = pixel_color;
    }
  }
}
