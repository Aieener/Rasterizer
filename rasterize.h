/*
 * rasterize.h
 * Yuding Ai
 * Penn ID: 31295008
 * CIS 560 HW3
 * Due Date: 2017 Mar 16
 */

#ifndef RASTERIZE_H
#define RASTERIZE_H
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include "tiny_obj_loader.h"
#include "mat4.h"
#include "vec4.h"
#include "ppm_filters.h"

//a little function that calc the slop
float calcm(float &qx1,float &qy1, float &qx2, float &qy2);

namespace rasterize {

    //read the camera
    bool Readcamera(const char *fname, std::array<float,6> &frustumpara,vec4 &eye, vec4 &center, vec4 &up);

    //calculate the projection matrix based on the camera data
    std::array<mat4,3> frustrum(std::array<float,6> &frustumpara, vec4 &eye, vec4 &center, vec4 &up);

    //------------------------------------------------------------
    //helper function for render()
    //------------------------------------------------------------

    //------------calc area---------------------------------------
    float calcarea(vec4 V1, vec4 V2, vec4 V3);
    float calcarea2d(vec4 V1, vec4 V2, vec4 V3);

    //------------calc z------------------------------------------
    float calcz(const std::array<float,2> &interz,const std::array<int,2>xintersec, const int &c);
    float correctcalcz(const std::array<float,2> &interz,const std::array<int,2>xintersec, const int &c);

    //--------------calc color------------------------------------

    vec4 calcnorgour(const std::array<std::array<float,2>,3 > &internorv, const std::array<int,2>&xintersec, const int &c);
    vec4 calcnorgourz(const std::array<std::array<float,2>,3 > &internorv,const float &ctargetz,const std::array<float,2> &correctinterz,const std::array<int,2> &xintersec, const int &c);

    vec4 calcnorbary(const std::vector<std::array<vec4,3> > &norvlist, const std::vector<std::array<vec4,3> > &norfacelist, const int &c, const int &y,int s);
    vec4 calcnorbaryz(const std::vector<std::array<vec4,3> > &norvlist, const std::vector<std::array<vec4,3> > &norfacelist, const int &c, const int &y,int s,const float &ctargetz);

    //--------------scan-line-------------------------------------
    std::array<int,2>edgetest(const std::array<vec4,3> &face,std::array<float,2> &interz,std::array<float,2> &correctinterz,const std::array<vec4,3> norv,std::array<std::array<float,2>,3 > &internorv,const std::array< int, 4> &bbox,const float &Y );


    //--------------drawing functions ----------------------------

    void ppmdraw(img_t *img, const std::vector<float> &zbuffer, const std::vector<float> &rawzlist,const std::vector<std::vector<std::array<int,2> > > &xinterseclist,std::vector<std::array<int,4> > bboxlist,const std::vector<std::array<float,3> >&colorlist);

    void ppmadvancedraw(img_t *img, const std::vector<float> &zbuffer, const std::vector<float> &rawzlist,const std::vector<std::vector<std::array<int,2> > > &xinterseclist,std::vector<std::array<int,4> > bboxlist,const std::vector<std::vector<std::array<float,3> > >&adcolorlist);

    
    //--------------The master rendering function ----------------
    img_t *render( int w,  int h,const std::vector<tinyobj::shape_t> &shapes, const std::array<mat4,3> frust,const std::vector<tinyobj::material_t> &materials,int trigger);
}


#endif // RASTERIZE_H
