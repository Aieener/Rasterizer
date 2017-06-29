/*
 * rasterize.cpp
 * Yuding Ai
 * Penn ID: 31295008
 * CIS 560 HW3
 * Due Date: 2017 Mar 16
 */

#include "rasterize.h"
#include <stdio.h>
#include <algorithm>

//------------------------------------------------------------
// directly from HW1 ppm_example.c
//------------------------------------------------------------
// Create a new image of specified size.
img_t *new_img(int w, int h) {
    /* Crash if width or height is 0 or negative.  This is pretty drastic,
     * but very convenient for debugging.  You normally use it as a sanity
     * check for errors that should not occur in production.
     *
     * In the debugger, your program will break at the line with the failed assertion,
     * so you can your variables and figure out what went wrong.  Outside the debugger
     * your program will generally crash with an error that states the failed assertion's
     * condition and its file and line number in the source code.
     */
    assert(w > 0);
    assert(h > 0);
    img_t *img = (img_t *) malloc(sizeof(img_t));

    // now initialize img appropriately
    img->w = w;
    img->h = h;

    // allocate memory for the image pixels
    img->data = (pixel_t *) malloc(w * h * sizeof(pixel_t));

    // zero out all the image pixels so they don't contain garbage
    // memset is quite handy for this
    memset(img->data, 0, w * h * sizeof(pixel_t));

    return img;
}

void destroy_img(img_t **img) {
    // step 1: free the image pixels
    free((*img)->data); // as long as you allocated img->data with malloc, free knows how big it is
    (*img)->data = NULL; // this is a sanity check to make sure we don't try and access img->data after deleting it
    free(*img); // now free the img_t structure
    *img = NULL; // finally, set the img pointer provided by the caller to NULL so the caller
    // can't accidentally access the freed memory
}

// Read in a PPM file
img_t *read_ppm(const char *fname) {
    int w, h;

    assert(fname != NULL); // crash if fname is NULL
    FILE *f = fopen(fname, "rb"); // open the ppm for reading in binary mode
    assert(f != NULL); // crash if the file didn't open

    fscanf(f, "P6 %d %d 255%*c", &w, &h); // read in the header and image width and height

    img_t *img = new_img(w, h); // create an empty image of the correct size

    fread(img->data, w * h, 3, f); // read the image data into img->data

    fclose(f); // close the file

    return img;
}

// Write out a PPM file
void write_ppm(const img_t *img, const char *fname) {
    assert(img != NULL); // crash if img is NULL
    assert(fname != NULL); // crash if fname is NULL

    FILE *f = fopen(fname, "wb"); // open fname for writing in binary mode; clobbers any existing file
    assert(f != NULL); // crash if file did not open

    fprintf(f, "P6\n%d %d 255\n", img->w, img->h); // write the image header
    fwrite(img->data, img->w * img->h, 3, f); // write the image data

    fclose(f); // close the file6
}

//------------------------------------------------------------
//  a little helper function to calculate slope
float calcm(float &qx1,float &qy1, float &qx2, float &qy2){return 1.0*(qy2 -qy1)/(qx2 - qx1);}

//------------------------------------------------------------

namespace rasterize{
    bool Readcamera(const char *fname, std::array<float,6> &frustumpara, vec4 &eye, vec4 &center, vec4 &up){
        //initially, before normalization, set the forth component to be 0;
        eye[3] = 0;
        center[3] = 0;
        up[3] = 0;
        assert(fname != NULL); // crash if fname is NULL;
        FILE *f =fopen(fname,"rb");// open the camera.txt for reading in binary mode
        assert(f != NULL);

        // read the first 6 frustrumparameters and assign then into frustrumpara
        fscanf(f,"%f %f %f %f %f %f",&frustumpara[0],&frustumpara[1],&frustumpara[2],&frustumpara[3],&frustumpara[4],&frustumpara[5]);
        fscanf(f,"%f %f %f", &eye[0],&eye[1],&eye[2]);
        fscanf(f,"%f %f %f", &center[0],&center[1],&center[2]);
        fscanf(f,"%f %f %f", &up[0],&up[1],&up[2]);
        fclose(f);

        //normalize the vec4 center (z) and vec4 up(y)
        center.norm();
        up.norm();

        // set the 4th component of vec4s to 1
        eye[3] = 1;
        center[3] = 1;
        up[3] = 1;
        //        right[3] = 1;
        return true;
    }

    std::array<mat4,3> frustrum(std::array<float,6> &frustumpara, vec4 &eye, vec4 &center, vec4 &up){
        //create an output mat4
        std::array<mat4,3> output;

        // The frustumparameter: left---> frustumpara[0]; right---> frustumpara[1]; top---> frustumpara[2];
        // bottom---> frustumpara[3]; near---> frustumpara[4]; far---> frustumpara[5]
        float l = frustumpara[0];
        float r = frustumpara[1];
        float b = frustumpara[2];
        float t = frustumpara[3];
        float n = frustumpara[4];
        float f = frustumpara[5];

        //float A = std::abs((r-l)/(t-b));
        float A = ((r-l)/(t-b));
        //float tanofhalf = std::abs((t-b)/(2*n));
        float tanofhalf = ((t-b)/(2*n));


        // Calculate the proj matrix using the one given in class
        /*       1/A*tanofhalf    0         (r+l)/(r-l)     0
         *       0          1/tanofhalf     (t+b)/(t-b)     0
         *       0                0             f/(f-n) -nf/(f-n)
         *       0                0              1          0
         *
         *   Same as
         *         2n/(r-l)    0        (r+l)/(r-l)     0
         *         0        2n/(t-b)    (t+b)/(t-b)     0
         *         0           0             f/(f-n) -nf/(f-n)
         *         0           0              1          0
         *
         */


        vec4 f0(1/(A*tanofhalf),0,0,0);
        vec4 f1(0,1/(tanofhalf),0,0);
        vec4 f2((r+l)/(r-l),(t+b)/(t-b),f/(f-n),1);
        vec4 f3(0,0,-n*f/(f-n),0);

        mat4 proj(f0,f1,f2,f3);
        output[0] = proj;

        // calculate the view matrix
        vec4 forward = center - eye;
        forward.norm();
        // compute the Right: x axis by cross product cross(center,up): z cross y since we have left-handed coordinate here
        vec4 right = cross(forward,up);
        right.norm();
        vec4 v0(right[0],up[0],forward[0],0);
        vec4 v1(right[1],up[1],forward[1],0);
        vec4 v2(right[2],up[2],forward[2],0);
        vec4 v3(0,0,0,1);

        mat4 view(v0,v1,v2,v3);
        output[1] = view;


        //calculate the view translate matrix
        mat4 translate;
        translate = translate.trans(-eye[0],-eye[1],-eye[2]);
        output[2] = translate;

        return output;
    }


    //------------------------------------------------------------
    //helper function for render()
    //------------------------------------------------------------
    
    
    //--------------calc area-------------------------------------
    float calcarea(vec4 V1, vec4 V2, vec4 V3){

        /*
         * V1 = (x1,y1,z1) V2 = (x2,y2,z2) V3 = (x3,y3,z3)
         * P = V2 - V1;  Q = V3 - V2
         * the area is then 1/2* |PxQ|
         * 
         *                    |   i        j       k   |
         * the Area =1/2 * det| x2-x1    y2-y1   z2-z1 | 
         *                    | x3-x1    y3-y1   z3-z1 |
         *          = 0.5 * (((y2-y1)*(z3-z1)-(z2-z1)*(y3-y1))**2 + ((z2-z1)*(x3-x1) - (x2-x1)*(z3-z1))**2 + ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1))**2)**0.5;
         * reference: http://www.had2know.com/academics/triangle-area-perimeter-angle-3-coordinates.html  
         */ 

        float x1,x2,x3,y1,y2,y3,z1,z2,z3;
        x1 = V1[0];
        y1 = V1[1];
        z1 = V1[2];

        x2 = V2[0];
        y2 = V2[1];
        z2 = V2[2];

        x3 = V3[0];
        y3 = V3[1];
        z3 = V3[2];
        float area;

        area= 0.5 * std::pow(std::pow((y2-y1)*(z3-z1)-(z2-z1)*(y3-y1),2) + std::pow((z2-z1)*(x3-x1) - (x2-x1)*(z3-z1),2) + std::pow((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1),2),0.5);
        return area;
    }

    float calcarea2d(vec4 V1, vec4 V2, vec4 V3){


        float x1,y1,x2,y2,x3,y3;
        x1 = V1[0];
        y1 = V1[1];

        x2 = V2[0];
        y2 = V2[1];

        x3 = V3[0];
        y3 = V3[1];
        float area;

        area= 0.5 * std::abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
        return area;
    }

    //--------------calc z----------------------------------------

    float calcz(const std::array<float,2> &interz,const std::array<int,2>xintersec, const int &c){
        // A helper function that calculate the z depth values for each triangle
        
        /* the input 'facelist' takes a collection of the vertexs for each triangle (each triangle with 3 vec4)
         * Using the three vertexs of each face:
         * 1.find t value : t = (X -xintersec[0])/(xintersec[1]-xintersec[0])
         * 2.linear interpolate z value by z = (z_1-z_0)*t + z_0;
         */ 
        
        float X = c*1.0 + 0.5;
        float t =(X-1.0*xintersec[0])/(1.0*xintersec[1]-1.0*xintersec[0]);

        float z = 1.0*t*(interz[1]-interz[0]) + 1.0*interz[0];

        return z;
    }

    float correctcalcz(const std::array<float,2> &interz,const std::array<int,2>xintersec, const int &c){
        // A helper function that calculate the perspective corrected z depth values for each triangle
        
        /* the input 'facelist' takes a collection of the vertexs for each triangle (each triangle with 3 vec4)
         * Using the three vertexs of each face:
         * 1.find t value : t = (c -xintersec[0])/(xintersec[1]-xintersec[0])
         * 2.linear interpolate z value by 1/z = (1/z_1-1/z_0)*t + 1/z_0;
         */ 
        float t =(1.0*c-1.0*xintersec[0])/(1.0*xintersec[1]-1.0*xintersec[0]);

        float z = pow(1.0*t*(1.0/interz[1]-1.0/interz[0]) + 1.0/interz[0],-1);

        return z;
    }


    //--------------calc color------------------------------------
    vec4 calcnorgour(const std::array<std::array<float,2>,3 > &internorv, const std::array<int,2> &xintersec, const int &c){
        // A helper function that calculate the norm_gouraud color values for each pixel in each triangle
        
        //reference:  https://www.stack.nl/~dimitri/3dsview/gouraud.html
        /* use the input internorv and xintersec to interpolate the intensity at the target pixel 
         * Using the fllowing equation
         * C = C_0*(1-t) + C_1*t = alpha*(C_1-C_0) + c_0
         */ 
        vec4 pixnormal(0,0,0,0);
        float delta = 0.0001;
        if(xintersec[1] - xintersec[0]> delta){
            float t =(1.0*c-1.0*xintersec[0])/(1.0*xintersec[1]-1.0*xintersec[0]);

            pixnormal[0] = 1.0*t*(internorv[0][1]-internorv[0][0]) + internorv[0][0];
            pixnormal[1] = 1.0*t*(internorv[1][1]-internorv[1][0]) + internorv[1][0];
            pixnormal[2] = 1.0*t*(internorv[2][1]-internorv[2][0]) + internorv[2][0];

            pixnormal.norm();
        }
        else{
            pixnormal[0] = internorv[0][0];
            pixnormal[1] = internorv[1][0];
            pixnormal[2] = internorv[2][0];
        }

        return pixnormal;

    }

    vec4 calcnorgourz(const std::array<std::array<float,2>,3 > &internorv,const float &ctargetz,const std::array<float,2> &correctinterz,const std::array<int,2> &xintersec, const int &c){
        //flowwing the instruction in class, using perspective correction to revise color at each pixel
        //also following the instructoin from www.scratchapixel.com

        vec4 pixnormal(0,0,0,0);
        float delta = 0.0001;
        if(xintersec[1] - xintersec[0]> delta){
            float t =(1.0*c-1.0*xintersec[0])/(1.0*xintersec[1]-1.0*xintersec[0]);

            pixnormal[0] = ctargetz*((1-t)*internorv[0][0]/(correctinterz[0])  + t*internorv[0][1]/(correctinterz[1]));
            pixnormal[1] = ctargetz*((1-t)*internorv[1][0]/(correctinterz[0])  + t*internorv[1][1]/(correctinterz[1]));
            pixnormal[2] = ctargetz*((1-t)*internorv[2][0]/(correctinterz[0])  + t*internorv[2][1]/(correctinterz[1]));

            pixnormal.norm();
        }
        else{
            pixnormal[0] = internorv[0][0];
            pixnormal[1] = internorv[1][0];
            pixnormal[2] = internorv[2][0];
        }

        pixnormal.norm();

        return pixnormal;
    }

    vec4 calcnorbary(const std::vector<std::array<vec4,3> > &norvlist, const std::vector<std::array<vec4,3> > &norfacelist, const int &c, const int &y,int s){
        // A helper function that calculate the norm_bary color values for each pixel in each triangle
        
        /* the input 'facelist' takes a collection of vertexs for each triangle (each triangle with 3 vec4)
         * Using the three vertexs of each face:
         * 1 . use the x,y value of the target pixel to calculate S1,S2,S3 (the area is calculate in 2d)
         * 2 . Then calculate color using C = C1*(S1/S) + C2*(S2/S) + C3*(S3/S)
         */ 

        //first get the target pixel a 2d coordinate; so even I use vec4, I only care about the first two component since we are in 2d
        //and not consider the distorsion 
        vec4 pix(c,y,1,1);

        //next the three vertex of the bounding 2d triangle would be:
        // vetex1 = norfacelist[s][0];
        // vetex2 = norfacelist[s][1];
        // vetex3 = norfacelist[s][2];
        // where s is the idx of triangle
         
        // now compute S, S1, S2 and S3
        float S,S1,S2,S3;
        // S:
        S = calcarea2d(norfacelist[s][0],norfacelist[s][1],norfacelist[s][2]);

        // S1:  center, facelist[s][1], facelist[s][2] corresponding to v1
        S1 = calcarea2d(pix,norfacelist[s][1],norfacelist[s][2]);
        
        // S2:  center, facelist[s][1], facelist[s][2] corresponding to v2
        S2 = calcarea2d(pix,norfacelist[s][2],norfacelist[s][0]);
        
        // S3:  center, facelist[s][2], facelist[s][0] corresponding to v3
        S3 = calcarea2d(pix,norfacelist[s][0],norfacelist[s][1]);

        // last compute norx, nory, norz using: C = C1*(S1/S) + C2*(S2/S) + C3*(S3/S)

        float nz,nz1,nz2,nz3;
        nz1 = norvlist[s][0][2];
        nz2 = norvlist[s][1][2];
        nz3 = norvlist[s][2][2];
        //nz = std::pow((1/nz1)*(S1/S) + (1/nz2)*(S2/S) + (1/nz3)*(S3/S),-1);
        nz = (nz1)*(S1/S) + (nz2)*(S2/S) + (nz3)*(S3/S);

        float ny,ny1,ny2,ny3;
        ny1 = norvlist[s][0][1];
        ny2 = norvlist[s][1][1];
        ny3 = norvlist[s][2][1];
        //ny = std::pow((1/ny1)*(S1/S) + (1/ny2)*(S2/S) + (1/ny3)*(S3/S),-1);
        ny = (ny1)*(S1/S) + (ny2)*(S2/S) + (ny3)*(S3/S);

        float nx,nx1,nx2,nx3;
        nx1 = norvlist[s][0][0];
        nx2 = norvlist[s][1][0];
        nx3 = norvlist[s][2][0];
        //nx = std::pow((1/nx1)*(S1/S) + (1/nx2)*(S2/S) + (1/nx3)*(S3/S),-1);
        nx = (nx1)*(S1/S) + (nx2)*(S2/S) + (nx3)*(S3/S);

        vec4 pixnormal(nx,ny,nz,0);

        pixnormal.norm();

        return pixnormal;
    }

    vec4 calcnorbaryz(const std::vector<std::array<vec4,3> > &norvlist, const std::vector<std::array<vec4,3> > &norfacelist, const int &c, const int &y,int s,const float &ctargetz){

        // reference the slide from piazza @175
         
        /* the input 'facelist' takes a collection of vertexs for each triangle (each triangle with 3 vec4)
         * Using the three vertexs of each face:
         * 1 . use the x,y,z value of the target pixel to calculate S1,S2,S3 (the area is calculate in 3d)
         * 2 . Then calculate color using C = Z*((C1/Z1)*(S1/S) + (C2/Z2)*(S2/S) + (C3/Z3)*(S3/S))
         */ 


        //first get the target pixel a 3d coordinate; 
        //and not consider the distorsion 
        vec4 pix(c,y,ctargetz,1);

        //next the three vertex of the bounding 2d triangle would be:
        // vetex1 = norfacelist[s][0];
        // vetex2 = norfacelist[s][1];
        // vetex3 = norfacelist[s][2];
        // where s is the idx of triangle
         
        // now compute S, S1, S2 and S3
        float S,S1,S2,S3;
        // S:
        S = calcarea(norfacelist[s][0],norfacelist[s][1],norfacelist[s][2]);

        // S1:  center, facelist[s][1], facelist[s][2] corresponding to v1
        S1 = calcarea(pix,norfacelist[s][1],norfacelist[s][2]);
        
        // S2:  center, facelist[s][1], facelist[s][2] corresponding to v2
        S2 = calcarea(pix,norfacelist[s][2],norfacelist[s][0]);
        
        // S3:  center, facelist[s][2], facelist[s][0] corresponding to v3
        S3 = calcarea(pix,norfacelist[s][0],norfacelist[s][1]);

        // last compute norx, nory, norz using: C = Ztarget*[C1/Z1*(S1/S) + C2/Z2*(S2/S) + C3/Z3*(S3/S)]
        float Z1,Z2,Z3;
        Z1 = norfacelist[s][0][2];
        Z2 = norfacelist[s][1][2];
        Z3 = norfacelist[s][2][2];

        float nz,nz1,nz2,nz3;
        nz1 = norvlist[s][0][2];
        nz2 = norvlist[s][1][2];
        nz3 = norvlist[s][2][2];

        nz = ctargetz*((nz1/Z1)*(S1/S) + (nz2/Z2)*(S2/S) + (nz3/Z3)*(S3/S));

        float ny,ny1,ny2,ny3;
        ny1 = norvlist[s][0][1];
        ny2 = norvlist[s][1][1];
        ny3 = norvlist[s][2][1];
        
        ny = ctargetz*((ny1/Z1)*(S1/S) + (ny2/Z2)*(S2/S) + (ny3/Z3)*(S3/S));

        float nx,nx1,nx2,nx3;
        nx1 = norvlist[s][0][0];
        nx2 = norvlist[s][1][0];
        nx3 = norvlist[s][2][0];

        nx = ctargetz*((nx1/Z1)*(S1/S) + (nx2/Z2)*(S2/S) + (nx3/Z3)*(S3/S));

        vec4 pixnormal(nx,ny,nz,0);

        pixnormal.norm();

        return pixnormal;
    }
    
    //--------------scan-line-------------------------------------
    std::array<int,2>edgetest(const std::array<vec4,3> &face,std::array<float,2> &interz,std::array<float,2> &correctinterz,const std::array<vec4,3> norv,std::array<std::array<float,2>,3 > &internorv,const std::array< int, 4> &bbox,const float &Y ){
        // return an array that contains the boundary of where line intersect on x-axis
        // calculation following the slide: Line segment intersection
        // the 3 vertex of the triangle are q1,q2,q3. here we take q1,q2 for example
        // p1.y = m2*X -m2q1.x + q1.y= Y
        // m2 = (q2.y - q1.y)/(q2.x-q1.x)
        // ==> Y ---> we will scan through a fix Y each time, so Y = 0.5,1.5,2.5 ....
        // ==> X = (Y-q1.y + m2q1.x)/m2
        // in such way, we find the all X values as the xintersection

        // then we scan though each row to see at where the line intersects
        
        /* In addition, this function also calculate the z value/perspective corrected z value on both side of the scanline and
         * the normal(color) on both side of the scanline. those values are loaded directly into the reference input: &interz, &correctinterz,
         * and &internorv
         */

        float delta = 0.0;

        int xmin = bbox[0];

        std::array<int,2> bound;
        bound[0] = xmin;
        bound[1] = xmin;

        
        //Now 'scan' through each line
        //It might not be an efficient way, but I have tried many method and 
        //this is the only one that works for me ...
        
        //first locate Y so see which two pairs of vertex to intersect 
        //(Notice, unless the there is a degenerate triangle, Y can't in between 
        //three pairs of vertex at the same time)
        
        //so for example: if Y is in between q1_y and q2_y and in between q2_y and q3_y
        //we then use X12 and X23 to check the intersection point
        
        std::array<vec4,3> XYlst = face;  //take a copy of face, then sort by Ycoor
        std::array<vec4,3> Noredge = norv;

        vec4 tempvec;
        vec4 nortempvec;

        //// Yeah, it's bubble sort, but the length is only 3.
        for (int c = 0 ; c < ( 3 - 1 ); c++)
        {
            for (int d = 0 ; d < 3 - c - 1; d++)
            {
                if (XYlst[d][1] > XYlst[d+1][1]) 
                {
                    tempvec  = XYlst[d];
                    nortempvec = Noredge[d];

                    XYlst[d]   = XYlst[d+1];
                    Noredge[d] = Noredge[d+1];

                    XYlst[d+1] = tempvec;
                    Noredge[d+1] = nortempvec;
                }
            }
        }
        // now the XYlst and Rawlist is sorted by y values
        // the "ymiddle" vetex has a y value of XYlst[1][1]
        // so does the Noredge is also sorted by y value on XYlst
        
        float q1_x, q1_y,q1_z,q2_x, q2_y,q2_z,q3_x, q3_y,q3_z;
        float nq1_x,nq2_x,nq3_x,nq1_y,nq2_y,nq3_y, nq1_z,nq2_z,nq3_z;

        // q1_y<=q2_y<=q3_y
        q1_x = XYlst[0][0];
        q1_y = XYlst[0][1];
        q1_z = XYlst[0][2];

        q2_x = XYlst[1][0];
        q2_y = XYlst[1][1];
        q2_z = XYlst[1][2];


        q3_x = XYlst[2][0];
        q3_y = XYlst[2][1];
        q3_z = XYlst[2][2];

        nq1_x = Noredge[0][0];
        nq1_y = Noredge[0][1];
        nq1_z = Noredge[0][2];

        nq2_x = Noredge[1][0];
        nq2_y = Noredge[1][1];
        nq2_z = Noredge[1][2];

        nq3_x = Noredge[2][0];
        nq3_y = Noredge[2][1];
        nq3_z = Noredge[2][2];

        //=================================== work on scanline to get xintersec and interz =================
        //=================================== internorv and perspective correction  ========================
        float X12,X13,X23;
        float m12,m13,m23;
        X12=X13=X23 = 0;
        
        if(Y-q2_y < delta ){
            // the "upper case"

            m12 =calcm(q2_x,q2_y,q1_x,q1_y);
            m13 =calcm(q3_x,q3_y,q1_x,q1_y);
            float xin12,xin13; 
            xin12 = (Y-q1_y)/m12 + q1_x;
            xin13 = (Y-q1_y)/m13 + q1_x;

            if(Y>=q1_y && xin12 >=bbox[0] &&xin12<=bbox[2] &&xin13>=bbox[0] &&xin13<=bbox[2]){
                for (int i = bbox[0]; i < bbox[2]; i++) { 
                    if (xin12 >= i &&xin12<=i+1){
                        X12 = i;
                    } 

                    if (xin13 >=i &&xin13<=i+1){
                        X13 = i;
                    }
                } 
                
                bound[0] = X12;
                bound[1] = X13;
                
                std::sort(bound.begin(),bound.end());
                

                //------------------------------------calc interz ------------------------------------------------
                // linear interpolation to find the z value for the side
                float t12 = 1.0*(Y-q1_y)/(q2_y - q1_y); 
                float t13 = 1.0*(Y-q1_y)/(q3_y - q1_y);

                if(X12>X13){
                    
                    interz[1] = 1.0*t12*(q2_z - q1_z) + 1.0*q1_z;
                    interz[0] = 1.0*t13*(q3_z - q1_z) + 1.0*q1_z;           

                    // ------------------------------------------------------------
                    // ++++++++++++++++++++perspective correstion +++++++++++++++++
                    // ------------------------------------------------------------

                    correctinterz[1] = pow(1.0*t12*(1.0/q2_z - 1.0/q1_z) + 1.0/q1_z,-1);
                    correctinterz[0] = pow(1.0*t13*(1.0/q3_z - 1.0/q1_z) + 1.0/q1_z,-1);           

                }
                else{
                    interz[0] = 1.0*t12*(q2_z - q1_z) + 1.0*q1_z;
                    interz[1] = 1.0*t13*(q3_z - q1_z) + 1.0*q1_z;           

                    // ------------------------------------------------------------
                    // ++++++++++++++++++++perspective correstion +++++++++++++++++
                    // ------------------------------------------------------------

                    correctinterz[0] = pow(1.0*t12*(1.0/q2_z - 1.0/q1_z) + 1.0/q1_z,-1);
                    correctinterz[1] = pow(1.0*t13*(1.0/q3_z - 1.0/q1_z) + 1.0/q1_z,-1);           

                }

                //------------------------------------calc internorv  --------------------------------------------
                
                if(X12>X13){
                    internorv[0][1] = 1.0*t12*(nq2_x - nq1_x)  + 1.0*nq1_x;
                    internorv[0][0] = 1.0*t13*(nq3_x - nq1_x)  + 1.0*nq1_x;

                    internorv[1][1] = 1.0*t12*(nq2_y - nq1_y)  + 1.0*nq1_y;
                    internorv[1][0] = 1.0*t13*(nq3_y - nq1_y)  + 1.0*nq1_y;

                    internorv[2][1] = 1.0*t12*(nq2_z - nq1_z)  + 1.0*nq1_z;
                    internorv[2][0] = 1.0*t13*(nq3_z - nq1_z)  + 1.0*nq1_z;

                }
                else{
                    internorv[0][0] = 1.0*t12*(nq2_x - nq1_x)  + 1.0*nq1_x;
                    internorv[0][1] = 1.0*t13*(nq3_x - nq1_x)  + 1.0*nq1_x;

                    internorv[1][0] = 1.0*t12*(nq2_y - nq1_y)  + 1.0*nq1_y;
                    internorv[1][1] = 1.0*t13*(nq3_y - nq1_y)  + 1.0*nq1_y;

                    internorv[2][0] = 1.0*t12*(nq2_z - nq1_z)  + 1.0*nq1_z;
                    internorv[2][1] = 1.0*t13*(nq3_z - nq1_z)  + 1.0*nq1_z;

                }
            }
        }

        else if(Y-q2_y >=delta ){
            // the "lower case"
            
            m23 =calcm(q2_x,q2_y,q3_x,q3_y);
            m13 =calcm(q3_x,q3_y,q1_x,q1_y);
            float xin23,xin13; 
            xin23 = (Y-q2_y)/m23 + q2_x;
            xin13 = (Y-q1_y)/m13 + q1_x;

            if(Y<=q3_y && xin23 >=bbox[0] &&xin23<=bbox[2] &&xin13>=bbox[0] &&xin13<=bbox[2]){
                for (int i = bbox[0]; i < bbox[2]; i++) { 
                    if (xin23 >= i &&xin23<=i+1){
                        X23 = i;
                    } 

                    if (xin13 >=i &&xin13<=i+1){
                        X13 = i;
                    }
                } 

                bound[0] = X23;
                bound[1] = X13;
                
                std::sort(bound.begin(),bound.end());


                //------------------------------------calc interz ------------------------------------------------
                // linear interpolation to find the z value for the side
                float t23 = 1.0*(Y-q2_y)/(q3_y - q2_y); 
                float t13 = 1.0*(Y-q1_y)/(q3_y - q1_y);

                if(X23>X13){
                    interz[1] = 1.0*t23*(q3_z - q2_z) + 1.0*q2_z;
                    interz[0] = 1.0*t13*(q3_z - q1_z) + 1.0*q1_z;           

                    // ------------------------------------------------------------
                    // ++++++++++++++++++++perspective correstion +++++++++++++++++
                    // ------------------------------------------------------------

                    correctinterz[1] = pow(1.0*t23*(1.0/q3_z - 1.0/q2_z) + 1.0/q2_z,-1);  
                    correctinterz[0] = pow(1.0*t13*(1.0/q3_z - 1.0/q1_z) + 1.0/q1_z,-1);           

                }
                else{
                    interz[0] = 1.0*t23*(q3_z - q2_z) + 1.0*q2_z;
                    interz[1] = 1.0*t13*(q3_z - q1_z) + 1.0*q1_z;           

                    // ------------------------------------------------------------
                    // ++++++++++++++++++++perspective correstion +++++++++++++++++
                    // ------------------------------------------------------------

                    correctinterz[0] = pow(1.0*t23*(1.0/q3_z - 1.0/q2_z) + 1.0/q2_z,-1);  
                    correctinterz[1] = pow(1.0*t13*(1.0/q3_z - 1.0/q1_z) + 1.0/q1_z,-1);           


                }

                //------------------------------------calc internorv  --------------------------------------------
                if(X23>X13){
                    internorv[0][1] = 1.0*t23*(nq3_x - nq2_x)  + 1.0*nq2_x;
                    internorv[0][0] = 1.0*t13*(nq3_x - nq1_x)  + 1.0*nq1_x;

                    internorv[1][1] = 1.0*t23*(nq3_y - nq2_y)  + 1.0*nq2_y;
                    internorv[1][0] = 1.0*t13*(nq3_y - nq1_y)  + 1.0*nq1_y;

                    internorv[2][1] = 1.0*t23*(nq3_z - nq2_z)  + 1.0*nq2_z;
                    internorv[2][0] = 1.0*t13*(nq3_z - nq1_z)  + 1.0*nq1_z;

                }
                else{
                    internorv[0][0] = 1.0*t23*(nq3_x - nq2_x)  + 1.0*nq2_x;
                    internorv[0][1] = 1.0*t13*(nq3_x - nq1_x)  + 1.0*nq1_x;

                    internorv[1][0] = 1.0*t23*(nq3_y - nq2_y)  + 1.0*nq2_y;
                    internorv[1][1] = 1.0*t13*(nq3_y - nq1_y)  + 1.0*nq1_y;

                    internorv[2][0] = 1.0*t23*(nq3_z - nq2_z)  + 1.0*nq2_z;
                    internorv[2][1] = 1.0*t13*(nq3_z - nq1_z)  + 1.0*nq1_z;

                }

            }

        }
        //=================================== All finish ==================================================
        return bound;

    }


    //--------------drawing functions ----------------------------
    void ppmadvancedraw(img_t *img, const std::vector<float> &zbuffer, const std::vector<float> &rawzlist,const std::vector<std::vector<std::array<int,2> > > &xinterseclist,std::vector<std::array<int,4> > bboxlist,const std::vector<std::vector<std::array<float,3> > >&adcolorlist){
        
        pixel_t *p = img->data ;// declare a pointer p to point to an array of pixel_t on the img
        int w = img->w;
        int h = img->h;

        float delta = 1e-200;
        int k = 0;
        for (unsigned int i = 0; i < bboxlist.size(); i++) { 
            int idx = 0;
            for(int y = bboxlist[i][1];y< bboxlist[i][3];y++ ){
                for (int c =xinterseclist[i][y-bboxlist[i][1]][0] ; c<xinterseclist[i][y-bboxlist[i][1]][1];c++){

                    //only draw the part inside the img
                    if(c>=0 &&c<w && y>=0 &&y<h){
                        if( rawzlist[k] - zbuffer[y*w+c] ==delta){

                            //draw the ppm pic
                            //================================================
                            p[y*w+ c].r = adcolorlist[i][idx][0];
                            p[y*w+ c].g = adcolorlist[i][idx][1];
                            p[y*w+ c].b = adcolorlist[i][idx][2];
                        }
                    }
                    k++;
                    idx++;
                }
            } 
        }
    }

    void ppmdraw(img_t *img, const std::vector<float> &zbuffer, const std::vector<float> &rawzlist,const std::vector<std::vector<std::array<int,2> > > &xinterseclist,std::vector<std::array<int,4> > bboxlist,const std::vector<std::array<float,3> >&colorlist){
        
        pixel_t *p = img->data ;// declare a pointer p to point to an array of pixel_t on the img
        int w = img->w;
        int h = img->h;

        float delta = 1e-200;
        int k = 0;
        for (unsigned int i = 0; i < bboxlist.size(); i++) { 

            for(int y = bboxlist[i][1];y< bboxlist[i][3];y++ ){

                for (int c =xinterseclist[i][y-bboxlist[i][1]][0] ; c<xinterseclist[i][y-bboxlist[i][1]][1];c++){

                    //only draw the part inside the img
                    if(c>=0 &&c<w && y>=0 &&y<h){
                        if( rawzlist[k] - zbuffer[y*w+c] ==delta){

                            //draw the ppm pic
                            //================================================
                            p[y*w+ c].r = colorlist[i][0];
                            p[y*w+ c].g = colorlist[i][1];
                            p[y*w+ c].b = colorlist[i][2];
                        }
                    }
                    k++;
                }
            } 
        }
    }

    //--------------The master rendering function ----------------

    img_t *render( int w,  int h,const std::vector<tinyobj::shape_t> &shapes, const std::array<mat4,3> frust,const std::vector<tinyobj::material_t> &materials,int trigger){

        img_t *img = new_img(w, h); // create an empty image of the correct size

        //vector lists
        std::vector<std::array<vec4,3> > norfacelist; // create a vector to hold all the normalizsed triangles
        std::vector<std::array<vec4,3> > norvlist; // create a vector to hold all the normal vectors for each vertex in triangles

        //useful lists for drawing
        std::vector<std::array<int,4> > bboxlist; // create a vector to hold all the bounding box;
        std::array<int,2> xintersec;              // the xintersec for each scanline
        std::vector<std::vector<std::array<int,2> > > xinterseclist; // create a vector to hold all the xintersections on each row 
        
        //zvalue lists
        std::vector<float> rawzlist; // create an array of float that hold zbuffer and initialize z to be 2
        std::vector<float> crawzlist; // create an array of float that hold corrected zbuffer and initialize z to be 2
        std::vector<float> zb;   // create an array of float that hold zbuffer and initialize z to be 2
        std::vector<float> czb;   // create an array of float that hold corrected zbuffer and initialize z to be 2

        //norv value list for the side of each scanline
        std::array<std::array<float,2>,3 > internorv; // the normals at for the edge of scanline 

        //z value for the side of each scanline
        std::array<float,2> interz;          // the z value for the side of triangle
        std::array<float,2> correctinterz;   // the z value for the side triangle using perspective correction

        float delta = 0;
        int size = w * h;
        for (int i = 0; i < size; i++) { 
            zb.push_back(2); 
            czb.push_back(2); 
        }

        //colors
        std::array<float,3>rcolor; //random color (not required, but I use it to check if my zbuffer are working correctly before implementing other color methods)
        //--------------------------------------------------
        std::array<float,3>dcolor; //diffuse color by default
        std::array<float,3>wcolor; //white color
        std::array<float,3>ncolor; //normal_flat color
        std::array<float,3>ngcolor; //normal_gouraud color
        std::array<float,3>nbcolor; //normal_bary color
        std::array<float,3>nbzcolor; //normal_bary_z color
        std::array<float,3>ngzcolor; //normal_gouraud_z color


        //color lists
        std::vector<std::array<float,3> >rcolorlist;   // a list to store random color
        std::vector<std::array<float,3> >dcolorlist;   // a list to store diffuse color
        std::vector<std::array<float,3> >wcolorlist;   // a list to store white color
        std::vector<std::array<float,3> >ncolorlist;   // a list to store normal_flat color

        //advance color list
        std::vector<std::vector<std::array<float,3> > >nbcolorlist;   // a list to store normal_bary color
        std::vector<std::vector<std::array<float,3> > >ngcolorlist;   // a list to store normal_gouraud color
        std::vector<std::vector<std::array<float,3> > >nbzcolorlist;   // a list to store normal_bary_z color
        std::vector<std::vector<std::array<float,3> > >ngzcolorlist;   // a list to store normal_gouraud color


        srand(time(NULL));
        //std::cout<<"shape# ="<<shapes.size()<<std::endl;

        //---------------- extract data from shapes------------

        // if all shapes are triangle: for example the cube.obj,
        for (unsigned int i = 0; i< shapes.size(); i++){

            unsigned int numofshape = shapes[i].mesh.indices.size()/3;
            //std::cout<<"#of triangle = "<<numofshape<<std::endl;
            for(unsigned int s = 0; s<numofshape;s++){
                // for each shape/triangle :
                // the mesh from shape
                
                tinyobj::mesh_t m= shapes[i].mesh;

                // the 3 indices of vertex for each shape/triangle
                std::array< int,3> idx;
                idx[0] = m.indices[3*s];
                idx[1] = m.indices[3*s + 1];
                idx[2] = m.indices[3*s + 2];

                //create a face data structure as an array of 3 vec4
                // to hold the pos of the 3 vertex of each triangle
                std::array<vec4,3> face;
                std::array<vec4,3> norv; // the normals at each vertex
                std::array< int,4> bbox = {{w,h,0,0}}; //bouding box:[xmin,ymin,xmax,ymax]

                // assign the 3 pos regarding to the 3 indices into face
                for (unsigned int j = 0; j < face.size(); j++) {

                    //--------------operate on face[j]---------------------- 
                    face[j][0] = m.positions[idx[j]*3];      //xpos
                    face[j][1] = m.positions[idx[j]*3+1];    //ypos
                    face[j][2] = m.positions[idx[j]*3+2];    //zpos
                    face[j][3] = 1;                          //set w = 1

                    //multiply all vertex coordinate by frustrum matrix(world to screen):
                    //thus to convert the triangle into raster space
                    //frus * coor --> the coor on raster space 
                    
                    vec4 buffer = frust[0]*frust[1]*frust[2]* face[j];
                    //divide it by w
                    float wcoor = buffer[3];
                    buffer /= wcoor;
                    face[j] = buffer;
                    //now face[j] stands for the coordinate on the raster space for jth vetex

                    //last convert it to pixel space:
                    face[j][0] = ((face[j][0]+1)*w/2);
                    face[j][1] = ((1-face[j][1])*h/2);
                    //---------------finish face[j]------------------------- 
                    
                   
                    //---------------operate on norv[j]--------------------- 
                    //load the normal values into vec4 norv
                    norv[j][0] = m.normals[idx[j]*3];     //xnor
                    norv[j][1] = m.normals[idx[j]*3+1];   //ynor
                    norv[j][2] = m.normals[idx[j]*3+2];   //znor
                    norv[j][3] = 0;                         //set w = 0
                    
                    //left multiply view metrix onto norv
                    vec4 buffer2 = frust[1]*frust[2]*norv[j];
                    norv[j] = buffer2;
                    
                    //then normalize the norv
                    norv[j].norm();

                    //---------------finish norv[j]------------------------- 

                }

                //check if all the z coor is in between 0 and 1 
                int counter = 0;
                    
                for (int i = 0; i < 3; i++) { 
                    if(face[i][2]<0|| face[i][2]>1.00001){
                        // as long as there is a signlevertex has a valid z value
                        // we should consider sucn triangle
                        counter++;
                    } 
                } 
                //if counter == 3 means all three vertex are outside the img, ignore such triangle

                if(counter!=3){
                    norvlist.push_back(norv);    // load the normal direction into norvlist
                    norfacelist.push_back(face); // load the norfacelist
                
                    std::array<float,3> Xlist;
                    Xlist[0] = face[0][0];
                    Xlist[1] = face[1][0];
                    Xlist[2] = face[2][0];
                    std::sort(Xlist.begin(),Xlist.end());
                    std::array<float,3> Ylist;
                    Ylist[0] = face[0][1];
                    Ylist[1] = face[1][1];
                    Ylist[2] = face[2][1];
                    std::sort(Ylist.begin(),Ylist.end());

                    bbox[0] = Xlist[0];
                    bbox[2] = Xlist[2]+1;
                    bbox[1] = Ylist[0];
                    bbox[3] = Ylist[2]+1;

                    bboxlist.push_back(bbox);
                    


                    //============================define the simple colors ===================================
                    //The simple colors: wcolor, rcolor, dcolor and ncolor
                    //----------------------------------------------------------------------------------------
                    //
                    //white
                    wcolor[0] = 255;
                    wcolor[1] = 255;
                    wcolor[2] = 255;

                    wcolorlist.push_back(wcolor);

                    //random
                    rcolor[0] = rand()%255;
                    rcolor[1] = rand()%255;
                    rcolor[2] = rand()%255;
                    rcolorlist.push_back(rcolor);

                    //diffuse color
                    if (materials.size()>= 1){
                        dcolor[0] = 255*materials[m.material_ids[s]].diffuse[0];
                        dcolor[1] = 255*materials[m.material_ids[s]].diffuse[1];
                        dcolor[2] = 255*materials[m.material_ids[s]].diffuse[2];
                    }
                    else{
                        dcolor[0] = 255; 
                        dcolor[1] = 255; 
                        dcolor[2] = 255; 
                    }
                    dcolorlist.push_back(dcolor);

                    //normal_flat color: base on the norm of the first vertex listed for the face
                    // -1 -> 0;  1 -> 255;
                    ncolor[0] = 255*((norv[0][0]+1)/2);
                    ncolor[1] = 255*((norv[0][1]+1)/2);
                    ncolor[2] = 255*((norv[0][2]+1)/2);
                    ncolorlist.push_back(ncolor);
                    

                    //============================finish the simple colors ===================================
                    
                    float targetz;
                    float ctargetz;

                    std::vector<std::array<int,2> > xintersecset; // create a vector to hold a set of xintersections on each box 
                    std::vector<std::array<float,3> > nbcolorset; // create a vector to hold a set of nbcolor for each pixel on each box 
                    std::vector<std::array<float,3> > ngcolorset; // create a vector to hold a set of ngbcolor for each pixel on each box 
                    std::vector<std::array<float,3> > nbzcolorset; // create a vector to hold a set of nbzcolor for each pixel on each box 
                    std::vector<std::array<float,3> > ngzcolorset; // create a vector to hold a set of ngzcolor for each pixel on each box 


                    for(int y = bbox[1];y< bbox[3];y++ ){
                        // the "ycoor" is actually y + 0.5
                        // so it scann through the middle of each
                        // pixel
                        float Y = y*1.0+0.5 ;
                        // create a helper function edgetest to do the scanline job
                        // furthermore, the edgetest also linear interpolate the z value/perspective correct z and normal/perspective correct 
                        // normal value on the side egde of each triangle
                        // so this 'edgetest' function has done a lot of job
                        xintersec = edgetest(face,interz,correctinterz,norv,internorv,bbox,Y);

                        xintersecset.push_back(xintersec);
                        
                        for (int c = xintersec[0]; c<xintersec[1];c++){
                            
                            //-------------------------------------------------------------------------------------- 
                            targetz = calcz(interz,xintersec,c);
                            ctargetz = correctcalcz(interz,xintersec,c);

                            rawzlist.push_back(targetz); //store the zvalue/depth for each triangle into a list
                            crawzlist.push_back(ctargetz); //store the perspective corrected zvalue/depth for each triangle into a list
                            
                            //============================define advanced colors ===================================
                            //nbcolor, ngcolor, nbzcolor, ngzcolor 
                            //-------------------------------------------------------------------------------------- 
                            //normal_bary color:calc color on 2d image space using barycentric interpolation
                            vec4 nbcol = calcnorbary(norvlist, norfacelist, c, y, s);
                            nbcolor[0] = 255*((nbcol[0]+1)/2);
                            nbcolor[1] = 255*((nbcol[1]+1)/2);
                            nbcolor[2] = 255*((nbcol[2]+1)/2);
                            nbcolorset.push_back(nbcolor);

                            //-------------------------------------------------------------------------------------- 
                            //normal_gouraud color: calc color on 2d image space using gouraud interpolation
                            
                            vec4 ngcol = calcnorgour(internorv, xintersec, c);
                            ngcolor[0] = 255*((ngcol[0]+1)/2);
                            ngcolor[1] = 255*((ngcol[1]+1)/2);
                            ngcolor[2] = 255*((ngcol[2]+1)/2);
                            ngcolorset.push_back(ngcolor);
                            

                            //-------------------------------------------------------------------------------------- 
                            //normal_gouraud_z color: calc color using gouraud interpolation with perspective correction

                            vec4 ngzcol= calcnorgourz(internorv,ctargetz,correctinterz,xintersec,c);
                            ngzcolor[0] = 255*((ngzcol[0]+1)/2);
                            ngzcolor[1] = 255*((ngzcol[1]+1)/2);
                            ngzcolor[2] = 255*((ngzcol[2]+1)/2);
                            ngzcolorset.push_back(ngzcolor);
                            

                            //-------------------------------------------------------------------------------------- 
                            //normal_bary_z color: calc color using barycentric interpolation with perspective correction
                            vec4 nbzcol = calcnorbaryz(norvlist,norfacelist,c,y,s,ctargetz);
                            nbzcolor[0] = 255*((nbzcol[0]+1)/2);
                            nbzcolor[1] = 255*((nbzcol[1]+1)/2);
                            nbzcolor[2] = 255*((nbzcol[2]+1)/2);
                            nbzcolorset.push_back(nbzcolor);
                            
                            //============================finish advanced colors ===================================
                            
                            //----------------------------record zbuffer values ------------------------------------
                            if(c>=0 && c<w &&y>=0 && y<h){
                                if (zb[y*w +c]-targetz > delta ){
                                    // update the zbufferlist
                                    zb[y*w + c] = targetz;
                                }
                            }
                            //----------------------------record perspective corrected zbuffer values --------------
                            if(c>=0 && c<w &&y>=0 && y<h){
                                if (czb[y*w +c]-ctargetz > delta ){
                                    // update the zbufferlist
                                    czb[y*w + c] = ctargetz;
                                }
                            }
                            //----------------------------finish zbuffer values ------------------------------------
                        }
                    }
                    
                    //----------- push xintersection value into xinterseclist -------------------
                    xinterseclist.push_back(xintersecset);

                    //----------- push advance color to advance colorlist -----------------------
                    nbcolorlist.push_back(nbcolorset);
                    ngcolorlist.push_back(ngcolorset);
                    ngzcolorlist.push_back(ngzcolorset);
                    nbzcolorlist.push_back(nbzcolorset);
                }
            }
        }

        // Last, color the img according to trigger values:
        // 0 -------- duffuse (defult)
        // 1 -------- white
        // 2 -------- norm_flat
        // 3 -------- norm_gouraud
        // 4 -------- norm_bary
        // 5 -------- norm_gouraud_z
        // 6 -------- norm_bary_z
        
        if(trigger == 0){
            std::cout<<"coloring the img using its diffuse color"<<std::endl;
            ppmdraw(img, czb,crawzlist,xinterseclist,bboxlist,dcolorlist);
        }
        else if(trigger == 1){
            std::cout<<"coloring the img using white color"<<std::endl;
            ppmdraw(img, czb,crawzlist,xinterseclist,bboxlist,wcolorlist);
        }
        else if(trigger == 2){
            std::cout<<"coloring the img using norm_flat color"<<std::endl;
            ppmdraw(img, czb,crawzlist,xinterseclist,bboxlist,ncolorlist);
        }
        else if(trigger == 3){
            std::cout<<"coloring the img using norm_gouraud color"<<std::endl;
            ppmadvancedraw(img, zb,rawzlist,xinterseclist, bboxlist,ngcolorlist);
        }
        else if(trigger == 4){
            std::cout<<"coloring the img using norm_bary color"<<std::endl;
            ppmadvancedraw(img, zb,rawzlist,xinterseclist, bboxlist,nbcolorlist);
        }
        else if(trigger == 5){
            std::cout<<"coloring the img using norm_gouraud_z color"<<std::endl;
            ppmadvancedraw(img, czb,crawzlist,xinterseclist, bboxlist,ngzcolorlist);
        }
        else if(trigger == 6){
            std::cout<<"coloring the img using norm_bary_z color"<<std::endl;
            ppmadvancedraw(img, czb,crawzlist,xinterseclist, bboxlist,nbzcolorlist);
        }

        return img;
    }
}
