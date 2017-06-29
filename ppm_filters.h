#ifndef PPM_FILTERS_H
#define PPM_FILTERS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
struct pixel_struct {
    unsigned char r, g, b; // a pixel contains three bytes, name r, g, and b
};


typedef struct pixel_struct pixel_t; // from now on "pixel_t p" is a synonym for "struct pixel_struct p"

typedef struct img_struct {
    pixel_t *data; // img is a pointer to a block of memory containing pixels, i.e. an array of pixels
    // see the discussion of argv below for an explanation of pointers
    int w, h; // image width and height
} img_t;

img_t *new_img(int w, int h);  // create a new image of specified width and height
void destroy_img(img_t **img); // delete img from memory

img_t *read_ppm(const char *fname); // read in an image in ppm format
void  write_ppm(const img_t *img, const char *fname); // write an image in ppm format
namespace filters {

// functions that modifies ppm img:
void grays_scale(img_t *img);       //done
void flipimage(img_t *img);         //done
void flopimage(img_t *img);         //done
void transpose(img_t *img);         //done
void boxblur(img_t *img,double n);  //done
void median(img_t *img,double n);   //done
void gaussian(img_t *img,double n,double s);         //done

// additional features
void sobel(img_t *img);             //done
void resize(img_t *img, int width,int height);       //done
void rotate(img_t *img,double angle);  // not implemented yet
void whirl(img_t *img, double angle);  // not implemented yet
}

#endif // PPM_FILTERS_H
