/*
 * ppm_example.c
 * HW_1  CIS 560
 * Author: Yuding Ai
 * Student ID: 31295008
 * Date: Jan 26 2017
 */ 

#include "ppm_filters.h"
/** int main(int argc, char **argv) { */
    /** // read in the input image */
    /** printf("Reading, %s!\n", argv[1]); */
    /** img_t *img = read_ppm(argv[1]); */
    /** for(int i = 1; i<argc;i++){ */
    /**     //check if there is a tag of "-grayscale" */
    /**     if(!strcmp(argv[i],"-grayscale")){  */
    /**         grays_scale(img); */
    /**         printf("grayscaled the img!\n"); */
    /**     }   */
    /**     //check if there is a tag of "-flipimage" */
    /**     if(!strcmp(argv[i],"-flip")){ */
    /**         flipimage(img); */
    /**         printf("flipped the img!\n"); */
    /**     }  */
    /**     //check if there is a tag of "-flopimage" */
    /**     if(!strcmp(argv[i],"-flop")){ */
    /**         flopimage(img); */
    /**         printf("flopped the img!\n"); */
    /**     } */
    /**     //check if there is a tag of "-transpose" */
    /**     if(!strcmp(argv[i],"-transpose")){ */
    /**         transpose(img); */
    /**         printf("transposed the img!\n"); */
    /**     } */
    /**     //check if there is a tag of "-boxblur" */
    /**     if(!strcmp(argv[i],"-boxblur")){ */
    /**         double n = 0; // if no input after -boxblur then n=0 */
    /**         n = atof (argv[i+1]); */
    /**         boxblur(img,n); */
    /**         printf("boxblurred the img by n = %f!\n",n); */
    /**     } */
    /**     //check if there is a tag of "-median" */
    /**     if(!strcmp(argv[i],"-median")){ */
    /**         double n = 0; */
    /**         n = atof (argv[i+1]); */
    /**         median(img,n); */
    /**         printf("medianblurred the img by n = %f!\n",n); */
    /**     } */
    /**     //check if there is a tag of "-gaussian" */
    /**     if(!strcmp(argv[i],"-gaussian")){ */
    /**         double n,s; */
    /**         n = s = 0; */
    /**         n = atof(argv[i+1]); */
    /**         s = atof(argv[i+2]); */
    /**         gaussian(img,n,s); */
    /**         printf("gaussianblurred the img by n = %f, s = %f!\n",n,s); */
    /**     } */
    /**     //check if there is a tag of "-gaussian" */
    /**     if(!strcmp(argv[i],"-sobel")){ */
    /**         sobel(img); */
    /**         printf("sobeled the img!\n"); */
    /**     } */
    /**  */
    /**     //check if there is a tag of "-scale" */
    /**     if(!strcmp(argv[i],"-scale")){ */
    /**         double factor = 1; // set default factor = 1 */
    /**         int width,height; */
    /**         factor = atof(argv [i+1]); */
    /**         width = (int)(factor * img->w);  */
    /**         height = (int)(factor * img->h);  */
    /**         resize(img,width,height); */
    /**         printf("scale the img by %f!\n",factor); */
    /**     } */
    /**     //check if there is a tag of "-size" */
    /**     if(!strcmp(argv[i],"-size")){ */
    /**         int width,height; */
    /**         width = (int)(img->w); // set default size to the original size */
    /**         height = (int)(img->h); */
    /**         width = atoi(argv[i+1]); */
    /**         height = atoi(argv[i+2]); */
    /**         resize(img,width,height); */
    /**         printf("size the img to %d by %d!\n",width,height); */
    /**     } */
    /** } */
    /**  */
    /** // now write the image */
    /** printf("Writing %s\n", argv[argc-1]); */
    /** write_ppm(img, argv[argc-1]); */
    /**  */
    /** // free up memory */
    /** destroy_img(&img); // &img is the address in memory where the img variable is stored */
    /** // Since img is of type (img *), &img is of type (img **) */
    /**  */
    /** return 0; */
/** } */

// Create a new image of specified size.
// Implementation of functions that modifies ppm img
 
namespace filters {
// make the img gray
void grays_scale(img_t *img){
    int size = img->w*img->h;
    for(pixel_t *p = img->data;p<img->data +size; p++){
        double y;
        y = p->b*0.114 + p->r*0.299 +p->g*0.587; // get the value for gray
        int value = abs((int) y); // round it back to pos interger
        p->r = p->g = p->b = (unsigned char)(value); // asign the gray value to the r g b component of each pixel 
    }
}

// flip the img
void flipimage(img_t *img){
    for (int row =0; row< img->h; row++){
        pixel_t *p = img ->data + row * img->w; // declare a pointer p to  point to an array of pixel_t on the img 
        for(int col = 0; col < img ->w/2; col++){
            //for each row
            //the location of left pixel is [col]
            //the location of the corresponding paired right pixel is [img->w -1 -col]
            pixel_t temp = p[col];              //store the left pixel at each row into temp
            p[col] = p[img->w -1 - col];        //assign the value of the right pixel into left ones 
            p[img->w -1 - col] = temp;          //assign the value of temp into the right ones
        }
    }
}

// flop the img
void flopimage(img_t *img){
    pixel_t *p = img->data;
    for(int col = 0; col<img->w; col++){
        for(int row = 0; row < img->h/2;row++){
            //for each col
            //the locatoin of upper pixel is [row*img->w + col]
            //the loc of the corresponding paired lower pixel is [img->w*(img->h -1-row) + col] (counting backwards)
            pixel_t temp = p[row*img->w + col ]; // store the upper pixel at each col into temp
            p[row*img->w + col] = p[img->w*(img->h -1-row) + col];// assign the value of lower pixel into upper pixel 
            p[img->w*(img->h -1-row) + col] = temp;  // assign the value of temp into the lower pixel
        } 
    }
}

// transpose the img
void transpose(img_t *img){
    //first to desize the w and h
    int bytesize, old_h,old_w;
    old_h = img->h;
    old_w = img->w;
    img->h = old_w;
    img->w = old_h;
    bytesize = img->w* img->h*sizeof(pixel_t);
    
    //allocate a block of memory to the buffer to temporary hold
    //the whole amount of pixels, like what img->data did, but in
    //transpose version
    pixel_t *buffer = (pixel_t *)malloc(bytesize);
    
    // declare a pointer points to the original data.
    pixel_t *p = img->data;
    int i = 0; // i serve as a counter.
    for(int col = 0;col<old_w;col++){
        for(int row = 0; row < old_h;row++){
            //each loop, assign pixels in cols 
            //into buffer hence store the tranposed
            //version of this img into buffer
            buffer[i] = p[old_w*row + col];
            i++;
        }
    }
    //Copy the data from buffer to img->data as to overwrite it
    memcpy(img->data,buffer,bytesize);
    //free the memory from buffer
    free(buffer);
}

// blut the img through boxblur
void boxblur(img_t *img,double n){
    int bytesize = img->w*img->h*sizeof(pixel_t);  // declare bytesize of the img
    pixel_t *buffer = (pixel_t *)malloc(bytesize); //allocte a block of mem to buffer to tempary hold the pixel value

    //declare a pointer points to the address of buffer
    pixel_t *pb =buffer;
    //declare a pointer points to the img->data
    pixel_t *p = img->data;
    //boxlength will be 2n+1
    double boxlength = 2*n+1;
    //boxsize 
    double boxsize = pow(boxlength,2);

    for(int row =0; row<img->h; row++){
        for(int col = 0; col < img->w; col ++){
            // if hits the boundary set to green
            if (col - n <0 || row -n <0 || col + n >img->w-1 || row + n >img->h-1){
                pb ->r = pb->b = 0;
                pb ->g = 255;
            }  
            // otherwise, calculate the average value and assign it to the target pixel
            else{
                double avr,avg,avb;
                avr = avg =avb = 0;
                for (int i = 0; i<boxlength; i++){
                    for(int j = 0; j < boxlength; j++){
                        // caculate sum of r,g,b over the box: 
                        // current loc: location of the target pixel p[row*img->w + col]
                        // left top box loc: current loc - n*img->w -(2*n)/2
                        // boxloc = left top + i*img->w + j
                        int boxloc;
                        boxloc  = row*(img->w) + col -n*(img->w) -n + i*img->w+ j;
                        // add to sum
                        avr +=p[boxloc].r ;
                        avb +=p[boxloc].b ;
                        avg +=p[boxloc].g ;
                    }
                }
                //normalize avr, avg, ave thus to get average value
                avr = avr/boxsize;
                avg = avg/boxsize;
                avb = avb/boxsize;
                //round to pos int
                int rvalue = abs((int) avr); // round it back to pos interger
                int gvalue = abs((int) avg); // round it back to pos interger
                int bvalue = abs((int) avb); // round it back to pos interger
                //assign the ave value to target pixel onto buffer
                pb->r = (unsigned char)(rvalue);
                pb->g = (unsigned char)(gvalue);
                pb->b = (unsigned char)(bvalue);
            }
            pb++; //increment pb to blur the next pixel
        }
    }
    //Copy the data from buffer to img->data as to overwrite it
    memcpy(img->data,buffer,bytesize);
    //free the memory from buffer
    free(buffer);
}

// blur the img through median blur
void median(img_t *img,double n){
    int bytesize = img->w*img->h*sizeof(pixel_t);  // declare bytesize of the img
    pixel_t *buffer = (pixel_t *)malloc(bytesize); //allocte a block of mem to buffer to tempary hold the pixel value

    //declare a pointer points to the address of buffer
    pixel_t *pb =buffer;
    //declare a pointer points to the img->data
    pixel_t *p = img->data;
    //boxlength will be 2n+1
    int boxlength = 2*n+1;
    //boxsize 
    int boxsize = pow(boxlength,2);

    for(int row =0; row<img->h; row++){
        for(int col = 0; col < img->w; col ++){
            // if hits the boundary set to green
            if (col - n <0 || row -n <0 || col + n >img->w-1 || row + n >img->h-1){
                pb ->r = pb->b = 0;
                pb ->g = 255;
            }  
            // otherwise, calculate the mean value and assign it to the target pixel
            else{
                int mr,mg,mb;
                int MR[boxsize]; // an array of int as to store the pixel value to find the mean
                int MG[boxsize];
                int MB[boxsize];
                int counter = 0;
                mr = mg =mb = 0;
                for (int i = 0; i<boxlength; i++){
                    for(int j = 0; j < boxlength; j++){
                        // caculate sum of r,g,b over the box: 
                        // current loc: location of the target pixel p[row*img->w + col]
                        // left top box loc: current loc - n*img->w -(2*n)/2
                        // boxloc = left top + i*img->w + j
                        int boxloc;
                        boxloc  = row*(img->w) + col -n*(img->w) -n + i*img->w+ j;
                        // assign to the int array
                        MR[counter] = p[boxloc].r;
                        MG[counter] = p[boxloc].g;
                        MB[counter] = p[boxloc].b;
                        counter++;
                    }
                }
                //find the mean value in MR, MG, MB and assign it to buffer
                
                double temp; // temp double use for swap
                // sort MR,MG and MB into ascending order
                for(int i = 0; i<boxsize -1; i++){
                    for (int j = i+1;j<boxsize;j++){
                        if (MR[j]<MR[i]){
                            //swap
                            temp = MR[i];
                            MR[i] = MR[j];
                            MR[j] = temp;
                        }
                        if (MG[j]<MG[i]){
                            //swap
                            temp = MG[i];
                            MG[i] = MG[j];
                            MG[j] = temp;
                        }
                         if (MB[j]<MB[i]){
                            //swap
                            temp = MB[i];
                            MB[i] = MB[j];
                            MB[j] = temp;
                        }                       
                    }
                }
                // find the mean value (locate at the middle of each array now)
                mr = MR[boxsize/2]; // since boxsize is always odd number thus boxsize/2 gives mean
                mb = MB[boxsize/2];
                mg = MG[boxsize/2];
                //round it back to int
                int rvalue = abs((int) mr); // round it back to pos interger
                int bvalue = abs((int) mb); // round it back to pos interger
                int gvalue = abs((int) mg); // round it back to pos interger

                //assign the ave value to target pixel onto buffer
                pb->r = (unsigned char)(rvalue);
                pb->g = (unsigned char)(gvalue);
                pb->b = (unsigned char)(bvalue);
            }
            pb++; //increment pb to blur the next pixel
        }
    }
    //Copy the data from buffer to img->data as to overwrite it
    memcpy(img->data,buffer,bytesize);
    //free the memory from buffer
    free(buffer);
}

// blur the img through gaussian blur
void gaussian(img_t *img,double n,double s){
    //assign weight exp(-x*x/(2*s*s)) to each pixel
    
    int bytesize = img->w*img->h*sizeof(pixel_t);  // declare bytesize of the img
    pixel_t *buffer = (pixel_t *)malloc(bytesize); //allocte a block of mem to buffer to tempary hold the pixel value
    int blurlength = 2*n +1; //essentially the boxlength

    //declare a pointer points to the address of buffer
    pixel_t *pb =buffer;
    //declare a pointer points to the img->data
    pixel_t *p = img->data;
    //First blur it by row
    for(int row = 0;row < img->h; row++){
        for (int col = 0; col < img->w;col++){
            // if hits the boundary set to green
            if (col - n <0 || row -n <0 || col + n >img->w-1 || row + n >img->h-1){
                pb ->r = pb->b = 0;
                pb ->g = 255;
            }  
            // otherwise, Gaussian blur it by row 
            else{
                double gr,gb,gg,weight;
                gr=gb=gg =weight= 0;
                for(int i = 0; i<blurlength; i++){
                    // evaluate the sum of gaussian r,g,b over the bulrrow: 
                    // current loc: location of the target pixel p[row*img->w + col]
                    // edgeloc = current loc - n
                    // target loc = edgeloc + i
                    int loc;
                    loc = row*img->w + col -n + i; 
                    gr += exp(-(n-i)*(n-i)/(2*s*s))*p[loc].r; 
                    gg += exp(-(n-i)*(n-i)/(2*s*s))*p[loc].g; 
                    gb += exp(-(n-i)*(n-i)/(2*s*s))*p[loc].b; 
                    weight += exp(-(n-i)*(n-i)/(2*s*s));
                }
                //now find the weighted ave and assign it to buffer
                gr = gr/weight;
                gg = gg/weight;
                gb = gb/weight;
                pb->r = (unsigned char)(gr);
                pb->g = (unsigned char)(gg);
                pb->b = (unsigned char)(gb);
            }
            pb++;
        } 
    }
    //Next blur it by col
    
    memcpy(img->data,buffer,bytesize);
    //Declare a new buffer and allocate mem 
    pixel_t *buffer2 = (pixel_t *)malloc(bytesize); //allocte a block of mem to buffer to tempary hold the pixel value

    for(int col = 0;col < img->w; col++){
        // declare a pointer points to buffer2
        pixel_t *pb2 =buffer2;
        pb2 +=col;
        for (int row = 0; row < img->h;row++){
            // if hits the boundary set to green
            if (col - n <0 || row -n <0 || col + n >img->w-1 || row + n >img->h-1){
                pb2 ->r = pb2->b = 0;
                pb2 ->g = 255;
            }  
            // otherwise, Gaussian blur the img by col 
            else{
                double gr,gb,gg,weight;
                gr=gb=gg =weight= 0;
                for(int i = 0; i<blurlength; i++){
                    // evaluate the sum of gaussian r,g,b over the blurcol: 
                    // current loc: location of the target pixel p[row*img->w + col]
                    // edgeloc (the upper edge) = current loc - n*img->w
                    // target loc = edgeloc + i*img->w
                    int loc;
                    loc = row*img->w + col -n*img->w + i*img->w; 
                    gr += exp(-(n-i)*(n-i)/(2*s*s))*p[loc].r; 
                    gg += exp(-(n-i)*(n-i)/(2*s*s))*p[loc].g; 
                    gb += exp(-(n-i)*(n-i)/(2*s*s))*p[loc].b; 
                    weight += exp(-(n-i)*(n-i)/(2*s*s));
                }
                //now find the weighted ave and assign it to buffer
                gr = gr/weight;
                gg = gg/weight;
                gb = gb/weight;
                //round it back to int
                int rvalue = abs((int) gr); // round it back to pos interger
                int gvalue = abs((int) gg); // round it back to pos interger
                int bvalue = abs((int) gb); // round it back to pos interger
                //assign the round value to buffer2
                pb2->r = (unsigned char)(rvalue);
                pb2->g = (unsigned char)(gvalue);
                pb2->b = (unsigned char)(bvalue);
            }
            pb2+=img->w;
        } 
    }
    //Copy the data from buffer2 to img->data as to overwrite it
    memcpy(img->data,buffer2,bytesize);
    //free the memory from buffer
    free(buffer);
    free(buffer2);
}


void sobel(img_t *img){
    //declare two kernels as sobel filter: xfilter and yfilter
    grays_scale(img);
    double xfilter[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
    double yfilter[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};
    // since the filter size is 3x3, we filter the img in a unit box of 3x3
    // so basicallg we are assuming "n=1" here
    // Furthermore, as we did before, here I set the boundary to bright green
    int bytesize = img->w*img->h*sizeof(pixel_t);  // declare bytesize of the img
    pixel_t *buffer = (pixel_t *)malloc(bytesize); //allocte a block of mem to buffer to tempary hold the pixel value

    //declare a pointer points to the address of buffer
    pixel_t *pb =buffer;
    //declare a pointer points to the img->data
    pixel_t *p = img->data;
    //boxlength will be 2n+1 = 3
    double boxlength = 3;
    
    // First scan through the img and find the max possible value the sobel filter could produce and assign it to variable "maxedge"
    double maxedge=0;
    for(int row =0; row<img->h; row++){
        for(int col = 0; col < img->w; col ++){
            // if hits boundary, do nothing
            if (col - 1 <0 || row -1 <0 || col + 1 >img->w-1 || row + 1 >img->h-1){
                // do nothing here
            }  
            else{           
                //declare double edge to hold the filtered value
                double  edg, xedg, yedg;
                edg = xedg = yedg = 0;
                for (int i = 0; i<boxlength; i++){
                    for(int j = 0; j < boxlength; j++){
                        // caculate sober filtered value of r,g,b over the box: 
                        // current loc: location of the target pixel p[row*img->w + col]
                        // left top box loc: current loc - img->w -1
                        // loc = left top + i*img->w + j
                        int loc;
                        loc  = row*(img->w) + col -img->w -1 + i*img->w+ j;
                        //operates onto a grayscale pic
                        xedg +=  xfilter[i][j]*p[loc].r;
                        yedg +=  yfilter[i][j]*p[loc].r;

                        edg = sqrt(xedg*xedg + yedg*yedg);
                    }
                }
                if (edg >=maxedge){
                    maxedge = edg;
                }    
            }
        }
    }
    /*printf("%f\n",maxedge);*/

    // Now do the filtering
    for(int row =0; row<img->h; row++){
        for(int col = 0; col < img->w; col ++){
            // if hits the boundary set to green
            if (col - 1 <0 || row -1 <0 || col + 1 >img->w-1 || row + 1 >img->h-1){
                pb ->r = pb->b = 0;
                pb ->g = 255;
            }  
            else{
                //sobel filter the img
                //declare double edge to hold the filtered value
                double  edge, xedge, yedge;
                edge = xedge = yedge = 0;
                for (int i = 0; i<boxlength; i++){
                    for(int j = 0; j < boxlength; j++){
                        // caculate sober filtered value of r,g,b over the box: 
                        // current loc: location of the target pixel p[row*img->w + col]
                        // left top box loc: current loc - img->w -1
                        // loc = left top + i*img->w + j
                        int loc;
                        loc  = row*(img->w) + col -img->w -1 + i*img->w+ j;
                        //operates onto a grayscale pic
                        xedge +=  xfilter[i][j]*p[loc].r;

                        yedge +=  yfilter[i][j]*p[loc].r;
                        
                        edge = sqrt(xedge*xedge + yedge*yedge);
                   }
                }

                // normalize it by dividing largest edge value and multiply 255
                edge = 255*(edge/maxedge); 

                int value = abs((int) edge); // round it back to pos interger

                /*if (value>=255){*/
                    /*//debug to check if value is out of scope*/
                    /*printf("%d\n",value);*/
                /*}*/

                //assigne the edge value to buffer
                pb->r = (unsigned char)(value);
                pb->g = (unsigned char)(value);
                pb->b = (unsigned char)(value);
            } 
            pb++;
        }
    }
    
    //Copy the data from buffer to img->data as to overwrite it
    memcpy(img->data,buffer,bytesize);
    //free the memory from buffer
    free(buffer);
}

// resize the img
// the following idea of implementation is imspired by an online resourse  http://tech-algorithm.com/articles/bilinear-image-scaling/ 
// It works almost fine only with one tiny issue, which is that when I resize an img into the its original size
// the img will became a little bit blurred, only a little bit.
void resize(img_t *img, int width,int height){
    double delta_x, delta_y;
    delta_x = delta_y = 0;
    double ratio_x = ((double)(img->w-1 ))/(width);
    double ratio_y = ((double)(img->h-1 ))/(height);

    int newbytesize = width*height*sizeof(pixel_t);  // declare newbytesize of the newimg
    pixel_t *buffer = (pixel_t *)malloc(newbytesize); //allocte a block of mem to buffer to tempary hold the pixel value

    //declare a pointer points to the address of buffer
    pixel_t *pb =buffer;
    //declare a pointer points to the img->data
    pixel_t *p = img->data;   
    /*int count = 0;*/
    for(int row =0; row<height; row++){
        for(int col = 0; col < width; col ++){
            // loc: location of the target pixel p[y*width + x]
            int loc ,Ar,Ag,Ab,Br,Bg,Bb,Cr,Cg,Cb,Dr,Dg,Db,r,g,b; // where A B C D stands for the 4 corner of unit box on the original img
            // and x,y maps to the coordingnate onto the original img for the cooresponding row and col
            int x,y;

            x = (int)(ratio_x * col);
            y = (int)(ratio_y * row);

            delta_y = (ratio_y *row) -y;
            delta_x = (ratio_x *col) -x; // this method is borrored from an online 
            //resourse: http://tech-algorithm.com/articles/bilinear-image-scaling/    and I DO NOT fully understand.
          
            loc = y*img->w + x;
            //assign rgb values to the four corners A B C D
            Ar=p[loc].r;
            Ag=p[loc].g;
            Ab=p[loc].b;
            Br=p[loc+1].r;
            Bg=p[loc+1].g;
            Bb=p[loc+1].b;
            Cr=p[loc+img->w].r;
            Cg=p[loc+img->w].g;
            Cb=p[loc+img->w].b;
            Dr=p[loc+img->w+1].r;
            Dg=p[loc+img->w+1].g;
            Db=p[loc+img->w+1].b;

            // Nest, Following the formula from wikipedia to bilinear extrapolate r g b values for the new img
            // f(x,y) = f(A_x,A_y)*(1-delta_x)*(1-delta_y) + f(B_x,B_y)*delta_x*(1-delta_y) + f(C_x,C_y)*(1-delta_x)*delta_y) + f(D_x,D_y)*delta_x*delta_y
            r = (int)(
                    Ar*(1-delta_x)*(1-delta_y) + Br*(delta_x)*(1-delta_y) + Cr*(delta_y)*(1-delta_x) + Dr*(delta_x*delta_y)
                    );
            g = (int)(
                    Ag*(1-delta_x)*(1-delta_y) + Bg*(delta_x)*(1-delta_y) + Cg*(delta_y)*(1-delta_x) + Dg*(delta_x*delta_y)
                    );
            b = (int)(
                    Ab*(1-delta_x)*(1-delta_y) + Bb*(delta_x)*(1-delta_y) + Cb*(delta_y)*(1-delta_x) + Db*(delta_x*delta_y)
                    );
            // assign the r g b value to the buffer
            pb->r = (unsigned char)(r);
            pb->g = (unsigned char)(g);
            pb->b = (unsigned char)(b);
            pb++;

        }
    }
    img->w = (int)(width);
    img->h = (int)(height);
    // free img->data 
    free(img->data);
    // and then re malloc it for the new size
    img->data = (pixel_t *)malloc(newbytesize); //allocte a block of mem to buffer to tempary hold the pixel value

    //Copy the data from buffer to img->data as to overwrite it
    memcpy(img->data,buffer,newbytesize);
    //free the memory from buffer
    free(buffer);
}

}









