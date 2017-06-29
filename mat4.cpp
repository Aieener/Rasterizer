//mat4.cpp
//matrix 4
//Author: Yuding Ai
//Date: 2017.02.08
#include "mat4.h"
#include <math.h>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
using namespace std;

// Global variables for exception strings to be thrown
//extern const string EXC_OUT_OF_BOUNDS = "Index-out-of-bounds exception!";
extern const string EXC_OUT_OF_BOUNDS;
///----------------------------------------------------------------------
/// Constructors
///----------------------------------------------------------------------

/// Default Constructor.  Initialize to identity matrix.
mat4::mat4(){
    // for a column-major matrix, it should be like (v0,v1,v2,v3);
    // where v# are column vectors
    vec4 v0(1,0,0,0);
    vec4 v1(0,1,0,0);
    vec4 v2(0,0,1,0);
    vec4 v3(0,0,0,1);

    data[0] =v0; 
    data[1] =v1; 
    data[2] =v2; 
    data[3] =v3; 
}

/// Initializes the diagonal values of the matrix to diag. All other values are 0.
mat4::mat4(float diag){
    vec4 v0(diag,0,0,0);
    vec4 v1(0,diag,0,0);
    vec4 v2(0,0,diag,0);
    vec4 v3(0,0,0,diag);

    data[0] =v0; 
    data[1] =v1; 
    data[2] =v2; 
    data[3] =v3; 
}

/// Initializes matrix with each vector representing a column in the matrix
mat4::mat4(const vec4 &col0, const vec4 &col1, const vec4 &col2, const vec4& col3){
    data[0] = col0;
    data[1] = col1;
    data[2] = col2;
    data[3] = col3;
}

/// copy constructor; copies values of m2 into this
mat4::mat4(const mat4 &m2){
    for(int i = 0; i<4; i++){
        data[i] = m2[i];
    }
}

///----------------------------------------------------------------------
/// Getters
///----------------------------------------------------------------------	

/// Returns the values of the column at the index
/// Does NOT check the array bound because doing so is slow
const vec4 &mat4::operator[](unsigned int index) const{
    //cout << "const []"<<endl;
    return data[index]; 
}

/// Returns a reference to the column at the index
/// Does NOT check the array bound because doing so is slow
vec4 &mat4::operator[](unsigned int index){
    //cout << "&[]" << endl;
    return data[index];
}

/// Returns the values of the column at the index
/// DOES check the array bound because doing so is slow
const vec4 &mat4::operator()(unsigned int index) const{
    //check bounds
    //unsigned int is always > 0
    if (index >3 ){
        throw EXC_OUT_OF_BOUNDS;
    }
    else{
        return data[index];
    }
}

/// Returns a reference to the column at the index
/// DOES check the array bound because doing so is slow
vec4 &mat4::operator()(unsigned int index){
    //check bounds
    if (index >3 ){
        throw EXC_OUT_OF_BOUNDS;
    }
    else{
        return data[index];
    }
}

///----------------------------------------------------------------------
/// Static Initializers
///----------------------------------------------------------------------	

/// Creates a 3-D rotation matrix.
/// Takes an angle in degrees and an axis represented by its xyz components, and outputs a 4x4 rotation matrix
mat4 mat4::rot(float angle, float x, float y, float z){
    // following the 460transforms3d.html#slide-23
    // here we apply the Rodrigue's Rotation formula
    // where the rotation axix is (x,y,z) a column vector
    float rangle = angle*M_PI/180; // since the c++ cos() and sin() takes the angle in radian, we convert degree to radian
    mat4 rotmatrix;
    //version 1: works
//    vec4 rotax(x,y,z,1);
//    for (int j =0; j<4;j++){
//        for (int i=0; i<4; i++){
//            //assign value in column major
//            //rotmatrix[j] will be the j'th col
//            //and rotmatrix[j][i] means the ith element in jth col
//            if(i == j){
//                //the diag element
//                rotmatrix[j][i] = cos(rangle) + rotax[i]*rotax[j]*(1-cos(rangle));
//            }
//            else if (i == 3 || j == 3){
//                //the 0 elements on the botton and right side
//                rotmatrix[j][i] = 0;
//            }
//            else if((((j-i)%3)+3)%3 == 1){
//                //cout <<"LOL"<<endl;
//                //the rotmatrix[2][1], rotmatrix[3][2] and rotmatrix[1][3],
//                rotmatrix[j][i] = -rotax[(j+1)%3]*sin(rangle) + rotax[(j+2)%3]*rotax[(j+3)%3]*(1-cos(rangle));
//            }
//            else if((((i-j)%3)+3)%3 == 1){
//                //the rotmatrix[1][2], rotmatrix[2][3] and rotmatrix[3][1],
//                rotmatrix[j][i] = rotax[(i+1)%3]*sin(rangle) + rotax[(i+2)%3]*rotax[(i+3)%3]*(1-cos(rangle));
//            }
//        }
//    }

    //version 2: easier to implement but basically identical to version1.
    // Both version works the same way;

    mat4 r1 = cos(rangle)*mat4(1.0);
    vec4 r21 = vec4(x*x,x*y,x*z,0.0);
    vec4 r22 = vec4(x*y,y*y,y*z,0.0);
    vec4 r23 = vec4(x*z,y*z,z*z,0.0);
    vec4 r24 = vec4(0,0,0,1);

    mat4 r2 = (1-cos(rangle))*mat4(r21,r22,r23,r24);

    vec4 r31 = vec4(0,z,-y,0.0);
    vec4 r32 = vec4(-z,0,x,0.0);
    vec4 r33 = vec4(y,-x,0,0.0);
    vec4 r34 = vec4(0,0,0,1);

    mat4 r3 = sin(rangle)*mat4(r31,r32,r33,r34);

    rotmatrix = r1+r2+r3;
    return rotmatrix;
}

/// Takes an xyz displacement and outputs a 4x4 translation matrix
mat4 mat4::trans(float x, float y, float z){
    // following the 460transforms3d.html#slide-25
    mat4 tran; // intialize a matrix as identity matrix
    //then modify it to the trans metrix by change the 4th col vector
    //to (x,y,z,1)
    vec4 fourthcol(x,y,z,1);
    tran[3] = fourthcol; 
    return tran;
}

/// Takes an xyz scale and outputs a 4x4 scale matrix
mat4 mat4::scale(float x, float y, float z){
    // following the 460transforms3d.html#slide-22
    mat4 scal; // intialize a matrix as identity matrix
    //then modify the first 3 diag elements into x,y,z
    scal[0][0] = x;
    scal[1][1] = y;
    scal[2][2] = z;
    return scal;
}

///----------------------------------------------------------------------
/// Operator Functions
///----------------------------------------------------------------------

/// Assign m2 to this and return this
mat4 &mat4::operator=(const mat4 &m2){
    //overwrite data
    for(int i =0; i<4;i++){
        data[i] = m2[i];
    }
    //return a reference to this
    return *this;
}

/// Test for equality
bool mat4::operator==(const mat4 &m2) const{
    // check if all vectors in the two matrixs are equal
    for (int i = 0; i < 4 ; i++) { 
        if (data[i]!=m2[i]) { 
            return false;
        } 
    } 
    return true;
}

/// Test for inequality
bool mat4::operator!=(const mat4 &m2) const{
    if (*this == m2) { 
        return false;
    } 
    else{
        return true;
    }
}

/// Element-wise arithmetic
/// e.g. += adds the elements of m2 to this and returns this (like regular +=)
///      +  returns a new matrix whose elements are the sums of this and v2
mat4& mat4::operator+=(const mat4 &m2){
    // concatencate all vectors in m2 to this
    for (int i = 0; i < 4; i++) { 
        data[i] += m2[i]; 
    } 
    return *this;
}

mat4& mat4::operator-=(const mat4 &m2){
    // substract all vectors in m2 to this
    for (int i = 0; i < 4; i++) { 
        data[i] -= m2[i]; 
    } 
    return *this;
}

// multiplication by a scalar
mat4& mat4::operator*=(float c){
    for (int i = 0; i < 4; i++) { 
        data[i] = data[i]*c; 
    } 
    return *this;
}                

// division by a scalar
mat4& mat4::operator/=(float c){
    for (int i = 0; i < 4; i++) { 
        data[i] = data[i]/c; 
    } 
    return *this;
}                

mat4  mat4::operator+(const mat4 &m2) const{
    // create a new mat4 and assign the value of this to it
    // so basically I'm making a copy of this here;
    mat4 newmat = *this; 
    for (int i = 0; i<4;i++){
        newmat[i] = newmat[i] + m2[i];
    }
    return newmat;
}
mat4  mat4::operator-(const mat4 &m2) const{
    mat4 newmat = *this; 
    for (int i = 0; i<4;i++){
        newmat[i] = newmat[i] - m2[i];
    }
    return newmat;
}

// multiplication by a scalar
mat4  mat4::operator*(float c) const{
    mat4 newmat = *this; 
    for (int i = 0; i<4;i++){
        newmat[i] = newmat[i]*c ;
    }
    return newmat;
}             

// division by a scalar
mat4  mat4::operator/(float c) const{
    mat4 newmat = *this; 
    for (int i = 0; i<4;i++){
        newmat[i] = newmat[i]/c ;
    }
    return newmat;
}            

/// Matrix multiplication (m1 * m2)
mat4 mat4::operator*(const mat4 &m2) const{
    mat4 result;
    for (int j = 0; j < 4;j++) { 
        for (int i = 0; i < 4; i++) { 
            // the row 'vector' from *this
            vec4 thisrow(data[0][i],data[1][i],data[2][i],data[3][i]);
            // now do the vector dot multiplication
            result[j][i] = dot(thisrow,m2[j]) ;
        }  
    } 
    return result;
}

/// Matrix/vector multiplication (m * v)
/// Assume v is a column vector (ie. a 4x1 matrix)
vec4 mat4::operator*(const vec4 &v) const{
    vec4 result;
    for (int i = 0; i < 4; i++) { 
        // the row 'vector' from *this
        vec4 thisrow(data[0][i],data[1][i],data[2][i],data[3][i]);
        result[i] = dot(thisrow,v); 
    } 
    return result;
}

///----------------------------------------------------------------------
/// Matrix Operations
///----------------------------------------------------------------------	

/// Returns the transpose of the input matrix (v_ij == v_ji)
mat4 mat4::transpose() const{
    mat4 result;
    for(int j =0; j<4;j++){
        for(int i= 0; i<4 ;i++){
            result[j][i] = data[i][j];
        }
    }
    return result;
}

/// Returns a column of the input matrix
// const version
const vec4 &mat4::col(unsigned int index) const{
    return data[index];
}

// non-const version
vec4 &mat4::col(unsigned int index){
    return data[index];
}             

///----------------------------------------------------------------------
/// other functions (not part of the mat4 class)
///----------------------------------------------------------------------		

/// Scalar multiplication (c * m)
mat4 operator*(float c, const mat4 &m){
    mat4 result(m);
    for (int i = 0; i < 4; i++) { 
        //each col vector in mat4 result muptiply by c
        //the order does not matter for scalar multiplication
        result[i] *=c;
    } 
    return result;
}

/// Vector/matrix multiplication (v * m)
/// Assume v is a row vector (ie. a 1x4 matrix)
vec4 operator*(const vec4 &v, const mat4 &m){
    vec4 result;
    for (int i = 0; i < 4; i++) { 
        result[i] = dot(v, m[i]);
    } 
    return result;
}

/// Prints the matrix to a stream in a nice format
std::ostream &operator<<(std::ostream &o, const mat4 &m){
    stringstream st;
    st<<std::setprecision(4);  // to make the matrix look good
    st<<std::showpoint;
    st<<"\n";// put a empty line above each matrix been printed
             // for better view
    // when print the matrix out, we do it by row
    for (int i = 0; i < 4; i++) { 
        st<<m[0][i]<<"    "<< m[1][i]<<"    "<<m[2][i]<<"    "<<m[3][i]<<"\n";
    } 
    string mg = st.str();
    o<<mg;
    return o;
}


