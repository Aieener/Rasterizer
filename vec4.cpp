//vect.cpp
//vectors
//Author: Yuding Ai
//Date: 2017.02.07
#include "vec4.h"
using namespace std;
#include <math.h>
#include <sstream>
// Global variables for exception strings to be thrown
extern const string EXC_OUT_OF_BOUNDS = "Index-out-of-bounds exception!";
///----------------------------------------------------------------------
/// Constructors
///----------------------------------------------------------------------
vec4::vec4(){
    //initialize the vector to (0,0,0,0)
    for (int i = 0; i<4; i++){
        data[i] = 0;
    }
}

vec4::vec4(float x, float y, float z, float w){
    data[0]=x;
    data[1]=y;
    data[2]=z;
    data[3]=w;
}

vec4::vec4(const vec4 &v2){
    //assign values from v2 to data
    data[0] = v2[0];
    data[1] = v2[1];
    data[2] = v2[2];
    data[3] = v2[3];
}

///----------------------------------------------------------------------
/// Getters/Setters
///----------------------------------------------------------------------		

float vec4::operator[](unsigned int index) const{
    float result = 0;
    result = data[index];
    // I've noticed that it would also work if we just return data[index] in one line of code
    // since this operator will only been called after a const vec4 object
    // so it won't altered anyway
    return result;
}

float& vec4::operator[](unsigned int index){
    return data[index];
}

float vec4::operator()(unsigned int index) const{
    //check bounds
    //unsigned int is always > 0
    if (index >3){
        throw EXC_OUT_OF_BOUNDS;
    }
    else{
        float result;
        result = data[index];
        return result;
    }
}

float& vec4::operator()(unsigned int index) {
    //check bounds
    if (index >3){
        throw EXC_OUT_OF_BOUNDS;
    }
    else{
        return data[index];
    }
}

///----------------------------------------------------------------------
/// Operator Methods
///----------------------------------------------------------------------		

/// Assign v2 to this and return a reference to this
vec4 &vec4::operator=(const vec4 &v2){
    //overwrite data
    for (int i = 0; i<4; i++){
        data[i]=v2[i];
    }
    //return a reference to this
    return *this;
}

/// Test for equality
bool vec4::operator==(const vec4 &v2) const{
    //check if all elements in the two vec are equal
    //are different: for example: (1,2,3,4) == (2,4,6,8)
    
    float eps = 10E-4; // SET PRECISION TO 10^-4
    for(int i = 0; i<4; i++){
//        if(abs(data[i]- v2[i])>eps){
        if(data[i]- v2[i]>eps || data[i]- v2[i]<-eps){
            return false;
        }
    }
    return true;
}

/// Test for inequality
bool vec4::operator!=(const vec4 &v2) const{
    if (*this == v2){
        return false;
    }
    else{
        return true;
    }
}
/// Arithmetic
vec4& vec4::operator+=(const vec4 &v2){
    // concatenate all values from v2 to this
    for(int i = 0; i<4; i++){
        data[i] += v2[i];
    }
    return *this;
}

vec4& vec4::operator-=(const vec4 &v2){
    for(int i = 0; i<4; i++){
        //substract the elements in v2
        //and overwrite the data
        data[i] -= v2[i];
    }
    return *this;
}

vec4& vec4::operator*=(float c){
    for(int i = 0; i<4; i++){
        //multiply all elements by c
        data[i] = data[i]*c;
    }
    return *this;
}                

vec4& vec4::operator/=(float c){
    for(int i = 0; i<4; i++){
        //divide all elements by c
        data[i] = data[i]/c;
    }
    return *this;
}                

vec4  vec4::operator+(const vec4 &v2) const{
    // create a new vec4 and assign the value of this to it
    // so basically I'm making a copy of this here;
    vec4 newvect = *this; 
    for (int i = 0; i<4;i++){
        newvect[i] = newvect[i] + v2[i];
    }
    return newvect;
}
vec4  vec4::operator-(const vec4 &v2) const{
    // create a new vec4 and assign the value of this to it
    // so basically I'm making a copy of this here;
    vec4 newvect = *this; 
    for (int i = 0; i<4;i++){
        newvect[i] = newvect[i] - v2[i];
    }
    return newvect;
}
vec4  vec4::operator*(float c) const{
    // create a new vec4 and assign the value of this to it
    // so basically I'm making a copy of this here;
    vec4 newvect = *this; 
    for (int i = 0; i<4;i++){
        newvect[i] = newvect[i]*c;
    }
    return newvect;
}          
vec4  vec4::operator/(float c) const {
    // create a new vec4 and assign the value of this to it
    // so basically I'm making a copy of this here;
    vec4 newvect = *this; 
    for (int i = 0; i<4;i++){
        newvect[i] = newvect[i]/c;
    }
    return newvect;   
}        

///----------------------------------------------------------------------
/// Other Methods
///----------------------------------------------------------------------		
/// Returns the geometric length of the input vector

float vec4::length() const{
    float len;
    // then calculate the geometric length
    len =sqrt(pow(data[0],2) + pow(data[1],2) + pow(data[2],2) + pow(data[3],2));
    return len;
}

/// return a new vec4 that is a normalized (unit-length) version of this one
vec4 vec4::normalize() const{
    // make a copy of 'this' vector
    vec4 newvector = *this;
    vec4 v0;
    if (newvector == v0){
        return v0;
    }
    else{
        float len = newvector.length();
        newvector/= len;
        return newvector;
    }
}

/// noralize this vector in place
void vec4::norm(){
    *this = normalize();
}

///----------------------------------------------------------------------
/// Other Functions (not part of the vec4 class)
///----------------------------------------------------------------------		

/// Dot Product
float dot(const vec4 &v1, const vec4 &v2){
    float result = 0;
    for(int i = 0; i< 4; i++){
        result += (v1[i])*(v2[i]);
    }
    return result;
}

/// Cross Product
// for this one, we use x,y,z elements onlu
vec4 cross(const vec4 &v1, const vec4 &v2) {
    vec4 crosvec;
    for(int i = 0; i<3;i++){
        crosvec[i] = v1[(i+1)%3]*v2[(i+2)%3] - (v2[(i+1)%3])*(v1[(i+2)%3]);
    }
    crosvec[3] = 0; //'the 4th element of the resultant vector should be 0'
    return crosvec;
}    
                                               
 /// Scalar Multiplication (c * v)
vec4 operator*(float c, const vec4 &v){
    vec4 result;
    for (int i = 0; i < 4; i++) { 
        result[i] = c*v[i]; 
    } 
    return result;
}

/// Prints the vector to a stream in a nice format for integration with cout
std::ostream &operator<<(std::ostream &o, const vec4 &v){
    stringstream st;
    st << "[ "<<v[0]<<", "<<v[1]<<", "<<v[2]<<", "<<v[3]<<" ]";
    string mg = st.str();
    o <<mg;
    return o;
}

                                             




