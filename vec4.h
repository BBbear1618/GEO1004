#ifndef GEO1004_VEC4_H
#define GEO1004_VEC4_H	// TODO -- why is this here? Be prepared to explain to the TA.


#include<iostream>
#include<cfloat>
#include<cmath>
#include<cassert>

#include "vec3.h"

/// four-dimensional vector class

class vec4 {
public:
	// default constructor: initialize all elements to be zero
	vec4();

	// initialized constructor
	vec4(float x, float y, float z, float w = 1.0f);

	// copy constructor
	vec4(const vec4& other);

	// construct a 4D vector from a 3D vector: w component set to 1 by default
	// useful for representing a 3D point in homogeneous coordinates
	vec4(const vec3& v);

	// destructor
	~vec4();

	// suqared length of the vector
	float length2() const;

	// length of the vector
	float length() const;

	// normalize to a unit vector
	void normalize();

	// compute the dotproduct between this and other
	float dot(const vec4& other) const;

	// compute the crossproduct between the first three components of this and other
	vec4 cross(const vec4& other) const;

	// operators -- vector-vector
	const vec4& operator+=(const vec4& other);		// cumulative addition
	const vec4& operator-=(const vec4& other);		// cumulative subtraction
	vec4 operator-() const;							// negation (unary operator)
	vec4 operator+(const vec4& other) const;		// addition
	vec4 operator-(const vec4& other) const;		// subtraction

	// operators -- vector-scalar
	const vec4& operator*=(float scalar);			// multiplication by scalar
	const vec4& operator/=(float scalar);			// division by scalar
	vec4 operator*(float scalar);					// vector times scalar
	vec4 operator/(float scalar);					// vector divided by scalar
	const vec4 operator=(const vec4& other);		// assignment operator

	// element access
	float& operator()(int index);       // RW access to element: return the index_th element
	float operator()(int index) const;	// RO access to element: return the index_th element

protected:
	float	m_data[4];	// data array
};


// TODO -- why is the following needed? Be prepared to explain to the TA
inline vec4 operator*(float scalar, const vec4& v) {
	// TODO -- multiply each component of v with scalar, in a new vector. return new vector

	return vec4(scalar * v(0), scalar * v(1), scalar * v(2), scalar * v(3)); // replace this line
}


inline std::ostream& operator<<(std::ostream& out, const vec4& v) {
	for (int i = 0; i < 4; ++i)
		out << v(i) << " ";
	return out;
}


inline std::istream& operator>>(std::istream& in, vec4& v) {
	// TODO: read a vector component-wise from the "in" stream
	in >> v(0) >> v(1) >> v(2) >> v(3);
	return in;
}


inline vec4::vec4() {
	// TODO -- initialize m_data with 0s
	m_data[0] = 0;
	m_data[1] = 0;
	m_data[2] = 0;
	m_data[3] = 0;
}


inline vec4::vec4(float x, float y, float z, float w) {
	// TODO -- initialize m_data with x,y,z,w. What does the "= float(0)" in the class definition do? Be prepared to explain to the TA.
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
	m_data[3] = w;
}


inline vec4::vec4(const vec3& v) {
	// TODO -- construct a 4D vector from a 3D vector: w component set to 1 by default
	// useful for representing a 3D point in homogeneous coordinates
	m_data[0] = v(0);
	m_data[1] = v(1);
	m_data[2] = v(2);
	m_data[3] = 1;
}


inline vec4::vec4(const vec4& other) {
	// TODO -- copy contents of other to this vector
	m_data[0] = other(0);
	m_data[1] = other(1);
	m_data[2] = other(2);
	m_data[3] = other(3);
}


inline vec4::~vec4() {
	// TODO -- is there anything to do here? Be prepared to explain to the TA.
}


inline float vec4::length2() const {
	// TODO -- compute squared length of the vector

	return (float)(m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2] + m_data[3] * m_data[3]); // replace this line
}


inline float vec4::length() const {
	// TODO -- compute the length of the vector

	return (float)sqrt(length2()); // replace this line
}


inline void vec4::normalize() {
	// TODO -- if length()==0, do nothing, otherwise normalize the vector
	float len = length();
	if (len != 0)
		*this /= len;
}


inline float vec4::dot(const vec4& other) const {
	// TODO -- compute dot product between this vector and other

	return m_data[0] * other(0) + m_data[1] * other(1) + m_data[2] * other(2) + m_data[3] * other(3); // replace this line.
}


inline vec4 vec4::cross(const vec4& other) const {
	// TODO -- compute the crossproduct between the first three components of this and other

	return vec4(m_data[1] * other(2) - m_data[2] * other(1), m_data[2] * other(0) - m_data[0] * other(2), m_data[0] * other(1) - m_data[1] * other(0)); // replace this line
}


inline const vec4& vec4::operator+=(const vec4& other) {
	// TODO -- add other to this vector component-wise, store in this vector
	m_data[0] += other(0);
	m_data[1] += other(1);
	m_data[2] += other(2);
	m_data[3] += other(3);
	return *this; // TODO -- why would you return a reference to *this? Be prepared to explain to the TA.
}


inline const vec4& vec4::operator-=(const vec4& other) {
	// TODO -- subtract other from this vector component-wise, store in this vector
	m_data[0] -= other(0);
	m_data[1] -= other(1);
	m_data[2] -= other(2);
	m_data[3] -= other(3);
	return *this;
}


inline vec4 vec4::operator-() const {
	// TODO -- why can't we return a reference, here? Be prepared to explain to the TA.
	// TODO -- return a new vector in which each component is the negated component from this vector

	return vec4(-m_data[0], -m_data[1], -m_data[2], -m_data[3]); // replace this line
}


inline vec4 vec4::operator+(const vec4& other) const {
	// TODO -- return a new vector in which each component is the sum of the components of this vector and other

	return vec4(m_data[0] + other(0), m_data[1] + other(1), m_data[2] + other(2), m_data[3] + other(3)); // replace this line
}


inline vec4 vec4::operator-(const vec4& other) const {
	// TODO -- return a new vector in which each component is the difference between the components of this vector and other

	return vec4(m_data[0] - other(0), m_data[1] - other(1), m_data[2] - other(2), m_data[3] - other(3)); // replace this line
}



inline const vec4& vec4::operator*=(float scalar) {
	// TODO -- replace each component of this with the matching component of this multiplied with the scalar
	//         Make sure to convert the scalar from S to float
	m_data[0] *= scalar;
	m_data[1] *= scalar;
	m_data[2] *= scalar;
	m_data[3] *= scalar;
	return *this;
}


inline const vec4& vec4::operator/=(float scalar) {
	// TODO -- why do we convert scalar and 0 to type float and not compare it as type S? Be prepared to explain to the TA.
	assert("vec4::operator/= -- invalid argument" && float(scalar) != float(0));
	// TODO -- replace each component of this with the matching component of this multiplied with the scalar
	//         Make sure to convert the scalar from S to float
	m_data[0] /= scalar;
	m_data[1] /= scalar;
	m_data[2] /= scalar;
	m_data[3] /= scalar;
	return *this;
}



inline vec4 vec4::operator*(float scalar) {
	// TODO -- return a new vector in which each component equals the matching component of this vector times the scalar
	//		   Make sure to convert the scalar from S to float

	return vec4(m_data[0] * scalar, m_data[1] * scalar, m_data[2] * scalar, m_data[3] * scalar); // replace this line
}



inline vec4 vec4::operator/(float scalar) {
	assert("vec4::operator/ -- invalid argument" && float(scalar) != float(0));
	// TODO -- return a new vector in which each component equals the matching component of this vector divided by the scalar

	return vec4(m_data[0] / scalar, m_data[1] / scalar, m_data[2] / scalar, m_data[3] / scalar); // replace this line
}


inline const vec4 vec4::operator=(const vec4& other) {
	// TODO -- overwrite each component in this vector with the matching component of other.
	m_data[0] = other(0);
	m_data[1] = other(1);
	m_data[2] = other(2);
	m_data[3] = other(3);
	return *this;
}


inline float& vec4::operator()(int n) {
	assert("vec4::operator() -- invalid argument" && n < 4);
	return m_data[n];
}


inline float vec4::operator()(int n) const {
	assert("vec4::operator() const -- invalid argument" && n < 4);
	return m_data[n];
}



#endif
