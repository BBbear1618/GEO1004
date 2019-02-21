#ifndef GEO1004_MAT3_H
#define GEO1004_MAT3_H

#include "vec3.h"

class mat3
{
public:
	// default constructor
	mat3();

	// initialized constructor
	mat3(
		float m00, float m01, float m02,
		float m10, float m11, float m12,
		float m20, float m21, float m22);

	// copy constructor
	mat3(const mat3& other);

	// destructor
	~mat3();

	// overwrite this mat3 with all zero entries
	void zeros();

	// overwrite this mat3 with identity
	void identity();

	// compute the transposed mat3
	mat3 transpose() const;

	// operators -- mat3-mat3
	mat3 operator+(const mat3& other) const; // addition
	mat3 operator-(const mat3& other) const; // subtraction
	mat3 operator*(const mat3& other) const; // multiplication
	mat3 operator-() const;			 // negation

	const mat3& operator+=(const mat3& other); // cumulative addition
	const mat3& operator-=(const mat3& other); // cumulative subtraction
	const mat3& operator*=(const mat3& other); // cumulative multiplication

	// operators -- mat3-vector
	vec3 operator*(const vec3& v);	// mat3-vector product

	// operators -- mat3-scalar
	const mat3& operator*=(float scalar);	// mat3-scalar product
	const mat3& operator/=(float scalar);	// mat3-scalar division
	mat3 operator*(float scalar);		// mat3-scalar product
	mat3 operator/(float scalar);		// mat3-scalar division

// assignment operator
	const mat3& operator=(const mat3& other);	// assignment

// access components
	float& operator()(int i, int j);	// RW access to components
	float operator()(int i, int j) const;	// RO access to components										// cast to float*

protected:
    // We changed this to make it more intuitive for a 3x3 matrix
    float	m_data[9];
};

inline std::ostream& operator<<(std::ostream& out, const mat3& M) {
    // output a mat3 row-by-row to the "out" stream
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            out << M(i, j) << " ";
        }
        out << std::endl;
    }
    return out;
}

inline mat3 operator*(float scalar, const mat3& M) {
    // DONE -- multiply each component of M with scalar, in a new mat3. return new mat3

    mat3 newmatrix;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++){
            newmatrix(i, j) = M(i, j) * scalar;
        }
    }

    return newmatrix;
}



inline std::istream& operator>>(std::istream& in, mat3& M) {
	// TODO: read a mat3 row-by-row from the "in" stream

    in >> M(0,0) >> M(0,1) >> M(0,2) >> M(1,0) >> M(1,1) >> M(1,2) >> M(2,0) >> M(2,1) >> M(2,2);

	return in;
}



inline mat3::mat3() {
    // DONE -- initialize m_data with 0s
    for (int i = 0; i < 9; i++) {
        m_data[i] = 0;
    }
}

inline mat3::mat3(float m00, float m01, float m02,
	float m10, float m11, float m12,
	float m20, float m21, float m22)
{
    m_data[0] = m00;
    m_data[1] = m01;
    m_data[2] = m02;
    m_data[3] = m10;
    m_data[4] = m11;
    m_data[5] = m12;
    m_data[6] = m20;
    m_data[7] = m21;
    m_data[8] = m22;
}


inline mat3::mat3(const mat3& other) {
    // DONE -- copy other to (*this) component by component
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            this->operator()(i,j) = other(i,j);
        }
    }
}


inline mat3::~mat3() {
}


inline void mat3::zeros() {
    // DONE -- overwrite this mat3 with an identity mat3.

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){

        this->operator()(i,j) = 0;
        }
    }
}


inline void mat3::identity() {
    // DONE -- overwrite this mat3 with an identity mat3.

    this->zeros();

    for (int i = 0; i < 3; i++)
    {
        this->operator()(i, i) = 1;
    }
}


inline mat3 mat3::transpose() const {
    // DONE -- compute the transpose of this mat3 in a new mat3 and return.

    mat3 newmatrix;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++){
            newmatrix(i,j) = this->operator()(j,i);
        }
    }
    return newmatrix;
}


inline mat3 mat3::operator+(const mat3& other) const {
    // DONE -- compute a new mat3 (*this)+other, return the new mat3

    mat3 newmatrix;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++){
            newmatrix(i,j) = this->operator()(i,j) + other(i,j);
        }
    }

    return newmatrix;
}


inline mat3 mat3::operator-(const mat3& other) const {
    // DONE -- compute a new mat3 (*this)-other, return the new mat3

    mat3 newmatrix;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++){
            newmatrix(i,j) = this->operator()(i,j) - other(i,j);
        }
    }

    return newmatrix;
}


inline mat3 mat3::operator*(const mat3& other) const {
    // DONE -- compute a new mat3 (*this) * other, return the new mat3

    mat3 newmatrix;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                newmatrix(i,j) += this->operator()(i,k) * other(k,j);
            }
        }
    }

    return newmatrix;
}


inline mat3 mat3::operator-() const {
	// TODO -- compute a new mat3 -(*this), return the new mat3

    mat3 newmatrix;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            newmatrix(i,j) = -this->operator()(i,j);
        }
    }

    return newmatrix;
}


inline const mat3& mat3::operator+=(const mat3& other) {
    // DONE -- add other to this mat3

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            this->operator()(i,j) += other(i,j);
        }
    }

	return *this;
}


inline const mat3& mat3::operator-=(const mat3& other) {
    // DONE -- subtract other from this mat3

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            this->operator()(i,j) -= other(i,j);
        }
    }

	return *this;
}


inline const mat3& mat3::operator*=(const mat3& other) {
    // DONE -- replace this mat3 by (*this) * other. Make sure you do not overwrite elements that you still need.
	//		   You may use mat3::operator*()

    *this = *this * other;

	return *this;
}


inline vec3 mat3::operator*(const vec3& v) {
    // DONE -- compute the mat3-vector product (*this) * v and return the result

    vec3 newvector;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            newvector(i) += v(j) * this->operator()(i,j);
        }
    }

    return newvector;
}


inline const mat3& mat3::operator*=(float scalar) {
    // DONE -- multiply each mat3 component by scalar.

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++){
            this->operator()(i,j) *= scalar;
        }
    }

	return *this;
}

inline const mat3& mat3::operator/=(float scalar) {
	assert("mat3::operator/= -- invalid argument" && scalar != 0);
    // DONE -- divide each mat3 component by scalar.

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++){
            this->operator()(i,j) /= scalar;
        }
    }
	return *this;
}


inline mat3 mat3::operator*(float scalar) {
    // DONE -- compute a new mat3 (*this) * scalar.

    mat3 newmatrix;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            newmatrix(i,j) = this->operator()(i,j) * scalar;
        }
    }

    return newmatrix;
}


inline mat3 mat3::operator/(float scalar) {
	assert("mat3::operator/ -- invalid argument" && scalar != 0);
    // DONE -- divide each mat3 component by scalar and store in a new mat3. return the new mat3.

    mat3 newmatrix;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            newmatrix(i,j) = this->operator()(i,j) / scalar;
        }
    }

    return newmatrix;
}


inline const mat3& mat3::operator=(const mat3& other) {
    // DONE -- overwrite each component in this mat3 by the matching component in other

    for (int i =0; i<3; i++){
        for (int j = 0; j<3; j++){
            this->operator()(i,j) = other(i,j);
        }
    }
	return *this;
}

// We changed this to make it more intuitive for a 3x3 matrix

inline float& mat3::operator()(int i, int j) {
    assert("mat3::operator() -- invalid arguments" && i < 3 && j < 3);
    return m_data[3 * i + j];
}


inline float mat3::operator()(int i, int j) const {
    assert("mat3::operator() const -- invalid arguments" && i < 3 && j < 3);
    return m_data[3 * i + j];
}



#endif
