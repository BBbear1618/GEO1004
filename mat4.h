#ifndef GEO1004_MAT_H
#define GEO1004_MAT_H

#include "vec4.h"

class mat4
{
public:
	// default constructor
	mat4();

	// initialized constructor
	mat4(
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33);

	// copy constructor
	mat4(const mat4& other);

	// destructor
	~mat4();

	// overwrite this mat4 with all zero entries
	void zeros();

	// overwrite this mat4 with identity
	void identity();

	// compute the transposed mat4
	mat4 transpose() const;

	// operators -- mat4-mat4
	mat4 operator+(const mat4& other) const; // addition
	mat4 operator-(const mat4& other) const; // subtraction
	mat4 operator*(const mat4& other) const; // multiplication
	mat4 operator-() const;			 // negation

	const mat4& operator+=(const mat4& other); // cumulative addition
	const mat4& operator-=(const mat4& other); // cumulative subtraction
	const mat4& operator*=(const mat4& other); // cumulative multiplication

	// operators -- mat4-vector
	vec4 operator*(const vec4& v);	// mat4-vector product

	// operators -- mat4-scalar
	const mat4& operator*=(float scalar);	// mat4-scalar product
	const mat4& operator/=(float scalar);	// mat4-scalar division
	mat4 operator*(float scalar);		// mat4-scalar product
	mat4 operator/(float scalar);		// mat4-scalar division

// assignment operator
	const mat4& operator=(const mat4& other);	// assignment

// access components
	float& operator()(int i, int j);	// RW access to components
	float operator()(int i, int j) const;	// RO access to components										// cast to float*

protected:
	float	m_data[16];
};


inline mat4 operator*(float scalar, const mat4& M) {
	// TODO -- multiply each component of M with scalar, in a new mat4. return new mat4

    mat4 newmatrix;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++){
            newmatrix(i, j) = M(i, j) * scalar;
        }
    }

    return newmatrix;
}

inline std::ostream& operator<<(std::ostream& out, const mat4& M) {
	// output a mat4 row-by-row to the "out" stream
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			out << M(i, j) << " ";
		}
		out << std::endl;
	}
	return out;
}

inline std::istream& operator>>(std::istream& in, mat4& M) {
	// TODO: read a mat4 row-by-row from the "in" stream

    in >> M(0,0) >> M(0,1) >> M(0,2) >> M(0,3) >> M(1,0) >> M(1,1) >> M(1,2) >> M(1,3) >>
    M(2,0) >> M(2,1) >> M(2,2) >> M(2,3) >> M(3,0) >> M(3,1) >> M(3,2) >> M(3,3);

	return in;
}



inline mat4::mat4() {
	// TODO -- initialize m_data with 0s

    for (int i = 0; i < 16; i++) {
        m_data[i] = 0;
    }
}

inline mat4::mat4(float m00, float m01, float m02, float m03,
	float m10, float m11, float m12, float m13,
	float m20, float m21, float m22, float m23,
	float m30, float m31, float m32, float m33)
{
	// TODO -- initialize m_data with the provided components.
    m_data[0] = m00;
    m_data[1] = m01;
    m_data[2] = m02;
    m_data[3] = m03;
    m_data[4] = m10;
    m_data[5] = m11;
    m_data[6] = m12;
    m_data[7] = m13;
    m_data[8] = m20;
    m_data[9] = m21;
    m_data[10] = m22;
    m_data[11] = m23;
    m_data[12] = m30;
    m_data[13] = m31;
    m_data[14] = m32;
    m_data[15] = m33;

}


inline mat4::mat4(const mat4& other) {
	// TODO -- copy other to (*this) component by component

    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            this->operator()(i,j) = other(i,j);
        }
    }
}


inline mat4::~mat4() {
}


inline void mat4::zeros() {
	// TODO -- overwrite this mat4 with an identity mat4.

    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){

        this->operator()(i,j) = 0;
        }
    }
}


inline void mat4::identity() {
	// TODO -- overwrite this mat4 with an identity mat4.

    this->zeros();

    for (int i = 0; i < 4; i++)
    {
        this->operator()(i, i) = 1;
    }
}


inline mat4 mat4::transpose() const {
	// TODO -- compute the transpose of this mat4 in a new mat4 and return.

    mat4 newmatrix;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++){
            newmatrix(i,j) = this->operator()(j,i);
        }
    }
    return newmatrix;
}


inline mat4 mat4::operator+(const mat4& other) const {
	// TODO -- compute a new mat4 (*this)+other, return the new mat4

    mat4 newmatrix;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++){
            newmatrix(i,j) = this->operator()(i,j) + other(i,j);
        }
    }

    return newmatrix;
}


inline mat4 mat4::operator-(const mat4& other) const {
	// TODO -- compute a new mat4 (*this)-other, return the new mat4

    mat4 newmatrix;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++){
            newmatrix(i,j) = this->operator()(i,j) - other(i,j);
        }
    }

    return newmatrix;
}


inline mat4 mat4::operator*(const mat4& other) const {
	// TODO -- compute a new mat4 (*this) * other, return the new mat4

    mat4 newmatrix;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                newmatrix(i,j) += this->operator()(i,k) * other(k,j);
            }
        }
    }

    return newmatrix;
}


inline mat4 mat4::operator-() const {
	// TODO -- compute a new mat4 -(*this), return the new mat4

    mat4 newmatrix;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            newmatrix(i,j) = -this->operator()(i,j);
        }
    }

    return newmatrix;
}


inline const mat4& mat4::operator+=(const mat4& other) {
	// TODO -- add other to this mat4

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            this->operator()(i,j) += other(i,j);
        }
    }

    return *this;
}


inline const mat4& mat4::operator-=(const mat4& other) {
	// TODO -- subtract other from this mat4

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            this->operator()(i,j) -= other(i,j);
        }
    }

    return *this;
}


inline const mat4& mat4::operator*=(const mat4& other) {
	// TODO -- replace this mat4 by (*this) * other. Make sure you do not overwrite elements that you still need.
	//		   You may use mat4::operator*()

    *this = *this * other;

    return *this;
}


inline vec4 mat4::operator*(const vec4& v) {
	// TODO -- compute the mat4-vector product (*this) * v and return the result

    vec4 newvector;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            newvector(i) += v(j) * this->operator()(i,j);
        }
    }

    return newvector;
}


inline const mat4& mat4::operator*=(float scalar) {
	// TODO -- multiply each mat4 component by scalar.

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++){
            this->operator()(i,j) *= scalar;
        }
    }

    return *this;
}


inline const mat4& mat4::operator/=(float scalar) {
	assert("mat4::operator/= -- invalid argument" && scalar != 0);
	// TODO -- divide each mat4 component by scalar.

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++){
            this->operator()(i,j) /= scalar;
        }
    }
    return *this;
}


inline mat4 mat4::operator*(float scalar) {
	// TODO -- compute a new mat4 (*this) * scalar.

    mat4 newmatrix;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            newmatrix(i,j) = this->operator()(i,j) * scalar;
        }
    }

    return newmatrix;
}


inline mat4 mat4::operator/(float scalar) {
	assert("mat4::operator/ -- invalid argument" && scalar != 0);
	// TODO -- divide each mat4 component by scalar and store in a new mat4. return the new mat4.

    mat4 newmatrix;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            newmatrix(i,j) = this->operator()(i,j) / scalar;
        }
    }

    return newmatrix;
}


inline const mat4& mat4::operator=(const mat4& other) {
	// TODO -- overwrite each component in this mat4 by the matching component in other

    for (int i =0; i<4; i++){
        for (int j = 0; j<4; j++){
            this->operator()(i,j) = other(i,j);
        }
    }
    return *this;
}



inline float& mat4::operator()(int i, int j) {
	assert("mat4::operator() -- invalid arguments" && i < 4 && j < 4);
	return m_data[4 * i + j];
}


inline float mat4::operator()(int i, int j) const {
	assert("mat4::operator() const -- invalid arguments" && i < 4 && j < 4);
	return m_data[4 * i + j];
}



#endif
