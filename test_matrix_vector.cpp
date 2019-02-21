#include "../lib/vec3.h"
#include "../lib/vec4.h"
#include "../lib/mat3.h"
#include "../lib/mat4.h"

#include <iostream>
using namespace std;


int main(int argc, char** argv) {
	std::cout << "Hello!!!\nThis is a test file..." << std::endl << std::endl;

	vec3 v0;
	cout << "v0 (constructed by default): " << v0 << endl;
    vec3 v1(1, 2, 3);
	cout << "v1 (constructed by initialization): " << v1 << endl;
	cout << "v1 * 2: " << v1 * 2 << endl;
	cout << "v1 / 2: " << v1 / 2 << endl;
	v1 *= 2;
	cout << "v1 *= 2: " << v1 << endl;
	v1 /= 2;
	cout << "v1 /= 2: " << v1 << endl;
	vec3 v2(v1);
	cout << "v2 (constructed as a copy): " << v2 << endl;
	cout << "v1 * v2 is: " << v1.dot(v2) << endl;	
	cout << "v1 x v2 is: " << v1.cross(v2) << endl;
	cout << "-v1: " << -v1 << endl;
	cout << "v1 + v2: " << v1 + v2 << endl;
	cout << "v1 - v2: " << v1 - v2 << endl;
	v1 += v2;
	cout << "v1 += v2: " << v1 << endl;
	v1 -= v2;
	cout << "v1 -= v2: " << v1 << endl;

	vec3 v3;
	cout << "Input v3 with each two of elements separated by a ws: " << endl;
    cin >> v3;
	cout << "The squared length of v3 is: " << v3.length2() << endl;
	cout << "The length of v3 is: " << v3.length() << endl;
	v3.normalize();
	cout << "The normalized v3 is: " << v3 << endl;

	vec4 v(1.2f, 2.3f, -1.7f);
	cout << "v4 (homogeneous coord.): "<< v << endl;

	mat3 matrix1(1, 2, 3,
		4, 5, 6,
		7, 8, 9);

	mat3 matrix2(9, 8, 7,
		6, 5, 4,
		3, 2, 1);

	cout << "matrix1: " << endl << matrix1 << endl;
	cout << "matrix2: " << endl << matrix2 << endl;
	cout << "-matrix1: " << endl << -matrix1 << endl;
	cout << "matrix1 * 2: " << endl << 2 * matrix1 << endl;
	matrix1.zeros();
	cout << "matrix1.zeros(): " << endl << matrix1 << endl;
	matrix1.identity();
	cout << "matrix1.identity(): " << endl << matrix1 << endl;
	cout << "matrix1.transpose(): " << endl << matrix1.transpose() << endl;
	cout << "matrix1 + matrix2: " << endl << matrix1 + matrix2 << endl;
	cout << "matrix1 - matrix2: " << endl << matrix1 - matrix2 << endl;
	cout << "matrix1 * matrix2: " << endl << matrix1 * matrix2 << endl;
	matrix1 += matrix2;
	cout << "matrix1 += matrix2: " << endl << matrix1 << endl;
	matrix1 -= matrix2;
	cout << "matrix1 -= matrix2: " << endl << matrix1 << endl;
	matrix1 *= matrix2;
	cout << "matrix1 *= matrix2: " << endl << matrix1 << endl;
	matrix1 /= 2;
	cout << "matrix1 /= 2: " << endl << matrix1 << endl;
	cout << "matrix1 * v1: " << endl << matrix1 * v1 << endl;


	return 0;
}

