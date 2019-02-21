#include "eigen_solver.h"
#include "eigen_symmetric.h"

using namespace eigen;


EigenSolver::EigenSolver()
{

}


EigenSolver::~EigenSolver()
{

}



// compute the eigen values and eigen vectors of matrix m.
// after computation:
//  - the eigen values are stored in the member "m_eigen_values"
//  - the eigen vectors are stored in the member "m_eigen_vectors"

void EigenSolver::solve(const mat3& m) {
	// TODO -> call the function defined in "eigen_symmetric.h"
	// Please read carefully the manual of the function.
	float mat[6], eigen_value[3], eigen_vector[9];
	// Only entries with respect to the upper triangle are used
	int k = 0;
	for (int j = 0; j < 3; j++)
		for (int i = 0; i <= j; i++)
		{
			mat[k++] = m(i, j);
		}
	eigen_symmetric<float>(mat, 3, eigen_vector, eigen_value);
	for (int i = 0; i < 3; i++)
	{
		this->m_eigen_values[i] = eigen_value[i];
	}
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{
			this->m_eigen_vectors[i](j) = eigen_vector[3 * i + j];
		}
}
