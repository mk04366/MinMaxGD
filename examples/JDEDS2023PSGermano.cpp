/**
MINMAXGD
Copyright (C) {{ 2007 }}  {{ Univeristy of Angers, L. Hardouin }}

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _WIN32
#include "../include/lminmaxgd.h"
#include "../src/gd.cpp"
#include "../src/poly.cpp"
#include "../src/serie.cpp"
#include "../src/smatrix.cpp"
#include "../src/tools.cpp"
#else
// Provide the two-argument minimum operation expected by the Windows build.
double min(double a, double b)
{
	// Return the smaller of the two floating-point arguments.
	if (a > b)
		// The second argument is smaller when the first is greater.
		return b;
	else
		// Otherwise the first argument is already the minimum.
		return a;
}
#include "..\include\lminmaxgd.h"
#include "..\src\gd.cpp"
#include "..\src\poly.cpp"
#include "..\src\serie.cpp"
#include "..\src\smatrix.cpp"
#include "..\src\tools.cpp"
using namespace std;
#endif

/////////////////////////////////////////////////////////////////////
/* Example :  Hardouin L., Cottenceau B., Lagrange S., Le Corronc E.

Performance analysis of linear systems over idempotent semiring with additive inputs
Wodes 2008, Goteborg, Sweden
	*/
/////////////////////////////////////////////////////////////////////
using namespace mmgd;

int main()
{
	// Allocate the 3-by-3 matrix of transfer-function series used as the system matrix.
	smatrix A(3, 3);
	// Reserve a matrix that will hold the Kleene star of A.
	smatrix A_star;
	// Allocate the 3-by-2 input matrix.
	smatrix B(3, 2);
	// Allocate the 1-by-3 output matrix.
	smatrix C(1, 3);
	// Reserve the state-to-input transfer matrix A_star otimes B.
	smatrix As_B(3, 2);
	// Reserve the complete input-to-output transfer matrix C otimes As_B.
	smatrix C_As_B(1, 2);
	// Hold the finite and ultimately periodic polynomial used for the reference signal.
	poly pho;
	// Hold the second polynomial used to build the reference matrix Z.
	poly z;
	// Hold the scalar reference signal as a 1-by-1 matrix.
	smatrix Z;
	// Hold the input obtained by left residuating Z by C_As_B.
	smatrix uopt;
	// Hold the output reconstructed from C_As_B and uopt.
	smatrix yopt;
	// Declared for an intermediate C_As_B series but not used by this example.
	serie CABpoly;

	// Reusable monomial whose values are changed by successive init calls.
	gd m;

	// Put the monomial (2,1) in the series at row 0, column 1 of A; epsilon and e are library identity/zero values.
	A(0, 1).init(epsilon, m.init(2, 1), e);
	// Put the monomial (0,3) in the series at row 1, column 0 of A.
	A(1, 0).init(epsilon, m.init(0, 3), e);
	// Put the monomial (1,2) in the series at row 1, column 2 of A.
	A(1, 2).init(epsilon, m.init(1, 2), e);
	// Put the monomial (0,4) in the series at row 2, column 1 of A.
	A(2, 1).init(epsilon, m.init(0, 4), e);

	// Set input 0 to state 0 to the zero-delay monomial (0,0).
	B(0, 0).init(epsilon, m.init(0, 0), e);
	// Set input 1 to state 1 to the zero-delay monomial (0,0).
	B(1, 1).init(epsilon, m.init(0, 0), e);
	// Connect state 2 to the single output with the zero-delay monomial (0,0).
	C(0, 2).init(epsilon, m.init(0, 0), e);

	// Compute the matrix Kleene star, which accumulates all finite paths represented by A.
	A_star = star(A);

	// Multiply A_star by B using semiring matrix multiplication to obtain state-to-input behavior.
	As_B = otimes(A_star, B);
	// Multiply C by As_B to obtain the input-to-output transfer matrix C_As_B.
	C_As_B = otimes(C, As_B);
	// Print the original system matrix.
	std::cout << " A :" << std::endl
			  << A << std::endl;
	// Print the Kleene-star matrix.
	std::cout << " A_star :" << std::endl
			  << A_star << std::endl;
	// Print the state-to-input transfer matrix.
	std::cout << " As_B :" << std::endl
			  << As_B << std::endl;
	// Print the complete input-to-output transfer matrix.
	std::cout << " C_As_B :" << std::endl
			  << C_As_B << std::endl;

	// Start pho with the monomial (0,3), then append the listed monomials with operator().
	pho.init(0, 3)(1, 4)(2, 5)(3, 9)(4, 10)(5, 11)(6, 17)(7, 18)(8, 23)(9, 24)(10, 25)(11, 26)(12, 30)(13, 31)(14, infinit);
	// Print the polynomial before the following semiring multiplication.
	std::cout << " pho :" << std::endl
			  << pho << std::endl;
	// Left-multiply pho by the monomial (0,1), shifting its exponents through the poly otimes overload.
	pho = otimes(m.init(0, 1), pho);

	// Start z with (0,14), append (1,23), (3,29), and terminate with (4,+infinity).
	z.init(0, 14)(1, 23)(3, 29)(4, infinit);
	// Build the 1-by-1 reference matrix Z from z, with epsilon as the transient part and e as the repeating monomial.
	Z(0, 0).init(epsilon, z, e);
	// Print the reference matrix used as the residuation target.
	std::cout << " Z :" << std::endl
			  << Z << std::endl;

	// Compute the greatest input matrix uopt such that C_As_B otimes uopt is below Z in the semiring order.
	uopt = lfrac(Z, C_As_B);
	// Print the residuated input.
	std::cout << " uopt :" << std::endl
			  << uopt << std::endl;

	// Reconstruct the output produced by C_As_B and the computed input uopt.
	yopt = otimes(C_As_B, uopt);
	// Print the reconstructed output for comparison with Z.
	std::cout << " yopt :" << std::endl
			  << yopt << std::endl;
}
