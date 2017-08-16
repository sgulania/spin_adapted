//
// //start of header guard
#ifndef LIBINT_INTERFACE
#define LIBINT_INTERFACE

#include <libint2.hpp>
#include "parser.h"
#include "ci_matrix.h"

void
get_libints(const libint2::BasisSet& basisset, const std::vector<libint2::Atom>& atoms,
            Matrix &S, Matrix &hV, Matrix &hT, vector<double>& AOInts);

Matrix
compute_1body_ints(libint2::Operator, 
		   const libint2::BasisSet, const std::vector<libint2::Atom>& atoms=std::vector<libint2::Atom>());

std::vector<double>
compute_2body_ints(const libint2::BasisSet bfs);

#endif
// //end of header guard
