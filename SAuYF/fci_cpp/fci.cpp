/******************************************************************************
 * 
 * fci.cpp
 *
 * This program loads basis functions and nuclear field information to 
 * form a Hamiltonian and compute its eigenvalues.  
 * 
 * JDWhitfield 
 * Dartmouth, 2016-2017
 *
 *
 *      This program is free software: you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation, either version 2 of the License, or
 *      (at your option) any later version.
 * 
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 * 
 *      You should have received a copy of the GNU General Public License
 *      along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 *****************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<assert.h>
#include<iomanip>
#include<random>
#include<chrono>
#include<algorithm>
#include"parser.h"
#include"ci_matrix.h"
#include"libint_interface.h"

int
main(int argc, char *argv[])
{ 
    using std::chrono::high_resolution_clock;

    std::cout.precision(10);
    std::cout << std::scientific;

    int debug=2;

    //timing variables
    //std::chrono::system_clock::time_point start_time;
    std::chrono::duration<double> time_elapsed;

    // ****************************
    // * Parse commandline inputs *
    // ****************************
    if (argc > 1)
    {	
	    if(!strcmp(argv[1],"-h"))
	    {
            std::cout << "fci [basis_func] [nuc_field]";
            std::cout << "\n\tDefault files are used when called with no parameters.\n\n";
		    return 0;
	    }
    }

    const auto basis_fname=(argc>1) ? argv[1] : "basis_funcs";
    const auto nuc_fname  =(argc>2) ? argv[2] : "nuc_field";


    // ****************************
    // * Output files             *
    // ****************************
    std::fstream fonebody;
    system("touch ham_ov.dat");
    system("rm ham_ov.dat");
    system("touch ham_ov.dat");
    fonebody.open("ham_ov.dat");

    


   /*
     Format ham_ov.dat-
     number of basis
     Ov matrix
     Core Hamiltonian matrix
     Nuclear Repulsion
    */

    std::fstream ftwobody;
    system("touch 2_ele.dat");
    system("rm 2_ele.dat");
    system("touch 2_ele.dat");
    ftwobody.open("2_ele.dat");
    /*
     Format 2_ele.dat-
     number of reduced integral
     i,j,k,l, ee(i,j,k,l)
    */

    //increase the printout precision
    fonebody.precision(10);
    ftwobody.precision(10);
    fonebody << std::scientific << std::showpos;
    ftwobody << std::scientific;

    // ****************************
    // * Parse data files         *
    // ****************************
    libint2::BasisSet shells = parse_basisfile(basis_fname);


    if(debug)
    {
        printf("Basis set\n");
        for( auto s : shells)
            std::cout << s << "\n";
    }

    //number of spin orbitals in the basis set
    int M=2*shells.nbf();

    fonebody << M << "\n";
    
    //number of electrons
    int N=-99;


    /*
    //number of 2e- integrals
    int m=(M/2) - 1;
    int nints=term4(m,m,m,m)+1;
    */

    auto atoms = parse_nucfile(nuc_fname,N);

    if(N==-99)
    {
        std::cout << "Nelec not set\n";
        throw;
    }
    
    Matrix S, h;
    std::vector<double> AOInts;

    //FCI matrix dimension
    int D=nchoosek(M,N);

    //libints
    auto start_time = high_resolution_clock::now();
    Matrix hT,hV;
    get_libints(shells,atoms,S,hV,hT,AOInts);
    h=hT+hV;
    time_elapsed=(high_resolution_clock::now()-start_time);
    std::cout << "Took " << time_elapsed.count() << " seconds to get AO Ints\n";

    if(debug)
    {
        std::cout << "S\n" << S << "\n\n";
        fonebody  << S << "\n";
        std::cout << "h\n" << h << "\n\n";
        fonebody  << h << "\n";
        std::cout << "AO 2body integrals (chemist's notation)\n";
	
	//temporary header to be overwritten by the number of integrals
        ftwobody  << "AO       \n";
	int nints=0;
        for(auto p=0; p!=M/2; p++) //unique integral labels, looping scheme from libint
            for(auto q=0; q<=p; q++)
                for(auto r=0; r<=p; r++)
                    for(auto s=0; s<= (p==r ? q : r) ; s++)
                        if(std::abs(AOInts[term4(p,q,r,s)]) > 1e-5)
			{
				std::cout << term4(p,q,r,s) << ": ["
					  << p << q << "|" << r << s << "] = " << AOInts[term4(p,q,r,s)] << std::endl;

 			        ftwobody << p << " " << q << " " << r << " " << s << " " << AOInts[term4(p,q,r,s)] << std::endl;
				nints++;
			}
	//go back to write the number of integrals at the top
        ftwobody.seekp(0);
	//header
	ftwobody  << nints;
    }




    Matrix h2;
    std::vector<double> MOInts;
    if(Matrix(S).isIdentity(1e-9))
    {
        //nothing changes if S is the identity matrix
        h2=h;
        MOInts=AOInts;
        std::cout<<"We don't need to orthogonalize\n";
    }
    else
    {// *     ORTHOGONALIZATION    *


        /***********************************************************/ 
        // Transformation matricies
        //
        //Basis transform to an orthogonal space
        //C = C(W) = XW for any WW^\dag =\id
        //

        Eigen::SelfAdjointEigenSolver<Matrix> Eigensystem(S);
        Matrix s=Matrix(Eigensystem.eigenvalues());
        //check which eigenvalues are non-trivial to avoid linear dependence
        //then compute s^{-1/2}
        Matrix shalf;
        shalf.resize(M/2,M/2);    // set size
        shalf.fill(0);            // make sure all elements are zero
        //make sure its a column vector
        if(s.rows()<s.cols())
            s.transposeInPlace();

        for(int j=0; j<s.rows(); j++)
        {
            if(std::abs(s(j,0))>1e-6)
                shalf(j,j)=1/sqrt(s(j,0));
            else
                shalf(j,j)=0;
        }

        //orthogonalization transform
        //Here we're just using W=\id
        Matrix C;
        bool CANONICAL=true;
        if(CANONICAL)
            C=Matrix(Eigensystem.eigenvectors()*shalf);
        else
            C=Matrix(Eigensystem.eigenvectors()*shalf*Eigensystem.eigenvectors().adjoint());

        //time the transformations
        start_time = high_resolution_clock::now();

        h2       =Matrix(C.adjoint())*h*C;   //One-body transform
        MOInts   = transform4(M,C,AOInts);//two-body transform

        time_elapsed = high_resolution_clock::now() - start_time;
        std::cout << "Transformation to (canonical) orthogonal basis took " << time_elapsed.count() << "ms \n";

    }

    if(debug)
    {
        std::cout << "h\n";
        std::cout << h2 << "\n\n";
        std::cout << "Orthogonalized Integrals\n";
        for(auto m : MOInts)
            std::cout << m << "\n";
    }

    //LAPACK
    start_time=high_resolution_clock::now();

    double* w=ci_eigenvals(M,N,MOInts,h2,debug);

    time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "Making matrix and getting eigenvalues took " << time_elapsed.count() 
              << "s to get eigenvalues.\n";

    std::cout << "\nEigenvalues( "<< D << " ) :\n --------- \n";

    std::cout.precision(10);
    std::cout << std::scientific;
    for(int i=0; i<D; i++)
        std::cout << w[i]  << "\n";

    
    // ****************************
    // * Close output files       *
    // ****************************

    fonebody.close();



    return 0;
}
