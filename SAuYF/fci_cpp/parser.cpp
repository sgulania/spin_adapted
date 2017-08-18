/******************************************************************************
 * 
 * parser.cpp
 *
 * This file has functions to parse inputs from various file formats for use in
 * other parts of the code base e.g. ci_matrix.cpp
 *
 * JDWhitfield 2017
 *
 ***
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


// standard C headers
#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>  /* for tolower() */

// standard C++ headers
#include <cmath>
#include <sstream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string> 

#include "parser.h"

//defined via extern.  Otherwise, I keep getting duplicates in the linker
const char *PeriodicTable[18]={
    "h",
    "he",
    "li","be",
    "b","c","n","o","f","ne",
    "na","mg",
    "al","si","p","s","cl","ar"};




const char * USAGE_STRING="Usage: integrals basis_functions_file nuc_field_filename\n  or   integrals basis_functions_file number_of_elecs";

//main function
/*
int 
main(int argc, char *argv[])
{

    using std::cout;
    using std::cerr;
    using std::endl;


    // ****************************
    // * Parse commandline inputs *
    // ****************************
    if (argc > 1)
    {	
	    if(!strcmp(argv[1],"-h"))
	    {
            cout << USAGE_STRING;
		    return 0;
	    }
    }

    const auto basis_fname=(argc>1) ? argv[1] : "basis_funcs";

    if(argc > 2)
    {
        //pull nuclear file or number of electrons
    }


    auto shells = parse_basisfile(basis_fname);

    for(auto s : shells)
        cout << s << '\n';

    


    libint2::initialize();

    // **** one electron integrals ****
    // compute overlap integrals
    auto S = compute_1body_ints(shells, libint2::OneBodyEngine::overlap);

    // dump to file




    // compute kinetic-energy integrals
    auto T = compute_1body_ints(shells, libint2::OneBodyEngine::kinetic);


    // dump to file
   
    
    // **** two electron integrals ****
    //setup storage array
    //dump to file
    //construct the electron repulsion integrals engine
    libint2::TwoBodyEngine<libint2::Coulomb> engine(max_nprim(shells), max_l(shells), 0);
    //compute integrals, chemist's notation
    double intval = *engine.compute(shells[s1],shells[s2],shells[s3],shells[s4]);

    cout << "\n\tOverlap Integrals:\n";
    cout << S << "\n\n";
    cout << "\n\tKinetic Integrals:\n";
    cout << T << "\n\n";
    for(auto a=0; a!=shells.size(); a++)
        for(auto b=0; b<=a; b++)
            for(auto c=0; c<=a; c++)
                for(auto d=0; d<= (a==c ? b : c) ; d++)
                    printf("%d: (%d %d | %d %d)= %lf\n",ijkl2intindex(a,b,c,d),a,b,c,d,
                        *engine.compute(shells[a],shells[b],shells[c],shells[d]));



    //compute integrals, chemist's notation
    //double intval = *engine.compute(shells[s1],shells[s2],shells[s3],shells[s4]);
        
    libint2::finalize();

    


    return 0;
}
*/






//My parsers 
std::vector<libint2::Atom>
parse_nucfile(const char* nuc_fname, int& N)
{
    std::vector<libint2::Atom> atoms;

    //variables for reading
    std::ifstream nuc_file(nuc_fname);
    std::string buffer;

    //variables for storing data read in
    int    Z;
    double x,y,z;
    libint2::Atom   atom;

    
    if(!nuc_file.is_open())
    {
        std::cout << "Unable to open file\n";
        std::cerr << "Error code: " << std::strerror(errno) << "\n";
        throw;
    }


    
    while(getline(nuc_file,buffer))
    {
        if(sscanf(buffer.c_str(),"electrons=%i",&N)==1)
        {
            continue;
        }

        if(sscanf(buffer.c_str(),"N %i",&N)==1)
        {
            std::cout << "read N";
            continue;
        }

        if(sscanf(buffer.c_str(),"%i %lf %lf %lf",&Z,&x,&y,&z)!=4)
        {
            std::cout << "\nError with formatting of nuclear file\n\n";
            throw;
        }
        else
        {
            //pull info
            atom.x=x;
            atom.y=y;
            atom.z=z;
            atom.atomic_number=Z;
            
            //save data
            atoms.push_back(atom);
        }
    }
    
    return atoms;
}

libint2::BasisSet 
parse_basisfile(const char * basis_fname)
{
    std::vector<libint2::Shell> shells;

    //variables for reading
    std::ifstream basis_file(basis_fname);
    double x=0;
    double y=0;
    double z=0;
    std::array<double,3> coords;
    double ai;
    double ci;
    double ci2;
    //int    nshells;
    int    nmatches;
    int    L,L2=-99;


    //Variables for making shell
    // note that libint can only handles one contraction at a time
    //      so we need to split up the contraction into separate 
    //      shells
    libint2::Shell::Contraction contr;
    libint2::Shell::Contraction contr2;
    std::vector<double> alphas;

    //for storing the contraction coeffs
    std::vector<double> d;
    std::vector<double> d2;

    //variables for reading in
    char   IType[25];
    int    NGauss;
    double Sc=0;


    std::string buffer;
    // start reading
    if(basis_file.is_open())
    {
        while(getline(basis_file,buffer))
        {
            //skip over empty lines
            if(buffer.empty()) continue;

            //if you see stars, the coordinate should follow
            if(buffer.find("****")!=std::string::npos)
            {

                //pull next line
                if(!getline(basis_file,buffer)) break;

                //load x,y,z coords
                nmatches=sscanf(buffer.data(),"%lf %lf %lf",&x,&y,&z);
                if(nmatches!=3)
                {
                    std::cout << "Error in file format (reading coords)" << "\n";
                    throw;
                }
                
                
                coords={{x,y,z}};
            
                //on to the next line 
                if(!getline(basis_file,buffer)) break;
            }

            //it should be a basis function header
            nmatches=sscanf(buffer.data(),"%s %i %lf",IType,&NGauss,&Sc);
            if(nmatches!=3)
            {
                std::cout << "Error in file format (reading basis header)" << "\n";
                throw;
            }
           
            //reset containers
            alphas.clear();
            d.clear();
            d2.clear();

            //What is the orbital type
            L2=-99;
            switch(strlen(IType))
            {
                case 2:
                    switch(tolower(IType[1]))
                    {
                    case 'p':
                        L2=1;
                        break;
                    case 'd':
                        L2=2;
                        break;
                    default:
                        std::cout << "Error basis type specification unknown (L2)"<<IType;
                        throw;
                    }
                    contr2.l=L2;
                    contr2.pure=false;
                case 1:
                    switch(tolower(IType[0]))
                    {
                    case 's':
                        L=0;
                        break;
                    case 'p':
                        L=1;
                        break;
                    case 'd':
                        L=2;
                        break;
                    case 'f':
                        L=3;
                        break;
                    default:
                        std::cout << "Error basis type specification unknown (L1)" << IType;
                        throw;
                    }
                    break;

                default:
                    std::cout << "Error basis type specification unknown" << IType;
                    throw;
            }
            contr.l=L;
            contr.pure=false;

            //read the coeffients
            for(int j=0; j<NGauss; j++)
            {

                //pull next line
                if(!getline(basis_file,buffer)){ std::cout<<"Empty"; break; }

                //parse line with one coefficent
                if(L2!=-99)
                {
                    nmatches=sscanf(buffer.data(),"%lf %lf %lf",&ai,&ci,&ci2);
                    if(nmatches!=3)
                    {
                        std::cout << "Error reading basis set data\n";
                        throw;
                    }
                    d2.push_back(ci2);

                }
                else
                {
                    nmatches=sscanf(buffer.data(),"%lf %lf",&ai,&ci);
                    if(nmatches!=2)
                    {
                        std::cout << "Error reading basis set data\n";
                        throw;
                    }
                }

                //store into appropriate arrays
                d.push_back(ci);
                alphas.push_back(ai);
            }


            //contraction coeffs
            std::vector<libint2::Shell::Contraction> contractions;

            //put in the coeffs
            contr.coeff=d;
            contractions.push_back(contr);

            //form the shell and push it into the list
            libint2::Shell shell(alphas,contractions,coords);

            //shell.alpha =alphas;
            //shell.O     =coords;
            //shell.contr =contractions;

            shells.push_back(shell);


            if(L2!=-99)
            {
                //clearn contractions so we can put in the next shell
                contractions.clear();

                contr2.coeff=d2;
                contractions.push_back(contr2);

                libint2::Shell shell(alphas,contractions,coords);

                shells.push_back(shell);
                
            }


            
        }
        /*
        for( auto s : shells)
            std::cout << s << "\n";
        */

    }
    else
    {
        std::cout << "\n\nUnable to open the file\n";
        throw ;//exception because file not open
    }

    //close the file
    basis_file.close();

    return libint2::BasisSet(shells);
}



