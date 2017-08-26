#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<assert.h>
#include<iomanip>
#include<random>
#include<chrono>
#include<algorithm>


int
PrimDiffs(const std::vector<int>& primRef, const std::vector<int>& primK,
          std::vector<int>&       virt  , std::vector<int>&        occ, 
	  std::vector<int>&         Pr  , std::vector<int>&         Pk ) ;
int 
main()
{
	using namespace std;


	vector<int> occ;
	vector<int> virt;
	vector<int> left  {1,5,3,7,2,2,6};
	vector<int> right {1,9,10,7,7,2,6};

	vector<int> Pr;
	vector<int> Pk;
	PrimDiffs(left,right,occ,virt,Pr,Pk);

	return 0;
}

int
PrimDiffs(const std::vector<int>& primRef, const std::vector<int>& primK,
          std::vector<int>&       virt  , std::vector<int>&        occ, 
	  std::vector<int>&         Pr  , std::vector<int>&         Pk ) 
{


    //empty containers
    virt.clear();
    occ.clear();
    Pk.clear();
    Pr.clear();

    int     n_ex_lvl=0;
    int     left_occupancy=0;
    int N = primRef.size();
    bool debug=true;
    bool found;


	//sort reference vector, save permutation P_r

    Pr.resize(N,0);

    for (int i = 0 ; i != Pr.size() ; i++) 
	    Pr[i] = i;

   sort(Pr.begin(),Pr.end(),[&](const int& a, const int& b)
  			                   {return (primRef[a] < primRef[b]);});


        //sort K vector, save permutation P_k

    Pk.resize(N,0);
	for (int i = 0 ; i != Pk.size() ; i++) 
	    Pk[i] = i;

    	sort(Pk.begin(),Pk.end(),[&](const int& a, const int& b)
  			                   {return (primK[a] < primK[b]);});

	//DEBUGGING/TESTING
	if(debug)
	{
		for (int i = 0 ; i != Pr.size() ; i++) 
		{
		   std::cout << Pr[i] << std::endl;
		}

		std::cout << "input ref: " ;
		for(int i=0; i!=primRef.size(); i++)
		   std::cout << primRef[i] << " ";
		std::cout << std::endl;
		std::cout << "sorted ref: ";
		for(int i=0; i!=Pr.size(); i++)
			std::cout << primRef[Pr[i]] << " ";

		std::cout << std::endl;

      	        for (int i = 0 ; i != Pk.size() ; i++) 
		{
		   std::cout << Pk[i] << std::endl;
		}

		std::cout << "input K: " ;
		for(int i=0; i!=primK.size(); i++)
		   std::cout << primK[i] << " ";
		std::cout << std::endl;
		std::cout << "sorted K: ";
		for(int i=0; i!=Pk.size(); i++)
			std::cout << primK[Pk[i]] << " ";

		std::cout << std::endl;

	}





    // %%%%%%%% find first orb in reference not in K %%%%%%%%
    for(int occ_orb_idx=0; occ_orb_idx<N; occ_orb_idx++) //loop through occ ref orbs
    {
        //pull occupied orbital label
        auto occ_orb=primRef[Pr[occ_orb_idx]];

        //flag to say if it was found
        found=false;

        for(int K_orb_idx=0; K_orb_idx<N; K_orb_idx++) //see if its in the K prim
        {
            if(occ_orb==primK[Pk[K_orb_idx]]) // found orb from ref in prim K
            {
                found=true;
                break;
            }
        }

        if(!found) // found orbital in ref not in K
        {
            //increase excitation level
            n_ex_lvl++;

            //L_occ1=left_occupancy=orb_index
            left_occupancy=occ_orb_idx;
            
            //% keep looking for second orb in ref not in K
            //%   L_occ2=left_occupancy=orb_index - 1 // extra sub for occ1
            if(n_ex_lvl==2)
                left_occupancy=occ_orb_idx-1;



            //% keep looking for third orb in ref not in K
            //%   if you find it, stop and return as H_JK=0
            if(n_ex_lvl==3)
                return 3;
            
            //save orb 
            occ.push_back(occ_orb);

        }
    }	

       //loop through occ ref orbs O(N) cost
         //check against orbs of K O(N)
	 //if not found
	   //increase excitation level
	   //include permutation of occ orbital to front in Pr (same info as orb_index marking it)
	   
    // if all the same return P_r, P_k
    //no need to go through K if all are orbs the same
    if(n_ex_lvl==0)
        return 0;


    // %%%%%%%% find orb in K not in reference %%%%%%%%
    for(int K_orb_idx=0; K_orb_idx<N; K_orb_idx++)
    {
        //pull orbital number
        auto K_orb=primK[Pk[K_orb_idx]];

        found=false;
        for(int occ_orb_idx=0; occ_orb_idx<N; occ_orb_idx++)
        {
            if(K_orb==primRef[Pr[occ_orb_idx]])
            {
                found=true;
                break;
            }
        }

        if(!found)
        {

            //save data
            virt.push_back(K_orb);
            
            //taking into account previously found
            left_occupancy=K_orb_idx - (virt.size()-1);

            //accumulate phase
            if(left_occupancy%2)
                phase*=-1;

            //stop looking if we only need one
            if(n_ex_lvl==1)
                return n_ex_lvl;
            
            //stop looking if we found two
            if(virt.size()==2)
                return n_ex_lvl;

        }

    }	
    
 
    printf("\n\nSomething went wrong in PrimDiffs, logic broken\n");
    throw;
    return -9;


       //loop through occ K orbs O(N) cost
         //check against orbs of ref, O(N)
	 //if not found
	    //include permutation of virt orb to front in Pk
            //do we need to find another? if yes, continue else return: Pr, Pk, virt_indicies, occ_indicies
}
