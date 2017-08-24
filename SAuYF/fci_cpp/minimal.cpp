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
main()
{
	using namespace std;


	vector<int> occ;
	vector<int> virt;
	vector<int> left  {1,5,3,7,2,2,6};
	vector<int> right {1,9,10,7,7,2,6};

	PrimDiffs(left,right,occ,virt);

	return 0;
}
