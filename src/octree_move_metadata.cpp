/*
Copyright (c)
	2015 by The University of Delaware
	Contributors: Boyu Zhang, Michela Taufer
	Affiliation: Global Computing Laboratory, Michela Taufer PI
	Url: http://gcl.cis.udel.edu/, https://github.com/TauferLab

All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

	1. Redistributions of source code must retain the above copyright notice, 
	this list of conditions and the following disclaimer.

	2. Redistributions in binary form must reproduce the above copyright notice,
	this list of conditions and the following disclaimer in the documentation
	and/or other materials provided with the distribution.

	3. If this code is used to create a published work, one of the following 
	papers must be cited.

		Boyu Zhang, Trilce Estrada, Pietro Cicotti, and Michela Taufer. On 
		Efficiently Capturing Scientific Properties in Distributed Big Data 
		without Moving the Data - A Case Study in Distributed Structural Biology 
		using MapReduce. In the Proceedings of the 16th IEEE International 
		Conferences on Computational Science and Engineering (CSE), December 
		2013, Sydney, Australia.

	4.  Permission of the PI must be obtained before this software is used
	for commercial purposes.  (Contact: taufer@acm.org)

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// MapReduce octree clustering, implementation 1: global move data once, count locally in reduce
// MPI Syntax:mpiexec -n 4 octree real|syn digits threshold save_file|none path_to_save_files file1 file2 dir1 ...
// map: (1) compute octkey for each ligand; (2) output key: octkey; value:id for each ligand
// collate(shuffle) : based on the number of reduce, distribute the key, values among processes
// reduce: based on the level to explroe, explore the level, decide to go down or up in the tree

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/stat.h"
#include "mapreduce.h"
#include "keyvalue.h"
#include "stdint.h"
#include "math.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
//#include "convert.h"


using namespace MAPREDUCE_NS;

/*---------------functions definition----------------------------------------------------------*/
void generate_octkey(int, char *, KeyValue *, void *);// the map, generate KeyValue object, key=octkey(int[15], value=id,rmsd,energy(char *))

double slope(double[], double[], int); //function inside generate_octkey
int key_partition(char *, int); //partition key,value pairs based on the first couple of digits of the octkey among the r reduce processes
void explore_level(int, int, MapReduce * ); //explore the int level of the tree
void sum(char *, int, char *, int, int *, KeyValue *, void *); //sum the keyvalue object locally to the process, to make sure the key is unique in the kv object
void save_to_file(bool, MapReduce *, char *, int, int); //save the kv object of a mr object to file or not based on the result signal
void gen_leveled_octkey(uint64_t, char *, int, char *, int, KeyValue *, void *); //the 2nd map, takes partially the entire octkey, and that is the key for the following reduce, e.g., keys[1...7]

/*for profiling and tracing*/
void data_shuffle(MapReduce *);


/* --------------------------------global vars-------------------------------------- */
int me,nprocs; //for mpi process, differenciate diff processes when print debug info
int digits=15; //width of the octkey, 15 default, set by main, used by map
int thresh=1; //number of points within the octant to be considered dense, and to be further subdivide, set by main, used by reduce
bool result=true; //save results to file, true= yes, set by main, used by map and reduce
bool branch=false;//true branch down, false branch up; used by main and sum functions
//char base_dir[] = "/home/bzhang/csd173/bzhang/mrmpi-20Jun11/mrmpi_clustering/scripts/";
//char base_dir[] = "/home/zhang/mrmpi-20Jun11/mrmpi_clustering/";
const char *base_dir; //save the out of core files, and the place to print kv, kmv info to files
int level; //level: explore this level of the oct-tree, used by main and gen_leveled_octkey
bool realdata=true; //use real dataset or synthetic dataset, true means using real, false means using synthetic, set by main, used by map generate octkey


/*---------------------------main execution----------------------------------------------*/

int main(int narg, char **args)
{
	/*MPI initialization*/
	MPI_Init(&narg,&args);
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	/*process input argument
	args[1]: real|syn - use real dataset of synthetic dataset
	args[2]: 15 - size of the octkey (the number of digits that the octkey contains)
	args[3]: 100 - the density threshold of the octant (an octant contains more
			  than 100 points is considered a dense octant)
	args[4]: save_file|none - whether to save the output of the results to file
	args[5]: ./ - the path to save spill files and results files
	args[6] and after: input files and directoreis*/

	if (narg <= 6){ //include the octree, args[0]
		if (me == 0){
			printf("Syntax: octree real|syn digits threshold save_file|none path_to_save_files file1 file2 dir1 ...\n");
		}
		MPI_Abort(MPI_COMM_WORLD,1);
	}else{
		const char *dataset = args[1];
		if (strcmp(dataset, "real") ==0){
			// using real dataset
			realdata=true;
		}else{
			realdata=false;
		}
		digits = atoi(args[2]);
		thresh = atoi(args[3]);
		const char *save_file = args[4];
		if (strcmp(save_file, "save_file") ==0){
			// save ligand results
			result=true;
		}else{
			result=false;
		}
		base_dir=args[5];
		if (me==0){
			printf("Proc %d: from command line: dataset is: %s, digits is: %d, thresh is: %d, base dir is: %s.\n",  me,  dataset, digits, thresh, base_dir);
		}
	}
	
	/*var initilization*/
	int min_limit, max_limit;
	min_limit=0;
	max_limit=digits+1;
	level=floor(max_limit+min_limit)/2;
	

	/*map: (1) compute octkey for each ligand; 
		   (2) key: octkey; value:id for each ligand*/	
	
	MapReduce *mr_convert = new MapReduce(MPI_COMM_WORLD);
	//mr_convert->memsize=1024;
	mr_convert->verbosity = 0;
	mr_convert->timer = 0;
	mr_convert->set_fpath(base_dir);
	
	/*compute octkey for each ligand, key: 15 digits char string; value: if,rmsd,energy char string*/

	int nwords = mr_convert->map(narg-6,&args[6],1,1,0,generate_octkey,NULL); //local

	data_shuffle(mr_convert);

	/*at this point, the key,value pairs are partitioned among the r reduce processes based on the first couple of digits of the key. 
	next the reduce perform local reduction, and decides which level to go next.
	no need to copy mr object? store kv and kmv */
	while ((min_limit+1) != max_limit){
		
		branch=false;// set branch to false, if exsit one dense octant, then branch=true, the next level is down
		
		explore_level(me, level, mr_convert);
		
		MPI_Barrier(MPI_COMM_WORLD); //enforce a barrier after each explore of each level, so branch is set correctly
		
		 /*compute next level to explore based on branch up or down*/
		if (branch) { //branch is true, branch down
			min_limit=level;
			level =  floor((max_limit+min_limit)/2);
		}else{//branch is flase, bracn up
			max_limit=level;
			level =  floor((max_limit+min_limit)/2);
		}
	}
		

	delete mr_convert;
	/*MPI finalization*/
	MPI_Finalize();
}

/*---------------function implementation-----------------------------*/

void data_shuffle (MapReduce *mr_convert){
        /*aggregate (kv-kv): the global communication among reduce processes and convert(kv-kmv)*/
        int nwords = mr_convert->aggregate(key_partition);

}

/*the map function
itask will have a value 0 <= itask < nfiles, where nfiles is the number
of filenames in the list of files that was generated. 
Your function is also passed a single filename, which it will presumably open and read. 
http://mapreduce.sandia.gov/doc/map.html
*/

void generate_octkey(int itask, char *fname, KeyValue *kv, void *ptr)
{
	std::ifstream ifile;
	ifile.open(fname);
	std::string line;
	std::string word;
	std::istringstream iss;
	std::string rmsd="0";
	std::string energy="0";
	double tmp=0, range_up=10.0, range_down=-10.0; //the upper and down limit of the range of slopes
	std::string real_key;
	char key[digits+1];
	std::string id;
	std::istringstream stm;



	while (std::getline(ifile, line)){
		std::vector<std::string> coor;
	
		iss.str(line);
		iss >> id;

		while (iss >> word){
			coor.push_back(word);
		}
		if (realdata==true){
			rmsd = coor.back();
			coor.pop_back();
			energy = coor.back();
			coor.pop_back();
		}else{
			//last one is real octkey
			real_key = coor.back();
			coor.pop_back();
			rmsd = coor.back();
			coor.pop_back();
			energy = coor.back();
			coor.pop_back();
			
		}


		const int num_atoms = coor.size()/3;
		double x[num_atoms], y[num_atoms], z[num_atoms];
		/*x,y,z double arries store the cooridnates of ligands */
		for (int i=0;i!=num_atoms; ++i){
			x[i] = atof(coor.at(3*i).c_str());
			y[i] = atof(coor.at(3*i+1).c_str());
			z[i] = atof(coor.at(3*i+2).c_str());

		}


		/*compute the b0, b1, b2 using linear regression*/
		double b0 = slope(x, y, num_atoms);
		double b1 = slope(y, z, num_atoms);
		double b2 = slope(x, z, num_atoms);


		/*compute octkey, "digit" many digits*/
		int count=0;//count how many digits are in the octkey
		double minx = range_down, miny = range_down, minz = range_down;
		double maxx = range_up, maxy = range_up, maxz = range_up;
		while (count < digits){
			int m0 = 0, m1 = 0, m2 = 0;
			double medx = minx + ((maxx - minx)/2);
			if (b0>medx){
				m0=1;
				minx=medx;
			}else{
			
				maxx=medx;
			}
			double medy = miny + ((maxy-miny)/2);
			if (b1>medy){
				m1=1;
				miny=medy;
			}else{
			
				maxy=medy;
			}
			double medz = minz + ((maxz-minz)/2);
			if (b2>medz){
				m2=1;
				minz=medz;
			}else{
			
				maxz=medz;
			}
			/*calculate the octant using the formula m0*2^0+m1*2^1+m2*2^2*/
			int bit=m0+(m1*2)+(m2*4);
			char bitc=(char)(((int)'0') + bit); //int 8 => char '8'
			key[count] = bitc;
		
			++count;
		}	
		if (realdata == false){
		//use the real octkey from reading the line of synthetic dataset
			strcpy(key, real_key.c_str());
		}

		/*put the keys, values to kv pair, key=keys, value=id,rmsd,energy*/
		std::string values;	
		values+=id.append(",").append(rmsd).append(",").append(energy);
		
		char *value = new char [values.length() + 1];
		strcpy(value, values.c_str());

		kv->add(key, count+1, value, (values.length() + 1)); //sizeof(key) does not work properly
		delete[] value;
		coor.clear();	
		iss.clear();
	}

}

double slope(double x[], double y[], int num_atoms){
	double slope=0.0;
	double sumx=0.0, sumy=0.0;
	for (int i=0; i!=num_atoms; ++i){
		sumx += x[i];
		sumy += y[i];
	}
	double xbar = sumx/num_atoms;
	double ybar = sumy/num_atoms;

	double xxbar =0.0, yybar =0.0, xybar =0.0;
	for (int i=0; i!=num_atoms; ++i){
		xxbar += (x[i] - xbar) * (x[i] - xbar);
		yybar += (y[i] - ybar) * (y[i] - ybar);
		xybar += (x[i] - xbar) * (y[i] - ybar);

	}
	slope = xybar / xxbar;
	return slope;
	
}


/*save key,value to file or not*/
void save_to_file(bool result, MapReduce *mr_convert, char *name, int kflag, int vflag){
	if (result){
		if (me==0)
			printf("Proc %d: need to print ligands selected to files.\n",me);
		char file_name[100];
		strcpy(file_name,base_dir);
		strcat(file_name,name);
		if (me==0)
			printf("Proc %d: the convert file is: %s.\n",me,file_name);
		mr_convert->print(file_name,1,-1,1,kflag,vflag);
	}else{
		if (me==0)
			printf("Proc %d: no need to print ligands selected to files.\n",me);
	}
}



/*partition the key value pairs among r reduce processes, 
based on the first multiple octkey digits
for example, if we have 64 reduce processes, we need to 
partition the octkeys based on the first 2 digits.
Octkeys with the same first 2 digits are communicated to the same reduce process,
then when the reduce process branch up and down through the octree,
it does not need to communciate the octkeys again. Assume that it will never
need to branch above the 2nd level of the octree.*/

int key_partition(char *key, int keybytes){
	/*compute how many digits of octkey to use to partition*/
	const int useful_bits = ceil (log(nprocs)/log(8.));
	int power_base=8;
	/*partition based on the useful_bits many bits of the key
	if using the first 3 digits, return (1st dig) + (2nd dig)x8 + (3 dit) x64*/
	int value=0;
	int key_int=0;
	double tmp=0;
	for (int i=0; i< useful_bits; i++){
		key_int=key[i] - '0'; //numberical value of key[i]
		
		tmp=pow(double(power_base),i);
		value+=(key_int * tmp);
		
	}
	return value;
}



/*based on the level to explore, count the points, lcoally*/
void explore_level(int me, int level, MapReduce *mr_convert){
	MapReduce *mr_level = new MapReduce(MPI_COMM_WORLD);
	//mr_level->memsize=1024;
	mr_level->verbosity = 0;
	mr_level->timer = 0;
	mr_level->set_fpath(base_dir);	
	int nkeys = mr_level->map(mr_convert, gen_leveled_octkey,NULL,0);

	int uniquekeys = mr_level->convert(); //duplicated kv pairs becomes kmv pair

	
	int unique = mr_level->reduce(sum,NULL);
	
	/*print to file or not*/
	char buffer[10];
	int tmp = sprintf(buffer, "level_%d",level);
	save_to_file(result, mr_level, buffer, 5, 1);

	/*if the reduce(sum) returns value >0, there is some octant across all processors that > thresh*/
	if (unique >0){
		branch=true;
	}
	delete mr_level;
}


/*-----read key, value pair from mr, take partially of the key, based on which level of the tree is searching, and emit key, value = null to the intenal KeyValue object*/


void gen_leveled_octkey(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue * kv, void *ptr){


	char *newkey = new char[level+1];
	char newvalue='1';
	memcpy(newkey, key, level); //newkey is the char[] stores the partial key that will be used to count in reduce
	newkey[level]='\0';
	kv->add(newkey, level+1, &newvalue,1);
	delete[] newkey;

}



void sum(char *key, int keybytes, char *multivalue,
	 int nvalues, int *valuebytes, KeyValue *kv, void *ptr) 
{
	uint64_t num_points=0;

	if ((multivalue!=NULL) && (nvalues!=0)){ //the kmv fits in one page memory

		num_points=nvalues;
		if (num_points >= thresh){
			kv->add(key,keybytes,(char *) &num_points,sizeof(int));
		}
	}else{ //the kmv does not fit in one page memory

		MapReduce *mr_tmp = (MapReduce *) valuebytes;
		int nblocks;
		uint64_t nvalues_total = mr_tmp->multivalue_blocks(nblocks);

		num_points=nvalues_total;
		if (num_points >= thresh){
            kv->add(key,keybytes,(char *) &num_points,sizeof(int));
        }
	}
}


