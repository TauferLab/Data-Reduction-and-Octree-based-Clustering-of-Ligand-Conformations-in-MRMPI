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

// MapReduce octree clustering, implementation 2: locally count, move octant id and local counts, count globaly in reduce
// MPI Syntax:mpirun -np 4 octree_gl digits threshold save_file|none file1 dir1 file2 dir2 ... (-np 4 specifies 4 processes to run concurrently, change 4 to the number of processes you would like to run; file1 dir1 ... specifies the input files and dirs to the program, if dir1 is specified, all the files within that dir is going to be considered as input files to the program)
// map: (1) compute octkey for each ligand; (2) based on the level to explore x, count partially the points in octants; (3) output octant id, partial count for each octant in level x
// collate(shuffle) : communicate octant id and partial count among processes
// reduce: sum all the partial counts associated with the same octant id, only output the ones that density >= threshold


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
void explore_level(int, int, MapReduce * ); //explore the int level of the tree
void gen_leveled_octkey(uint64_t, char *, int, char *, int, KeyValue *, void *); //the 2nd map, takes partially the entire octkey, and that is the key for the following reduce, e.g., keys[1...7]
void sum(char *, int, char *, int, int *, KeyValue *, void *); //sum the keyvalue object locally to the process, to make sure the key is unique in the kv object

void save_to_file(bool, MapReduce *, char *,int, int); //save the kv object of a mr object to file or not based on the result signal
void compress_leveled_octkey(char *, int, char *, int, int *, KeyValue *, void *); //compress the keyvalue object locally to the process, to make sure the key is unique in the kv object


/* --------------------------------global vars-------------------------------------- */
int me,nprocs; //for mpi process, differenciate diff processes when print debug info
int digits=15; //width of the octkey, 15 default, set by main, used by map
int thresh=1; //number of points within the octant to be considered dense, and to be further subdivide, set by main, used by reduce
bool result=true; //save results to file, true= yes, set by main, used by map and reduce
bool branch=false;//true branch down, false branch up; used by main and sum functions
//char base_dir[] = "/home/bzhang/csd173/bzhang/mrmpi-20Jun11/mrmpi_clustering/scripts/";
//char base_dir[] = "/home/zhang/mrmpi-20Jun11/mrmpi_clustering/";
const char *base_dir;
int level; //level: explore this level of the oct-tree, used by main and gen_leveled_octkey
bool realdata=true; //use real dataset or synthetic dataset, true means using real, false means using synthetic, set by main, used by map generate octkey

/*-------------------------------------------------------------------------*/

int main(int narg, char **args)
{
	/*MPI initialization*/
	MPI_Init(&narg,&args);
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	/*process input argument*/
	if (narg <= 6){
		if (me == 0){
			printf("Syntax: octree real|syn digits threshold save_file|none file1 file2 dir1 ...\n");
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
	int min_limit, max_limit; //level: explore this level of the oct-tree
	min_limit=0;
	max_limit=digits+1;
	level=floor((max_limit+min_limit)/2);
	
	
	/*map: (1) compute octey for each ligand; (2) based on the level to explore x,  partial count of each ocatant (3) output octant id, and partial counts*/
	MPI_Barrier(MPI_COMM_WORLD);

	
	MapReduce *mr_convert = new MapReduce(MPI_COMM_WORLD);
	mr_convert->memsize=1024;
	mr_convert->verbosity = 0;
	mr_convert->timer = 0;
	mr_convert->set_fpath(base_dir);
	
	/*compute octkey for each ligand, key: 15 digits char string; value: if,rmsd,energy char string*/
	int nwords = mr_convert->map(narg-6,&args[6],1,1,0,generate_octkey,NULL);
	
	/*save key,value to file or not*/
	//save_to_file(result, mr_convert, "convert", 5, 5);
	
	MPI_Barrier(MPI_COMM_WORLD); //barreir after map compute octkey
	
	
	while ((min_limit+1) != max_limit){
		
		branch=false;// set branch to false, if exsit one dense octant, then branch=true, the next level is down
		
		explore_level(me, level, mr_convert); //map: partial count level x, reduce, global count level x
		
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


//	printf("in the function: genearte_octkey.\n");

	while (std::getline(ifile, line)){
		std::vector<std::string> coor;
	
		//std::cout << "data line:" << line <<std::endl;
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
			//printf("from process: %d, the real octkey is: %s.\n",itask,real_key.c_str());
			
		}


		const int num_atoms = coor.size()/3;
		double x[num_atoms], y[num_atoms], z[num_atoms];
		/*x,y,z double arries store the cooridnates of ligands */
		for (int i=0;i!=num_atoms; ++i){
			x[i] = atof(coor.at(3*i).c_str());
			y[i] = atof(coor.at(3*i+1).c_str());
			z[i] = atof(coor.at(3*i+2).c_str());

		}
//		for (int i=0;i<num_atoms;++i){
//			printf("Proc %d: the x,y,z coordinates are: %f, %f, %f.\n",me,x[i],y[i],z[i]);
//		}

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
			//keys[count] =  bit;
			//printf("Proc %d: the %d level octkey is: %c.\n",me, count,bitc);
			key[count] = bitc;
		
			++count;
		}	
		if (realdata == false){
		//use the real octkey from reading the line of synthetic dataset
			strcpy(key, real_key.c_str());
			//printf("from process: %d, the real octkey is: %s, the copied real octkey is: %s.\n", itask, real_key.c_str(), key);
		}

		/*put the keys, values to kv pair, key=keys, value=id,rmsd,energy*/
		//char *key = reinterpret_cast<char*>(keys);
		std::string values;	
		values+=id.append(",").append(rmsd).append(",").append(energy);
		
		char *value = new char [values.length() + 1];
		strcpy(value, values.c_str());

		kv->add(key, count+1, value, (values.length() + 1)); //sizeof(key) does not work properly
		//printf("from process %d, insert %s, %s to kv.\n", itask, key, value);
		delete[] value;
		//delete[] key;
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
//	printf("Proc %d: the xybar, xxbar, and slope are: %f, %f, %f.\n",me, xybar, xxbar, slope);
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


/*-------explore the x level of the tree, each map counts partially points----------------*/
void explore_level(int me, int level,MapReduce *mr_convert){

	/*count partially density of level x octants*/
	MapReduce *mr_level = new MapReduce(MPI_COMM_WORLD);
	mr_level->memsize = 1024;
	mr_level->verbosity = 0;
	mr_level->timer = 0;
	mr_level->set_fpath(base_dir);

	//in mr_level, each key is octant id of level x
	int nkeys = mr_level->map(mr_convert, gen_leveled_octkey,NULL,0);
	
	//partially count for each octant
	int unique_local = mr_level->compress(compress_leveled_octkey,NULL);

	//printf("Proc %d: unique keys locally to the process are %d.\n",me,unique_local);

	/*gather kv pairs across the processes, make kv to kmv object*/
	mr_level->collate(NULL);
	
	int unique = mr_level->reduce(sum,NULL);
	//printf("Proc %d: the total number of etries in the kmv is: %d.\n", me, unique);


	/*print to file or not*/
	char buffer[10];
	int tmp = sprintf(buffer, "level_%d",level);
	save_to_file(result, mr_level, buffer, 5, 1);

	/*if the reduce(sum) returns value >0, there is some octant across all processors that > thresh*/
	if (unique >0){
		branch=true;
//		printf("the key with max num of points is: %s, the max num of points is: %d.\n", max_key, max_points);  
	}

	delete mr_level;

	
}




/*-----read key, value pair from mr, take partially of the key, based on which level of the tree is searching, and emit key, value = null to the intenal KeyValue object*/


void gen_leveled_octkey(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue * kv, void *ptr){


	//printf("the key to gen_leveled_octkey is: %s.\n",key);
	char *newkey = new char[level+1];
	memcpy(newkey, key, level); //newkey is the char[] stores the partial key that will be used to count in reduce
	newkey[level]='\0';
	kv->add(newkey, level+1, NULL,0);
	//printf("the newkey to gen_leveld_octkey is :%s.\n",newkey);
	delete[] newkey;

}
/*---compress the values with the same key locally per process to make the kv object has unique keys--------------*/

void compress_leveled_octkey(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr){

	kv->add(key,keybytes,(char *) &nvalues,sizeof(int));
	//printf("the added key is: %s, the added values is: %d.\n",key, nvalues);
}


void sum(char *key, int keybytes, char *multivalue,
	 int nvalues, int *valuebytes, KeyValue *kv, void *ptr) 
{
	int num_points=0;

        if ((multivalue!=NULL) && (nvalues!=0)){ //the kmv fits in one page memory
              int *multivalueptr = (int *) multivalue;
              for (int i=0; i!=nvalues; ++i){
                      num_points += *multivalueptr;
                      multivalueptr++;
              }
                if (num_points >= thresh){
                        kv->add(key,keybytes,(char *) &num_points,sizeof(int));
                //      printf("Proc %d: from sum function, the added key is: %s, the added val is: %d.\n", me, key, num_points);
                }
        }else{ //the kmv does not fit in one page memory

                MapReduce *mr_tmp = (MapReduce *) valuebytes;
                int nblocks;
                uint64_t nvalues_total = mr_tmp->multivalue_blocks(nblocks);

              for (int iblock = 0; iblock < nblocks; iblock++){
                      int nv = mr_tmp->multivalue_block(iblock,&multivalue,&valuebytes);
                       int *multivalueptr = (int *) multivalue;
                      for (int i = 0; i < nv; i++){
                               num_points += *multivalueptr;
                               multivalueptr++;
                      }
              }
                if (num_points >= thresh){
                        kv->add(key,keybytes,(char *) &num_points,sizeof(int));
                //      printf("Proc %d: from sum function, the added key is: %s, the added val is: %d.\n", me, key, num_points);
                }
        }


}

