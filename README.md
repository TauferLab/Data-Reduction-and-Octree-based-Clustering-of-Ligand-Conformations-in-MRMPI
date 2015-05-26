# Octree-based Clustering in MapReduce-MPI 

This repository contains two variations of a general analysis algorithm for 
large datasets residing in distributed memory systems. Both variations avoid 
the need to move data among nodes by extracting relevant data properties 
locally and concurrently and transforming the analysis problem (e.g., clustering 
or classification) into a search for property aggregates. The programs are 
adapted for MapReduce and implemented in MapReduce-MPI.

For detailed design and implementation aspects, please refer to our following papers:

Boyu Zhang, Trilce Estrada, Pietro Cicotti, and Michela Taufer. On Efficiently 
Capturing Scientific Properties in Distributed Big Data without Moving the 
Data - A Case Study in Distributed Structural Biology using MapReduce. In the 
Proceedings of the 16th IEEE International Conferences on Computational Science 
and Engineering (CSE), December 2013, Sydney, Australia.



## Directories and files 

### data
Contains two text files serve as sample input to the Hadoop program.

* ligand_small.txt 

Contains 10 ligand conformations, one conformation per line. Each line is in the 
format of [ligand_id x1 y1 z1 x2 y2 z2 … xN yN ZN energy rmsd], in which 
ligand_id is the ligand conformation ID, x1 y1 z1 x2 y2 z2 … xN yN ZN are the 
coordinates in the Cartesian space of the N ligand atoms, energy is the 
conformation energy, and rmsd is the RMSD from the crystal structure in LPDB. 
From a representation point of view, ligand_id is represented as the string 
type and the others values in each line are represented as the double type.

* ligand_1a9u.txt 

Contains 76889 ligand conformations with the same format.

### src
Contains source files for the two programs. Both variations are based on the 
general theme of moving computation to the data. 

* octree_move_metadata.cpp
Contains source code for the first variation of the clustering. In the first 
variation, called GlobalToLocal, the extracted properties are moved across nodes 
to build a global view of the dataset; these properties are iteratively and 
locally analyzed while searching for some class or cluster convergence. 

* octree_move_density.cpp
Contains source code for the second variation of the clustering. In the second 
variation, called LocalToGlobal, partial properties are analyzed locally on each 
node to generate a property aggregate, which is a scalar number that summarizes 
local properties of the data. In this variation, instead of exchanging 
extracted properties (i.e., metadata), LocalToGlobal only exchanges a few 
scalar property aggregates (i.e., local densities) across the nodes. 

* Makefile
The Makefile to compile both programs.

## Compilation instructions 

### Supported platform
GNU/Linux

### Required software
* C++ compiler must be installed, gcc/g++ or icc/i++.
* MPI library must be installed, MPICH, MVAPICH, etc.
* MapReduce-MPI library must be installed. This library is available here:
http://mapreduce.sandia.gov

### Compilation instructions using command line
The following instructions assume that you have the MPI library and the MapReduce-MPI
library correctly installed and configured.

* Edit Makefile to include your own MapReduce-MPI library 

```
CFLAGS = -W -O -I/home/bzhang/mrmpi-7Apr14/src
USERLIBS = /home/bzhang/mrmpi-7Apr14/src/libmrmpi_mpicc.a
```

* Run Makefile to generate executables

```
[bzhang@geronimo src]$ make all
mpicxx -W -O -I/home/bzhang/mrmpi-7Apr14/src -c octree_move_density.cpp
octree_move_density.cpp: In function ‘int main(int, char**)’:
octree_move_density.cpp:93: warning: converting to ‘int’ from ‘double’
octree_move_density.cpp:125: warning: converting to ‘int’ from ‘double’
octree_move_density.cpp:128: warning: converting to ‘int’ from ‘double’
mpicxx -g -O  octree_move_density.o /home/bzhang/mrmpi-7Apr14/src/libmrmpi_mpicc.a -o octree_move_density
mpicxx -W -O -I/home/bzhang/mrmpi-7Apr14/src -c octree_move_metadata.cpp
octree_move_metadata.cpp: In function ‘int main(int, char**)’:
octree_move_metadata.cpp:103: warning: converting to ‘int’ from ‘double’
octree_move_metadata.cpp:135: warning: converting to ‘int’ from ‘double’
octree_move_metadata.cpp:138: warning: converting to ‘int’ from ‘double’
octree_move_metadata.cpp: In function ‘int key_partition(char*, int)’:
octree_move_metadata.cpp:333: warning: converting to ‘int’ from ‘double’
octree_move_metadata.cpp:344: warning: converting to ‘int’ from ‘double’
octree_move_metadata.cpp: In function ‘void sum(char*, int, char*, int, int*, MAPREDUCE_NS::KeyValue*, void*)’:
octree_move_metadata.cpp:404: warning: comparison between signed and unsigned integer expressions
octree_move_metadata.cpp:414: warning: comparison between signed and unsigned integer expressions
mpicxx -g -O  octree_move_metadata.o /home/bzhang/mrmpi-7Apr14/src/libmrmpi_mpicc.a -o octree_move_metadata
```

### Execution instructions 

* Run with data ligand_1a9u.txt:
``` 
[bzhang@geronimo src]$ mpiexec -n 2 ./octree_move_metadata real 15 500 save_file ./ ../data/ligand_1a9u.txt 
Proc 0: from command line: dataset is: real, digits is: 15, thresh is: 500, base dir is: ./.
Proc 0: need to print ligands selected to files.
Proc 0: the convert file is: ./level_8.
Proc 0: need to print ligands selected to files.
Proc 0: the convert file is: ./level_12.
Proc 0: need to print ligands selected to files.
Proc 0: the convert file is: ./level_10.
Proc 0: need to print ligands selected to files.
Proc 0: the convert file is: ./level_9.
[bzhang@geronimo src]$ ll
total 784
-rw-r--r-- 1 bzhang users      0 May 26 09:57 level_10.0
-rw-r--r-- 1 bzhang users      0 May 26 09:57 level_10.1
-rw-r--r-- 1 bzhang users      0 May 26 09:57 level_12.0
-rw-r--r-- 1 bzhang users      0 May 26 09:57 level_12.1
-rw-r--r-- 1 bzhang users   1994 May 26 09:57 level_8.0
-rw-r--r-- 1 bzhang users    365 May 26 09:57 level_8.1
-rw-r--r-- 1 bzhang users   1633 May 26 09:57 level_9.0
-rw-r--r-- 1 bzhang users      0 May 26 09:57 level_9.1
-rw-r--r-- 1 bzhang users    644 May 25 21:02 Makefile
-rwxr-xr-x 1 bzhang users 202158 May 25 21:02 octree_move_density
-rw-r--r-- 1 bzhang users  14309 May 25 15:59 octree_move_density.cpp
-rw-r--r-- 1 bzhang users 161496 May 25 21:02 octree_move_density.o
-rwxr-xr-x 1 bzhang users 202620 May 25 21:02 octree_move_metadata
-rw-r--r-- 1 bzhang users  13247 May 25 15:59 octree_move_metadata.cpp
-rw-r--r-- 1 bzhang users 161736 May 25 21:02 octree_move_metadata.o
```

* Command line arguments:
```
[real|syn]: real ligand conformation datasets or synthetic datasets with specific distributions.
15: length of the octkey.
500: threshold of the density.
[save_file|none]: whether to save the octant id and density to files.
./: the directory to save the octant id and density files.
../data/ligand_1a9u.txt: input file.
```

### Examine output
```
[bzhang@geronimo src]$ more level_9.0 
KV pair: proc 0, sizes 10 4, key 433323467, value 564
KV pair: proc 0, sizes 10 4, key 433323657, value 842
KV pair: proc 0, sizes 10 4, key 433323643, value 1352
KV pair: proc 0, sizes 10 4, key 433323654, value 1880
KV pair: proc 0, sizes 10 4, key 433323655, value 958
KV pair: proc 0, sizes 10 4, key 433323645, value 1754
KV pair: proc 0, sizes 10 4, key 433323646, value 1350
KV pair: proc 0, sizes 10 4, key 433323656, value 2136
KV pair: proc 0, sizes 10 4, key 433323652, value 1638
KV pair: proc 0, sizes 10 4, key 433323647, value 2226
KV pair: proc 0, sizes 10 4, key 433323640, value 806
KV pair: proc 0, sizes 10 4, key 433327201, value 914
KV pair: proc 0, sizes 10 4, key 433332625, value 708
KV pair: proc 0, sizes 10 4, key 433323642, value 736
KV pair: proc 0, sizes 10 4, key 433327210, value 910
KV pair: proc 0, sizes 10 4, key 433332630, value 808
KV pair: proc 0, sizes 10 4, key 433323607, value 920
KV pair: proc 0, sizes 10 4, key 433323650, value 2004
KV pair: proc 0, sizes 10 4, key 433332616, value 850
KV pair: proc 0, sizes 10 4, key 433323644, value 840
KV pair: proc 0, sizes 10 4, key 433332661, value 548
KV pair: proc 0, sizes 10 4, key 433327200, value 574
KV pair: proc 0, sizes 10 4, key 433323653, value 1038
KV pair: proc 0, sizes 10 4, key 433323651, value 1030
KV pair: proc 0, sizes 10 4, key 433323616, value 1280
KV pair: proc 0, sizes 10 4, key 433323641, value 1762
KV pair: proc 0, sizes 10 4, key 433323617, value 744
KV pair: proc 0, sizes 10 4, key 433332634, value 1020
KV pair: proc 0, sizes 10 4, key 433323605, value 662
KV pair: proc 0, sizes 10 4, key 433332670, value 600
```
Note: Key is the octant id, value is the density of this octant.



