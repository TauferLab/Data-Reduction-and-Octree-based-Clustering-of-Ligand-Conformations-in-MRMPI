CC = mpicc
CPP = mpicxx
CLINK = mpicc
CPPLINK = mpicxx
DEBUG = -g
CFLAGS = -W -O -I/home/bzhang/mrmpi-7Apr14/src
LFLAGS = -g -O 
USERLIBS = /home/bzhang/mrmpi-7Apr14/src/libmrmpi_mpicc.a


all:    octree_move_density octree_move_metadata

octree_move_density: octree_move_density.o $(USERLIBS)
	$(CPPLINK) $(LFLAGS) octree_move_density.o $(USERLIBS) -o octree_move_density
octree_move_metadata: octree_move_metadata.o $(USERLIBS)
	$(CPPLINK) $(LFLAGS) octree_move_metadata.o $(USERLIBS) -o octree_move_metadata

clean:
	rm *.o  octree_move_density octree_move_metadata

# rules
%.o:%.cpp
	$(CPP) $(CFLAGS) -c $<

%.o:%.c
	$(CC) $(CFLAGS) -c $<
