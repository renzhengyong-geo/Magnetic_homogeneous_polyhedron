
CXX              :=  g++
CXXFLAGS         :=  -O3 -fopenmp 
M3DINCLUDE       := -I../src

INCS:= $(wildcard ../src/*.h) 
SRCS:= $(wildcard ../src/*.cpp) 
OBJS:= $(patsubst %.cpp, %.o, $(SRCS))

	
# How to compile C++
%.o : %.cpp 
	$(CXX) $(CXXFLAGS) $(M3DINCLUDE) -c $< -o $@


#--------------------------------------------------------------
m3d: $(OBJS) $(INCS) test_m3d.o
	$(CXX) $(CXXFLAGS) -o m3d test_m3d.o $(OBJS) 
clean:	
	@rm -rf *.o *~ $(OBJS)  *.vtk m3d
		
