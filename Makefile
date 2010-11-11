CXXFLAGS = -Wall -O3
#CXXFLAGS = -Wall -g

OBJS = ikrobot.o robot.o tree.o

LIBS = -lGL -lglut -lGLU -larmadillo
LD_FLAGS = -L/usr/sww/lib # find glut for inst machines

default: compile

compile: ikrobot

ikrobot: $(OBJS)
	$(CXX) $(LD_FLAGS) $(LIBS) -o $@ $(OBJS)

clean:
	rm -f *.o ikrobot

ikrobot.o: ikrobot.cpp
	$(CXX) -c $(CXXFLAGS) -o $@ $<

robot.o: robot/robot.cpp robot/robot.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<

tree.o: robot/tree.cpp robot/tree.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<
