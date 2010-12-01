CXXFLAGS = -Wall -O3
#CXXFLAGS = -Wall -g

OBJS = context.o ikrobot.o link.o robot.o tree.o

LIBS = -lGL -lglut -lGLU -larmadillo
LD_FLAGS = -L/usr/sww/lib # find glut for inst machines

default: compile

compile: ikrobot

ikrobot: $(OBJS)
	$(CXX) $(LD_FLAGS) $(LIBS) -o $@ $(OBJS)

clean:
	rm -f *.o ikrobot

context.o: robot/context.cpp robot/context.h util/util.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<

ikrobot.o: ikrobot.cpp util/util.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<

link.o: robot/link.cpp robot/link.h util/util.h util/buffer.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<

robot.o: robot/robot.cpp robot/robot.h util/util.h util/buffer.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<

tree.o: robot/tree.cpp robot/tree.h util/util.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<
