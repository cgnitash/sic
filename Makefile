
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -Werror 

sicrun: a2m.o ensemble.o pwms.o clap.o
	 $(CXX) $(CXXFLAGS) a2m.o ensemble.o pwms.o clap.o -o sicrun

a2m.o: src/a2m.cpp src/ensemble.hpp src/pwms.hpp
	$(CXX) -c $(CXXFLAGS) src/a2m.cpp 

pwms.o: src/pwms.cpp src/pwms.hpp src/ensemble.hpp
	$(CXX) -c $(CXXFLAGS) src/pwms.cpp 

ensemble.o: src/ensemble.cpp src/ensemble.hpp
	$(CXX) -c $(CXXFLAGS) src/ensemble.cpp 

clap.o: src/clap.cpp src/clap.hpp
	$(CXX) -c $(CXXFLAGS) src/clap.cpp 

