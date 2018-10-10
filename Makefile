CXX = g++ 
CXXFLAGS = -std=c++11
SOURCE = main.cpp simulation.cpp getparameters.cpp competition.cpp
EXEC = ecolattice

ecolat: $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(SOURCE)

clean:
	rm -f $(TARGET) *.o

