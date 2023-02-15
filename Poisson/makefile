CXX = g++
CXXFLAGS = -O3 -I /home/magu/library/eigen
TARGET = poisson
SOURCES = poisson.cpp
OBJECTS = $(SOURCES:.cpp=.o)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

.PHONY: clean
clean:
	rm $(TARGET) $(OBJECTS)

.PHONY: distclean
distclean: clean
	rm *.dat

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
