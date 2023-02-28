CXX = g++
CXXFLAGS = -O3 -I /home/magu/library/eigen
TARGET = ex3
SOURCES = ex3.cpp
OBJECTS = $(SOURCES:.cpp=.o)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

.PHONY: clean
clean:
	@rm -f $(TARGET) $(OBJECTS)

.PHONY: distclean
distclean: clean
	@rm -f *.dat

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
