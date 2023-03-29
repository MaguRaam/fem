CXX = g++
CXXFLAGS = -O2 -I /home/magu/library/eigen
TARGET = main
SOURCES = main.cpp mesh.cpp fe.cpp
OBJECTS = $(SOURCES:.cpp=.o)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

.PHONY: clean distclean
clean:
	@rm -f $(TARGET) $(OBJECTS)

distclean: clean
	@rm -f *.dat 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
