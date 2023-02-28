CXX = g++
CXXFLAGS = -O2 -I /home/magu/library/eigen -I../
TARGET = klocal
SOURCES = klocal.cpp
OBJECTS = $(SOURCES:.cpp=.o)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

.PHONY: clean distclean
clean:
	@rm -f $(TARGET) $(OBJECTS)

distclean: clean
	@rm -f *.vtk *.png *.mtx

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
