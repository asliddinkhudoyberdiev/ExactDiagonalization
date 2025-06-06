CXX = g++
CXXFLAGS = -O3 -march=native
OMPFLAGS = -fopenmp

# Executable name
TARGET = zeitentwicklung

# Source and object files
OBJS = einlesen.o operator.o initial.o main.o

# Default target: build, run, clean
all: $(TARGET)
	./$(TARGET)
	@$(MAKE) clean

# Link objects into final executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $(TARGET) $(OBJS)

# Compile source files
einlesen.o: einlesen.cpp
	$(CXX) -c $(CXXFLAGS) einlesen.cpp

operator.o: operator.cpp
	$(CXX) -c $(CXXFLAGS) $(OMPFLAGS) operator.cpp

initial.o: initial.cpp
	$(CXX) -c $(CXXFLAGS) initial.cpp

main.o: main.cpp
	$(CXX) -c $(CXXFLAGS) main.cpp

# Clean up build artifacts
clean:
	rm -f $(TARGET) *.o
