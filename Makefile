CXX = g++
CXXFLAGS = -std=c++20 -Wall -O3 -I./inc -g -march=native
LDFLAGS = -lhts

SRCDIR = src
OBJDIR = obj
BINDIR = bin

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))
DEPS = $(OBJECTS:.o=.d)

EXEC = main

all: $(EXEC)

$(EXEC): $(OBJECTS)
	@echo "Linking $@..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(OBJDIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

-include $(DEPS)

clean:
	rm -rf $(OBJDIR) *.o *.d $(EXEC)

test: $(EXEC)
	@echo "Running test case:"
	./$(EXEC)
