CXX = g++
CXXFLAGS = -Wall -O2 -std=c++17 -I./src -I/Users/Hashiren1/local/ntl/include
LDFLAGS = -L/Users/Hashiren1/local/ntl/lib -L/opt/homebrew/lib -lntl -lgmp

OBJDIR = obj
SRCDIR = src

SRCS = $(SRCDIR)/main.cpp $(SRCDIR)/tool.cpp $(SRCDIR)/kleinSamplar.cpp
OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))
TARGET = $(OBJDIR)/main

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(OBJDIR)/*.o $(TARGET)