HOME_BREW ?= $(HOME)/homebrew/opt
NTL_ROOT ?= $(HOME_BREW)/ntl
GMP_ROOT ?= $(HOME_BREW)/gmp

CXX = g++
CXXFLAGS = -Wall -O2 -std=c++17 -I./src -I$(NTL_ROOT)/include -I$(GMP_ROOT)/include
NTL_LIB_DIR = $(NTL_ROOT)/lib
GMP_LIB_DIR = $(GMP_ROOT)/lib
RPATH_FLAGS = -Wl,-rpath,$(NTL_LIB_DIR) -Wl,-rpath,$(GMP_LIB_DIR)
LDFLAGS = -L$(NTL_LIB_DIR) -L$(GMP_LIB_DIR) -lntl -lgmp $(RPATH_FLAGS)

OBJDIR = obj
SRCDIR = src

SRCS = $(SRCDIR)/main.cpp $(SRCDIR)/tool.cpp $(SRCDIR)/kleinSamplar.cpp $(SRCDIR)/gaussSieve.cpp
OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))
TARGET = $(OBJDIR)/gsieve

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(OBJDIR)/*.o $(TARGET)
