src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)
exe = $(src:.cpp=)

LDFLAGS = -lm -pthread
CXXFLAGS = -O3 -g -std=c++17 -Wall

clique_sampler: clique_sampler.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	
all: $(exe)
.PHONY: clean
clean:
	rm -f $(obj) $(exe)
