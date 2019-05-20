src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)
exe = $(src:.cpp=)

LDFLAGS = -lm
CXXFLAGS = -O3 -g -std=c++11 -Wall

clique_sampler: clique_sampler.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

adaptive_enumerator: adaptive_enumerator.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

maximal_indep_set: maximal_indep_set.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	
all: $(exe)
.PHONY: clean
clean:
	rm -f $(obj) $(exe)
