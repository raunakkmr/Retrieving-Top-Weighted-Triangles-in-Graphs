src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)
exe = $(src:.cpp=)

LDFLAGS = -fopenmp
CXXFLAGS = -fopenmp -Ofast -march=native -std=c++11 -Wall #-D_GLIBCXX_PARALLEL

test_reading: test_reading.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

convert_data: convert_data.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

brute_force_enumerator: brute_force_enumerator.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clique_sampler: clique_sampler.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

compare_deterministic: compare_deterministic.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

compare_sampling: compare_sampling.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	
all: $(exe)
.PHONY: clean
clean:
	rm -f $(obj) $(exe)
