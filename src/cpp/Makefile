src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)
exe = $(src:.cpp=)

LDFLAGS = -lm -pthread
CXXFLAGS = -O3 -g -std=c++17 -Wall

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
