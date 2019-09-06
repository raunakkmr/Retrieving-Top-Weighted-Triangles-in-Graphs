src = $(wildcard *.cpp)
headers = $(wildcard *.h)
obj = $(src:.cpp=.o)
exe = $(src:.cpp=)

LDFLAGS = -fopenmp -fPIC
CXXFLAGS = -fopenmp -Ofast -march=native -std=c++11 -Wall -fPIC -I./ #-D_GLIBCXX_PARALLEL
OBJ_GFLAGS = $(addprefix gflags/, gflags.o gflags_nc.o gflags_completions.o gflags_reporting.o)

gflags/%.o: include/gflags/%.cc
	$(CXX) -c $(DEBUG) $(CXXFLAGS) $< -o $@

test_reading: $(OBJ_GFLAGS) $(headers) test_reading.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

convert_data: $(OBJ_GFLAGS) $(headers) convert_data.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

brute_force_enumerator: $(OBJ_GFLAGS) $(headers) brute_force_enumerator.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clique_sampler: $(OBJ_GFLAGS) $(headers) clique_sampler.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

compare_deterministic: $(OBJ_GFLAGS) $(headers) compare_deterministic.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

compare_sampling: $(OBJ_GFLAGS) $(headers) compare_sampling.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

static_heavy_light: $(OBJ_GFLAGS) $(headers) static_heavy_light.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

compare_parallel_sampling: $(OBJ_GFLAGS) $(headers) compare_parallel_sampling.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

print_edge_weights: $(OBJ_GFLAGS) $(headers) print_edge_weights.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	
all: $(exe)
.PHONY: clean
clean:
	rm -f $(obj) $(exe)
