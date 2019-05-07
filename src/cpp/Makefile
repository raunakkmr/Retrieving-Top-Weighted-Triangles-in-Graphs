src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)

LDFLAGS = -lm
CXXFLAGS = -O3 -g -std=c++11 -Wall

clique_sampler: $(obj)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) clique_sampler
