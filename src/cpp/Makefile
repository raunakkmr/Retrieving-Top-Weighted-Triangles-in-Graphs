src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)

LDFLAGS = -lm
CXXFLAGS = -O3 -std=c++11 -Wall

clique_sampler: $(obj)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) clique_sampler
