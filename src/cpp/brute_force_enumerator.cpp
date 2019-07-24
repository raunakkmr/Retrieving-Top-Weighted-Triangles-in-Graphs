#include <bits/stdc++.h>

#include "graph.h"
#include "triangle_sampler.h"

using namespace std;
using namespace wsdm_2019_graph;

int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    srand(0);

    // To use binary data files, add a character after dataset.
    // So run as ./brute_force_enumerator filename dataset {something here if we want to use binary otherwise nothing}.
    auto G = read_graph(argv[1], argc > 3);
    string dataset = argv[2];

    auto all_tris = brute_force_sampler(G);
    int num_all_tris = all_tris.size();

    auto get_bytes = [](int x) {
      if (x <= numeric_limits<char>::max()) {
        return 0;
      } else if (x <= numeric_limits<short>::max()) {
        return 1;
      } else if (x <= numeric_limits<int>::max()) {
        return 3;
      }
      return -1;
    };

    // ofstream triangle_file(dataset+"-triangles.txt");
    ofstream triangle_file(dataset+"-triangles.binary", ios::binary);
    binary_write(triangle_file, num_all_tris);
    int est_bytes = 0;
    map<long long, long long> num_tris;
    for (auto t : all_tris) {
      int w = t.weight;  // TODO: change to long long.
      int a = get<0>(t.vertices);
      int b = get<1>(t.vertices);
      int c = get<2>(t.vertices);
      // binary_write(triangle_file, w);
      // binary_write(triangle_file, a);
      // binary_write(triangle_file, b);
      // binary_write(triangle_file, c);
      int bytes = 0;
      // uses 2 as a special code that the weight is 1.
      // call it premature optimization but it saves a lot
      bytes |= (w == 3 ? 2 : get_bytes(w));
      bytes <<= 2;
      bytes |= get_bytes(a);
      bytes <<= 2;
      bytes |= get_bytes(b);
      bytes <<= 2;
      bytes |= get_bytes(c);

      binary_compressed_write(triangle_file, bytes);
      binary_compressed_write(triangle_file, a);
      binary_compressed_write(triangle_file, b);
      binary_compressed_write(triangle_file, c);
      est_bytes += 4 + get_bytes(a) + get_bytes(b) + get_bytes(c);
      if (w > 3) {
        binary_compressed_write(triangle_file, w);
        est_bytes += 1 + get_bytes(w);
      }

      num_tris[get<0>(t.vertices)]++;
      num_tris[get<1>(t.vertices)]++;
      num_tris[get<2>(t.vertices)]++;
    }
    triangle_file.close();

    // ofstream tris_to_weight_file(dataset+"-ntris-to-weight.txt");
    ofstream tris_to_weight_file(dataset+"-ntris-to-weight.binary", ios::binary);
    int g_sz = G.size();
    binary_write(tris_to_weight_file, g_sz);
    for (int i = 0; i < (int) G.size(); i++) {
      double weight_sum = 0;
      for (auto e : G[i]) {
        weight_sum += e.wt;
      }
      // tris_to_weight_file << weight_sum << " " << num_tris[i] << '\n';
      binary_write(tris_to_weight_file, weight_sum);
      binary_write(tris_to_weight_file, num_tris[i]);
    }
    tris_to_weight_file.close();

    // ofstream degree_file(dataset+"-degree-to-weight.txt");
    ofstream degree_file(dataset+"-degree-to-weight.binary", ios::binary);
    binary_write(degree_file, g_sz);
    for (int i = 0; i < (int) G.size(); i++) {
      double weight_sum = 0;
      for (auto e : G[i]) {
        weight_sum += e.wt;
      }
      // degree_file << G[i].size() << " " << weight_sum << '\n';
      int gi_sz = G[i].size();
      binary_write(degree_file, gi_sz);
      binary_write(degree_file, weight_sum);
    }
    degree_file.close();

    return 0;
}