#include <Rcpp.h>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <string>

using namespace Rcpp;

// Helper function to create a unique key for each coordinate
std::string coord_to_key(int x, int y, int z) {
  return std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z);
}

// [[Rcpp::export]]
List build_neighbor_list(IntegerMatrix coords) {
  int V = coords.nrow();
  List neighbor_list(V);
  
  // Create a hash set for fast lookup of existing voxels
  std::unordered_set<std::string> voxel_set;
  std::unordered_map<std::string, int> coord_to_index;
  
  for (int i = 0; i < V; ++i) {
    int x = coords(i, 0);
    int y = coords(i, 1);
    int z = coords(i, 2);
    std::string key = coord_to_key(x, y, z);
    voxel_set.insert(key);
    coord_to_index[key] = i;
  }
  
  // Offsets for immediate neighbors along x, y, z axes
  int offsets[6][3] = {
    { -1,  0,  0 },
    {  1,  0,  0 },
    {  0, -1,  0 },
    {  0,  1,  0 },
    {  0,  0, -1 },
    {  0,  0,  1 }
  };
  
  // For each voxel, find its neighbors
  for (int i = 0; i < V; ++i) {
    int x = coords(i, 0);
    int y = coords(i, 1);
    int z = coords(i, 2);
    
    std::vector<int> neighbors;
    
    // Check each of the six possible neighbors
    for (int o = 0; o < 6; ++o) {
      int nx = x + offsets[o][0];
      int ny = y + offsets[o][1];
      int nz = z + offsets[o][2];
      
      // Create key for neighbor coordinate
      std::string neighbor_key = coord_to_key(nx, ny, nz);
      
      // Check if neighbor exists in the voxel set
      if (voxel_set.find(neighbor_key) != voxel_set.end()) {
        int neighbor_idx = coord_to_index[neighbor_key];
        neighbors.push_back(neighbor_idx + 1); // Convert to 1-based indexing for R
      }
    }
    
    neighbor_list[i] = wrap(neighbors);
  }
  
  return neighbor_list;
}

// [[Rcpp::export]]
NumericVector compute_average_inclusion(NumericVector delta_vec, List neighbor_list) {
  int V = delta_vec.size();
  NumericVector pj_vec(V);
  
  for (int i = 0; i < V; ++i) {
    double sum = 0.0;
    
    IntegerVector neighbors = neighbor_list[i];
    for (int j = 0; j < neighbors.size(); ++j) {
      int idx = neighbors[j] - 1; // Convert to 0-based indexing
      sum += delta_vec[idx];
    }
    
    pj_vec[i] = sum / neighbors.size();
  }
  
  return pj_vec;
}