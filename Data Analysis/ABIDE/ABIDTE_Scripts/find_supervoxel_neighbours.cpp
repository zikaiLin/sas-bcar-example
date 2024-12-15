#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// [[Rcpp::export]]
List findSupervoxelNeighborsOptimized(NumericVector Acube, IntegerVector dims, IntegerVector unique_labels) {
  int nx = dims[0], ny = dims[1], nz = dims[2];
  int n_voxels = nx * ny * nz;
  
  // Convert to a flat vector for easier indexing
  std::vector<int> Acube_flat(Acube.begin(), Acube.end());
  
  // Initialize the list to store neighbors
  std::unordered_map<int, std::unordered_set<int>> neighbors_map;
  
  // Define 6-connected neighborhood offsets
  int offsets[6][3] = {{-1, 0, 0}, {1, 0, 0},
                       {0, -1, 0}, {0, 1, 0},
                       {0, 0, -1}, {0, 0, 1}};
  
  // Iterate over unique labels
  for (int label : unique_labels) {
    Rcpp::Rcout << "Processing label: " << label << "\n";
    std::unordered_set<int> neighbors;
    if (label == 0){
      continue;
    }
    int voxel_count = 0;
    
    // Find indices of all voxels with the current label
    for (int idx = 0; idx < n_voxels; ++idx) {
      if (Acube_flat[idx] == label) {
        voxel_count++;
        int z = idx / (nx * ny);
        int y = (idx % (nx * ny)) / nx;
        int x = idx % nx;
        
        // Check all neighbors
        for (int k = 0; k < 6; ++k) {
          int nx = x + offsets[k][0];
          int ny = y + offsets[k][1];
          int nz = z + offsets[k][2];
          
          // Ensure neighbor is within bounds
          if (nx >= 0 && nx < dims[0] && ny >= 0 && ny < dims[1] && nz >= 0 && nz < dims[2]) {
            int neighbor_idx = nx + dims[0] * (ny + dims[1] * nz);
            int neighbor_label = Acube_flat[neighbor_idx];
            
            // Add neighbor label if valid and not the same
            if (neighbor_label > 0 && neighbor_label != label) {
              neighbors.insert(neighbor_label);
            }
          }
        }
      }
    }
    
    // Debug print: neighbors found
    // Rcpp::Rcout << "  Total voxels for label " << label << ": " << voxel_count << "\n";
    // Rcpp::Rcout << "  Neighbors for label " << label << ": ";
    // for (const auto &neighbor : neighbors) {
    //   Rcpp::Rcout << neighbor << " ";
    // }
    // Rcpp::Rcout << "\n";
    
    // Store unique neighbors for the current label
    neighbors_map[label] = neighbors;
  }
  
  // Convert the map to an R list
  List result(unique_labels.size());
  for (size_t i = 0; i < unique_labels.size(); ++i) {
    int label = unique_labels[i];
    result[i] = IntegerVector(neighbors_map[label].begin(), neighbors_map[label].end());
  }
  result.attr("names") = unique_labels;
  return result;
}