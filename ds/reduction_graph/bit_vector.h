#include <iostream>
#include <vector>

template <typename NT>
class BitVector {
 private:
  std::vector<unsigned long long> data;

 public:
  using NodeType = NT;

  BitVector(NodeType num_nodes) {
    // Calculate the number of unsigned long longs needed to store the bits for all nodes
    NodeType num_ulongs = (num_nodes + 63) / 64;
    // Resize the data vector accordingly
    data.resize(num_ulongs);
  }

  // Constructor
  BitVector() {}

  void resize(NodeType num_nodes) {
    NodeType num_ulongs = (num_nodes + 63) / 64;
    // Resize the data vector accordingly
    data.resize(num_ulongs);
  }

  // Function to find the positions of all set bits in the bit vector
  std::vector<NodeType> findSetBits() const {
    std::vector<NodeType> positions;
    for (size_t i = 0; i < data.size(); ++i) {
      unsigned long long num = data[i];
      while (num != 0) {
        int bitPos = __builtin_ctzll(num) + i * 64;  // Adjust position by index
        positions.push_back(bitPos);
        num &= (num - 1);  // Clear the least significant set bit
      }
    }
    return positions;
  }

  // Function to set a specific bit at position 'pos' to 1
  void insert(NodeType pos) {
    NodeType index = pos / 64;   // Get index of the unsigned long long
    NodeType offset = pos % 64;  // Get bit offset within the unsigned long long
    if (index >= data.size()) data.resize(index + 1);
    data[index] |= (1ULL << offset);  // Set the bit at the specified position
  }

  // Function to clear a specific bit at position 'pos' to 0
  void erase(NodeType pos) {
    NodeType index = pos / 64;         // Get index of the unsigned long long
    NodeType offset = pos % 64;        // Get bit offset within the unsigned long long
    if (index >= data.size()) return;  // Do nothing if index is out of bounds
    data[index] &= ~(1ULL << offset);  // Clear the bit at the specified position
  }

  // Function to check if the bit at position 'pos' is set
  bool find(NodeType pos) const {
    NodeType index = pos / 64;                     // Get index of the unsigned long long
    NodeType offset = pos % 64;                    // Get bit offset within the unsigned long long
    if (index >= data.size()) return false;        // Bit is not set if index is out of bounds
    return ((data[index] >> offset) & 1ULL) != 0;  // Check if the bit is set
  }

  const NodeType getBitVectorSize() { return data.size(); }
};