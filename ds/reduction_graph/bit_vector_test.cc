#include "bit_vector.h"

#include <vector>

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(BitVectorTest, insert) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  BitVector vec = BitVector<int>(freeNodes.size());
  vec.insert(0);
  EXPECT_EQ(vec.getBitVectorSize(), 1);

  EXPECT_EQ(vec.find(0), true);
  EXPECT_EQ(vec.find(1), false);
  EXPECT_EQ(vec.find(2), false);
  vec.insert(1);
  EXPECT_EQ(vec.find(0), true);
  EXPECT_EQ(vec.find(1), true);
  EXPECT_EQ(vec.find(2), false);
}

TEST(BitVectorTest, erase) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  BitVector vec = BitVector<int>(freeNodes.size());
  vec.insert(0);
  EXPECT_EQ(vec.getBitVectorSize(), 1);
  EXPECT_EQ(vec.find(0), true);
  EXPECT_EQ(vec.find(1), false);
  EXPECT_EQ(vec.find(2), false);
  vec.erase(0);
  vec.erase(1);
  EXPECT_EQ(vec.find(0), false);
  EXPECT_EQ(vec.find(1), false);
  EXPECT_EQ(vec.find(2), false);
}

TEST(BitVectorTest, findSetBits) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  BitVector vec = BitVector<int>(freeNodes.size());
  std::vector<int> setNodes = vec.findSetBits();
  EXPECT_EQ(setNodes.size(), 0);
  vec.insert(0);
  EXPECT_EQ(vec.getBitVectorSize(), 1);
  EXPECT_EQ(vec.find(0), true);
  EXPECT_EQ(vec.find(1), false);
  EXPECT_EQ(vec.find(2), false);
  setNodes = vec.findSetBits();
  EXPECT_EQ(setNodes.size(), 1);
  EXPECT_EQ(setNodes[0], 0);
  vec.insert(2);
  setNodes = vec.findSetBits();
  EXPECT_EQ(setNodes.size(), 2);
  EXPECT_EQ(setNodes[0], 0);
  EXPECT_EQ(setNodes[1], 2);
  vec.erase(0);
  setNodes = vec.findSetBits();
  EXPECT_EQ(setNodes.size(), 1);
  EXPECT_EQ(setNodes[0], 2);
}

TEST(BitVectorTest, resuize) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  BitVector vec = BitVector<int>(freeNodes.size());
  EXPECT_EQ(vec.getBitVectorSize(), 1);
  vec.resize(100);
  EXPECT_EQ(vec.getBitVectorSize(), 2);
  vec.resize(1000);
  EXPECT_EQ(vec.getBitVectorSize(), (1000 / 64) + 1);
  vec.resize(64);
  EXPECT_EQ(vec.getBitVectorSize(), 1);
}