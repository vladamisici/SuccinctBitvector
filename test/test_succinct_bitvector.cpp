#include "succinct_bitvector.hpp"
#include <iostream>
#include <random>
#include <chrono>
#include <cassert>
#include <iomanip>

using namespace succinct;

// Test helper: brute force rank
uint64_t brute_force_rank1(const std::vector<bool>& bits, uint64_t i) {
    uint64_t count = 0;
    for (uint64_t j = 0; j < i && j < bits.size(); ++j) {
        if (bits[j]) ++count;
    }
    return count;
}

// Test helper: brute force select
uint64_t brute_force_select1(const std::vector<bool>& bits, uint64_t k) {
    uint64_t count = 0;
    for (uint64_t i = 0; i < bits.size(); ++i) {
        if (bits[i]) {
            if (count == k) return i;
            ++count;
        }
    }
    return bits.size();
}

// Basic functionality tests
void test_basic_operations() {
    std::cout << "Testing basic operations..." << std::endl;
    
    // Test 1: Small bitvector
    {
        std::vector<bool> bits = {1, 0, 1, 1, 0, 0, 1, 0, 1, 1};
        SuccinctBitvector bv(bits);
        
        assert(bv.size() == 10);
        assert(bv.count_ones() == 6);
        assert(bv.count_zeros() == 4);
        
        // Test rank1
        assert(bv.rank1(0) == 0);
        assert(bv.rank1(1) == 1);
        assert(bv.rank1(5) == 3);
        assert(bv.rank1(10) == 6);
        
        // Test rank0
        assert(bv.rank0(0) == 0);
        assert(bv.rank0(2) == 1);
        assert(bv.rank0(6) == 3);
        assert(bv.rank0(10) == 4);
        
        // Test select1
        assert(bv.select1(0) == 0);
        assert(bv.select1(1) == 2);
        assert(bv.select1(2) == 3);
        assert(bv.select1(5) == 9);
        
        // Test select0
        assert(bv.select0(0) == 1);
        assert(bv.select0(1) == 4);
        assert(bv.select0(3) == 7);
        
        // Test access
        for (uint64_t i = 0; i < bits.size(); ++i) {
            assert(bv.access(i) == bits[i]);
            assert(bv[i] == bits[i]);
        }
    }
    
    // Test 2: String constructor
    {
        SuccinctBitvector bv("10110010110");
        assert(bv.size() == 11);
        assert(bv.count_ones() == 6);
        assert(bv.rank1(5) == 3);
        assert(bv.select1(3) == 6);
    }
    
    // Test 3: Edge cases
    {
        // Empty bitvector
        SuccinctBitvector empty(std::vector<bool>{});
        assert(empty.size() == 0);
        assert(empty.count_ones() == 0);
        assert(empty.empty());
        
        // All zeros
        std::vector<bool> zeros(100, false);
        SuccinctBitvector bv_zeros(zeros);
        assert(bv_zeros.count_ones() == 0);
        assert(bv_zeros.count_zeros() == 100);
        assert(bv_zeros.rank1(50) == 0);
        assert(bv_zeros.rank0(50) == 50);
        
        // All ones
        std::vector<bool> ones(100, true);
        SuccinctBitvector bv_ones(ones);
        assert(bv_ones.count_ones() == 100);
        assert(bv_ones.count_zeros() == 0);
        assert(bv_ones.rank1(50) == 50);
        assert(bv_ones.rank0(50) == 0);
    }
    
    std::cout << "✓ Basic operations test passed" << std::endl;
}

// Correctness test against brute force
void test_correctness() {
    std::cout << "Testing correctness with random data..." << std::endl;
    
    std::mt19937 gen(42);
    std::bernoulli_distribution dist(0.3);  // 30% density
    
    for (int test = 0; test < 10; ++test) {
        // Generate random bitvector
        uint64_t n = 1000 + test * 1000;
        std::vector<bool> bits(n);
        for (uint64_t i = 0; i < n; ++i) {
            bits[i] = dist(gen);
        }
        
        SuccinctBitvector bv(bits);
        
        // Test rank1 at random positions
        for (int i = 0; i < 100; ++i) {
            uint64_t pos = std::uniform_int_distribution<uint64_t>(0, n)(gen);
            assert(bv.rank1(pos) == brute_force_rank1(bits, pos));
        }
        
        // Test select1 for all valid k
        uint64_t ones = bv.count_ones();
        if (ones > 0) {
            for (int i = 0; i < std::min(100ULL, ones); ++i) {
                uint64_t k = std::uniform_int_distribution<uint64_t>(0, ones - 1)(gen);
                assert(bv.select1(k) == brute_force_select1(bits, k));
            }
        }
        
        // Test access
        for (int i = 0; i < 100; ++i) {
            uint64_t pos = std::uniform_int_distribution<uint64_t>(0, n - 1)(gen);
            assert(bv.access(pos) == bits[pos]);
        }
    }
    
    std::cout << "✓ Correctness test passed" << std::endl;
}

// Test large bitvectors
void test_large_scale() {
    std::cout << "Testing large bitvectors..." << std::endl;
    
    // Test with 100M bits
    uint64_t n = 100'000'000;
    std::mt19937 gen(12345);
    std::bernoulli_distribution dist(0.5);
    
    // Generate in chunks to avoid excessive memory for vector<bool>
    std::vector<uint64_t> raw_words((n + 63) / 64);
    for (auto& word : raw_words) {
        word = std::uniform_int_distribution<uint64_t>()(gen);
    }
    
    SuccinctBitvector bv(raw_words.data(), n);
    
    // Memory usage
    auto mem = bv.memory_usage();
    std::cout << "  Size: " << n / 1'000'000.0 << "M bits" << std::endl;
    std::cout << "  Memory: " << mem.total_bytes / (1024.0 * 1024.0) << " MB" << std::endl;
    std::cout << "  Overhead: " << mem.bits_per_element - 1.0 << " bits/element" << std::endl;
    std::cout << "  Compression: " << (mem.bits_per_element / 8.0) * 100.0 << "%" << std::endl;
    
    // Performance test
    auto start = std::chrono::high_resolution_clock::now();
    
    // Random rank queries
    uint64_t checksum = 0;
    for (int i = 0; i < 1'000'000; ++i) {
        uint64_t pos = std::uniform_int_distribution<uint64_t>(0, n)(gen);
        checksum += bv.rank1(pos);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    
    std::cout << "  Rank performance: " << duration / 1'000'000.0 << " ns/query" << std::endl;
    std::cout << "  Throughput: " << 1'000'000'000.0 / duration << " M queries/sec" << std::endl;
    
    // Select performance
    uint64_t ones = bv.count_ones();
    start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < 1'000'000; ++i) {
        uint64_t k = std::uniform_int_distribution<uint64_t>(0, ones - 1)(gen);
        checksum += bv.select1(k);
    }
    
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    
    std::cout << "  Select performance: " << duration / 1'000'000.0 << " ns/query" << std::endl;
    std::cout << "  Throughput: " << 1'000'000'000.0 / duration << " M queries/sec" << std::endl;
    
    // Prevent optimization
    volatile uint64_t dummy = checksum;
    (void)dummy;
    
    std::cout << "✓ Large scale test passed" << std::endl;
}

// Test serialization
void test_serialization() {
    std::cout << "Testing serialization..." << std::endl;
    
    std::mt19937 gen(999);
    std::bernoulli_distribution dist(0.4);
    
    // Create original bitvector
    uint64_t n = 10000;
    std::vector<bool> bits(n);
    for (uint64_t i = 0; i < n; ++i) {
        bits[i] = dist(gen);
    }
    
    SuccinctBitvector bv1(bits);
    
    // Save and load
    bv1.save("test_bitvector.bin");
    SuccinctBitvector bv2;
    bv2.load("test_bitvector.bin");
    
    // Verify
    assert(bv1.size() == bv2.size());
    assert(bv1.count_ones() == bv2.count_ones());
    
    for (uint64_t i = 0; i < n; ++i) {
        assert(bv1.access(i) == bv2.access(i));
    }
    
    for (uint64_t i = 0; i <= n; ++i) {
        assert(bv1.rank1(i) == bv2.rank1(i));
    }
    
    uint64_t ones = bv1.count_ones();
    for (uint64_t k = 0; k < ones; ++k) {
        assert(bv1.select1(k) == bv2.select1(k));
    }
    
    std::cout << "✓ Serialization test passed" << std::endl;
}

// Test iterator interface
void test_iterator() {
    std::cout << "Testing iterator interface..." << std::endl;
    
    std::vector<bool> bits = {0, 1, 0, 0, 1, 1, 0, 1, 0, 1};
    SuccinctBitvector bv(bits);
    
    // Collect positions using iterator
    std::vector<uint64_t> positions;
    for (uint64_t pos : bv) {
        positions.push_back(pos);
    }
    
    // Verify
    std::vector<uint64_t> expected = {1, 4, 5, 7, 9};
    assert(positions == expected);
    
    std::cout << "✓ Iterator test passed" << std::endl;
}

int main() {
    std::cout << "Running Succinct Bitvector Tests" << std::endl;
    std::cout << "=================================" << std::endl;
    
    try {
        test_basic_operations();
        test_correctness();
        test_large_scale();
        test_serialization();
        test_iterator();
        
        std::cout << "\n✅ All tests passed!" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "\n❌ Test failed: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}