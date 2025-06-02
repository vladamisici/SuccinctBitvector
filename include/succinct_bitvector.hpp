#pragma once

#include <vector>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <memory>
#include <cassert>

namespace succinct {

/**
 * @brief A succinct bitvector supporting constant-time rank and select operations
 * 
 * This data structure stores a bitvector in n + o(n) bits while supporting:
 * - rank₁(i): number of 1-bits in positions [0, i)
 * - select₁(k): position of the k-th 1-bit (0-indexed)
 * - rank₀(i): number of 0-bits in positions [0, i)
 * - select₀(k): position of the k-th 0-bit (0-indexed)
 * 
 * All operations run in O(1) time after O(n) preprocessing.
 */
class SuccinctBitvector {
public:
    // Constructors
    SuccinctBitvector() = default;
    explicit SuccinctBitvector(const std::vector<bool>& bits);
    explicit SuccinctBitvector(const std::string& bit_string);
    SuccinctBitvector(const uint64_t* raw_bits, uint64_t num_bits);
    
    // Copy and move semantics
    SuccinctBitvector(const SuccinctBitvector&) = default;
    SuccinctBitvector(SuccinctBitvector&&) noexcept = default;
    SuccinctBitvector& operator=(const SuccinctBitvector&) = default;
    SuccinctBitvector& operator=(SuccinctBitvector&&) noexcept = default;
    
    // Core operations
    /**
     * @brief Count of 1-bits in positions [0, i)
     * @param i Position (0 ≤ i ≤ size())
     * @return Number of 1-bits before position i
     */
    uint64_t rank1(uint64_t i) const;
    
    /**
     * @brief Position of the k-th 1-bit (0-indexed)
     * @param k Which 1-bit to find (0 ≤ k < rank1(size()))
     * @return Position of the k-th 1-bit
     */
    uint64_t select1(uint64_t k) const;
    
    /**
     * @brief Count of 0-bits in positions [0, i)
     * @param i Position (0 ≤ i ≤ size())
     * @return Number of 0-bits before position i
     */
    uint64_t rank0(uint64_t i) const { return i - rank1(i); }
    
    /**
     * @brief Position of the k-th 0-bit (0-indexed)
     * @param k Which 0-bit to find (0 ≤ k < rank0(size()))
     * @return Position of the k-th 0-bit
     */
    uint64_t select0(uint64_t k) const;
    
    // Utility operations
    /**
     * @brief Access bit at position i
     * @param i Position (0 ≤ i < size())
     * @return Value of bit at position i
     */
    bool access(uint64_t i) const;
    bool operator[](uint64_t i) const { return access(i); }
    
    /**
     * @brief Total number of bits
     */
    uint64_t size() const { return n_; }
    
    /**
     * @brief Total number of 1-bits
     */
    uint64_t count_ones() const { return rank1(n_); }
    
    /**
     * @brief Total number of 0-bits
     */
    uint64_t count_zeros() const { return n_ - count_ones(); }
    
    /**
     * @brief Check if bitvector is empty
     */
    bool empty() const { return n_ == 0; }
    
    // Memory usage statistics
    struct MemoryUsage {
        uint64_t bits_capacity;      // Capacity for storing bits
        uint64_t superblock_bytes;   // Bytes for superblock index
        uint64_t block_bytes;        // Bytes for block index
        uint64_t total_bytes;        // Total memory usage
        double bits_per_element;     // Average bits per element
    };
    
    MemoryUsage memory_usage() const;
    
    // Serialization
    void save(const std::string& filename) const;
    void load(const std::string& filename);
    
private:
    // Constants for the two-level index structure
    static constexpr uint64_t WORD_BITS = 64;
    static constexpr uint64_t SUPERBLOCK_BITS = 512;
    static constexpr uint64_t WORDS_PER_SUPERBLOCK = SUPERBLOCK_BITS / WORD_BITS;
    
    // Bit storage
    std::vector<uint64_t> words_;  // Packed bits
    uint64_t n_ = 0;              // Total number of bits
    
    // Index structures
    std::vector<uint64_t> superblocks_;  // Cumulative popcount every 512 bits
    std::vector<uint16_t> blocks_;       // Relative popcount within superblock
    
    // Select acceleration table
    struct SelectTable {
        static constexpr size_t TABLE_SIZE = 1 << 16;
        alignas(64) uint8_t positions[TABLE_SIZE][16];
        
        SelectTable();
    };
    static const SelectTable select_table_;
    
    // Helper functions
    static uint64_t popcount64(uint64_t x) noexcept;
    static uint64_t select64(uint64_t x, uint64_t k) noexcept;
    
    // Build index structures
    void build_index();
    
    // Select implementation helpers
    uint64_t select_in_word(uint64_t word, uint64_t k) const noexcept;
};

// Iterator for traversing set bits
class SuccinctBitvectorIterator {
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = uint64_t;
    using difference_type = std::ptrdiff_t;
    using pointer = const uint64_t*;
    using reference = const uint64_t&;
    
    SuccinctBitvectorIterator(const SuccinctBitvector* bv, uint64_t pos)
        : bv_(bv), pos_(pos) {}
    
    uint64_t operator*() const { return pos_; }
    
    SuccinctBitvectorIterator& operator++() {
        if (pos_ < bv_->size()) {
            uint64_t rank = bv_->rank1(pos_ + 1);
            if (rank < bv_->count_ones()) {
                pos_ = bv_->select1(rank);
            } else {
                pos_ = bv_->size();
            }
        }
        return *this;
    }
    
    bool operator==(const SuccinctBitvectorIterator& other) const {
        return bv_ == other.bv_ && pos_ == other.pos_;
    }
    
    bool operator!=(const SuccinctBitvectorIterator& other) const {
        return !(*this == other);
    }
    
private:
    const SuccinctBitvector* bv_;
    uint64_t pos_;
};

// Range-based for loop support
inline SuccinctBitvectorIterator begin(const SuccinctBitvector& bv) {
    return SuccinctBitvectorIterator(&bv, bv.count_ones() > 0 ? bv.select1(0) : bv.size());
}

inline SuccinctBitvectorIterator end(const SuccinctBitvector& bv) {
    return SuccinctBitvectorIterator(&bv, bv.size());
}

} // namespace succinct