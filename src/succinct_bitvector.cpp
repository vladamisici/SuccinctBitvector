#include <algorithm>
#include <cstdint>
#include <fstream>
#include <cstring>

#include "succinct_bitvector.hpp"

namespace succinct {

// Initialize the static select acceleration table
const SuccinctBitvector::SelectTable SuccinctBitvector::select_table_;

SuccinctBitvector::SelectTable::SelectTable() {
    // Build lookup table for all 16-bit values
    for (std::uint32_t v = 0; v < TABLE_SIZE; ++v) {
        uint16_t bits = static_cast<uint16_t>(v);
        uint8_t pos = 0;
        
        // Find position of each set bit
        for (uint8_t i = 0; i < 16; ++i) {
            if ((bits >> i) & 1) {
                if (pos < 16) {
                    positions[v][pos] = i;
                }
                ++pos;
            }
        }
        
        // Fill remaining positions with sentinel value
        for (; pos < 16; ++pos) {
            positions[v][pos] = 255;
        }
    }
}

// Constructor from vector<bool>
SuccinctBitvector::SuccinctBitvector(const std::vector<bool>& bits) : n_(bits.size()) {
    if (n_ == 0) return;
    
    // Allocate storage
    uint64_t num_words = (n_ + WORD_BITS - 1) / WORD_BITS;
    words_.resize(num_words, 0);
    
    // Pack bits into words
    for (uint64_t i = 0; i < n_; ++i) {
        if (bits[i]) {
            uint64_t word_idx = i / WORD_BITS;
            uint64_t bit_idx = i % WORD_BITS;
            words_[word_idx] |= (1ULL << bit_idx);
        }
    }
    
    build_index();
}

// Constructor from string
SuccinctBitvector::SuccinctBitvector(const std::string& bit_string) : n_(bit_string.size()) {
    if (n_ == 0) return;
    
    // Allocate storage
    uint64_t num_words = (n_ + WORD_BITS - 1) / WORD_BITS;
    words_.resize(num_words, 0);
    
    // Pack bits from string
    for (uint64_t i = 0; i < n_; ++i) {
        if (bit_string[i] == '1') {
            uint64_t word_idx = i / WORD_BITS;
            uint64_t bit_idx = i % WORD_BITS;
            words_[word_idx] |= (1ULL << bit_idx);
        } else if (bit_string[i] != '0') {
            throw std::invalid_argument("Bit string must contain only '0' and '1' characters");
        }
    }
    
    build_index();
}

// Constructor from raw bits
SuccinctBitvector::SuccinctBitvector(const uint64_t* raw_bits, uint64_t num_bits) : n_(num_bits) {
    if (n_ == 0) return;
    
    // Copy raw bits
    uint64_t num_words = (n_ + WORD_BITS - 1) / WORD_BITS;
    words_.resize(num_words);
    std::memcpy(words_.data(), raw_bits, num_words * sizeof(uint64_t));
    
    // Clear unused bits in last word
    uint64_t used_bits_in_last_word = n_ % WORD_BITS;
    if (used_bits_in_last_word != 0) {
        uint64_t mask = (1ULL << used_bits_in_last_word) - 1;
        words_.back() &= mask;
    }
    
    build_index();
}

// Bit-parallel popcount implementation
uint64_t SuccinctBitvector::popcount64(uint64_t x) noexcept {
    // Use compiler intrinsic if available
    #ifdef __GNUC__
        return __builtin_popcountll(x);
    #elif defined(_MSC_VER)
        return __popcnt64(x);
    #else
        // Fallback to bit-parallel algorithm
        x = x - ((x >> 1) & 0x5555555555555555ULL);
        x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
        x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
        x = x + (x >> 8);
        x = x + (x >> 16);
        x = x + (x >> 32);
        return x & 0x7F;
    #endif
}

// Build the two-level index structure
void SuccinctBitvector::build_index() {
    uint64_t num_words = words_.size();
    if (num_words == 0) return;
    
    // Allocate index structures
    uint64_t num_superblocks = (num_words + WORDS_PER_SUPERBLOCK - 1) / WORDS_PER_SUPERBLOCK;
    superblocks_.resize(num_superblocks);
    blocks_.resize(num_words);
    
    // Build index in a single pass
    uint64_t total_ones = 0;
    uint64_t superblock_ones = 0;
    
    for (uint64_t w = 0; w < num_words; ++w) {
        // Start of new superblock?
        if (w % WORDS_PER_SUPERBLOCK == 0) {
            uint64_t sb_idx = w / WORDS_PER_SUPERBLOCK;
            superblocks_[sb_idx] = total_ones;
            superblock_ones = 0;
        }
        
        // Store relative popcount within superblock
        blocks_[w] = static_cast<uint16_t>(superblock_ones);
        
        // Update counters
        uint64_t word_ones = popcount64(words_[w]);
        superblock_ones += word_ones;
        total_ones += word_ones;
    }
}

// Rank operation: count 1-bits in [0, i)
uint64_t SuccinctBitvector::rank1(uint64_t i) const {
    if (i == 0) return 0;
    if (i > n_) {
        throw std::out_of_range("rank1: position out of range");
    }
    
    // Decompose position
    uint64_t pos = i - 1;  // We want bits [0, i), so look at position i-1
    uint64_t word_idx = pos / WORD_BITS;
    uint64_t bit_idx = pos % WORD_BITS;
    
    // Sum contributions from three levels
    uint64_t sb_idx = word_idx / WORDS_PER_SUPERBLOCK;
    uint64_t superblock_sum = superblocks_[sb_idx];
    uint64_t block_sum = blocks_[word_idx];
    
    // Count bits within the word up to bit_idx (inclusive)
    uint64_t mask = (bit_idx == 63) ? ~0ULL : ((1ULL << (bit_idx + 1)) - 1);
    uint64_t word_sum = popcount64(words_[word_idx] & mask);
    
    return superblock_sum + block_sum + word_sum;
}

// Select operation: find position of k-th 1-bit
uint64_t SuccinctBitvector::select1(uint64_t k) const {
    uint64_t total_ones = count_ones();
    if (k >= total_ones) {
        throw std::out_of_range("select1: k out of range");
    }
    
    // Binary search on superblocks to find containing superblock
    uint64_t left = 0;
    uint64_t right = superblocks_.size() - 1;
    
    while (left < right) {
        uint64_t mid = left + (right - left + 1) / 2;
        if (superblocks_[mid] <= k) {
            left = mid;
        } else {
            right = mid - 1;
        }
    }
    
    uint64_t sb_idx = left;
    uint64_t sb_start = sb_idx * WORDS_PER_SUPERBLOCK;
    uint64_t k_in_sb = k - superblocks_[sb_idx];
    
    // Linear search within superblock to find containing word
    uint64_t word_idx = sb_start;
    uint64_t word_end = std::min(sb_start + WORDS_PER_SUPERBLOCK, words_.size());
    
    while (word_idx + 1 < word_end) {
        uint64_t next_block_sum = blocks_[word_idx + 1];
        if (next_block_sum > k_in_sb) {
            break;
        }
        ++word_idx;
    }
    
    // Find bit within word
    uint64_t k_in_word = k_in_sb - blocks_[word_idx];
    uint64_t bit_pos = select_in_word(words_[word_idx], k_in_word);
    
    return word_idx * WORD_BITS + bit_pos;
}

// Helper: find k-th set bit within a word using lookup table
uint64_t SuccinctBitvector::select_in_word(uint64_t word, uint64_t k) const noexcept {
    // Process word in 16-bit chunks
    for (int chunk = 0; chunk < 4; ++chunk) {
        uint16_t bits = (word >> (chunk * 16)) & 0xFFFF;
        uint64_t chunk_ones = popcount64(bits);
        
        if (k < chunk_ones) {
            // k-th bit is in this chunk
            return chunk * 16 + select_table_.positions[bits][k];
        }
        k -= chunk_ones;
    }
    
    // Should never reach here if k is valid
    assert(false);
    return 0;
}

// Select operation for 0-bits
uint64_t SuccinctBitvector::select0(uint64_t k) const {
    uint64_t total_zeros = count_zeros();
    if (k >= total_zeros) {
        throw std::out_of_range("select0: k out of range");
    }
    
    // Binary search to find position where rank0(pos) = k + 1
    uint64_t left = 0;
    uint64_t right = n_;
    
    while (left < right) {
        uint64_t mid = left + (right - left) / 2;
        if (rank0(mid) <= k) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    
    return left - 1;
}

// Access bit at position
bool SuccinctBitvector::access(uint64_t i) const {
    if (i >= n_) {
        throw std::out_of_range("access: position out of range");
    }
    
    uint64_t word_idx = i / WORD_BITS;
    uint64_t bit_idx = i % WORD_BITS;
    return (words_[word_idx] >> bit_idx) & 1;
}

// Memory usage statistics
SuccinctBitvector::MemoryUsage SuccinctBitvector::memory_usage() const {
    MemoryUsage usage;
    usage.bits_capacity = words_.capacity() * sizeof(uint64_t);
    usage.superblock_bytes = superblocks_.capacity() * sizeof(uint64_t);
    usage.block_bytes = blocks_.capacity() * sizeof(uint16_t);
    usage.total_bytes = usage.bits_capacity + usage.superblock_bytes + usage.block_bytes;
    usage.bits_per_element = n_ > 0 ? (usage.total_bytes * 8.0) / n_ : 0.0;
    return usage;
}

// Save to file
void SuccinctBitvector::save(const std::string& filename) const {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    // Write header
    ofs.write(reinterpret_cast<const char*>(&n_), sizeof(n_));
    
    // Write data
    uint64_t words_size = words_.size();
    ofs.write(reinterpret_cast<const char*>(&words_size), sizeof(words_size));
    ofs.write(reinterpret_cast<const char*>(words_.data()), words_size * sizeof(uint64_t));
    
    uint64_t sb_size = superblocks_.size();
    ofs.write(reinterpret_cast<const char*>(&sb_size), sizeof(sb_size));
    ofs.write(reinterpret_cast<const char*>(superblocks_.data()), sb_size * sizeof(uint64_t));
    
    uint64_t blocks_size = blocks_.size();
    ofs.write(reinterpret_cast<const char*>(&blocks_size), sizeof(blocks_size));
    ofs.write(reinterpret_cast<const char*>(blocks_.data()), blocks_size * sizeof(uint16_t));
}

// Load from file
void SuccinctBitvector::load(const std::string& filename) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Cannot open file for reading: " + filename);
    }
    
    // Read header
    ifs.read(reinterpret_cast<char*>(&n_), sizeof(n_));
    
    // Read data
    uint64_t words_size;
    ifs.read(reinterpret_cast<char*>(&words_size), sizeof(words_size));
    words_.resize(words_size);
    ifs.read(reinterpret_cast<char*>(words_.data()), words_size * sizeof(uint64_t));
    
    uint64_t sb_size;
    ifs.read(reinterpret_cast<char*>(&sb_size), sizeof(sb_size));
    superblocks_.resize(sb_size);
    ifs.read(reinterpret_cast<char*>(superblocks_.data()), sb_size * sizeof(uint64_t));
    
    uint64_t blocks_size;
    ifs.read(reinterpret_cast<char*>(&blocks_size), sizeof(blocks_size));
    blocks_.resize(blocks_size);
    ifs.read(reinterpret_cast<char*>(blocks_.data()), blocks_size * sizeof(uint16_t));
}

} // namespace succinct