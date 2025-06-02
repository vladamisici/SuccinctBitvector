# succinct-bitvector

Fast rank/select operations on compressed bitstrings.

## What is this?

A C++17 implementation of a succinct bitvector - basically a compressed bitstring that still allows O(1) access to:
- `rank(i)`: count of 1s up to position i
- `select(k)`: position of the k-th 1

Uses ~37.5% extra space on top of the raw bits. Not the theoretical minimum but fast as hell.

## Performance

On my i7-8700K:

```
rank:   52ns  (19M ops/sec)
select: 382ns (2.6M ops/sec)  
access: 26ns  (38M ops/sec)
```

Cache effects are significant - sequential access is 3.5x faster than random.

## Build

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```

Run tests:
```bash
./test/test_succinct_bitvector
```

Run benchmarks:
```bash
./bench/bench_succinct_bitvector
```

## Usage

```cpp
#include "succinct_bitvector.hpp"

// Build from string
SuccinctBitvector bv("10110010");

// Operations
uint64_t ones_before_5 = bv.rank1(5);      // 3
uint64_t third_one_pos = bv.select1(2);    // 3 (0-indexed)
bool bit_at_4 = bv[4];                     // false

// Build from large data
std::vector<bool> bits(1'000'000);
// ... fill bits ...
SuccinctBitvector large(bits);
```

## How it works

Two-level index structure:
1. **Superblocks** (every 512 bits): Store absolute popcount
2. **Blocks** (every 64 bits): Store relative popcount within superblock
3. **Words**: Raw 64-bit chunks

For select, I use a 16-bit lookup table to find bit positions in O(1) within words.

```
Memory layout:
[word0|word1|...|word7] [word8|word9|...|word15] ...
        block0                  block1
|--------- superblock0 ---------|
```

## Space breakdown

For n bits:
- Raw bits: n bits
- Superblocks: n/512 * 64 bits = n/8 bits
- Blocks: n/64 * 9 bits ≈ n/7 bits
- Total overhead: ~0.375n bits (37.5%)

Could get this down to ~12% with fancier encoding but the speed tradeoff isn't worth it IMO.

## CLI tool

```bash
# Build bitvector
./bvtool build data.bv 101100101

# Query
./bvtool rank data.bv 5
./bvtool select data.bv 2

# Analyze
./bvtool analyze data.bv
```

## Implementation notes

- Uses `__builtin_popcountll` when available, falls back to bit-parallel algorithm
- 16-bit select table is built once at startup (1MB static memory)
- All query ops are thread-safe (read-only after construction)
- Tried SIMD but single-word operations are already fast enough

## Benchmarks vs other libraries

Tested on 100M random bits:

| Library | rank (ns) | select (ns) | space overhead |
|---------|-----------|-------------|----------------|
| This | 52 | 382 | 37.5% |
| sdsl-lite | 45 | 410 | 25% |
| folly::BitVector | 78 | n/a | 50% |

## TODO / Future work

- [ ] Compressed representations for sparse/dense vectors (RRR encoding)
- [ ] SIMD batch operations
- [ ] Better select0 (currently uses binary search)
- [ ] Memory-mapped file support

## Why?

Needed this for a text indexing project. Existing libraries were either too heavy (sdsl-lite pulls in a ton of stuff) or didn't support select (folly).

Plus it's a fun bit-twiddling exercise.

## References

- Jacobson, G. (1989). Space-efficient static trees and graphs
- Raman, Raman, and Rao (2002). Succinct indexable dictionaries

## License

MIT