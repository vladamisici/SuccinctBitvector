#include "../include/succinct_bitvector.hpp"
#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <numeric>

using namespace succinct;
using namespace std::chrono;

class Benchmark {
    std::string name_;
    std::vector<double> times_;
    high_resolution_clock::time_point start_;

public:
    explicit Benchmark(const std::string& name) : name_(name) {}

    void start() {
        start_ = high_resolution_clock::now();
    }

    void stop() {
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<nanoseconds>(end - start_).count();
        times_.push_back(duration);
    }

    void report(uint64_t ops) const {
        if (times_.empty()) return;

        // Calculate statistics
        std::vector<double> sorted_times = times_;
        std::sort(sorted_times.begin(), sorted_times.end());

        double mean = std::accumulate(sorted_times.begin(), sorted_times.end(), 0.0) / sorted_times.size();
        double median = sorted_times[sorted_times.size() / 2];
        double p95 = sorted_times[static_cast<size_t>(sorted_times.size() * 0.95)];
        double p99 = sorted_times[static_cast<size_t>(sorted_times.size() * 0.99)];

        std::cout << std::left << std::setw(25) << name_ << " | ";
        std::cout << std::right << std::setw(7) << std::fixed << std::setprecision(1)
                  << mean / ops << " ns/op | ";
        std::cout << std::setw(6) << std::setprecision(1)
                  << (1e9 * ops) / mean << " Mops/s | ";
        std::cout << "median: " << std::setw(5) << std::setprecision(1) << median / ops << " ns | ";
        std::cout << "p95: " << std::setw(5) << std::setprecision(1) << p95 / ops << " ns | ";
        std::cout << "p99: " << std::setw(5) << std::setprecision(1) << p99 / ops << " ns";
        std::cout << std::endl;
    }
};

// Benchmark configuration
struct Config {
    uint64_t bitvector_size = 100'000'000;  // 100M bits
    uint64_t num_queries = 10'000'000;      // 10M queries
    uint64_t warmup_queries = 100'000;      // 100K warmup
    double density = 0.5;                   // 50% ones
    uint64_t num_runs = 5;                  // Number of benchmark runs
};

// Generate random positions for queries
std::vector<uint64_t> generate_random_positions(uint64_t n, uint64_t count, std::mt19937& gen) {
    std::uniform_int_distribution<uint64_t> dist(0, n - 1);
    std::vector<uint64_t> positions(count);
    for (auto& pos : positions) {
        pos = dist(gen);
    }
    return positions;
}

// Generate random k values for select queries
std::vector<uint64_t> generate_random_ks(uint64_t max_k, uint64_t count, std::mt19937& gen) {
    std::uniform_int_distribution<uint64_t> dist(0, max_k - 1);
    std::vector<uint64_t> ks(count);
    for (auto& k : ks) {
        k = dist(gen);
    }
    return ks;
}

void benchmark_construction(const Config& config) {
    std::cout << "\n=== Construction Benchmark ===" << std::endl;
    std::cout << "Bitvector size: " << config.bitvector_size / 1'000'000.0 << "M bits" << std::endl;

    std::mt19937 gen(42);
    std::bernoulli_distribution bit_dist(config.density);

    // Generate random bits
    std::vector<bool> bits(config.bitvector_size);
    for (size_t i = 0; i < bits.size(); ++i) {
        bits[i] = bit_dist(gen);
    }

    Benchmark bench("Construction");

    for (uint64_t run = 0; run < config.num_runs; ++run) {
        bench.start();
        SuccinctBitvector bv(bits);
        bench.stop();

        if (run == 0) {
            auto mem = bv.memory_usage();
            std::cout << "Memory usage: " << mem.total_bytes / (1024.0 * 1024.0) << " MB" << std::endl;
            std::cout << "Bits per element: " << mem.bits_per_element << std::endl;
            std::cout << "Space overhead: " << (mem.bits_per_element - 1.0) * 100.0 << "%" << std::endl;
        }
    }

    bench.report(config.bitvector_size);
}

void benchmark_queries(const Config& config) {
    std::cout << "\n=== Query Performance Benchmark ===" << std::endl;
    std::cout << "Queries per test: " << config.num_queries / 1'000'000.0 << "M" << std::endl;

    std::mt19937 gen(12345);

    // Create bitvector with raw data for speed
    std::vector<uint64_t> raw_words(config.bitvector_size / 64 + 1);
    for (auto& word : raw_words) {
        word = std::uniform_int_distribution<uint64_t>()(gen);
    }

    SuccinctBitvector bv(raw_words.data(), config.bitvector_size);
    uint64_t ones = bv.count_ones();
    uint64_t zeros = bv.count_zeros();

    std::cout << "Density: " << (100.0 * ones) / config.bitvector_size << "%" << std::endl;

    // Generate query data
    auto rank_positions = generate_random_positions(config.bitvector_size, config.num_queries, gen);
    auto select1_ks = generate_random_ks(ones, config.num_queries, gen);
    auto select0_ks = generate_random_ks(zeros, config.num_queries, gen);
    auto access_positions = generate_random_positions(config.bitvector_size, config.num_queries, gen);

    // Warmup
    volatile uint64_t warmup_sum = 0;
    for (uint64_t i = 0; i < config.warmup_queries; ++i) {
        warmup_sum += bv.rank1(rank_positions[i % rank_positions.size()]);
    }

    std::cout << "\nOperation                 | Time/op  | Throughput | Percentiles" << std::endl;
    std::cout << "--------------------------|----------|------------|-----------------------------" << std::endl;

    // Benchmark rank1
    {
        Benchmark bench("rank1");
        volatile uint64_t checksum = 0;

        for (uint64_t run = 0; run < config.num_runs; ++run) {
            bench.start();
            for (uint64_t i = 0; i < config.num_queries; ++i) {
                checksum += bv.rank1(rank_positions[i]);
            }
            bench.stop();
        }

        bench.report(config.num_queries);
    }

    // Benchmark rank0
    {
        Benchmark bench("rank0");
        volatile uint64_t checksum = 0;

        for (uint64_t run = 0; run < config.num_runs; ++run) {
            bench.start();
            for (uint64_t i = 0; i < config.num_queries; ++i) {
                checksum += bv.rank0(rank_positions[i]);
            }
            bench.stop();
        }

        bench.report(config.num_queries);
    }

    // Benchmark select1
    {
        Benchmark bench("select1");
        volatile uint64_t checksum = 0;

        for (uint64_t run = 0; run < config.num_runs; ++run) {
            bench.start();
            for (uint64_t i = 0; i < config.num_queries; ++i) {
                checksum += bv.select1(select1_ks[i]);
            }
            bench.stop();
        }

        bench.report(config.num_queries);
    }

    // Benchmark select0
    {
        Benchmark bench("select0");
        volatile uint64_t checksum = 0;

        for (uint64_t run = 0; run < config.num_runs; ++run) {
            bench.start();
            for (uint64_t i = 0; i < config.num_queries; ++i) {
                checksum += bv.select0(select0_ks[i]);
            }
            bench.stop();
        }

        bench.report(config.num_queries);
    }

    // Benchmark access
    {
        Benchmark bench("access");
        volatile bool checksum = false;

        for (uint64_t run = 0; run < config.num_runs; ++run) {
            bench.start();
            for (uint64_t i = 0; i < config.num_queries; ++i) {
                checksum ^= bv.access(access_positions[i]);
            }
            bench.stop();
        }

        bench.report(config.num_queries);
    }
}

void benchmark_density_impact() {
    std::cout << "\n=== Density Impact on Performance ===" << std::endl;
    std::cout << "Testing different bit densities..." << std::endl;

    const uint64_t size = 10'000'000;  // 10M bits
    const uint64_t queries = 1'000'000; // 1M queries
    std::mt19937 gen(999);

    std::cout << "\nDensity | rank1    | select1  | rank0    | select0" << std::endl;
    std::cout << "--------|----------|----------|----------|----------" << std::endl;

    for (double density : {0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99}) {
        // Create bitvector with specific density
        std::bernoulli_distribution dist(density);
        std::vector<bool> bits(size);
        for (size_t i = 0; i < bits.size(); ++i) {
            bits[i] = dist(gen);
        }

        SuccinctBitvector bv(bits);
        uint64_t ones = bv.count_ones();
        uint64_t zeros = bv.count_zeros();

        // Generate queries
        auto positions = generate_random_positions(size, queries, gen);
        auto select1_ks = ones > 0 ? generate_random_ks(ones, queries, gen) : std::vector<uint64_t>();
        auto select0_ks = zeros > 0 ? generate_random_ks(zeros, queries, gen) : std::vector<uint64_t>();

        std::cout << std::fixed << std::setprecision(2) << std::setw(6) << density * 100 << "% | ";

        // Benchmark each operation
        volatile uint64_t checksum = 0;

        // rank1
        auto start = high_resolution_clock::now();
        for (uint64_t i = 0; i < queries; ++i) {
            checksum += bv.rank1(positions[i]);
        }
        auto end = high_resolution_clock::now();
        double rank1_ns = duration_cast<nanoseconds>(end - start).count() / static_cast<double>(queries);
        std::cout << std::setw(7) << std::setprecision(1) << rank1_ns << " ns | ";

        // select1
        if (ones > 0) {
            start = high_resolution_clock::now();
            for (uint64_t i = 0; i < queries; ++i) {
                checksum += bv.select1(select1_ks[i]);
            }
            end = high_resolution_clock::now();
            double select1_ns = duration_cast<nanoseconds>(end - start).count() / static_cast<double>(queries);
            std::cout << std::setw(7) << std::setprecision(1) << select1_ns << " ns | ";
        } else {
            std::cout << "    N/A | ";
        }

        // rank0
        start = high_resolution_clock::now();
        for (uint64_t i = 0; i < queries; ++i) {
            checksum += bv.rank0(positions[i]);
        }
        end = high_resolution_clock::now();
        double rank0_ns = duration_cast<nanoseconds>(end - start).count() / static_cast<double>(queries);
        std::cout << std::setw(7) << std::setprecision(1) << rank0_ns << " ns | ";

        // select0
        if (zeros > 0) {
            start = high_resolution_clock::now();
            for (uint64_t i = 0; i < queries; ++i) {
                checksum += bv.select0(select0_ks[i]);
            }
            end = high_resolution_clock::now();
            double select0_ns = duration_cast<nanoseconds>(end - start).count() / static_cast<double>(queries);
            std::cout << std::setw(7) << std::setprecision(1) << select0_ns << " ns";
        } else {
            std::cout << "    N/A";
        }

        std::cout << std::endl;
    }
}

void benchmark_scaling() {
    std::cout << "\n=== Scaling Benchmark ===" << std::endl;
    std::cout << "Testing performance vs bitvector size..." << std::endl;

    std::mt19937 gen(777);
    const uint64_t queries_per_size = 1'000'000;

    std::cout << "\nSize (M) | Memory (MB) | Overhead | rank1 (ns) | select1 (ns)" << std::endl;
    std::cout << "---------|-------------|----------|------------|-------------" << std::endl;

    for (uint64_t size_mb : {1, 10, 50, 100, 500}) {
        uint64_t size = size_mb * 1'000'000;

        // Create random bitvector
        std::vector<uint64_t> raw_words((size + 63) / 64);
        for (auto& word : raw_words) {
            word = std::uniform_int_distribution<uint64_t>()(gen);
        }

        SuccinctBitvector bv(raw_words.data(), size);
        auto mem = bv.memory_usage();
        uint64_t ones = bv.count_ones();

        // Generate queries
        auto positions = generate_random_positions(size, queries_per_size, gen);
        auto select_ks = generate_random_ks(ones, queries_per_size, gen);

        std::cout << std::setw(8) << size_mb << " | ";
        std::cout << std::setw(11) << std::fixed << std::setprecision(1)
                  << mem.total_bytes / (1024.0 * 1024.0) << " | ";
        std::cout << std::setw(7) << std::setprecision(2)
                  << (mem.bits_per_element - 1.0) * 100.0 << "% | ";

        volatile uint64_t checksum = 0;

        // Benchmark rank1
        auto start = high_resolution_clock::now();
        for (uint64_t i = 0; i < queries_per_size; ++i) {
            checksum += bv.rank1(positions[i]);
        }
        auto end = high_resolution_clock::now();
        double rank_ns = duration_cast<nanoseconds>(end - start).count() / static_cast<double>(queries_per_size);
        std::cout << std::setw(10) << std::setprecision(1) << rank_ns << " | ";

        // Benchmark select1
        start = high_resolution_clock::now();
        for (uint64_t i = 0; i < queries_per_size; ++i) {
            checksum += bv.select1(select_ks[i]);
        }
        end = high_resolution_clock::now();
        double select_ns = duration_cast<nanoseconds>(end - start).count() / static_cast<double>(queries_per_size);
        std::cout << std::setw(12) << std::setprecision(1) << select_ns << std::endl;
    }
}

void benchmark_cache_effects() {
    std::cout << "\n=== Cache Effects Benchmark ===" << std::endl;
    std::cout << "Testing sequential vs random access patterns..." << std::endl;

    const uint64_t size = 100'000'000;  // 100M bits
    const uint64_t queries = 10'000'000; // 10M queries
    std::mt19937 gen(555);

    // Create bitvector
    std::vector<uint64_t> raw_words((size + 63) / 64);
    for (auto& word : raw_words) {
        word = std::uniform_int_distribution<uint64_t>()(gen);
    }
    SuccinctBitvector bv(raw_words.data(), size);

    // Sequential positions
    std::vector<uint64_t> seq_positions(queries);
    for (uint64_t i = 0; i < queries; ++i) {
        seq_positions[i] = (i * size) / queries;
    }

    // Random positions
    auto rand_positions = generate_random_positions(size, queries, gen);

    // Stride positions (cache-unfriendly)
    std::vector<uint64_t> stride_positions(queries);
    const uint64_t stride = 4096 * 8;  // Jump by page size in bits
    for (uint64_t i = 0; i < queries; ++i) {
        stride_positions[i] = (i * stride) % size;
    }

    std::cout << "\nAccess Pattern | rank1 (ns) | Throughput (M/s)" << std::endl;
    std::cout << "---------------|------------|------------------" << std::endl;

    volatile uint64_t checksum = 0;

    // Sequential
    auto start = high_resolution_clock::now();
    for (uint64_t i = 0; i < queries; ++i) {
        checksum += bv.rank1(seq_positions[i]);
    }
    auto end = high_resolution_clock::now();
    double seq_ns = duration_cast<nanoseconds>(end - start).count() / static_cast<double>(queries);
    std::cout << "Sequential     | " << std::setw(10) << std::fixed << std::setprecision(1) << seq_ns
              << " | " << std::setw(16) << (1000.0 / seq_ns) << std::endl;

    // Random
    start = high_resolution_clock::now();
    for (uint64_t i = 0; i < queries; ++i) {
        checksum += bv.rank1(rand_positions[i]);
    }
    end = high_resolution_clock::now();
    double rand_ns = duration_cast<nanoseconds>(end - start).count() / static_cast<double>(queries);
    std::cout << "Random         | " << std::setw(10) << std::fixed << std::setprecision(1) << rand_ns
              << " | " << std::setw(16) << (1000.0 / rand_ns) << std::endl;

    // Stride
    start = high_resolution_clock::now();
    for (uint64_t i = 0; i < queries; ++i) {
        checksum += bv.rank1(stride_positions[i]);
    }
    end = high_resolution_clock::now();
    double stride_ns = duration_cast<nanoseconds>(end - start).count() / static_cast<double>(queries);
    std::cout << "Stride         | " << std::setw(10) << std::fixed << std::setprecision(1) << stride_ns
              << " | " << std::setw(16) << (1000.0 / stride_ns) << std::endl;

    std::cout << "\nCache efficiency: Sequential is " << std::setprecision(1)
              << (rand_ns / seq_ns) << "x faster than random access" << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "Succinct Bitvector Performance Benchmarks" << std::endl;
    std::cout << "=========================================" << std::endl;

    Config config;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--size" && i + 1 < argc) {
            config.bitvector_size = std::stoull(argv[++i]) * 1'000'000;
        } else if (arg == "--queries" && i + 1 < argc) {
            config.num_queries = std::stoull(argv[++i]) * 1'000'000;
        } else if (arg == "--density" && i + 1 < argc) {
            config.density = std::stod(argv[++i]);
        } else if (arg == "--runs" && i + 1 < argc) {
            config.num_runs = std::stoull(argv[++i]);
        } else if (arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]\n"
                      << "Options:\n"
                      << "  --size N      Bitvector size in millions (default: 100)\n"
                      << "  --queries N   Number of queries in millions (default: 10)\n"
                      << "  --density D   Bit density 0.0-1.0 (default: 0.5)\n"
                      << "  --runs N      Number of benchmark runs (default: 5)\n"
                      << "  --help        Show this help message\n";
            return 0;
        }
    }

    // Run benchmarks
    benchmark_construction(config);
    benchmark_queries(config);
    benchmark_density_impact();
    benchmark_scaling();
    benchmark_cache_effects();

    std::cout << "\n✅ Benchmark completed successfully!" << std::endl;
    return 0;
}