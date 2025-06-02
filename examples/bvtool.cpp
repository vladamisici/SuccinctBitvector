#include "succinct_bitvector.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <chrono>
#include <random>

using namespace succinct;

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " <command> [arguments]\n\n"
              << "Commands:\n"
              << "  build <output_file> <bit_string>     Build bitvector from bit string\n"
              << "  build <output_file> -f <input_file>  Build bitvector from file\n"
              << "  rank <input_file> <position>         Query rank1(position)\n"
              << "  rank0 <input_file> <position>        Query rank0(position)\n"
              << "  select <input_file> <k>              Query select1(k)\n"
              << "  select0 <input_file> <k>             Query select0(k)\n"
              << "  access <input_file> <position>       Get bit at position\n"
              << "  info <input_file>                    Show bitvector information\n"
              << "  dump <input_file> [start] [end]      Dump bits in range\n"
              << "  analyze <input_file>                 Analyze bitvector structure\n"
              << "  help                                 Show this help message\n\n"
              << "Examples:\n"
              << "  " << program_name << " build my_bits.bv 101100101\n"
              << "  " << program_name << " rank my_bits.bv 5\n"
              << "  " << program_name << " select my_bits.bv 2\n";
}

void build_command(const std::string& output_file, const std::string& input, bool from_file) {
    try {
        SuccinctBitvector bv;
        
        if (from_file) {
            // Read bits from file
            std::ifstream ifs(input);
            if (!ifs) {
                std::cerr << "Error: Cannot open input file: " << input << std::endl;
                return;
            }
            
            std::string bit_string;
            char c;
            while (ifs.get(c)) {
                if (c == '0' || c == '1') {
                    bit_string += c;
                }
            }
            
            std::cout << "Building bitvector from " << bit_string.size() << " bits..." << std::endl;
            bv = SuccinctBitvector(bit_string);
        } else {
            // Build from command line string
            std::cout << "Building bitvector from " << input.size() << " bits..." << std::endl;
            bv = SuccinctBitvector(input);
        }
        
        // Save to file
        bv.save(output_file);
        
        // Print statistics
        auto mem = bv.memory_usage();
        std::cout << "✓ Bitvector built successfully!\n"
                  << "  Size: " << bv.size() << " bits\n"
                  << "  Ones: " << bv.count_ones() << " (" 
                  << std::fixed << std::setprecision(2) 
                  << (100.0 * bv.count_ones() / bv.size()) << "%)\n"
                  << "  Memory: " << mem.total_bytes << " bytes ("
                  << mem.bits_per_element << " bits/element)\n"
                  << "  Saved to: " << output_file << std::endl;
                  
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void query_command(const std::string& input_file, const std::string& query_type, uint64_t arg) {
    try {
        SuccinctBitvector bv;
        bv.load(input_file);
        
        auto start = std::chrono::high_resolution_clock::now();
        uint64_t result;
        
        if (query_type == "rank") {
            result = bv.rank1(arg);
            std::cout << "rank1(" << arg << ") = " << result << std::endl;
        } else if (query_type == "rank0") {
            result = bv.rank0(arg);
            std::cout << "rank0(" << arg << ") = " << result << std::endl;
        } else if (query_type == "select") {
            result = bv.select1(arg);
            std::cout << "select1(" << arg << ") = " << result << std::endl;
        } else if (query_type == "select0") {
            result = bv.select0(arg);
            std::cout << "select0(" << arg << ") = " << result << std::endl;
        } else if (query_type == "access") {
            bool bit = bv.access(arg);
            std::cout << "access(" << arg << ") = " << (bit ? '1' : '0') << std::endl;
            return;
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        
        std::cout << "Query time: " << duration << " ns" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void info_command(const std::string& input_file) {
    try {
        SuccinctBitvector bv;
        bv.load(input_file);
        
        auto mem = bv.memory_usage();
        
        std::cout << "Bitvector Information\n"
                  << "====================\n"
                  << "File: " << input_file << "\n"
                  << "Size: " << bv.size() << " bits\n"
                  << "Ones: " << bv.count_ones() << " ("
                  << std::fixed << std::setprecision(2)
                  << (100.0 * bv.count_ones() / bv.size()) << "%)\n"
                  << "Zeros: " << bv.count_zeros() << " ("
                  << (100.0 * bv.count_zeros() / bv.size()) << "%)\n\n"
                  << "Memory Usage\n"
                  << "------------\n"
                  << "Bits storage: " << mem.bits_capacity << " bytes\n"
                  << "Superblocks: " << mem.superblock_bytes << " bytes\n"
                  << "Blocks: " << mem.block_bytes << " bytes\n"
                  << "Total: " << mem.total_bytes << " bytes\n"
                  << "Overhead: " << std::setprecision(3) 
                  << (mem.bits_per_element - 1.0) * 100.0 << "%\n"
                  << "Bits/element: " << mem.bits_per_element << std::endl;
                  
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void dump_command(const std::string& input_file, uint64_t start, uint64_t end) {
    try {
        SuccinctBitvector bv;
        bv.load(input_file);
        
        if (end == 0 || end > bv.size()) {
            end = bv.size();
        }
        if (start >= end) {
            std::cerr << "Error: Invalid range" << std::endl;
            return;
        }
        
        std::cout << "Bits [" << start << ".." << end << "):\n";
        
        const uint64_t bits_per_line = 64;
        for (uint64_t i = start; i < end; i += bits_per_line) {
            std::cout << std::setw(8) << i << ": ";
            
            for (uint64_t j = i; j < std::min(i + bits_per_line, end); ++j) {
                if (j > i && (j - i) % 8 == 0) std::cout << ' ';
                std::cout << (bv.access(j) ? '1' : '0');
            }
            
            // Show rank info at end of line
            uint64_t line_end = std::min(i + bits_per_line, end);
            std::cout << "  [rank1=" << bv.rank1(line_end) << "]";
            std::cout << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void analyze_command(const std::string& input_file) {
    try {
        SuccinctBitvector bv;
        bv.load(input_file);
        
        std::cout << "Bitvector Structure Analysis\n"
                  << "===========================\n";
        
        // Analyze bit distribution
        const uint64_t chunks = 100;
        uint64_t chunk_size = bv.size() / chunks;
        
        std::cout << "\nBit Density Distribution (100 chunks):\n";
        std::cout << "Chunk | Density | Visualization\n";
        std::cout << "------|---------|" << std::string(50, '-') << "\n";
        
        for (uint64_t i = 0; i < chunks; ++i) {
            uint64_t start = i * chunk_size;
            uint64_t end = (i == chunks - 1) ? bv.size() : (i + 1) * chunk_size;
            
            uint64_t ones = bv.rank1(end) - bv.rank1(start);
            double density = static_cast<double>(ones) / (end - start);
            
            std::cout << std::setw(5) << i << " | "
                      << std::setw(6) << std::fixed << std::setprecision(2) 
                      << density * 100 << "% | ";
            
            // Visualization bar
            int bar_length = static_cast<int>(density * 50);
            std::cout << std::string(bar_length, '█') << std::endl;
        }
        
        // Find runs of consecutive bits
        std::cout << "\nConsecutive Bit Runs Analysis:\n";
        uint64_t max_ones_run = 0, max_zeros_run = 0;
        uint64_t current_ones = 0, current_zeros = 0;
        uint64_t total_runs = 0;
        
        bool last_bit = false;
        for (uint64_t i = 0; i < bv.size(); ++i) {
            bool bit = bv.access(i);
            if (i == 0 || bit != last_bit) {
                total_runs++;
                max_ones_run = std::max(max_ones_run, current_ones);
                max_zeros_run = std::max(max_zeros_run, current_zeros);
                current_ones = current_zeros = 0;
            }
            if (bit) current_ones++;
            else current_zeros++;
            last_bit = bit;
        }
        max_ones_run = std::max(max_ones_run, current_ones);
        max_zeros_run = std::max(max_zeros_run, current_zeros);
        
        std::cout << "Total runs: " << total_runs << "\n"
                  << "Average run length: " << std::fixed << std::setprecision(2)
                  << static_cast<double>(bv.size()) / total_runs << "\n"
                  << "Longest 1-run: " << max_ones_run << "\n"
                  << "Longest 0-run: " << max_zeros_run << "\n";
        
        // Sample queries for performance indication
        std::cout << "\nSample Query Performance:\n";
        std::mt19937 gen(42);
        std::uniform_int_distribution<uint64_t> pos_dist(0, bv.size() - 1);
        std::uniform_int_distribution<uint64_t> k1_dist(0, bv.count_ones() - 1);
        
        auto benchmark_op = [](const auto& op, const std::string& name) {
            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < 1000; ++i) {
                volatile auto result = op();
            }
            auto end = std::chrono::high_resolution_clock::now();
            auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000;
            std::cout << name << ": " << ns << " ns/query\n";
        };
        
        benchmark_op([&]() { return bv.rank1(pos_dist(gen)); }, "rank1");
        benchmark_op([&]() { return bv.select1(k1_dist(gen)); }, "select1");
        benchmark_op([&]() { return bv.access(pos_dist(gen)); }, "access");
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string command = argv[1];
    
    if (command == "help" || command == "-h" || command == "--help") {
        print_usage(argv[0]);
        return 0;
    }
    
    if (command == "build") {
        if (argc == 4) {
            build_command(argv[2], argv[3], false);
        } else if (argc == 5 && std::string(argv[3]) == "-f") {
            build_command(argv[2], argv[4], true);
        } else {
            std::cerr << "Error: Invalid arguments for build command\n";
            print_usage(argv[0]);
            return 1;
        }
    } else if (command == "rank" || command == "rank0" || 
               command == "select" || command == "select0" || 
               command == "access") {
        if (argc != 4) {
            std::cerr << "Error: Invalid arguments for " << command << " command\n";
            print_usage(argv[0]);
            return 1;
        }
        query_command(argv[2], command, std::stoull(argv[3]));
    } else if (command == "info") {
        if (argc != 3) {
            std::cerr << "Error: Invalid arguments for info command\n";
            print_usage(argv[0]);
            return 1;
        }
        info_command(argv[2]);
    } else if (command == "dump") {
        if (argc < 3 || argc > 5) {
            std::cerr << "Error: Invalid arguments for dump command\n";
            print_usage(argv[0]);
            return 1;
        }
        uint64_t start = (argc >= 4) ? std::stoull(argv[3]) : 0;
        uint64_t end = (argc == 5) ? std::stoull(argv[4]) : 0;
        dump_command(argv[2], start, end);
    } else if (command == "analyze") {
        if (argc != 3) {
            std::cerr << "Error: Invalid arguments for analyze command\n";
            print_usage(argv[0]);
            return 1;
        }
        analyze_command(argv[2]);
    } else {
        std::cerr << "Error: Unknown command: " << command << std::endl;
        print_usage(argv[0]);
        return 1;
    }
    
    return 0;
}