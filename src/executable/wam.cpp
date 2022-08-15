#include "Filesystem.hpp"
#include "CLI11.hpp"
#include "Bam.hpp"

using ghc::filesystem::path;
using ghc::filesystem::exists;
using ghc::filesystem::create_directories;
using gfase::SamElement;
using gfase::Bam;

#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <utility>
#include <string>

using std::unordered_map;
using std::sort;
using std::runtime_error;
using std::ofstream;
using std::cerr;
using std::pair;
using std::string;


template<class T1, class T2> void write_sorted_distribution_to_file(const unordered_map<T1,T2>& distribution, path output_path){
    vector <pair <T1, T2> > sorted_distribution(distribution.size());

    size_t i=0;
    for (auto& [key,count]: distribution){
        sorted_distribution[i] = {key, count};
        i++;
    }

    sort(sorted_distribution.begin(), sorted_distribution.end(),[](pair<T1,T2>& a, pair<T1,T2>& b){
        return a.first < b.first;
    });

    ofstream file(output_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: file could not be written: " + output_path.string());
    }

    for (auto& [key,count]: sorted_distribution){
        file << key << ',' << count << '\n';
    }
}


void get_identity_from_bam(path bam_path, path output_dir){
    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    Bam bam_reader(bam_path);

    unordered_map<double, int64_t> identity_distribution;
    unordered_map<size_t, int64_t> length_distribution;

    bam_reader.for_alignment_in_bam(true, [&](SamElement& e){
//        cerr << e.ref_name << ' ' << e.query_name << ' ' << int(e.mapq) << ' ' << e.flag << '\n';
        if (e.is_not_primary()){
            return;
        }

        if (not e.is_supplementary()){
            length_distribution[e.query_length]++;
        }

        if (e.mapq < 1){
            return;
        }

        int64_t matches = 0;
        int64_t nonmatches = 0;

        e.for_each_cigar([&](auto type, auto length){
            if (type == '='){
                matches += length;
            }
            else if (type == 'X' or type == 'I' or type == 'D'){
                nonmatches += length;
            }
            else if (type == 'M'){
                throw runtime_error("ERROR: alignment contains ambiguous M operations, cannot determine mismatches "
                                    "without = or X operations");
            }
        });

        double numerator = double(matches);
        double denominator = double(nonmatches) + double(matches);

        double identity = 0;
        if (denominator > 0) {
            // 7 decimals of precision is probably enough?
            identity = round(10000000*numerator / denominator)/10000000;
        }

        identity_distribution[identity]++;
    });

    write_sorted_distribution_to_file(identity_distribution, output_dir / "identity_distribution.csv");
    write_sorted_distribution_to_file(length_distribution, output_dir / "length_distribution.csv");
}


int main (int argc, char* argv[]){
    path bam_path;
    path output_dir;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input_bam",
            bam_path,
            "Path to BAM")
            ->required();

    app.add_option(
            "-o,--output_dir",
            output_dir,
            "Path to directory which will be created for output (must not exist already)")
            ->required();

    CLI11_PARSE(app, argc, argv);

    get_identity_from_bam(bam_path, output_dir);

    return 0;
}
