#include "Filesystem.hpp"
#include "CLI11.hpp"
#include "Bam.hpp"

using ghc::filesystem::path;
using ghc::filesystem::exists;
using ghc::filesystem::create_directories;
using gfase::SamElement;
using gfase::Bam;

#include <unordered_map>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>

using std::unordered_map;
using std::runtime_error;
using std::ofstream;
using std::cerr;
using std::string;


void get_identity_from_bam(path bam_path, path output_dir){
    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    Bam bam_reader(bam_path);

    unordered_map<int64_t, int64_t> identity_ppm_distribution;

    bam_reader.for_alignment_in_bam(true, [&](SamElement& e){
//        cerr << e.ref_name << ' ' << e.query_name << ' ' << int(e.mapq) << ' ' << e.flag << '\n';

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
                throw runtime_error("ERROR: alignment contains ambiguous M operations, cannot determine mismatches without = or X operations");
            }
        });

        int64_t identity_ppm = (matches * 1000000) / (nonmatches + matches);

        identity_ppm_distribution[identity_ppm] += 1;
    });

    ofstream file(output_dir / "identity_distribution.csv");
    for (auto& [ppm,count]: identity_ppm_distribution){
        file << double(ppm)/10000 << ',' << count << '\n';
    }
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
