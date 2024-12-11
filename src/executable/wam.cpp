#include "Filesystem.hpp"
#include "CLI11.hpp"
#include "Bam.hpp"

using ghc::filesystem::path;
using ghc::filesystem::exists;
using ghc::filesystem::create_directories;
using gfase::SamElement;
using gfase::Bam;
using AlignmentSummary = gfase::Bam::AlignmentSummary;

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

// bool compareRecords(const Record& a, const Record& b) {
//     if (a.chromosome == b.chromosome) {
//         return a.startPosition < b.startPosition;
//     }
//     return a.chromosome < b.chromosome;
// }

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

void sort_bedgraph(const std::string& input_file, const std::string& output_file) {

    if (std::system("command -v bedtools > /dev/null 2>&1") != 0) {
        std::cerr << "Error: bedtools is not installed or not in PATH." << std::endl;
    }
    // Make the bedtools sort command
    std::string command = "bedtools sort -i " + input_file + " > " + output_file;

    // Run the command
    int ret_code = std::system(command.c_str());

    // Check the return code
    if (ret_code != 0) {
        std::cerr << "Error: bedtools sort failed with return code " << ret_code << std::endl;
    } else {
        std::cout << "Successfully sorted the BEDGraph file: " << output_file << std::endl;
    }
}

void write_sorted_alignment_summary_to_file(const unordered_map<string, AlignmentSummary>& distribution, path output_path) {
    // Vector of pairs to sort the unordered_map
    vector<pair<string, AlignmentSummary>> sorted_distribution(distribution.size());

    // size_t i = 0;
    // for (auto& [key, summary] : distribution) {
    //     sorted_distribution[i] = {key, summary};
    //     i++;
    // }

    ofstream file(output_path);
    if (!(file.is_open() && file.good())) {
        throw runtime_error("ERROR: file could not be written: " + output_path.string());
    }

    // write the header to the csv
    file << "#chr" << "\tstart_pos" << "\tend_pos" << "\tidentity" << "\tmatches" << "\tnonmatches"
                                                                << "\tinferred_len" << "\talignmentName" << "\n";
      
    // add bedGraph header 
    file << "track type=bedGraph name=\"identity\" autoScale=on\n";

                                                 

    for (auto& [key, summary] : distribution) {
        file << summary.ref_name << '\t'
             << summary.start << '\t'
             << summary.end << '\t'
             << summary.identity << '\t'
             << summary.matches << '\t'
             << summary.nonmatches << '\t'
             << summary.inferred_length << '\t'
             << key << '\n';
    }
    file.close();

    sort_bedgraph(output_path, output_path.string()+".sorted.bed");

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
    // Unordered map to store alignment summaries by unique key
    unordered_map<string, AlignmentSummary> alignment_summaries;


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
        int64_t inferred_query_length = 0;
        int64_t alignment_end = e.start_pos;

        e.for_each_cigar([&](auto type, auto length){
            if (type == '='){
                matches += length;
                inferred_query_length += length;
                alignment_end += length;
            }
            else if (type == 'X'){
                nonmatches += length;
                inferred_query_length += length;
                alignment_end += length;
            }
            else if (type == 'I'){
                nonmatches += length;
                inferred_query_length += length;
            }
            else if (type == 'D'){
                nonmatches += length;
                alignment_end += length;
            }
            else if (type == 'S' or type == 'H'){
                inferred_query_length += length;
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

        // make a unique name for each alignment and insert the summary data into the map
        string uniqueName = bam_reader.createUniqueKey(e.ref_name, e.start_pos,
                                                         alignment_end, matches, nonmatches, e.query_name);
        AlignmentSummary summary = {e.ref_name, e.start_pos, alignment_end, matches, nonmatches, inferred_query_length, identity};
        alignment_summaries[uniqueName] = summary;


    });

    write_sorted_distribution_to_file(identity_distribution, output_dir / "identity_distribution.csv");
    write_sorted_distribution_to_file(length_distribution, output_dir / "length_distribution.csv");
    write_sorted_alignment_summary_to_file(alignment_summaries, output_dir / "alignment_summary.tsv");
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
