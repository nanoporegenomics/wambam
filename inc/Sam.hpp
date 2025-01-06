#pragma once

#include "Filesystem.hpp"

using ghc::filesystem::path;

#include <functional>
#include <iostream>
#include <bitset>
#include <string>
#include <vector>

using std::function;
using std::ostream;
using std::bitset;
using std::string;
using std::vector;


namespace gfase {


class SamElement {
public:
    string query_name;
    string ref_name;
    vector<uint32_t> cigars;
    int32_t query_length;
    uint16_t flag;
    uint8_t mapq;
    int32_t start_pos;

    SamElement();
    SamElement(string& read_name, string& ref_name, uint16_t flag, uint8_t mapq, int32_t start_pos);
    bool is_first_mate() const;
    bool is_second_mate() const;
    bool is_not_primary() const;
    bool is_primary() const;
    bool is_supplementary() const;
    void for_each_cigar(const function<void(char type, uint32_t length)>& f) const;
};


void for_element_in_sam_file(path sam_path, const function<void(SamElement& e)>& f);


}

ostream& operator<<(ostream& o, const gfase::SamElement& a);

