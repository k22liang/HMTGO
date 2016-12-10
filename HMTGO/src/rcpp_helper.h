
#ifndef HMM_RCPP_HELPER_H
#define HMM_RCPP_HELPER_H

#include <string>
#include <set>
#include <unordered_map>
#include <vector>
#include <utility>
#include <fstream>

/**
 * Write a string to an output file stream.
*/
void write_binary_string(std::ofstream & ofs, std::string const & s){
    int size = s.size();
    ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
    ofs.write(s.data(), sizeof(char) * size);
}

/**
 * Write a set of strings into the output file stream.
 */
void write_binary_strset(std::ofstream & ofs, std::set<std::string> const & strset){
    for(std::set<std::string>::const_iterator it=strset.cbegin(); it!=strset.cend(); ++it){
        write_binary_string(ofs, *it);
    }
    int size = -1;
    ofs.write(reinterpret_cast<char*>(&size), sizeof(size));
}

/**
 * Convert a Rcpp::List to an unordered map (hash table).
*/
std::unordered_map< std::string , std::set< std::string > > list2map(Rcpp::List & list){
    // get names of the list
    std::vector<std::string> names(Rcpp::as< std::vector<std::string> >(list.names()));
    std::unordered_map< std::string , std::set< std::string > > m;
    int size = list.size();
    for(int i = 0; i < size; ++i){
        std::vector<std::string> content = Rcpp::as< std::vector<std::string> >(list[i]);
        m.insert(std::make_pair(names[i], std::set<std::string>(content.cbegin(), content.cend())));
    }
    return m;
}

void write_binary_map(char const * file, std::unordered_map< std::string, std::set<std::string> > const & m){
    std::ofstream ofs(file, std::ios::binary|std::ios::out);
    for(std::unordered_map< std::string, std::set<std::string> >::const_iterator it = m.cbegin(); it != m.cend(); ++it){
        write_binary_string(ofs, it -> first);
        write_binary_strset(ofs, it -> second);
    }
}

#endif // HMM_RCPP_HELPER_H
