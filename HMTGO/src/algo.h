
#ifndef HMM_ALGO_H
#define HMM_ALGO_H

#include <fstream>
#include <deque>
#include <string>
#include <set>
#include <unordered_map>
#include <algorithm>
#ifndef NoRcpp 
#include <Rcpp.h>
#endif

namespace hmm {
    /*
     * Convert a set of integers to 1 line of string.
     */
    std::string to_string(std::set<int> const & genes){
        std::string s;  
        for(std::set<int>::const_iterator it=genes.cbegin(); it!=genes.cend(); ++it){
            s += std::to_string(*it) + " ";
        }
        return s;
    }
    /*
     * Check whether a set<int> s1 contains another set<int> s2.
     * This function tries to be smart and check the size first.
     */
    bool genes_include(std::set<int> const & gs1, std::set<int> const & gs2)
    {
        return gs1.size()>=gs2.size() && std::includes(gs1.cbegin(), gs1.cend(), gs2.cbegin(), gs2.cend());
    }

    /*
     * Calculate the difference between 2 sets s1 and s2 and returns a set.
    */
    std::set<int> genes_difference(std::set<int> const & gs1, std::set<int> const & gs2)
    {
        std::set<int> diff;
        std::set_difference(gs1.cbegin(), gs1.cend(), gs2.cbegin(), gs2.cend(), 
                std::inserter(diff, diff.end()));
        return diff;
    }

    void erase_genes_from(std::set<int> & gs1, std::set<int> const & gs2)
    {
        for(auto it=gs2.cbegin(); it!=gs2.cend(); ++it){
            gs1.erase(*it);
        }
    }
#ifndef NoRcpp 
    /**
     * Convert a Rcpp::List to an unordered map (hash table).
    */
    std::unordered_map<std::string, std::set<std::string>> list_to_map(Rcpp::List & list){
        // get names of the list
        std::vector<std::string> names(Rcpp::as<std::vector<std::string> >(list.names()));
        std::unordered_map<std::string, std::set<std::string>> m;
        int size = list.size();
        for(int i=0; i<size; ++i){
            std::vector<std::string> content = Rcpp::as<std::vector<std::string>>(list[i]);
            m.insert(std::make_pair(names[i], std::set<std::string>(content.cbegin(), content.cend())));
        }
        return m;
    }
#endif
	/**
	 * Read in map from a binary file.
	 */
	std::unordered_map<std::string, std::set<std::string>> read_binary_map(const char * file){
		std::ifstream ifs(file, std::ios::binary | std::ios::in);
		std::unordered_map< std::string, std::set<std::string> >  m;
		std::deque<std::string> contents;
		std::string temp;
		int size;
		while(ifs){
		    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
		    if(size < 0){
		        m.insert(std::make_pair(contents.front(), std::set<std::string>(++contents.begin(), contents.end())));
		        contents.clear();
		        ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
		    }
		    if(!ifs || size < 0){
		        break;
		    }
		    temp.resize(size);
		    ifs.read(&temp[0], sizeof(char) * size);
		    contents.push_back(temp);
		}
		return m;
	}
}

#endif // HMM_ALGO_H
