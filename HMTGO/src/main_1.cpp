// [[Rcpp::plugins(cpp0x)]]
#include "tree.h"
#include <stdlib.h>
#include <Rcpp.h>
#include <iostream>

// [[Rcpp::export]]
Rcpp::List tree_transformation(Rcpp::List children_db, Rcpp::List genes_db, std::string root){
int log = 0;
bool debug = 0;
// get the children database
std::unordered_map<std::string, std::set<std::string>> d_children = hmm::list_to_map(children_db);
// get the genes database
std::unordered_map<std::string, std::set<std::string>> d_genes = hmm::list_to_map(genes_db);
hmm::Tree tree;
tree.build_graph(root, d_children, d_genes, log, debug);
tree.build_tree(log, debug);
return tree.r_tree();
}
