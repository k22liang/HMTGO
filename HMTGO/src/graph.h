
#ifndef HMM_GRAPH_H
#define HMM_GRAPH_H

#include "node.h"
#include <vector>
#include <unordered_set>
#include <forward_list>
#include <map>

namespace hmm {

class Graph
{
private:
    void create_nodes(std::string const & root,
            std::unordered_map<std::string, std::set<std::string>> const & d_children) 
    {
        std::unordered_set<std::string> nodes_created;
        create_nodes_1(root, d_children, nodes_created);
    }

    void create_nodes_1(std::string const & node,
            std::unordered_map<std::string, std::set<std::string>> const & d_children, 
            std::unordered_set<std::string> & nodes_created)
    {
        if (node.substr(0, 3) != "GO:") {
            throw std::runtime_error("Invalid node name: " + node + "\n");
        }
        nodes.insert(new Node);
        nodes_created.insert(node);
        if(d_children.count(node)){
            std::set<std::string> const & children = d_children.at(node);
            for(auto it=children.cbegin(); it!=children.cend(); ++it){
                if(nodes_created.count(*it)==0){
                    hmm::Graph::create_nodes_1(*it, d_children, nodes_created);
                }
            }
        }
    }
    /*
     * Set up a node all its offsprings.
     * This is a recursive helper function for Graph::build.
     */
    void build_graph_1(std::set<hmm::Node*>::const_iterator & it_node,
            std::unordered_map<std::string, std::set<std::string>> const & d_children, 
            std::unordered_map<std::string, std::set<int>> const & d_gene_index, 
            std::unordered_map<std::string, hmm::Node*> & nodes_created, int log)
    {
        hmm::Node* node = *it_node;  
        nodes_created[node->get_name()] = node;
        // assign genes to node
        if(d_gene_index.count(node->get_name())){
            node->set_genes(d_gene_index.at(node->get_name()));
        }
        // add relationship
        if(d_children.count(node->get_name())){
            std::set<std::string> const & child_names = d_children.at(node->get_name());
            for(auto it=child_names.cbegin(); it!=child_names.cend(); ++it){
                if(nodes_created.count(*it)){
					// an existing node
                    hmm::Node* child = nodes_created[*it];
                    node->add_child(child, log);
                    child->add_parent(node, log);
                }else{
					// add a new node
                    ++it_node;
                    hmm::Node* child = *it_node;
                    child->set_name(*it);
                    node->add_child(child, log);
                    child->add_parent(node, log);
					// recursive calling, depth first
                    hmm::Graph::build_graph_1(it_node, d_children, d_gene_index, nodes_created, log);
                }
            }
        }
    }

    /*
     * Convert the database of genes in the format of unordered_map<string, std::set<string>> 
     * to the format of unordered_map<string, std::set<int>>.
     * All unique genes together is stored in all_genes.
     */
    std::unordered_map<std::string, std::set<int>> 
    convert_gene_database(std::unordered_map<std::string, std::set<std::string >> const & d_genes)
    {
        // merge all genes 
        std::set<std::string> genes;
        for(auto it=d_genes.cbegin(); it!=d_genes.cend(); ++it){
            std::set<std::string> const & g = it->second;
            genes.insert(g.cbegin(), g.cend());
        }
        all_genes.insert(all_genes.begin(), genes.cbegin(), genes.cend());
        // map each gene to its index
        std::unordered_map<std::string, int> index;
        for(int i=all_genes.size()-1; i>=0; --i){
            index[all_genes[i]] = i;
        }
        // create a gene database using genes' indices
        std::unordered_map<std::string, std::set<int>> d_gene_index;
        for(auto it=d_genes.cbegin(); it!=d_genes.cend(); ++it){
            d_gene_index[it->first] = gene_index(index, it->second);
        }
        return d_gene_index;
    }

    /*
     * Find the set of indices (integer) corresponding to a set of genes (string).
     */
    std::set<int> gene_index(std::unordered_map<std::string, int> const & index, 
            std::set<std::string> const & genes)
    {
        std::set<int> gi;
        for(auto it=genes.cbegin(); it!=genes.cend(); ++it){
            gi.insert(index.at(*it));
        }
        return gi;
    }
    
    //void fill_parents(hmm::Node* node);

    /*
     * Check whether the "root" is indeed root.
     *     a. it has no parent
     *     b. all other node have at least 1 parent
     *
     * @param active_only If true, only active parents are counted while checking the graph;
     * otherwise, all parents are counted.
     */
    void check_root() const
    {
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            int const np = (*it)->get_parents().size();
            if((*it)->get_name() == root->get_name()){
                if(np > 0){
                    throw std::runtime_error((*it)->summary() + " is root but has "  + std::to_string(np) 
                            + " parent(s)!\nThis typically means that the children database is broken.\n");
                }
            }else if(np == 0){
                throw std::runtime_error((*it)->summary() + " is not root but has no parents!\nThis typically means that the children database is broken.\n");
            }
        }
    }

    /*
     * For each node, check whether its probes is a super set of the probes of its children.
     */
    void check_probes() const
    {
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            std::set<int> const & genes = (*it)->get_genes();
            std::set<hmm::Node*> children = (*it)->get_children();
            for(auto cit=children.cbegin(); cit!=children.cend(); ++cit){
                std::set<int> const & child_genes = (*cit)->get_genes();
                if(!(*it)->contains_node(*cit)){
                    throw std::runtime_error("Genes associated with " + (*it)->summary() 
                            + " is not a super set of the genes associated with its child " 
                            + (*cit)->summary() + ".\nThis typically means that the children database is broken.\n\n" 
                            + (*it)->summary() + ":\n" + hmm::to_string(genes) + "\n\n" 
                            + (*cit)->summary() + ":\n" + hmm::to_string(child_genes) 
                            + "\n\nGenes in " + (*cit)->summary() + " but not in " + (*it)->summary() + ":\n" 
                            + hmm::to_string(hmm::genes_difference(child_genes, genes)) + "\n\n");
                }
            }
        }
    }

    /**
     * Summarize the graph.
     *
     * @param active_only If true, only active children/parents of each node are counted;
     * otherwise, all children/parents of each node are counted.
     */
    void summarize_graph(int log) const
    {
        // number of parents for each node
        std::vector<int> np(size());
        // number of children for each node
        std::vector<int> nc(size());
        // size of each node
        std::vector<int> nsize(size());
        // number of nodes having multiple parents
        int nknot = 0;
        if(log>=1){
            int i = 0;
            for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it, ++i){
                np[i] = (*it)->get_parents().size();
                nc[i] = (*it)->get_children().size();
                nsize[i] = (*it)->get_genes().size();
                if(np[i] > 1){
                    ++nknot;
                }
            }
            //-----------------------------------------------------------------------
            int max_np = *std::max_element(np.cbegin(), np.cend());
            int max_nc = *std::max_element(nc.cbegin(), nc.cend());
            int max_size = *std::max_element(nsize.cbegin(), nsize.cend());
            int min_size = *std::min_element(nsize.cbegin(), nsize.cend());
            // report statistics about the graph
            std::cout << "Total number of nodes: " << size() << "\n";
            std::cout << "Max number of parents: " << max_np << "\n";
            std::cout << "Max number of children: " << max_nc << "\n";
            std::cout << "Max number of genes: " << max_size << "\n";
            std::cout << "Min number of genes: " << min_size << "\n";
            std::cout << "Number of knots: " << nknot << "\n";
        }
        if(log>=3){
            // size: freq
            std::map<int, int> freq;
            for(int i=size()-1; i>=0; --i){
                ++freq[nsize[i]];
            }
            // freq: list of sizes
            std::map<int, std::set<int>> freq_compact;
            for(auto it=freq.cbegin(); it!=freq.cend(); ++it){
                freq_compact[it->second].insert(it->first);
            }
            std::cout << "Size table:\n";
            for(auto it=freq_compact.cbegin(); it!=freq_compact.cend(); ++it){
                std::set<int> const & sizes = it->second;
                for(auto jt=sizes.cbegin(); jt!=sizes.cend(); ++jt){
                    std::cout << *jt << " ";
                }
                std::cout << ": " << it->first << "\n";
            }
            std::cout << "\n";
        }
    }

    /*
     * Verify that if node A has B as one of its children, 
     * then B has A as one of its parents.
     */
    void check_pair() const
    {
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            std::set<hmm::Node*> const & children = (*it)->get_children();
            for(auto cit=children.cbegin(); cit!=children.cend(); ++cit){
                if(!(*cit)->get_parents().count(*it)){
                    throw std::runtime_error((*it)->summary() + " has " + (*cit)->summary() 
                            + " as one of its children but " + (*cit)->summary() + " does not have " 
                            + (*it)->summary() + " as one of its parents.");
                }
            }
            std::set<hmm::Node*> const & parents = (*it)->get_parents();
            for(auto pit=parents.cbegin(); pit!=parents.cend(); ++pit){
                if(!(*pit)->get_children().count(*it)){
                    throw std::runtime_error((*it)->summary() + " has " + (*pit)->summary() 
                            + " as one of its parents but " + (*pit)->summary() + " does not have " 
                            + (*it)->summary() + " as one of its children.");
                }
            }
        }
    }

    /*
     * Check whether there are any empty nodes in the graph.
     */
    void check_empty() const
    {
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            if(!(*it)->size()){
                throw std::runtime_error((*it)->summary() 
                        + " is empty, but there shouldn't be any empty nodes at this stage.\n");
            }
        }
    }

    void check_clones() const
    {
        std::unordered_map<int, std::forward_list<hmm::Node const *>> group_by_size;
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            group_by_size[(*it)->size()].push_front(*it);
        }
        for(auto it=group_by_size.cbegin(); it!=group_by_size.cend(); ++it){
            check_clones_1(it->second);
        }
    }

    void check_clones_1(std::forward_list<hmm::Node const *> const & nodes) const
    {
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            for(auto jt=std::next(it); jt!=nodes.cend(); ++jt){
                if((*it)->contains_node(*jt)){
                    // *it and *jt are clones
                    throw std::runtime_error((*it)->summary() + " is a clone of " + (*jt)->summary() 
                            + ", but there shouldn't be any clones in the graph at this stage.");
                }
            }
        }
    }

    void check_relationship() const
    {
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            check_relationship_1(*it);
        }    
    }

    void check_relationship_1(hmm::Node const * node) const
    {
        std::unordered_set<hmm::Node*> uid = get_larger_unidentified(node);
        for(auto it=uid.cbegin(); it!=uid.cend(); ++it){
            if((*it)->contains_node(node)){
                // *it is a parent of node
                throw std::runtime_error("\nA missing relationship between " + (*it)->summary() 
                        + " and " + node->summary() 
                        + " is found, but there shouldn't be any missing relationship in the graph at this stage.\n");
            }
        }
    }

    void check_shortcuts() const
    {
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            check_shortcuts_1(*it);
        }
    }

    void check_shortcuts_1(hmm::Node const * node) const
    {
        std::set<hmm::Node*> const & children = node->get_children();
        std::vector<hmm::Node const *> children_sorted(children.begin(), children.end());
        // sort children_sorted ascendingly according to size
        std::sort(children_sorted.begin(), children_sorted.end(),
                [](hmm::Node const * n1, hmm::Node const * n2){
                    return n1->less_than(n2);
                } ); 
        for(auto it=children_sorted.cbegin(); it!=children_sorted.cend(); ++it){
            int size = (*it)->size();
            for(auto jt=children_sorted.crbegin(); jt!=children_sorted.crend(); ++jt){
                if((*jt)->size()<=size){
                    break;
                }
                // if(parents.count(*jt)){ // problematic!!! the correct way is as follows
                if((*jt)->contains_node(*it)){ 
                    // *jt is a parent of *it, which means that node->*it is a shortcut
                    throw std::runtime_error(node->summary() + " --> " + (*it)->summary() 
                            + " is a shortcut (" + node->summary() + "). However, there shouldn't be any shortcut in the graph at this stage.\n"); 
                }
            }
        }
    }
protected:

    hmm::Node* root;
    /* 
     * All nodes in the graph. 
     * It is updated if a node is removed from the graph.
     */
    std::set<hmm::Node*> nodes;
    /*
     * All UNIQUE genes together from all nodes.
     * It acts as the pool of all genes such that genes for each node can be represented as set<int>,
     * which greatly reduces the memory usage.
     */
    std::vector<std::string> all_genes;
    // database of genes: node - genes (index)
    std::unordered_map<std::string, std::set<int>> d_gene_index;
    /*
     * Write all nodes into the specified file.
     * This is for debugging purpose.
     * 
     * @param file a file to write the nodes into.
     */
    void write_nodes(std::string const file) const
    {
        std::ofstream ofs(file);
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            ofs << (*it)->summary() << "\n";
        }
    }

    void write_clones(std::string const file) const
    {
        std::ofstream ofs(file);
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            if((*it)->has_clones()){
                ofs << (*it)->get_name() << "\t";
                std::set<std::string> const & clones = (*it)->get_clones();
                for(auto jt=clones.cbegin(); jt!=clones.cend(); ++jt){
                    ofs << *jt << "\t";
                }
                ofs << "\n";
            }
        }
    }

    /**
     * Check whether the graph is valid.
     * To make sure your operations on the graph are correct,
     * you can run this function after the graph is changed.
     * The following conditions are checked.
     *
     * 1. the "root" node is indeed a root
     *     a. it has no parent
     *     b. each other node has at least 1 parent
     * 2. a node must indeed a parent for all its children
     *
     * Whenever a condition above is violated, a runtime error exception is thrown.
     * If the graph is valid, summary statistics about it are printed.
     *
     * @param active_only If true, only active parents are counted while checking the graph;
     * otherwise, all parents are counted.
     *
     * @param log If true, then extensive log information are printed.
     */
    void validate_graph(int log) const
    {
        if(log >= 1){
            std::cout << "\nValidating the structure of graph.\n";
        }
        hmm::Graph::check_root();
        hmm::Graph::check_probes();
        hmm::Graph::check_pair();
        summarize_graph(log);
    }

    /**
     * Validate more on the graph structure.
     * 1. check for empty nodes.
     * 2.
     */
    void validate_graph_more(int log) const
    {
        if(log >= 1){
            std::cout << "\nValidating more on the structure of graph.\n";
        }
        check_empty();
        check_clones();
        check_relationship();
        check_shortcuts();
    }
    
    /*
     * Best to return a vector as you will need to sort the nodes according to their size.
     */
    std::unordered_set<hmm::Node*> get_larger_unidentified(hmm::Node const * node) const
    { 
        std::set<hmm::Node*> ancestors = node->get_ancestors(); 
        std::set<hmm::Node*> larger_nodes;
        int size = node->size();
        std::copy_if(nodes.cbegin(), nodes.cend(), 
                std::inserter(larger_nodes, larger_nodes.end()), [size](hmm::Node const * node){
            return node->size() > size;
        });
        std::unordered_set<hmm::Node*> unidentified;
        std::set_difference(larger_nodes.cbegin(), larger_nodes.cend(), 
                ancestors.cbegin(), ancestors.cend(), std::inserter(unidentified, unidentified.end()));
        return unidentified;
    } 

public:
    /*
     * Get the number of nodes/vertices in the graph.
     */
    int size() const
    {
        return nodes.size();
    }
    /*
     * Build the graph from data exported from R.
     * You should make sure that d_children and d_genes are compatible.
     * Don't do it here in your code. It makes your complicated. 
     * Do it when you prepare the data.
     */
    void build_graph(std::string const & name_root, 
            std::unordered_map<std::string, std::set<std::string>> const & d_children, 
            std::unordered_map<std::string, std::set<std::string>> const & d_genes, int log, bool debug)
    {
        // create nodes 
        create_nodes(name_root, d_children); 
        // build the graph
        std::set<hmm::Node*>::const_iterator it_node = nodes.cbegin();
        // record the root node 
        root = *it_node;
        root -> set_name(name_root);
        d_gene_index = convert_gene_database(d_genes);
        std::unordered_map<std::string, hmm::Node*> nodes_created;
        build_graph_1(it_node, d_children, d_gene_index, nodes_created, log);
        if(debug){
            validate_graph(log);
        }
    }

protected:
    virtual ~Graph() 
    {
        for(auto it=nodes.begin(); it!=nodes.end(); ++it){
            delete (*it);
        }
    }

    //void fill_parents();
    
    // int number_parents(Node const * node, bool active_only) const;

    // int number_children(Node const * node, bool active_only) const;

};

}

#endif // HMM_GRAPH_H


/*
 * Fill in the parents information for each node in the graph.
void hmm::Graph::fill_parents(){
    for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
        fill_parents(*it);
    }
}
 */
/*
 * Mark a node as a parent for each of its child.
void hmm::Graph::fill_parents(hmm::Node const * node){
    std::set<Node const *> const & children = (node->children).first;
    for(auto it=children.cbegin(); it!=children.cend(); ++it){
        ((*it)->parents).first.insert(node);
    }
}
 */

/*
 * Force each node's genes to be a subset of all its parents nodes' genes.
void force_graph(){
    reset_skip();
}

 */
/*
 * For the node's offsprings to be true parent - children relationship.
void force_graph(Node const * node){
    std::unordered_map<Node const *, bool> children = node->children;
    for(auto it=children.cbegin(); it!=children.cend(); ++it){
        force_relationship(node, it->first);
    }
}
void force_relationship(Node const * parent, Node * child){
}
 */
