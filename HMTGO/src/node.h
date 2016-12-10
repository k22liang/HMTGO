
#ifndef HMM_NODE_H
#define HMM_NODE_H

#include <string>
#include <set>
#include <unordered_map>
#include <tuple>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "algo.h"


namespace hmm {

class Node
{
public:

    bool less_than(hmm::Node const * other) const
    {
        int s1 = size();
        int s2 = other->size();
        if (s1 < s2) {
            return true;
        }
        if (s2 < s1) {
            return false;
        }
        return (this < other);
    }
    /*
     * Check whether this node contains another node,
     * i.e., whether the genes of this node is a super set of the genes of another node.
     */
    bool contains_node(hmm::Node const * node) const
    {
        return hmm::genes_include(get_genes(), node->get_genes());
    }

    /*
     * Check whether this node contains any node in nodes.
     */
    bool contains_any_node(std::set<hmm::Node*> const & nodes) const
    {
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            if(contains_node(*it)){
                return true;
            }
        }
        return false;
    }

    void add_clone(hmm::Node const * clone)
    {
        clones.insert(clone->get_name());
        clones.insert(clone->clones.cbegin(), clone->clones.cend());
    }

    void set_genes(std::set<int> const & node_genes)
    {
        genes = node_genes;
    }   

    std::set<int> const & get_genes() const
    {
        return genes;
    }

    std::set<hmm::Node*> const & get_children() const
    {
        return children;
    }
    
    // void set_children(std::set<Node*> const & _children, bool use_second, bool changed);
    
    std::set<hmm::Node*> const & get_parents() const
    {
        return parents;
    }

    // void set_parents(std::set<Node*> const & _parents, bool use_second, bool changed);

    void remove_child(Node* child, int log)
    {
        if(log >= 4){ 
            // extensive logging information
            std::cout << "Removing " << child->summary() << " from the children of " << summary() << ".\n";
        }
        children.erase(child);
    }

    void remove_parent(Node* parent, int log)
    {
        if(log >= 4){
            // extensive logging information
            std::cout << "Removing " << parent->summary() << " from the parents of " << summary() << ".\n";
        }
        parents.erase(parent);
    }

    void add_child(Node* child, int log)
    {
        if(child != this){
            children.insert(child);
            if(log >= 4){
                // extensive logging information
                std::cout << child->summary() << " is added as a child of " << summary() << ".\n";
            }
        }
    }

    void add_parent(Node* parent, int log)
    {
        if(parent != this){
            parents.insert(parent);
            if(log >= 4){
                // extensive logging information
                std::cout << parent->summary() << " is added as a parent of " << summary() << ".\n";
            }
        }
    }
    /**
     * Get the size of the node, i.e., number of genes associated with the node.
     */
    int size() const
    {
        return get_genes().size();
    }
    
    std::string summary() const
    {
        std::stringstream ss;
        ss << this;
        return name + "<" + ss.str() + ", P: " + std::to_string(get_parents().size()) + 
            ", C: " + std::to_string(get_children().size()) + ", S: " + std::to_string(size()) + ">";
    }

    std::set<hmm::Node*> get_ancestors() const
    {
        std::set<hmm::Node*> ancestors;
        get_ancestors_1(ancestors);
        // bool has_root = std::count_if(ancestors.cbegin(), ancestors.cend(), [](hmm::Node* node){
        //     return node->name=="GO:0008150";
        // }) + (name=="GO:0008150");
        // if(!has_root){
        //     throw std::runtime_error(summary()+" does not have the root as one of its ancestors.\n");
        // }
        return ancestors;
    }

    std::set<hmm::Node*> get_offsprings() const
    {
        std::set<hmm::Node*> offsprings;
        get_offsprings_1(offsprings);
        return offsprings;
    }

    bool is_knot() const
    {
        return get_parents().size() > 1;
    }
    
    void set_name(std::string node_name)
    {
        name = node_name;
    }

    std::string get_name() const
    {
        return name;
    }

    bool is_changed() const;

    bool has_parent(hmm::Node* node) const
    {
        return get_parents().count(node);
    }

    bool has_child(hmm::Node* node) const
    {
        return get_children().count(node);
    }

    bool has_clones() const
    {
        return clones.size() > 0;
    }

    std::set<std::string> const & get_clones() const
    {
        return clones;
    }

    hmm::Node* only_parent() const
    {
        std::set<hmm::Node*> const & parents = get_parents();
        int size = parents.size();
        if(size != 1){
            throw std::runtime_error("The number of parents of " + summary() + " is not 1.");
        }
        return *parents.cbegin();
    }

    void set_index(int i)
    {
        index = i;
    }

    int get_index() const
    {
        return index;
    }

    // bool linked_with(hmm::Node* node) const;

    // int get_depth(Node* node) const;

    // int get_depth(std::unordered_map<hmm::Node*, int> & depth, hmm::Node* node) const;

    /*
     * The default constructor.
     * The memeber skip is set to false by default.
     */
    /**
     * Implement the default constructor.
     */
    Node() {}
private:
   
    std::string name;

    std::set<int> genes;
  
    std::set<Node*> children;
   
    std::set<Node*> parents;
    
    std::set<std::string> clones;
    // index has 2 usages:
    // 1. uses as a sorting criteria to get deterministic results 
    // 2. index nodes when converting the tree to data structure in R
    int index;
    
    // bool skip;
   
    void get_ancestors_1(std::set<hmm::Node*> & ancestors) const
    {
        std::set<hmm::Node*> const & parents = get_parents();
        for(auto it=parents.cbegin(); it!=parents.cend(); ++it){
            if(!ancestors.count(*it)){
                ancestors.insert(*it);
                (*it)->get_ancestors_1(ancestors);
            }
        }
    }

    void get_offsprings_1(std::set<hmm::Node*> & offsprings) const
    {
        std::set<hmm::Node*> const & children = get_children();
        for(auto it=children.cbegin(); it!=children.cend(); ++it){
            if(!offsprings.count(*it)){
                offsprings.insert(*it);
                (*it)->get_offsprings_1(offsprings);
            }
        }
    }
    /**
     * Get the size and index of this node as a pair.
     */
    /**
	std::pair<int, hmm::Node*> get_size_augumented() const
    {
        std::make_pair(size(), this);
    }
	*/

};
}

#endif // HMM_NODE_H

/*
void hmm::Node::set_children(std::set<Node*> const & _children, bool use_second, bool changed)
{
    set_set(children, _children, use_second, changed);
}
*/

/*
void hmm::Node::set_parents(std::set<Node*> const & _parents, bool use_second, bool changed)
{
    set_set(parents, _parents, use_second, changed);
}
*/

/*
int hmm::Node::depth() const
{
    std::unordered_map<hmm::Node*, int> node_depths;
    return depth(node_depths);
}

int hmm::Graph::get_depth(std::unordered_map<hmm::Node*, int> & depth) const
{
    if(depth.count(node)){
        return depth[node];
    }
    std::set<hmm::Node*> const & children = node -> get_children();
    std::vector<int> child_depths(children.size());
    int i = 0;
    for(auto it=children.cbegin(); it!=children.cend(); ++it, ++i){
        child_depths[i] = get_depth(depth, *it);
    }
    int d = *std::max_element(child_depths.cbegin(), child_depths.cend()) + 1;
    depth[node] = d;
    return d;
}
*/

/*
 * Get the number of parents of a node.
 * @param active_only If true, only active parents are counted;
 * otherwise, all parents are counted.
 */
/*
int hmm::Graph::number_parents(Node const * node, bool active_only) const {
    if(active_only){
        int n = 0;
        std::unordered_map<hmm::Node*, bool> const & parents = node -> parents;
        for(auto it=parents.cbegin(); it!=parents.cend(); ++it){
            n += it -> second;
        }
        return n;
    }
    return (node -> parents).size();
}
*/

/*
 * Get the number of children of a node.
 * @param active_only If true, only active children are counted;
 * otherwise, all children are counted.
 */
/*
int hmm::Graph::number_children(Node const * node) const {
    if(active_only){
        int n = 0;
        std::unordered_map<hmm::Node*, bool> const & children = node -> children;
        for(auto it=children.cbegin(); it!=children.cend(); ++it){
            n += it -> second;
        }
        return n;
    }
    return (node -> children).size();
}
*/

// bool hmm::Node::linked_with(hmm::Node* node) const
// {
//     return has_child(node) || has_parent(node);
// }

