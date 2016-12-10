
#ifndef HMM_TREE_H
#define HMM_TREE_H

#include "algo.h"
#include "graph.h"
#include <iostream>
#include <forward_list>
#include <list>
#include <queue>
#include <algorithm>
#include <unordered_map>
#ifndef NoRcpp
#include <Rcpp.h>
#endif

namespace hmm {

struct node_size_comparator {
    bool operator() (hmm::Node const * lhs, hmm::Node const * rhs) const 
    {
        return lhs->size() < rhs->size();
    }
};

class Tree : public Graph
{
private:
    long counter;

#ifndef NoRcpp
    Rcpp::List clones;
#endif

    // void debug_check_clones() const;

    /**
    * Remove node if it is a clone of one of its parents.
    */
    void clear_children_clones_1(Node* node, int log)
    {
        int size = node->size();
        // get parents
        std::set<hmm::Node*> const & parents = node->get_parents();
        for(auto it=parents.cbegin(); it!=parents.cend(); ++it){
            if((*it)->size() == size){
                // node is a clone of its parent *it
                if(log >= 2){
                    std::cout << node->summary() + " is a clone of its parent " 
                        + (*it)->summary() + " and will be removed.\n\n";
                }
                remove_clone(*it, node, log);
                return;
            }
        }
    }

    /**
     * Remove empty nodes from the graph.
     */
    void remove_empty_nodes(int log)
    {
        if(log >= 1){
            std::cout << "Removing empty nodes from the graph. If a node has no genes, it will be removed.\n"; 
        }
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ){
            auto cit = it;
            ++it;
            if((*cit)->size()==0){
                if(log >= 2){
                    std::cout << "\n" + (*cit)->summary() + " is empty and will be removed.\n";
                }
                remove_node(*cit, log);
            }
        }
    }

    /**
     * Clear children of a node which are clones of the node.
     */
    void clear_children_clones(int log)
    {
        if(log >= 1){
            std::cout << "\nClearing clones in the graph.\n" 
                << "If a node is a clone of one of its parents, then it is removed.\n\n";
        }
        // the iterator order does not matter here
        // make a copy of nodes to avoid problems
        for(auto it=nodes.begin(); it!=nodes.end(); ){
            auto jt = it;
            ++it;
            clear_children_clones_1(*jt, log);
        }
    }

    /*
     * Remove a node from the graph.
     *
     * @param node the node to be removed.
     * 
     */
    void remove_node(Node* node, int log)
    {
        if(log >= 4){
            std::cout << "Removing " + node->summary() + " from the graph.\n\n";
        }
        if(node==root){
            throw std::runtime_error(root->summary() + " is root and should never be removed.");
        }
        // parents
        std::set<hmm::Node*> const & parents = node->get_parents();
        for(auto it=parents.cbegin(); it!=parents.cend(); ++it){
            (*it)->remove_child(node, log);
            if(log >= 2){
                std::cout << "Removed " << node->summary() 
                    << " from the children of its parent " << (*it)->summary() << ".\n";
            }
        }
        // children
        std::set<hmm::Node*> const & children = node->get_children();
        for(auto it=children.cbegin(); it!=children.cend(); ++it){
            (*it)->remove_parent(node, log);
            if(log >= 2){
                std::cout << "Removed " << node->summary() 
                    << " from the parents of its its child " << (*it)->summary() << ".\n";
            }
        }
        // parents <-> children
        for(auto cit=children.cbegin(); cit!=children.cend(); ++cit){
            for(auto pit=parents.cbegin(); pit!=parents.cend(); ++pit){
                (*pit)->add_child(*cit, log);
                (*cit)->add_parent(*pit, log);
            }
        }
        // remove the node from the pool of nodes
        nodes.erase(node);
        // delete the node
        delete node;
        node = nullptr;
    }

    /*
     * Update the structure of the graph. 
     * The following actions are performed for each node in the graph.
     * 1. clearing clones
     * 2. identifying missing relationship
     * 3. removing shortcuts to children
     */
    void update_graph_structure(int log, bool debug)
    {
        if(log >= 1){
            std::cout << "Updating the structure of the graph ...\n"
                << "The following actions will be done for each node in the graph.\n" 
                << "1. removing empty nodes\n" 
                << "2. clearing clones\n"
                << "3. identifying missing relationship\n" 
                << "4. removing shortcuts to children\n";
        }
        // remove empty nodes
        remove_empty_nodes(log);
        if(debug){
            validate_graph(log);
        }
        if(debug){
            write_nodes("nodes/1/nodes_ben.txt");
        }
        // clear clones
        clear_children_clones(log);
        if(debug){
            validate_graph(log);
        }
        identify_remove_clones(log);
        if(debug){
            validate_graph(log);
        }
        if(debug){
            write_nodes("nodes/2/nodes_ben.txt");
            write_clones("nodes/2/clones_ben.txt");
        }
        // identify missing relationships
        identify_missing_parents(log, debug);
        if(debug){
            validate_graph(log);
        }
        if(debug){
            write_nodes("nodes/3/nodes_ben.txt");
            write_clones("nodes/3/clones_ben.txt");
        }
        // remove shortcuts
        remove_shortcuts_to_children(log);
        if(debug){
            validate_graph(log);
        }
        if(debug){
            write_nodes("nodes/4/nodes_ben.txt");
            write_clones("nodes/4/clones_ben.txt");
        }
    }

    // std::set<hmm::Node*> get_unidentified(Node* node) const;

    void identify_remove_clones_1(std::forward_list<hmm::Node*> const & nodes, int log)
    {
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            for(auto jt=std::next(it); jt!=nodes.cend(); ++jt){
                if((*it)->contains_node(*jt)){
                    // *it and *jt are clones, remove clone *it
                    if(log >= 2){
                        std::cout << (*it)->summary() << "is a clone of "
                            << (*jt)->summary() << " and will be removed.\n\n";
                    }
                    remove_clone(*jt, *it, log);
                    break;
                }
            }
        }
    }

    void identify_remove_clones_1(hmm::Node* node, int log)
    {
        int size = node->size();
        for(auto it=nodes.begin(); it!=nodes.end(); ){
            auto jt = it;
            ++it;
            if((*jt)->size()==size && (*jt)->contains_node(node) && *jt!=node){
                // *jt is a clone of node
                remove_clone(node, *jt, log);
            }
        }
    }

    /**
     * Further remove identify and remove clones 
     * after removing clones in children/parents.
     * Nodes are first grouped by size and then clones are searched in each group.
     */
    void identify_remove_clones(int log)
    {
        if(log >= 1){
            std::cout << "Further identify and remove clones.\n";
        }
        std::unordered_map<int, std::forward_list<hmm::Node*>> group_by_size;
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            group_by_size[(*it)->size()].push_front(*it);
        }
        for(auto it=group_by_size.cbegin(); it!=group_by_size.cend(); ++it){
            identify_remove_clones_1(it->second, log);
        }
    }

    /**
     * For each node in the graph,
     * identify the relationship between it and other nodes.
     */
    void identify_missing_parents(int log, bool debug)
    {
        counter = 0;
        if(log >= 1){
            std::cout << "\nIdentifying relationship betweeen each node and other unidentified nodes.\n";
        }
        // nodes sorted ascendingly according to size
        std::vector<hmm::Node*> nodes_sorted(nodes.cbegin(), nodes.cend());
        std::sort(nodes_sorted.begin(), nodes_sorted.end(), 
                [](hmm::Node const * n1, hmm::Node const * n2){
                    return n1->less_than(n2);
                });
        for(auto it=nodes_sorted.cbegin(); it!=nodes_sorted.cend(); ++it){
            identify_missing_parents_1(*it, log, debug);
        }
        if(log >= 1){
            std::cout << counter << " missing relationships are identified.\n";
        }
    }

    /**
     * If A is a missing child, then its offsprings are not;
     * if A.genes is not a subset of node.genes, 
     * then its ancesters are not children of node.genes.
     */
    std::forward_list<hmm::Node*> identify_missing_parents_1(hmm::Node* node, int log, bool debug)
    {
        std::forward_list<hmm::Node*> missing_parents;
        std::unordered_set<hmm::Node*> uid = get_larger_unidentified(node);
        if(log >= 2){
            std::cout << "\nIdentifying missing parents for " 
                << node->summary() << ".\n"
                << "Number of unidentified nodes: " << uid.size() << "\n";
        }
        // sort uid according to size in ascending order
        std::vector<hmm::Node*> uid_sorted(uid.begin(), uid.end());
        std::sort(uid_sorted.begin(), uid_sorted.end(),
                [](hmm::Node const * n1, hmm::Node const * n2){
                    return n1->less_than(n2);
                } ); 
        // check if node is a parent of any element in uid
        for(auto it=uid_sorted.cbegin(); it!=uid_sorted.cend(); ++it){
            if(uid.count(*it) && (*it)->contains_node(node)){
                // *it is a parent of node
                if(log >= 3){
                    std::cout << "\nA missing relationship between " 
                        << (*it)->summary() << " and " << node->summary() << " is found.\n";
                }
                // ------------------------- debugging begins -----------------------------------
                if(debug){
                    if(node->has_parent(*it) || (*it)->has_child(node)){
                        throw std::runtime_error(node->summary() 
                                + " is already in the parents of " + (*it)->summary());
                    }
                }
                // ------------------------- debugging ends -----------------------------------
                ++counter;
                (*it)->add_child(node, log);
                node->add_parent(*it, log);
                missing_parents.push_front(*it);
                // remove ancestors of *it from uid
                // @TODO: this part can be optimized as you can just loop through ancestors, 
                // there's no need to create a set of ancestors
                uid.erase(*it);
                // @TODO: hmm::erase_from?
                std::set<hmm::Node*> ancestors = (*it)->get_ancestors();
                for(auto it=ancestors.cbegin(); it!=ancestors.cend(); ++it){
                    uid.erase(*it);
                }
            }
        }
        return missing_parents;
    }

    /*
     * Check whether there is a (missing) relationship between the node parent 
     * and the offsprings of the node child. 
     */ 
    void identify_missing_children_2(hmm::Node* parent, hmm::Node const * child, 
            std::unordered_set<hmm::Node const *> & identified, int log)
    {
        if(log>=3){
            std::cout << "Identifying missing children for " << parent->summary() 
                << " in the offsprings of " << child->summary() << ".\n";
        }
        std::set<hmm::Node*> const & grandchildren = child->get_children();
        for(auto it=grandchildren.cbegin(); it!=grandchildren.cend(); ++it){
            if(identified.count(*it)){
                // offsprings of *it has already been identified
                continue;
            }
            // check if there is a (missing) relationship parent -> *it
            if(parent->contains_node(*it)){
                if(log>=4){
                    std::cout << "Missing relationship " << parent->summary() << " --> "
                        << (*it)->summary() << " identified.\n";
                }
                parent->add_child(*it, log);
                (*it)->add_parent(parent, log);
                identified.insert(*it);
                continue;
            }
            identify_missing_children_2(parent, *it, identified, log);
            identified.insert(*it);
        }
    }

    /*
     * Identify relationship between parent and these breaking nodes (including the knot).
     */
    void identify_missing_children_1(hmm::Node* parent, hmm::Node const * knot, 
            std::vector<bool> const & flag, std::vector<hmm::Node*> const & children_sorted, int log)
    {
        if(log>=2){
            std::cout << "Identifying missing relationship for " << parent->summary() << ".\n";
        }
        std::unordered_set<hmm::Node const *> identified;
        identify_missing_children_2(parent, knot, identified, log);
        int size = flag.size();
        for(int i=0; i<size; ++i){
            if(flag[i]){
                identify_missing_children_2(parent, children_sorted[i], identified, log);
            }
        }
    }

    /**
     * Remove the specified clone of a node.
     * @param node a node whose clone is to be removed.
     * @param clone a clone of node.
     */
    void remove_clone(hmm::Node* node, hmm::Node* clone, int log)
    {
        if(log >= 3){
            std::cout << "Removing " << clone->summary() 
                << " which is a clone of " << node->summary() << ".\n\n";
        }
        if (clone == node) {
            // avoid remove a node itself
            return;
        }
        // get clone's parents
        std::set<hmm::Node*> const & parents = clone->get_parents();
        for(auto it=parents.cbegin(); it!=parents.cend(); ++it){
            (*it)->add_child(node, log);
            node->add_parent(*it, log);
            (*it)->remove_child(clone, log);
        }   
        // get clone's children
        std::set<hmm::Node*> const & children = clone->get_children();
        for(auto it=children.cbegin(); it!=children.cend(); ++it){
            node->add_child(*it, log);
            (*it)->add_parent(node, log);
            (*it)->remove_parent(clone, log);
        }
        // record aliases
        node->add_clone(clone);
        // remove the clone node
        nodes.erase(clone);
        // delete clone
        delete clone;
        clone = nullptr;
    }
    
    /**
     * Remove shortcuts to children.
     * If A is a grand child of B but there exists an edge from B to A,
     * removes the edge from B to A.
     */
    void remove_shortcuts_to_children(int log)
    {
        counter = 0;
        if(log >= 1){
            std::cout << "\nRemoving shortcuts to children for each node in the graph ...\n";
        }
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            remove_shortcuts_to_children_1(*it, log);
        }
        if(log >= 1){
            std::cout << counter << " shorcuts are removed.\n";
        }
    }

    /**
     * Remove shortcuts from node to its children.
     * If A is a grand parents of B but there exists an edge from A to B,
     * then remove the edge from A to B.
     * @param node 
     */
    void remove_shortcuts_to_children_1(hmm::Node* node, int log)
    {
        if(log >= 2){
            std::cout << "\nRemoving shorcuts from " << node->summary() << " to its children.\n";
        }
        std::set<hmm::Node*> const & children = node->get_children();
        std::vector<hmm::Node*> children_sorted(children.begin(), children.end());
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
                    // *jt is an ancestor of *it, which means that node->*it is a shortcut
                    if(log >= 3){
                        std::cout << node->summary() << " --> " << (*it)->summary() << " is a shortcut ("
                            << (*jt)->summary() << ").\n"; 
                    }
                    ++counter;
                    node->remove_child(*it, log);
                    (*it)->remove_parent(node, log);
                    break;
                }
            }
        }
    }

    /**
     * Find the largest knot (a vertex with multiple parents) if any in the graph.
     * Note since a knot might not exist, 
     * you always have to check whether the returned node is a knot or not.
     *
     * @return A (potential) knot with the largest size.
     */
    hmm::Node* find_largest_knot(int log) const 
    {
        if(log >= 2){
            std::cout << "\n\nFinding the largest node with multiple parents ...\n\n";
        }
        // no need to check bound as nodes is not empty
        return *std::max_element(nodes.begin(), nodes.end(), 
                [](hmm::Node const * n1, hmm::Node const * n2){
                    return std::make_pair(n1->is_knot(), n1->size()) 
                        < std::make_pair(n2->is_knot(), n2->size());
                });
    }
    
    // std::list<hmm::Node*> path_to_remove(std::set<hmm::Node*> const & ancestors, 
            // hmm::Node * remove) const;

    /*
     * Get the path to be removed.
     */
    std::list<hmm::Node*> path_to_remove(hmm::Node* knot, hmm::Node* remove) const
    {
        std::list<hmm::Node*> path;
        std::set<hmm::Node*> const & parents = knot->get_parents();
        while(!remove->contains_any_node(parents)){
            path.push_back(remove);
            remove = remove->only_parent();
        }
        return path;
    }

    /*
     * Get the path to be removed.
     * Unlike path_to_remove, 
     * this function is for calculating the removing cost and is slightly different.
     */
    std::list<hmm::Node*> path_to_remove_cost(hmm::Node const * keep, hmm::Node* remove) const
    {
        std::list<hmm::Node*> path; 
        while(!remove->contains_node(keep)){
            path.push_back(remove);
            remove = remove->only_parent();
        }
        return path;
    }

    /*
     * Keep unique genes in "genes",
     * i.e., get rid of genes in children (excluding exclude) of parent.
     */
    void unique_genes(std::set<int> & genes, hmm::Node const * parent, hmm::Node* exclude) const
    {
        std::set<hmm::Node*> children = parent->get_children();
        for(auto it=children.cbegin(); it!=children.cend(); ++it){
            if(*it != exclude){
                hmm::erase_genes_from(genes, (*it)->get_genes());
            }
        }
    }

    int cost_remove_parent(hmm::Node* knot, hmm::Node const * keep, hmm::Node* remove) const
    {
        auto path = path_to_remove_cost(keep, remove);
        int size = path.size();
        int num_safe = 0;
        std::set<int> ug = knot->get_genes();
        hmm::Node* exclude = knot;
        for(auto it=path.cbegin(); it!=path.cend(); ++it){
            unique_genes(ug, *it, exclude);
            if(ug.size()==0){
                break;
            }
            ++num_safe;
            exclude = *it;
        }
        if(num_safe==0){
            return 3 * size;
        }
        return 2*size - num_safe;
    }

    /*
     * keep is a parent which is in parents.
     */
    int cost_remove_parents(hmm::Node* knot, hmm::Node const * keep, 
            std::set<hmm::Node*> const & parents) const
    {
        int total_cost = 0;
        for(auto it=parents.cbegin(); it!=parents.cend(); ++it){
            if(*it != keep){
                total_cost += cost_remove_parent(knot, keep, *it);
            }
        }
        return total_cost;
    }

    hmm::Node* parent_to_keep(hmm::Node* knot, std::set<hmm::Node*> const & parents) const
    {
        int min_cost = std::numeric_limits<int>::max();
        hmm::Node* keep = nullptr;
        for(auto it=parents.cbegin(); it!=parents.cend(); ++it){
            int cost = cost_remove_parents(knot, *it, parents);
            if(cost<min_cost){
                min_cost = cost;
                keep = *it;
            }
        }
        return keep;
    }

    /**
     * Remove multiparent relationship from the graph.
     */
    void remove_multiparents(int log, bool debug)
    {
        if(log >= 1){
            std::cout << "\nRemoving multiparent relationship for each node in the graph ...\n";
        }
        hmm::Node* knot = find_largest_knot(log);
        while(knot->is_knot()){
            remove_multiparents_1(knot, log, debug);
            knot = find_largest_knot(log);
        }
    }

    /**
     * Remove multiparent relationship of knot.
     *
     * @param knot a node that has multiple parents.
     *
     * @TODO: don't forget to update the parent information,
     * this can be trick as when to update ...
     */
    void remove_multiparents_1(hmm::Node* knot, int log, bool debug)
    {
        /**
         * make a copy of the parents of knots
         * strictly speaking this is not necessary but might resulting confusing and error-prone code
         */
        std::set<hmm::Node*> const & ps = knot->get_parents();
        std::vector<hmm::Node*> parents(ps.cbegin(), ps.cend());
        if(log >= 2){
            std::cout << "Removing multiparent relationship of " << knot->summary() << " with parents:\n";
            for(auto it=parents.cbegin(); it!=parents.cend(); ++it){
                std::cout << (*it)->summary() << " ";
            }
            std::cout << "\n";
        }
        // sort parents from large to small to remove parents from large to small
        std::sort(parents.begin(), parents.end(), 
            [](hmm::Node const * n1, hmm::Node const * n2){
                return n1->less_than(n2);
            });
        hmm::Node const * keep = parent_to_keep(knot, ps);
        if(log >= 2){
            std::cout << "\nKeeping " << keep->summary() 
                << " as the only parent of " << knot->summary() << ".\n"; 
        }
        for(auto it=parents.cbegin(); it!=parents.cend(); ++it){
            if(*it != keep){
                remove_multiparents_2(knot, *it, log, debug);
            }
        }
    }

    /**
     * Remove a parent-children relationship remove->knot. 
     */
    void remove_multiparents_2(hmm::Node* knot, hmm::Node* remove, int log, bool debug)
    {
        if(log >= 2){
            std::cout << "Relationship " << remove->summary() 
                << " --> " << knot->summary() << " is to be removed.\n";
        }
        // break the relatinship remove --> knot for path_to_remove to work well
        remove->remove_child(knot, log);
        knot->remove_parent(remove, log);
        // get the path to be removed
        std::list<hmm::Node*> path = path_to_remove(knot, remove);
        if(log >= 2){
            std::cout << "The path (bottom-up) to be remove is:\n";
            for(auto it=path.cbegin(); it!=path.cend(); ++it){
                std::cout << (*it)->summary() << " ";
            }   
            std::cout << "\n";
        }
        // remove the path
        for(auto it=path.cbegin(); it!=path.cend(); ++it){
            // remove relationship: *it->knot
            remove_multiparents_3(knot, *it, log, debug);
        }
    }

    /**
     * Remove a parent in a path.
     */
    void remove_multiparents_3(hmm::Node* knot, hmm::Node* parent, int log, bool debug)
    {
        if(log >= 3){
            std::cout << "\nBreaking parent-child relationship between " 
                << parent->summary() << " and " << knot->summary() << ".\n";
        } 
        // break the relationship parent --> knot so that it's easier to find breaking set of children
        parent->remove_child(knot, log);
        knot->remove_parent(parent, log);
        // sorted children of the node parent according to size (from large to small)  
        std::set<hmm::Node*> const & cs = parent->get_children();
        std::vector<hmm::Node*> children_sorted(cs.cbegin(), cs.cend());
        std::sort(children_sorted.begin(), children_sorted.end(), 
                [](hmm::Node const * n1, hmm::Node const * n2){
                    return n1->less_than(n2);
                });
        // find breaking nodes to remove
        int size = children_sorted.size();
        std::vector<bool> flag(size, false);
        for(int i=0; i<size; ++i){
            if(safe_to_remove(knot, flag, children_sorted)){
                break;
            }
            flag[i] = true;
        }
        for(int i=0; i<size; ++i){
            if(flag[i]){
                flag[i] = false;
                if(!safe_to_remove(knot, flag, children_sorted)){
                    flag[i] = true;
                }
            }else{
                break;
            }
        }
        if(log >= 3){
            std::cout << "The following nodes\n";
            for(int i=0; i<size; ++i){
                if(flag[i]){
                    std::cout << children_sorted[i]->summary() << " ";
                }
            }
            std::cout << "\ntogether with the knot " << knot->summary() 
                << " are identified as breaking nodes.\n";
        }
        // update genes of parent
        parent->set_genes(reluctant_diff(parent, knot, flag, children_sorted));
        // remove breaking set of children
        // hmm::Node const * grandpa = *parent->get_parents().cbegin();
        hmm::Node* grandpa = parent->only_parent();
        for(int i=0; i<size; ++i){
            if(flag[i]){
                // remove relationship between parent and children_sorted[i]
                parent->remove_child(children_sorted[i], log);
                children_sorted[i]->remove_parent(parent, log);
                // add relationship between grandpa and children_sorted[i]
                grandpa->add_child(children_sorted[i], log);
                children_sorted[i]->add_parent(grandpa, log);
            }
        }
        // add knot as a child of grandpa
        // stricktly speaking, this is not required to generate the right tree structure 
        // but you'd better keep it for passing debugging tests
        grandpa->add_child(knot, log);
        knot->add_parent(grandpa, log);
        // remove clones of the "parent" node
        identify_remove_clones_1(parent, log);
        // identified possible missing children for the "parent" node
        identify_missing_children_1(parent, knot, flag, children_sorted, log);
        // identified possible missing parents for the "parent" node
        std::forward_list<hmm::Node*> const missing_parents = identify_missing_parents_1(parent, log, debug);
        // remove possible shorcuts from grandpa
        remove_shortcuts_to_children_1(grandpa, log);
        // remove possible shorcuts from parents
        remove_shortcuts_to_children_1(parent, log);
        // remove possible shorcuts from missing parents of the "parent" node
        for(auto it=missing_parents.cbegin(); it!=missing_parents.cend(); ++it){
            remove_shortcuts_to_children_1(*it, log);
        }
        // validate graph
        if(debug){
            validate_graph(log);
            if(log){
                std::cout << "Debugging after remove multiparent relationship " 
                    << parent->summary() << " --> " << knot->summary() << ".\n";
            }
            validate_graph_more(log);
        }
    }

    /**
     */
    bool safe_to_remove(hmm::Node const * knot, std::vector<bool> const & flag, 
            std::vector<hmm::Node*> const & children_sorted) const
    {
        std::set<int> other_genes;
        int size = children_sorted.size();
        for(int i=0; i<size; ++i){
            if(!flag[i]){
                std::set<int> const & genes = children_sorted[i]->get_genes();
                other_genes.insert(genes.cbegin(), genes.cend());
            }
        }
        if(hmm::genes_include(other_genes, knot->get_genes())){
            return false;
        }
        for(int i = 0; i<size; ++i){
            if(flag[i]){
                if(hmm::genes_include(other_genes, children_sorted[i]->get_genes())){
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * It seems to me that this can be optimized a little bit. 
     * You can operate on the parent node directly.
     */
    std::set<int> reluctant_diff(hmm::Node const * parent, hmm::Node const * knot, 
            std::vector<bool> const & flag, std::vector<hmm::Node*> const & children_sorted) const
    {
        std::set<int> rdiff = parent->get_genes();
        hmm::erase_genes_from(rdiff, knot->get_genes());
        int size = flag.size();
        // remove genes of children to be removed
        for(int i=0; i<size; ++i){
            if(flag[i]){
                hmm::erase_genes_from(rdiff, children_sorted[i]->get_genes());
            }
        }
        // merge genes of keeping children
        for(int i=0; i<size; ++i){
            if(!flag[i]){
                std::set<int> const & genes = children_sorted[i]->get_genes();
                rdiff.insert(genes.cbegin(), genes.cend());
            }
        }
        return rdiff;
    }

    void index_nodes(int base)
    {
        root->set_index(base);
        ++base;
        index_nodes_1(root, base);
    }

    void index_nodes_1(hmm::Node * const node, int & base)
    {
        // get children and index them
        std::set<hmm::Node*> const & children = node->get_children();
        for(auto it=children.cbegin(); it!=children.cend(); ++it){
            (*it)->set_index(base);
            ++base;
        }
        // recur
        for(auto it=children.cbegin(); it!=children.cend(); ++it){
            index_nodes_1(*it, base);
        }
    }

    int max_number_children() const
    {
        int max = 0;
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            int number_children = (*it)->get_children().size();
            if(number_children > max){
                max = number_children;
            }
        }
        return max;
    }

    int max_number_genes() const
    {
        int max = 0;
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            int number_genes = (*it)->get_genes().size();
            if(number_genes > max){
                max = number_genes;
            }
        }
        return max;
    }

#ifndef NoRcpp
    void record_clones() 
    {
        clones = Rcpp::List(size());
        int i = 0;
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            std::set<std::string> cs = (*it)->get_clones();
            cs.insert((*it)->get_name());
            clones[i] = Rcpp::CharacterVector(cs.begin(), cs.end());
            ++i;
        }
    }

    Rcpp::List get_mapping() const
    {
        int size = clones.size();
        Rcpp::List mapping(size);
        for(int i=0; i<size; ++i){
            Rcpp::CharacterVector cs = clones[i];
            std::list<int> mapping_it = get_mapping_1(Rcpp::as<std::string>(cs[0]));
            mapping[i] = Rcpp::NumericVector(mapping_it.begin(), mapping_it.end());
        }
        return mapping;
    }

    std::list<int> get_mapping_1(std::string const & node) const
    {
        // a priority queue containing nodes to be checked
        std::priority_queue<hmm::Node*, std::vector<hmm::Node*>, hmm::node_size_comparator> nodes;
        nodes.push(root);
        std::set<int> const & genes_node = d_gene_index.at(node);
        std::list<int> mapping;
        get_mapping_2(nodes, genes_node, mapping);
        return mapping;
    }

    void get_mapping_2(std::priority_queue<hmm::Node*, std::vector<hmm::Node*>, hmm::node_size_comparator> & nodes,
            std::set<int> const & genes, std::list<int> & mapping) const
    {
        int size = genes.size();
        while(nodes.size()){
            // check if the genes of the top element is a subset of genes
            hmm::Node* top = nodes.top();
            nodes.pop();
            std::set<int> const & genes_top = top->get_genes();
            if(hmm::genes_include(genes, genes_top)){
                mapping.push_back(top->get_index()+1);
                // no need to continue if a node with identical genes is found
                if(size==genes_top.size()){
                    return;
                }
                continue;
            }
            // push the children of top into nodes
            std::set<hmm::Node*> const & children = top->get_children();
            for(auto it=children.cbegin(); it!=children.cend(); ++it){
                nodes.push(*it);
            }
        }
    }

    Rcpp::NumericMatrix get_children_index() const 
    {
        Rcpp::NumericMatrix children_index(max_number_children(), size());
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            int column = (*it)->get_index();
            std::set<hmm::Node*> const & children = (*it)->get_children();
            int row = 0;
            for(auto jt=children.cbegin(); jt!=children.cend(); ++jt){
                children_index(row, column) = (*jt)->get_index() + 1;
                ++row;
            }
        }
        return children_index;
    }

    Rcpp::NumericMatrix get_genes_index() const
    {
        Rcpp::NumericMatrix genes_index(max_number_genes(), size());
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            int column = (*it)->get_index();
            std::set<int> const & genes = (*it)->get_genes();
            int row = 0;
            for(auto jt=genes.cbegin(); jt!=genes.cend(); ++jt){
                genes_index(row, column) = *jt + 1;
                ++row;
            }
        }
        return genes_index;
    }

    Rcpp::NumericVector get_parents_index() const
    {
        Rcpp::NumericVector parents_index(size());
        for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
            std::set<hmm::Node*> const & parents = (*it)->get_parents();
            if(parents.size()){
                parents_index[(*it)->get_index()] = (*parents.cbegin())->get_index() + 1;
            }
        }
        return parents_index;
    }
#endif

public:
#ifndef NoRcpp
    Rcpp::List r_tree() 
    {
        index_nodes(0);
        return Rcpp::List::create(
                Rcpp::Named("all_genes") = Rcpp::wrap(all_genes),
                Rcpp::Named("genes_index") = get_genes_index(),
                Rcpp::Named("children_index") = get_children_index(),
                Rcpp::Named("parent_index") = get_parents_index(),
                Rcpp::Named("clones") = clones,
                Rcpp::Named("mapping") = get_mapping()
            );
    }
#endif
    /**
     */
    void build_tree(int log, bool debug)
    {
        if(log >= 1){
            std::cout << "\nTransforming the graph to a tree ...\n\n";
        }
        update_graph_structure(log, debug);
#ifndef NoRcpp
        record_clones();
#endif
        if(debug){
            validate_graph_more(log);
        }
        // transform the graph to a tree
        remove_multiparents(log, debug);
    }
    // std::set<hmm::Node*> root_offsprings() const; 

    // std::set<hmm::Node*> root_unidentified() const;
};

}

#endif // HMM_TREE_H

// std::list<hmm::Node const *> hmm::Tree::path_to_remove_2(std::set<hmm::Node const *> const & ancestors, 
//         hmm::Node const * remove) const
// {
//     std::list<hmm::Node const *> path;
//     while(!ancestors.count(remove)){
//         path.push_back(remove);
//         remove = remove->only_parent();
//     }
//     return path;
// }
//

/**
 * Remove the offsprings of node from nodes.
 * @param nodes
 * @param node
 */
// void remove_offsprings_from(std::set<hmm::Node*> & nodes, hmm::Node* node){
//     
// }

/*
 * Testing whether the get_offsprings method works well.
std::set<hmm::Node*> hmm::Tree::root_offsprings() const 
{
    return get_offsprings(root);
}
*/
/*
std::set<hmm::Node*> hmm::Tree::root_unidentified() const
{
    return hmm::Tree::get_smaller_unidentified(root);
}
*/

// void hmm::Tree::debug_check_clones() const
// {
//     for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
//         int size = (*it)->size();
//         // get parents
//         std::set<hmm::Node*> const & parents = (*it)->get_parents();
//         for(auto pit=parents.cbegin(); pit!=parents.cend(); ++pit){
//             if((*pit)->size()==size){
//                 throw std::runtime_error((*it)->summary() + " is a clone of its parent " 
//                 + (*pit)->summary() + ", but it was not removed.");
//             }
//         }
//     }
// }

/*
 * Reset the skip indicator (to false) for all nodes.
void reset_skip(){
    for(auto it=nodes.cbegin(); it!=nodes.cend(); ++it){
        (*it)->skip = false;
    }
}

 */
