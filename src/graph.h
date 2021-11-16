#ifndef GRAPH_H
#define GRAPH_H

#include <assert.h>
#include <stdlib.h>
#include <vector>
#include <utility>
#include <map>

template <class N, class A>
class graph{
    public:
    std::vector<std::pair<N,A>> nodes;
    std::vector< std::map<size_t, A> > arcs;
    std::map< N, size_t> reverse_index;

    graph(){}

    void add( N node){
        reverse_index[node] = nodes.size();
        nodes.push_back(std::make_pair(node,A{}));
        arcs.emplace_back();
    }

    A &arc(size_t i, size_t j){
        assert(i < nodes.size());
        assert(j < nodes.size());
        return arcs[i][j];
    }
    A &arc(const N &a, const N &b){
        size_t i = reverse_index[a];
        size_t j = reverse_index[b];
        return arc(i,j);
    }
    auto begin(){
        return nodes.begin();
    }
    auto end(){
        return nodes.end();
    }
    class neighbour{
        graph &owner;
        size_t target;
        public:
        neighbour(graph &owner, size_t target) : owner(owner), target(target) {}
        class nei_iter{
            graph &owner;
            size_t target;
            decltype(owner.arcs[0].begin()) index;
            public:

            nei_iter(graph &owner, size_t target) : owner(owner), target(target), index(owner.arcs[target].begin()) {}
            nei_iter(graph &owner, size_t target, decltype(owner.arcs[0].begin()) index) : owner(owner), target(target), index(index) {}
            auto &operator *(){
                return owner.nodes[index->first].first;
            }
            nei_iter &operator++(){
                ++index;
                return *this;
            }
            nei_iter operator++(int){
                nei_iter tmp = *this;
                ++(*this);
                return tmp;
            }
            friend bool operator== (const nei_iter& a, const nei_iter& b) { return a.target == b.target && a.index == b.index; };
            friend bool operator!= (const nei_iter& a, const nei_iter& b) { return !(a==b); };
        };
        nei_iter begin(){
            return nei_iter(owner,target);
        }

        nei_iter end(){
            return nei_iter(owner,target,owner.arcs[target].end());
        }
    };
    neighbour neighbours(size_t index){
        assert(index < nodes.size());
        assert(index >= 0);
        return neighbour(*this, index);
    }
    neighbour neighbours(const N &key){
        size_t index = reverse_index[key];
        return neighbour(*this, index);
    }

    bool in( const N &key){
        return reverse_index.find(key) != reverse_index.end();
    }
    A& value(const N &key){
        size_t index = reverse_index[key];
        assert(in(key));

        return nodes[index].second;
    }
};

#endif
