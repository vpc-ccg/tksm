
#ifndef TREE_H
#define TREE_H

#include <functional>
#include <utility>
#include <map>

template <class N, class A>
class tree{

    public:
    std::map<N,std::pair<A,tree<N,A>>> children;
    tree<N,A> *parent;
    N data;


    tree() : parent(nullptr), data(N{}){}
    tree(const N &data) :parent(nullptr), data(data){}
    void add_child( const N &n, const A &a){
        children.emplace(n, std::make_pair(a,tree<N,A>{n}));
        children.at(n).second.parent = this;
    }

    tree<N,A> &operator[](const N &n){
        return children.at(n).second;
    }

    tree<N,A> &try_get(const N &n, const A &a){
        if(children.find(n) == children.end()){
            add_child(n,a);
        }
        return children.at(n).second;
    }
    A& value(const N &n) {
        return children.at(n).first;
    }
    A value(const N &n) const {
        return children.at(n).first;
    }


    template< class Func>
    void df_execute2(int depth, Func foo) const{
        foo(depth, this);
        for( auto &p : children){
            p.second.second.df_execute2(depth+1,foo);
        }
    }

    template< class Func>
    void df_execute2(Func foo)const {
        foo(0,this);
        for( auto &p : children){
            p.second.second.df_execute2(1,foo);
        }
    }
    template< class Func>
    void df_execute(int depth, A& arc_val, Func foo){
        foo(depth,data, arc_val);
        for( auto &p : children){
            p.second.second.df_execute(depth+1,p.second.first,foo);
        }
    }
    template< class Func>
    void df_execute(Func foo){
        foo(0,data,A{});
        for( auto &p : children){
            p.second.second.df_execute(1,p.second.first,foo);
        }
    }
};

#endif
