
#pragma once
#ifndef TREE_H
#define TREE_H

#include <functional>
#include <map>
#include <utility>

template <class N, class A>
class tree {
public:
    std::map<N, tree<N, A>> children;
    tree<N, A> *parent;
    A data;
    N identity;

    tree() : parent(nullptr), data(A{}), identity(N{}) {}
    tree(const A &data) : parent(nullptr), data(data), identity(N{}) {}
    void add_child(const N &n, const A &a) {
        children.emplace(n, tree<N, A>{a});
        children.at(n).parent   = this;
        children.at(n).identity = n;
    }

    tree<N, A> &operator[](const N &n) { return children.at(n); }

    tree<N, A> &try_get(const N &n, const A &a) {
        if (children.find(n) == children.end()) {
            add_child(n, a);
        }
        return children.at(n);
    }
    A &value(const N &n) { return children.at(n).data; }
    A value(const N &n) const { return children.at(n).data; }

    template <class Func>
    void df_execute2(int depth, Func foo) const {
        foo(depth, this);
        for (auto &p : children) {
            p.second.df_execute2(depth + 1, foo);
        }
    }

    template <class Func>
    void df_execute2(Func foo) const {
        foo(0, this);
        for (auto &p : children) {
            p.second.df_execute2(1, foo);
        }
    }
    template <class Func>
    void df_execute(int depth, A &arc_val, Func foo) {
        foo(depth, data, arc_val);
        for (auto &p : children) {
            p.second.df_execute(depth + 1, p.second.first, foo);
        }
    }
    template <class Func>
    void df_execute(Func foo) {
        foo(0, data, A{});
        for (auto &p : children) {
            p.second.df_execute(1, p.second.first, foo);
        }
    }
};

#endif
