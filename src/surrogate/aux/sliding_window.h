#pragma once

#include <vector>

template <class T>
struct SlidingWindow {
    using size_type = typename std::vector<T>::size_type;
    using difference_type = typename std::vector<T>::difference_type;
    using value_type = T;

private:
    std::vector<T> data;
    int curr_begin;

public:
    template <class Itr, class ValT>
    class iterator {
        Itr begin_itr;
        Itr end_itr;
        Itr curr;
        bool done;

    public:
        typedef iterator self_type;
        typedef T& reference;
        typedef ValT value_type;
        typedef std::forward_iterator_tag iterator_category;
        typedef int difference_type;

        iterator(Itr begin, Itr end, int curr_pos, bool finish)
            : begin_itr(begin)
            , end_itr(end)
            , curr(begin + curr_pos)
            , done(finish)
        {
        }

        self_type begin() { return begin_itr; }
        self_type end() { return end_itr; }

        self_type operator++()
        {
            curr++;
            if (curr == end_itr) {
                curr = begin_itr;
                done = true;
            }
            return *this;
        }

        self_type operator++(int dummy)
        {
            return operator++();
        }

        bool operator==(const self_type& rhs) const
        {
            return curr == rhs.curr && done == rhs.done;
        }

        bool operator!=(const self_type& rhs) const
        {
            return curr != rhs.curr || done != rhs.done;
        }

        Itr operator->()
        {
            return curr;
        }
        reference operator*() { return *curr; }
    };

    SlidingWindow(const int size, const T& init)
        : data(size, init)
        , curr_begin(0)
    {
    }

    SlidingWindow(const std::initializer_list<T>& init)
        : data(init)
        , curr_begin(0)
    {
    }

    typedef typename std::vector<T>::iterator it_type;

    iterator<it_type, T> begin()
    {
        return iterator<it_type, T>(data.begin(), data.end(), curr_begin, false);
    }
    iterator<it_type, T> end()
    {
        return iterator<it_type, T>(data.begin(), data.end(), curr_begin, true);
    }

    iterator<it_type, const T> begin() const
    {
        return iterator<it_type, const T>(data.begin(), data.end(), curr_begin, false);
    }
    iterator<it_type, const T> end() const
    {
        return iterator<it_type, const T>(data.begin(), data.end(), curr_begin, true);
    }

    std::size_t size() const { return data.size(); }

    T& operator[](const size_type& index) { return data[index]; }
    const T& operator[](const size_type& index) const { return data[index]; }

    void append(T el)
    {
        data[curr_begin] = el;
        curr_begin = (curr_begin + 1) % data.size();
    }

    void clear(T el) { std::fill(data.begin(), data.end(), el); }

    SlidingWindow& operator<<(T& el)
    {
        append(el);
        return *this;
    }
};
