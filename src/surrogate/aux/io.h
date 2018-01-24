#pragma once

#include <iostream>
#include <string>
#include <memory>
#include <cassert>
#include <vector>
#include <sstream>

template<class T>
inline std::vector<T> str2Vector(const std::string& str) {
    std::stringstream ss(str);
    std::vector<T> v;
    T el;
    while (!ss.eof())
    {
        ss >> el;
        v.emplace_back(el);
    }
    return v;
}

template <typename Itr1, typename Itr2>
class printableSeq {
public:
    printableSeq(Itr1 beg, Itr2 end, const std::string& sep=" ")
        : beg(beg)
        , end(end)
        , sep(sep)
    {
    }

    std::ostream& operator()(std::ostream& os=std::cout) const
    {
        return os << *this;
    }

    friend std::ostream& operator<<(std::ostream& os,
                                    const printableSeq<Itr1, Itr2>& pa)
    {
        if (pa.beg != pa.end)
            os << *pa.beg;
        Itr1 p = pa.beg;
        while (++p != pa.end)
            os << pa.sep << *p;
        return os;
    }

    const std::string sep;
    Itr1 beg;
    Itr2 end;
};

template <typename Itr1, typename Itr2>
printableSeq<Itr1, Itr2> printSeq(Itr1 beg, Itr2 end,
                                  const std::string& sep = " ")
{
    return printableSeq<Itr1, Itr2>(beg, end, sep);
}

template <typename Itr1>
printableSeq<Itr1, Itr1> printSeq_n(Itr1 beg, const size_t n,
                                    const std::string& sep = " ")
{
    return printableSeq<Itr1, Itr1>(beg, beg + n, sep);
}

template <typename T>
struct SerializedWriter {
    T& dt;

    SerializedWriter(T&& dt) : dt(dt) {}

    friend std::ostream& operator<<(std::ostream& out,
                                    const SerializedWriter<T>& s)
    {
        out.write(reinterpret_cast<const char*>(&s.dt), sizeof(T));
        return out;
    }
};

template <typename T>
struct SerializedReader {
    SerializedReader(T& dt) : dt(dt) {}

    friend std::istream& operator>>(std::istream& in, SerializedReader<T>&& s)
    {
        in.read(reinterpret_cast<char*>(&s.dt), sizeof(T));
        return in;
    }

    T& dt;
};

template <typename T>
SerializedWriter<const T> serialize(T&& dt)
{
    return SerializedWriter<const T>(std::forward<T>(dt));
}

template <typename T>
SerializedReader<T> read_serialized(T& dt)
{
    return SerializedReader<T>(dt);
}

char* streamToString(FILE *f);
char* fileToString(FILE* in);
std::string exec2Str(const char* cmd);
