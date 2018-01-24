#pragma once

#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <memory>

class Plotter {
    FILE* proc;
    static int count;
    std::string file_name;

public:
    Plotter()
        : proc(popen("gnuplot -noraise -persist", "w"))
        , file_name("plot" + std::to_string(count) + ".dat")
    {
        if (proc == nullptr)
            std::cerr << "error: could not open gnuplot proc" << std::endl;
        count++;
    }

    Plotter(const Plotter& p)
        : Plotter()
    {}

    Plotter(Plotter&& p)
        : proc(p.proc)
        , file_name(std::move(p.fileName()))
    {
        p.proc = nullptr;
    }

    Plotter& operator=(Plotter p)
    {
        std::swap(proc, p.proc);
        return *this;
    }

    Plotter& operator=(Plotter&& p)
    {
        proc = p.proc;
        file_name = p.fileName();
        p.proc = nullptr;
    }

    virtual ~Plotter() { close(); }

    template <typename T>
    Plotter(const T data[], const int size, const int dim)
        : Plotter()
    {
        plot(data, size, dim);
    }

    std::string fileName() const { return file_name; }

    void command(const std::string& cmd);
    void close();
    void plot(const std::string& file_name, const unsigned dim = 1);

    template <typename T>
    void plot(const T data[], const int size, const unsigned dim = 1)
    {
        std::ofstream df(file_name);
        if (dim <= 3) {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < dim; j++)
                    df << data[i * dim + j] << " ";
                df << "\n";
            }
        } else {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < dim; j++)
                    df << j << " " << data[i * dim + j] << " " << i << '\n';
                df << "\n";
            }
        }
        df.close();
        plot(file_name, dim);
    }
};
