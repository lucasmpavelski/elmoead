#include "plotter.h"

int Plotter::count = 0;

void Plotter::command(const std::string& cmd)
{
    if (proc != nullptr) {
        fwrite(cmd.c_str(), cmd.length(), 1, proc);
        fflush(proc);
    }
}

void Plotter::close()
{
    command("exit\n");
    if (proc != nullptr)
        pclose(proc);
    proc = nullptr;
}

void Plotter::plot(const std::string& file_name, const unsigned dim)
{
    if (dim <= 2)
        command("plot '" + file_name + "'\n");
    else if (dim == 3)
        command("splot '" + file_name + "'\n");
    else if (dim > 3)
        command("plot '" + file_name + "' u 1:2:3 with lines palette\n");
}
