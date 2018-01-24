#include "io.h"

char* streamToString(FILE *f)
{
    if (f == NULL)
        return NULL;
    const size_t curr_pos = (size_t)ftell(f); // current position
    fseek(f, 0, SEEK_END);
    const size_t last_pos = (size_t)ftell(f); // last position
    const size_t diff = last_pos - curr_pos;
    char* str = new char[diff + 1];
    fseek(f, -diff, SEEK_END);
    const size_t read = fread(str, sizeof(char), diff, f);
    assert((int)read == diff);
    str[diff] = '\0';
    fseek(f, -diff, SEEK_END);
    return str;
}

char* fileToString(FILE *f)
{
    if (f == NULL)
        return NULL;
    fseek(f, 0, SEEK_END);
    const size_t pos = ftell(f); // last position or file size in chars
    rewind(f);
    char* str = new char[pos + 1];
    const size_t read = fread(str, sizeof(char), pos, f);
    if (pos != read)
        return NULL;
    str[pos] = '\0';
    return str;
}

std::string exec2Str(const char *cmd)
{
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe))
    {
        if(fgets(buffer, 128, pipe) != NULL)
            result += buffer;
    }
    pclose(pipe);
    return result;
}
