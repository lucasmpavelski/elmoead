#include "error.h"

void warning(const bool condition, const char* message)
{
    if (condition)
        warning(message);
}

void warning(const char* message)
{
    std::cerr << "WARNING: " << message << std::endl;
}

AssertionFailureException::AssertionFailureException(const char* expression,
                                                     const char* file,
                                                     const int line,
                                                     const std::string& message)
    : expression(expression)
    , file(file)
    , line(line)
    , message(message)
{
    std::ostringstream outputStream;
    if (!message.empty())
        outputStream << message << ": ";

    if (expression == "false" || expression == "0" || expression == "FALSE")
        outputStream << "Unreachable code assertion";
    else
        outputStream << "Assertion '" << expression << "'";

    outputStream << " failed in file '" << file << "' line " << line;
    report = outputStream.str();
}
