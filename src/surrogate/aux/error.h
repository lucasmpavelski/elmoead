#pragma once

#include <iostream>
#include <sstream>
#include <exception>
#include <string>
#include <cassert>

void warning(const bool condition, const char* message);
void warning(const char* message);

struct AssertionFailureException : public std::exception
{
    using string = std::string;

private:
    const char* expression;
    const char* file;
    int line;
    string message;
    string report;

public:
    AssertionFailureException(const char* expression,
                              const char* file,
                              const int line,
                              const std::string& message);

    virtual const char* what() const noexcept { return report.c_str(); }
    const char* getExpression() const noexcept { return expression; }
    const char* getFile() const noexcept { return file; }
    int getLine() const noexcept { return line; }
    const char* Message() const noexcept { return message.c_str(); }
};

/// Assert that EXPRESSION evaluates to true, otherwise raise AssertionFailureException with associated MESSAGE (which may use C++ stream-style message formatting)
#define throw_assert(EXPRESSION, MESSAGE) if(EXPRESSION) {} else throw AssertionFailureException(#EXPRESSION, __FILE__, __LINE__, static_cast<std::ostringstream&>(std::ostringstream() << MESSAGE).str())

