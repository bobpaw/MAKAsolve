# Style guide

Use `clang-format` to format files before committing. It can be installed via
the system package manager or pip. Then use `git clang-format` before
committing. Using a tool allows developers to write as they please and commit
as recommended.

## Line length

Line length should be kept to 80 characters at most.

## Whitespace

Use tabs for indentation and continuation. Do put spaces between flow-control
statements (like `if`, `while`, or `switch`) and parentheses. Do not put
spaces between function names and parentheses.

## Short if

# Functions 

## Function declaration

Do not put whitespace between the function name and the opening parenthesis.
Put the arguments on a new line if there are a lot of them. Put the semicolon
or opening curly brace on the same line as the closing parenthesis for the
argument list.

For example:

```c++
int func(int x);

ReturnType longfunctionname(
    Namespace1::AnotherLongType variable1,
    Namespace1::ALongerTypeName variable2
);
```

## Short functions

Put short class functions (e.g. class getters/setters) on a single line.
Non-class functions should always use multiple lines. The exception is that
empty functions can always be put on a single line.

## Long function calls

Ideally function calls that would go over 80 characters should be continued
onto multiple lines. Arguments which are the result of some other function
should be kept together to make them easier to understand.

```c++
func(
    argument1, argument2, argument3, meshname,
    functionforargument(meshname, elementid)
);
```

# Documentation

## Comments

Prefer `/* */` for multiline comments.

## Doxygen

Document public headers. Use Doxygen `\brief` style with backslashes instead
of at symbols. Prefer to introduce blocks with `/**` and member documentation
with `/**< */` or `///`. Try to document all function parameters.

For example:

```c++
/**
 * \brief Perform a function.
 * \param x First value.
 * \param y Second value.
 * \return A value based on x and y.
 */
int foo(int x, int y);
```

