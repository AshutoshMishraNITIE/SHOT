# ifndef CPPAD_CONFIGURE_HPP
# define CPPAD_CONFIGURE_HPP
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-19 Bradley M. Bell

CppAD is distributed under the terms of the
             Eclipse Public License Version 2.0.

This Source Code may also be made available under the following
Secondary License when the conditions for such availability set forth
in the Eclipse Public License, Version 2.0 are satisfied:
      GNU General Public License, Version 2.0 or later.
---------------------------------------------------------------------------- */

/*!
$begin configure.hpp$$
$spell
    noexcept
    pragmas
    unreferenced
    CppAD
    cppad
    yyyymmdd
    yyyy
    mm
    dd
    adolc
    cmake
    colpack
    eigen
    ipopt
    gettimeofday
    namespace
    mkstemp
    tmpnam
    nullptr
    sizeof
    std
    hpp
    addr
$$

$section Preprocessor Symbols Set By CMake Command$$

$head CPPAD_COMPILER_HAS_CONVERSION_WARN$$
is the compiler a variant of g++ and has conversion warnings
$srccode%hpp% */
# define CPPAD_COMPILER_HAS_CONVERSION_WARN 1
/* %$$

$head CPPAD_DISABLE_SOME_MICROSOFT_COMPILER_WARNINGS$$
This macro is only used to document the pragmas that disables the
follow warnings:

$subhead C4100$$
unreferenced formal parameter.

$subhead C4127$$
conditional expression is constant.

$srccode%hpp% */
# define CPPAD_DISABLE_SOME_MICROSOFT_COMPILER_WARNINGS 1
# if _MSC_VER
# pragma warning( disable : 4100 )
# pragma warning( disable : 4127 )
# endif
# undef CPPAD_DISABLE_SOME_MICROSOFT_COMPILER_WARNINGS
/* %$$

$head CPPAD_USE_CPLUSPLUS_2011$$
Should CppAD use C++11 features. This will be true if the current
compiler flags request C++11 features and the install procedure
determined that all the necessary features are available.
$srccode%hpp% */
# if     _MSC_VER
# define    CPPAD_USE_CPLUSPLUS_2011 1
# else   //
# if         __cplusplus >= 201100
# define         CPPAD_USE_CPLUSPLUS_2011 1
# else       //
# define         CPPAD_USE_CPLUSPLUS_2011 0
# endif      //
# endif //
/* %$$

$head CPPAD_PACKAGE_STRING$$
cppad-yyyymmdd as a C string where yyyy is year, mm is month, and dd is day.
$srccode%hpp% */
# define CPPAD_PACKAGE_STRING "cppad-20190831"
/* %$$

$head CPPAD_HAS_ADOLC$$
Was a adolc_prefix specified on the cmake command line.
$srccode%hpp% */
# define CPPAD_HAS_ADOLC 0
/* %$$

$head CPPAD_HAS_COLPACK$$
Was a colpack_prefix specified on the cmake command line.
$srccode%hpp% */
# define CPPAD_HAS_COLPACK 0
/* %$$

$head CPPAD_HAS_EIGEN$$
Was a eigen_prefix specified on the cmake command line.
$srccode%hpp% */
# define CPPAD_HAS_EIGEN 0
/* %$$

$head CPPAD_HAS_IPOPT$$
Was a ipopt_prefix specified on the cmake command line.
$srccode%hpp% */
# define CPPAD_HAS_IPOPT 0
/* %$$

$head CPPAD_DEPRECATED$$
This symbol is not currently being used.
$srccode%hpp% */
# define CPPAD_DEPRECATED 
/* %$$

$head CPPAD_BOOSTVECTOR$$
If this symbol is one, and _MSC_VER is not defined,
we are using boost vector for CPPAD_TESTVECTOR.
It this symbol is zero,
we are not using boost vector for CPPAD_TESTVECTOR.
$srccode%hpp% */
# define CPPAD_BOOSTVECTOR 0
/* %$$

$head CPPAD_CPPADVECTOR$$
If this symbol is one,
we are using CppAD vector for CPPAD_TESTVECTOR.
It this symbol is zero,
we are not using CppAD vector for CPPAD_TESTVECTOR.
$srccode%hpp% */
# define CPPAD_CPPADVECTOR 1
/* %$$

$head CPPAD_STDVECTOR$$
If this symbol is one,
we are using standard vector for CPPAD_TESTVECTOR.
It this symbol is zero,
we are not using standard vector for CPPAD_TESTVECTOR.
$srccode%hpp% */
# define CPPAD_STDVECTOR 0
/* %$$

$head CPPAD_EIGENVECTOR$$
If this symbol is one,
we are using Eigen vector for CPPAD_TESTVECTOR.
If this symbol is zero,
we are not using Eigen vector for CPPAD_TESTVECTOR.
$srccode%hpp% */
# define CPPAD_EIGENVECTOR 0
/* %$$

$head CPPAD_HAS_GETTIMEOFDAY$$
If this symbol is one, and _MSC_VER is not defined,
this system supports the gettimeofday function.
Otherwise, this symbol should be zero.
$srccode%hpp% */
# define CPPAD_HAS_GETTIMEOFDAY 1
/* %$$

$head CPPAD_TAPE_ADDR_TYPE$$
Is the type used to store address on the tape. If not size_t, then
<code>sizeof(CPPAD_TAPE_ADDR_TYPE) <= sizeof( size_t )</code>
to conserve memory.
This type must support std::numeric_limits,
the <= operator,
and conversion to size_t.
Make sure that the type chosen returns true for is_pod<CPPAD_TAPE_ADDR_TYPE>
in pod_vector.hpp.
This type is later defined as addr_t in the CppAD namespace.
$srccode%hpp% */
# define CPPAD_TAPE_ADDR_TYPE unsigned int
/* %$$

$head CPPAD_TAPE_ID_TYPE$$
Is the type used to store tape identifiers. If not size_t, then
<code>sizeof(CPPAD_TAPE_ID_TYPE) <= sizeof( size_t )</code>
to conserve memory.
This type must support std::numeric_limits,
the <= operator,
and conversion to size_t.
Make sure that the type chosen returns true for is_pod<CPPAD_TAPE_ID_TYPE>
in pod_vector.hpp.
This type is later defined as tape_id_t in the CppAD namespace.
$srccode%hpp% */
# define CPPAD_TAPE_ID_TYPE unsigned int
/* %$$

$head CPPAD_MAX_NUM_THREADS$$
Specifies the maximum number of threads that CppAD can support
(must be greater than or equal four).

The user may define CPPAD_MAX_NUM_THREADS before including any of the CppAD
header files.  If it is not yet defined,
$srccode%hpp% */
# ifndef CPPAD_MAX_NUM_THREADS
# define CPPAD_MAX_NUM_THREADS 48
# endif
/* %$$

$head CPPAD_HAS_MKSTEMP$$
It true, mkstemp works in C++ on this system.
$srccode%hpp% */
# define CPPAD_HAS_MKSTEMP 1
/* %$$

$head CPPAD_HAS_TMPNAM_S$$
It true, tmpnam_s works in C++ on this system.
$srccode%hpp% */
# define CPPAD_HAS_TMPNAM_S 0
/* %$$

$head Symbols Conditional on C++11$$
The following symbols has two definitions, one when
$cref/C++11/configure.hpp/CPPAD_USE_CPLUSPLUS_2011/$$
is available and another when it is not.

$subhead CPPAD_NULL$$
This preprocessor symbol is used for a null pointer.
It is $code nullptr$$ when C++11 is available and $code 0$$ otherwise.

$subhead CPPAD_NOEXCEPT$$
This preprocessor symbol is
$code noexcept$$ when C++11 is available and empty otherwise.

$subhead CPPAD_NDEBUG_NOEXCEPT$$
This preprocessor symbol is
$code noexcept$$ when C++11 is available and $code NDEBUG$$ is defined.
Otherwise it is empty.


$end
*/
// -------------------------------------------------
# if CPPAD_USE_CPLUSPLUS_2011
# define CPPAD_NULL                nullptr
# define CPPAD_NOEXCEPT            noexcept
//
# ifdef NDEBUG
# define CPPAD_NDEBUG_NOEXCEPT     noexcept
# else
# define CPPAD_NDEBUG_NOEXCEPT
# endif
// -------------------------------------------------
# else
# define CPPAD_NULL                 0
# define CPPAD_NOEXCEPT
# define CPPAD_NDEBUG_NOEXCEPT
# endif
// -------------------------------------------------

# endif
