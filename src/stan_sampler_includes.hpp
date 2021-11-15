#if (defined(__clang__) && (__clang_major__ > 3 || (__clang_major__ == 3 && __clang_minor__ >= 7))) || \
    (defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)))
#  define SUPPRESS_DIAGNOSTIC 1
#endif

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS 1
#ifdef SUPPRESS_DIAGNOSTIC
#  ifdef __clang__
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wunknown-pragmas"
#    pragma clang diagnostic ignored "-Wunused-variable"
#    pragma clang diagnostic ignored "-Wunused-parameter"
#    pragma clang diagnostic ignored "-Wunused-local-typedef"
#    pragma clang diagnostic ignored "-Wunused-function"
#    pragma clang diagnostic ignored "-Wsign-compare"
#    pragma clang diagnostic ignored "-Wlanguage-extension-token"
#    pragma clang diagnostic ignored "-Winfinite-recursion"
#    pragma clang diagnostic ignored "-Wignored-qualifiers"
#    pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#    pragma clang diagnostic ignored "-Wdeprecated-copy"
#    pragma clang diagnostic ignored "-Wshorten-64-to-32"
#    pragma clang diagnostic ignored "-Wfloat-conversion"
#    if __clang_major__ >= 8
#      pragma clang diagnostic ignored "-Wimplicit-float-conversion"
#    endif
#    if __clang_major__ >= 10
#      pragma clang diagnostic ignored "-Wimplicit-int-float-conversion"
#    endif
#  else
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wunknown-pragmas"
#    pragma GCC diagnostic ignored "-Wunused-variable"
#    pragma GCC diagnostic ignored "-Wunused-parameter"
#    pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#    pragma GCC diagnostic ignored "-Wunused-function"
#    pragma GCC diagnostic ignored "-Wsign-compare"
#    pragma GCC diagnostic ignored "-Wignored-qualifiers"
#    if __GNUC__ >= 6
#      pragma GCC diagnostic ignored "-Wignored-attributes"
#    endif
#    pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#    pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#    if defined(__MINGW32__) || defined(__MINGW64__)
#      pragma GCC diagnostic ignored "-Wattributes"
#    endif
// This is for gcc under -Wpedantic, since some some warnings can't
// be silenced. It is highly undesirably, as it can suppress other,
// useful warnings with later code.
#    include <boost/math/tools/config.hpp>
#    ifdef BOOST_MATH_USE_FLOAT128
#      pragma GCC system_header
#    endif
#  endif
#endif

#include <stan/callbacks/interrupt.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <stan/callbacks/stream_writer.hpp>

#if defined(_WIN32)
#  define BOOST_MATH_DISABLE_DEPRECATED_03_WARNING 1
#endif
#include <stan/io/dump.hpp>
#include <stan/io/empty_var_context.hpp>
#include <stan/io/var_context.hpp>

#ifdef SUPPRESS_DIAGNOSTIC
#  ifdef __clang__
#    pragma clang diagnostic pop
#  else
#    pragma GCC diagnostic pop
#  endif
#endif

