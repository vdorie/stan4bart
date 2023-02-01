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
#    pragma clang diagnostic ignored "-Wshorten-64-to-32"
#    if __has_warning("-Wdeprecated-copy")
#      pragma clang diagnostic ignored "-Wdeprecated-copy"
#    endif
#    if __has_warning("-Wfloat-conversion")
#      pragma clang diagnostic ignored "-Wfloat-conversion"
#    endif
#    if __has_warning("-Wimplicit-float-conversion")
#      pragma clang diagnostic ignored "-Wimplicit-float-conversion"
#    endif
#    if __has_warning("-Wimplicit-int-float-conversion")
#      pragma clang diagnostic ignored "-Wimplicit-int-float-conversion"
#    endif
#    if __has_warning("-Wunused-but-set-variable")
#      pragma clang diagnostic ignored "-Wunused-but-set-variable"
#    endif
#    if __has_warning("-Wdeprecated-declarations")
#      pragma clang diagnostic ignored "-Wdeprecated-declarations"
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
// be silenced. It is highly undesirable, as it can suppress other,
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

#if __cplusplus >= 201703L
#  define _HAS_AUTO_PTR_ETC 0
#endif

#include <stan/io/dump.hpp>
#include <stan/io/empty_var_context.hpp>
#include <stan/io/var_context.hpp>

#undef _HAS_AUTO_PTR_ETC

#ifdef SUPPRESS_DIAGNOSTIC
#  ifdef __clang__
#    pragma clang diagnostic pop
#  else
#    pragma GCC diagnostic pop
#  endif
#endif

