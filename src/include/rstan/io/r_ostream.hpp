#ifndef R__IO_R_OSTREAM_HPP
#define R__IO_R_OSTREAM_HPP

#include <streambuf>
#include <ostream>
// #include <Rinternals.h>
#include <ext/R.h> // R_FlushConsole, R_ProcessEvents
#include <R_ext/Print.h>

/**
 * Similar version of both std::cout and std::cerr are implemented for
 * RStan to write to cout and cerr of R.
 *
 * See
 * http://gcc.gnu.org/onlinedocs/libstdc++/manual/bk01pt11ch25.html#io.streambuf.derived
 * http://goo.gl/mKmeP
 * and http://www.cplusplus.com/reference/iostream/streambuf/overflow/
 *
 *
 */

namespace rstan {

  namespace io {

    class r_cout_streambuf : public std::streambuf {
    public:
      r_cout_streambuf() {}

    protected:
      /**
       * @param  c  An additional character to consume.
       * @return  EOF to indicate failure, something else (usually
       *          @a c, or not_eof())
       */
      virtual int_type overflow(int_type  c) {
        if (c != EOF)
          Rprintf("%c", c);
        return c;
      }

      virtual int sync() {
         R_FlushConsole();
         R_ProcessEvents();
         return 0;
      }

      virtual std::streamsize xsputn(const char_type* s, std::streamsize n) {
        Rprintf("%.*s", static_cast<int>(n), s);
        return n;
      }
    };

    class r_cerr_streambuf : public std::streambuf {
    public:
      r_cerr_streambuf() {}

    protected:
      virtual int_type overflow(int_type c) {
        if (c != EOF)
          REprintf("%c", c);
        return c;
      }

      virtual int sync() {
         R_FlushConsole();
         R_ProcessEvents();
         return 0;
      }

      virtual std::streamsize xsputn(const char_type* s, std::streamsize n) {
        REprintf("%.*s", static_cast<int>(n), s);
        return n;
      }
    };


    template <class T>
    class r_ostream : public std::ostream {
    protected:
      T buf;
    public:
      r_ostream(bool u) : std::ostream(&buf) {
        if (u)
          setf(std::ios_base::unitbuf);
      }
    };

    /**
     * Define global rstan::io::rcout and rstan::io::rcerr,
     * which can be used similarly as std::cout and std::cerr.
     *
     */
    // extern
    static r_ostream<r_cout_streambuf> rcout(false);
    // extern
    static r_ostream<r_cerr_streambuf> rcerr(true);
  }

}

#endif
