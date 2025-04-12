#ifndef CM_SCANNER_H
#define CM_SCANNER_H

#include <cctype>
#include <cstdio>

namespace cm
{

template <class _Type, size_t _buf_size>
class buffer_reader
{
protected:
  FILE        *src;
  _Type *const buff;
  _Type       *buff_end;
  _Type       *buff_pos;

  void _flush()
  {
    buff_end = buff + fread(buff, sizeof(_Type), _buf_size, src);
    buff_pos = buff;
  }

public:
  explicit buffer_reader(FILE *_src) : src(_src), buff(new _Type[_buf_size])
  {
    _flush();
  }

  buffer_reader(const buffer_reader &) = delete;
  buffer_reader(buffer_reader &&)      = delete;
  buffer_reader &operator=(const buffer_reader &) = delete;
  buffer_reader &operator=(buffer_reader &&) = delete;

  _Type get() const
  {
    return *buff_pos;
  }
  _Type next()
  {
    _Type result = get();
    buff_pos++;
    if (buff_pos == buff_end)
      _flush();
    return result;
  }

  bool end() {
    if (buff_pos == buff_end)
      _flush();
    return buff_pos == buff_end;
  }

  ~buffer_reader()
  {
    if (src != stdin)
      fclose(src);
    delete[] buff;
  }
};

template <class _BufferReader>
class scanner : protected _BufferReader
{
private:
  using _BufferReader::get;
  using _BufferReader::next;
  inline bool _is_ws(char c)
  {
    return c <= ' ';
  }
  inline bool _is_cr(char c)
  {
    return c == '\n' || c == '\r';
  }

  int _get_sign()
  {
    while (!isdigit(get()) && get() != '-')
      next();
    if (get() == '-')
      return next(), -1;
    return 1;
  }

public:
  scanner() = delete;
  using _BufferReader::_BufferReader;

  char next_char()
  {
    while (_is_ws(get()))
      next();
    return next();
  }

  char *next_token(char *s)
  {
    while (_is_ws(get()))
      next();
    while (!_is_ws(get()))
      *s++ = next();
    *s = '\0';
    return s;
  }

  char *next_line(char *s)
  {
    while (_is_ws(get()))
      next();
    while (!_is_cr(get()))
      *s++ = next();
    *s = '\0';
    return s;
  }

  template <class _Integer>
  _Integer next_integer()
  {
    _Integer sign = _get_sign();
    _Integer result(0);
    while (isdigit(get()))
      result = result * _Integer(10) + _Integer(next() - '0');
    return sign * result;
  }

  int next_int()
  {
    return next_integer<int>();
  }

  long long next_long()
  {
    return next_integer<long long>();
  }

  double next_double()
  {
    int    sign   = _get_sign();
    double result = 0;
    while (isdigit(get()))
      result = result * 10 + (next() - '0');
    if (get() == '.')
    {
      next();
      double cur_ep = 0.1;
      while (isdigit(get()))
        result += cur_ep * (next() - '0'), cur_ep *= 0.1;
    }
    return sign * result;
  }
};

} // namespace cm

#ifdef CM_DEBUG

namespace cm
{
using optimal_reader = buffer_reader<char, 1>;
}
#define DEFINE_SCANNER(sc, filename) cm::scanner<cm::optimal_reader> sc(stdin);

#else

namespace cm
{
using optimal_reader = buffer_reader<char, 1 << 16>;
}

#define DEFINE_SCANNER(sc, filename)                                           \
  cm::scanner<cm::optimal_reader> sc(fopen(filename, "r"));

#endif

#endif
