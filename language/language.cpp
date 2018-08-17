// language.cpp

#undef NDEBUG
#include "../common.inc"

#include "language.hpp"



namespace Language_sp
{


// Utf8

Utf8::Utf8 (const string &fName)
: is (fName, ios_base::binary | ios_base::in)
{ 
  static_assert (sizeof (Codepoint) >= 4, "Too small Codepoint");
  ASSERT (is. good ());
  char c;
  EXEC_ASSERT (getChar (is, c));
  ASSERT (c == '\xEF');  // ??
  EXEC_ASSERT (getChar (is, c));
  ASSERT (c == '\xBB');  // ??
  EXEC_ASSERT (getChar (is, c));
  ASSERT (c == '\xBF');  // ??
}



bool Utf8::get (Codepoint &v)
{
  char c;
  if (! getChar (is, c))
    return false;
  const size_t next = nextBytes (c);
  ASSERT (next < sizeof (Codepoint));
  v = static_cast <unsigned char> (c);
  // Continuation bytes
  FOR (size_t, i, next)
  {
    EXEC_ASSERT (getChar (is, c));
    ASSERT ((c & '\xC0') == '\x80');
    v *= 64;
    v += static_cast <unsigned char> (c & '\x3F');
  }
  return true;
}



size_t Utf8::nextBytes (char &c)
{
  if (! (c & '\x80'))
    return 0;
  size_t n = 0;
  while (c & '\x80')
  {
    c = c << 1;
    n++;
  }
  ASSERT (n > 1);
  c = c >> n;
  return n - 1;
}




// Alphabet

const Vector<Codepoint> Alphabet::delimiters 
  ({ 0x0D , 0x0A, ' ', 0x09, '.', ',', ';', ':', '!', '?', '&', '(', ')' 
  , '-', '+', '*', '/', '=', '[', ']', '{', '}', '|', '\\', '<', '>'
  , '"', '\'', 0xAB/*<<*/, 0xBB/*>>*/, 0xA9/*Copyright*/, 0x2014/*dash*/
  , '#', '$', '%', '@'
  });
                 


void Alphabet::qc () const
{ 
  if (! qc_on)
    return;

  ASSERT (size);
  ASSERT (startCapital > 32);
  ASSERT (startSmall >= startCapital + (int) size);
  ASSERT (extraCapital. size () == extraSmall.size ());
  if (! extraCapital. empty ())
  { 
    ASSERT (extraSmall [0] > extraCapital [0]);
    FOR (size_t, i, extraCapital. size ())
    {
      ASSERT (extraCapital [i]);
      ASSERT (! toSmall_ (extraCapital [i]));
      ASSERT (! toSmall_ (extraSmall   [i]));
      ASSERT (extraSmall [i] == extraCapital [i] + extraCapitalShift ());
    }
  }
  for (const Codepoint v : vowels)
  {
    ASSERT (v);
    ASSERT (toSmall (v) == v);
  }
  // delimiters
  Set<Codepoint> delimiterSet;
  for (const Codepoint d : delimiters)
  {
    ASSERT (d);
    ASSERT (! toSmall (d));
    delimiterSet << d;
  }
  ASSERT (delimiters. size () == delimiterSet. size ());
  // names
  ASSERT (names. size () == size + extraCapital. size ());
  Set<string> nameSet;
  for (const string& name : names)
  {
    ASSERT (! name. empty ());
    nameSet << name;
  }
  ASSERT (nameSet. size () == names. size());
}



Codepoint Alphabet::toSmall_ (Codepoint c) const
{ 
  ASSERT (c >= 0);
  if (   c >= startSmall
      && c < startSmall + (int) size
      )
    return c;
  if (   c >= startCapital
      && c < startCapital + (int) size
      )
    return c + capitalShift ();
  return 0;
}



Codepoint Alphabet::toSmall (Codepoint c) const
{ 
  if (const Codepoint sc = toSmall_ (c))
    return sc;
  if (extraSmall. contains (c))
    return c;
  if (extraCapital. contains (c))
    return c + extraCapitalShift ();
  return 0;
}



string Alphabet::small2name (Codepoint c) const
{ 
  if (c < 128)  // ASCII
    return string (1, c);
  if (c >= startSmall && c < startSmall + size)
    return names [c - startSmall];
  size_t index;
  EXEC_ASSERT (extraSmall. find (c, index));
  return names [size + index];
}



}  // namespace

