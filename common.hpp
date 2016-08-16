// common.hpp

#ifndef COMMON_HPP_64052  // random number
#define COMMON_HPP_64052


#ifdef _MSC_VER
  #pragma warning (disable : 4514)  // '...': unreferenced inline function has been removed
  #pragma warning (disable : 4592)  // symbol will be dynamically initialized (implementation limitation)
  #pragma warning (disable : 4625)  // copy constructor was implicitly defined as deleted
  #pragma warning (disable : 4626)  // assignment operator was implicitly defined as deleted
  #pragma warning (disable : 4710)  // function not inlined
  #pragma warning (disable : 4800)  // 'const char *': forcing value to bool 'true' or 'false' (performance warning)
  #pragma warning (disable : 4820)  // '...' bytes padding added after data member '...'
  #pragma warning (disable : 5026)  // move constructor was implicitly defined as deleted
  #pragma warning (disable : 5027)  // move assignment operator was implicitly defined as deleted
#endif

#include <time.h>
#include <cstring>
#include <string>
#include <stdexcept>
#include <limits>
#include <list>
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>



namespace Common_sp
{
using namespace std;



bool initCommon ();
  // Module initialization
  // Invoked automaticallly



extern vector<string> programArgs;
extern string programName;
extern ostream* logPtr;

void errorExit (const char* msg,
                bool segmFault = false);
  // Update: *logPtr
  // Invokes: if segmFault then abort() else exit(1)



template <typename T, typename S>
  T const_static_cast (const S* s)
    { return static_cast <T> (const_cast <S*> (s)); }
template <typename T, typename S>
  T const_static_cast (const S &s)
    { return static_cast <T> (const_cast <S&> (s)); }



typedef  unsigned int   uint; 
typedef  unsigned long  ulong; 
typedef  long int       lint;



extern const vector<bool> Bool;  // {false, true}

enum ebool {EFALSE = false, 
            ETRUE = true, 
            UBOOL = true + 1};

inline bool operator<= (ebool a, ebool b)
  { static const char rank [3/*ebool*/] = {0, 2, 1};
  	return rank [a] <= rank [b];
  }


inline void toggle (bool &b)
  { b = ! b; }


inline void advance (size_t &index, 
                     size_t size)
  // Update: index: < size
  { index++; if (index == size) index = 0; }
  

template <typename T> bool maximize (T &a, T b)
  { if (a < b) { a = b; return true; } return false; }

template <typename T> bool minimize (T &a, T b)
  { if (a > b) { a = b; return true; } return false; }

template <typename T /*:number*/> bool between (T x, T low, T high)
  { return x >= low && x < high; }

template <typename T/*:number*/> bool betweenEqual (T x, T low, T high)
  { return x >= low && x <= high; }

inline bool divisible (uint n,
                       uint divisor)
  { return ! (n % divisor); }
  
inline uint remainder (int n, uint div)
  { if (const int rem = n % (int) div) 
      return (uint) (rem > 0 ? rem : (rem + (int) div)); 
    return 0; 
  }
  
inline uint gcd (uint a,
                 uint b)
	// Greatest common divisor
	// Euclid's algorithm
	{ if (! b)
			return a;
		if (a < b)
			return gcd (b, a);
		return gcd (b, a % b);
	}

uint powInt (uint a,
             uint b);
  // Return: a^b
  // Time: O(log(b))


const size_t NO_INDEX = SIZE_MAX;


// char

inline bool isAlpha (char c)
  { return strchr ("abcdefghijklmnopqrstuvwxyz", tolower (c)); }

inline bool isDigit (char c)
  { return strchr ("0123456789", c); }
  
inline bool isLetter (char c)
  { return isAlpha (c) || isDigit (c) || c == '_'; }

inline char toUpper (char c)
  { return (char) toupper (c); }

inline char toLower (char c)
  { return (char) tolower (c); }

inline bool isUpper (char c)
  { return toUpper (c) == c; }

inline bool isLower (char c)
  { return toLower (c) == c; }

inline bool printable (char c)
  { return between (c, ' ', (char) 127); }
 


/* Usage:
  for (Iter<T> iter (t); iter. next (); )
    ...
*/
template <typename T>
struct Iter
{
private:
  T& t;
  typename T::iterator itNext;
public:
  typename T::iterator it;
    
  Iter (T &t_arg)
    : t (t_arg)
    , itNext (t_arg. begin ())
    {}
    
  bool next () 
    { it = itNext;
      if (it == t. end ())
        return false;
      itNext = it;
      itNext++;
      return true;
    }
  size_t getIndex () const
    { return (size_t) (it - t. begin ()); }
  typename T::value_type& operator* () const
    { return const_cast <typename T::value_type&> (*it); }  // *set::iterator = const set::value_type
  typename T::value_type* operator-> () const
    { return & const_cast <typename T::value_type&> (*it); }  // *set::iterator = const set::value_type
  typename T::value_type erase ()
    { typename T::value_type val = *it;
      itNext = t. erase (it); 
      return val;
    }
  void insert (const typename T::value_type &val)
    { itNext = t. insert (itNext, val) + 1; }
    // Opposite to erase()
};



template <typename T>
struct List : list<T>
{
private:
	typedef  list<T>  P;
public:

	List ()
	  {}
	  
	T at (size_t index) const
	  { size_t i = 0;
	  	CONST_ITER (typename List<T>, it, *this)
	  	  if (i == index)
	  	  	return *it;
	  	  else
	  	  	i++;
	  	throw runtime_error ("List index is out of range");
	  }
	size_t find (const T &t) const
	  { size_t i = 0;
	  	CONST_ITER (typename List<T>, it, *this)
	  	  if (*it == t)
	  	  	return i;
	  	  else
	  	  	i++;
	  	return NO_INDEX;
	  }
  bool isPrefix (const List<T> &prefix) const
    { typename List<T>::const_iterator wholeIt  =      P::begin ();
    	typename List<T>::const_iterator prefixIt = prefix. begin ();
    	for (;;)
    	{ if (prefixIt == prefix. end ())
	    		return true;
	    	if (wholeIt == P::end ())
	    		return false;
	      if (*prefixIt != *wholeIt)
	      	return false;
	      wholeIt ++;
	      prefixIt++;
	    }
    }
  List<T>& operator<< (T t) 
    { P::push_back (t); 
    	return *this;
    }    
  T popFront ()
    { const T t = P::front ();
      P::pop_front ();
      return t;
    }
  T popBack ()
    { const T t = P::back ();
      P::pop_back ();
      return t;
    }
};



// string

extern const string noString;

inline string nvl (const string& s,
                   const string& nullS = "-")
  { return s. empty () ? nullS : s; }
  	
inline string prepend (const string &prefix,
                    	 const string &s)
  { if (s. empty ())
  	  return string ();
  	return prefix + s;
  }

inline bool isQuoted (const string &s)
  { return ! s. empty () && s [0] == '\"' && s [s. size () - 1] == '\"'; }

inline string strQuote (const string &s)
  { return "\"" + s + "\""; }

inline string unQuote (const string &s)
  { return s. substr (1, s. size () - 2); }

inline bool isSpace (char c)
  { return c > '\0' && c <= ' ' && isspace (c); }

bool strBlank (const string &s);

template <typename T>
  string toString (const T t)
    { ostringstream oss;
      oss << t;
      return oss. str ();
    }

template <typename T>
  T str2 (const string &s)
    { T i;
      istringstream iss (s);
      iss >> i;
      if (   Common_sp::strBlank (s)
          || ! iss. eof ()
         )
        throw runtime_error ("Converting \"" + s + "\"");    
      return i;
    }

template <typename T>
  bool str2 (const string &s,
             T &t)
    { try { t = str2<T> (s); return true; } 
        catch (...) { return false; } 
    }

inline bool isLeft (const string &s,
                    const string &left)
  { return s. substr (0, left. size ()) == left; }

bool isRight (const string &s,
              const string &right);

bool trimPrefix (string &s,
                 const string &prefix);
  // Return: success

bool trimSuffix (string &s,
                 const string &suffix);
  // Return: success

void trimSuffixNonAlphaNum (string &s);

bool trimTailAt (string &s,
                 const string &tailStart);
  // Return: trimmed
  // Update: s

bool goodName (const string &name);

void strUpper (string &s);

void strLower (string &s);

bool isUpper (const string &s);

bool isLower (const string &s);

inline bool charInSet (char c,
		                   const string &charSet)
  { return charSet. find (c) != string::npos; }

string::const_iterator stringInSet (const string &s,
                    	           	  const string &charSet);
  // Return: != s.end() => *Return is not in charSet

size_t strCountSet (const string &s,
		                const string &charSet);

void strDeleteSet (string &s,
		               const string &charSet);

void trimLeading (string &s);

void trimTrailing (string &s);

inline bool contains (const string &hay,
                      const string &needle)
  { return hay. find (needle) != string::npos; }

inline bool contains (const string &hay,
                      char needle)
  { return hay. find (needle) != string::npos; }

size_t containsWord (const string& hay,
                     const string& needle);
  // Return: position of needle in hay; string::npos if needle does not exist

inline void trim (string &s)
  { trimTrailing (s);
    trimLeading (s); 
  }

void replace (string &s,
              char from,
              char to);
  
void replace (string &s,
              const string &fromChars,
              char to);
  
void replaceStr (string &s,
                 const string &from,
                 const string &to);

string to_c (const string &s);
  // " --> \", etc.

void collapseSpace (string &s);

  
string str2streamWord (const string &s,
                       size_t wordNum);
  // Return: May be string()
  // Input: wordNum: 0-based  

string str2sql (const string &s);
  // Return: ' s ' (replacing ' by '')

string sql2escaped (const string &s);


string findSplit (string &s,
                  char c = ' ');
	// Return: prefix of s+c before c
	// Update: s

string rfindSplit (string &s,
                   char c = ' ');
	// Return: suffix of c+s after c
	// Update: s

List<string> str2list (const string &s,
                       char c = ' ');
  // Invokes: findSplit(s,c)

string list2str (const List<string> &strList,
                 const string &sep = " ");


const char fileSlash = 
  #ifdef _MSC_VER
    '\\'
  #else  // UNIX
    '/'
  #endif
  ;


inline string getFileName (const string &path)  
  { const size_t pos = path. rfind ('/');
  	if (pos == string::npos)
  		return path;
  	return path. substr (pos + 1);
  }

bool fileExists (const string &fName);


#ifndef _MSC_VER
bool directoryExists (const string &dirName);
#endif


size_t strMonth2num (const string& month);
// Input: month: "Jan", "Feb", ... (3 characters)


// istream

bool getChar (istream &is,
              char &c);
  // Output: s if (bool)Return

void skipLine (istream &is);

void readLine (istream &is,
               string &s);
  // Output: s

string getToken (istream &is,
                 const string &skip,
                 const string &delimeters);
  // Return: empty() <=> eof

inline void pressAnyKey ()
  { cout << "Press any key..."; char c; cin >> c; }



struct Rand
// Numerical Recipes in C, p. 279
{
	static const long /*actually ulong*/ max_;
private:
	long seed;
	  // 0 < seed < max_
public:
	
	explicit Rand (ulong seed_arg = 1)
  	{ setSeed (seed_arg); }
	void setSeed (ulong seed_arg)
    {	seed = (long) (seed_arg % (ulong) max_);
	    qc ();
    }
	  
	ulong get (ulong max);
    // Return: in 0 .. max - 1
    // Input: 0 < max <= max_
  double getProb ()
    { return (double) get ((ulong) max_) / ((double) max_ - 1); }
    // Return: [0,1]
private:
	void qc () const;
	void run ();	
};




bool verbose (uint inc = 0);

class Verbose
{
	int verbose_old;
public:
	explicit Verbose (int verbose_arg);
	Verbose ();
	  // Increase verbosity
 ~Verbose ();
};

struct Unverbose 
{
	Unverbose ();
 ~Unverbose ();
};
  


struct Nocopy
{
protected:
	Nocopy ()
	  {}
private:
  Nocopy (const Nocopy &);
  Nocopy& operator= (const Nocopy &);
};



template <typename T>
struct Singleton : Nocopy
{
private:
	static bool beingRun;
protected:
	Singleton ()
	  { if (beingRun)
	  	  throw runtime_error ("Singleton");
	  	beingRun = true;
	  }
 ~Singleton ()  // virtual: segmentation fault
    { beingRun = false; }
};
template <typename T> bool Singleton<T>::beingRun = false;



struct Json;
struct JsonContainer;



struct Root
{
protected:
  Root () 
    {}
public:
  virtual ~Root () noexcept
    {}
    // A desrtructor should be virtual to be automatically invoked by a descendent class destructor
  virtual Root* copy () const
    { NOT_IMPLEMENTED; }
    // Return: the same type
    
  virtual void qc () const
    {}
  virtual void saveText (ostream& /*os*/) const 
    { NOT_IMPLEMENTED; }
    // Parsable output
  void saveFile (const string &fName) const;
    // if fName.empty() then do nothing
    // Invokes: saveText()
  virtual void print (ostream& os) const
    { saveText (os); }
    // Human-friendly
  virtual Json* toJson (JsonContainer* /*parent_arg*/,
                        const string& /*name_arg*/) const
    { NOT_IMPLEMENTED; }
	virtual bool empty () const
	  { NOT_IMPLEMENTED; }
  virtual void clear ()
    { NOT_IMPLEMENTED; }
    // Postcondition: empty()
  virtual void read (istream &/*is*/)
	  { NOT_IMPLEMENTED; }
	  // Input: a line of is
};

 

inline ostream& operator<< (ostream &os,
                            const Root &r) 
  { r. saveText (os);
    return os;
  }
  


template <typename T /*Root*/> 
struct AutoPtr : unique_ptr<T>
{
private:
	typedef  unique_ptr<T>  P;
public:

	explicit AutoPtr (T* t = nullptr) throw ()
	  : P (t)
	  {}
	AutoPtr (const AutoPtr<T> &t) 
	  : P (t->copy ())
	  {}
	AutoPtr<T>& operator= (T* t)
	  { P::reset (t);
	  	return *this;
	  }
	AutoPtr<T>& operator= (const AutoPtr<T> &t)
	  { P::reset (t->copy ());
	  	return *this;
	  }
};



struct Named : Root
{
  string name;
    // !empty(), no spaces at the ends, printable ASCII characeters

  Named ()
    {}  
  explicit Named (const string &name_arg); 
  Named* copy () const
    { return new Named (*this); }

  void saveText (ostream& os) const
    { os << name; }
	bool empty () const
	  { return name. empty (); }
  void clear ()
    { name. clear (); }
  void read (istream &is)
	  { is >> name; }
};


inline string named2name (const Named* n)
  { return n ? n->name : "<anonymous>"; }



typedef int (*CompareInt) (const void*, const void*);



// STL algorithms

template <typename To, typename From>
  void insertAll (To &to,
                  const From &from)
    { to. insert (to. begin (), from. begin (), from. end ()); }

template <typename To, typename From>
  void insertIter (To &to,
                   const From &from)
    { CONST_ITER (typename From, it, from)
        to << *it;
    }

template <typename T>
  void sort (T &t)
    { std::sort (t. begin (), t. end ()); }

template <typename T, typename C>
  void sort (T &t,
             C compare)
    { std::sort (t. begin (), t. end (), compare); }

template <typename Key, typename Value>
  bool contains (const map <Key, Value> &m,
                 const Key& key)
    { return m. find (key) != m. end (); }

template <typename Key, typename Value>
  bool find (const map <Key, Value> &m,
             const Key& key,
             Value &value)
    // Return: success
    // Output: value, if Return
    { const auto it = m. find (key);
    	if (it == m. end ())
    		return false;
    	value = it->second; 
    	return true;
    }

template <typename Key, typename Value>
  const Value* findPtr (const map <Key, Value> &m,
                        const Key& key)
    { const auto it = m. find (key);
    	if (it == m. end ())
    		return nullptr;
    	return & it->second; 
    }

template <typename Key, typename Value>
  const Value* findPtr (const map <Key, const Value* /*!nullptr*/> &m,
                        const Key& key)
    { const Value* value;
      if (find (m, key, value))
        return value;
      return nullptr;
    }

template <typename Key, typename Value>
  const Value& findMake (map <Key, const Value* /*!nullptr*/> &m,
                         const Key& key)
    { if (! contains (m, key))
        m [key] = new Value ();
      return * m [key];
    }



template <typename T>
struct Vector : vector<T>
{
private:
	typedef  vector<T>  P;
public:
	

	Vector ()
	  {}
	Vector (const vector<T> &x)
	  { *this = x; }
	Vector<T>& operator= (const vector<T> &x)
	  { return P::operator= (x); }
	explicit Vector (size_t n, 
	                 const T &value = T ())
	  : P (n, value)
	  {}
  Vector (initializer_list<T> init)
    : P (init)
    {}
/*static Vector<T> make (T a)
    { Vector<T> v (1);
      v [0] = a;
      return v; 
    }
  static Vector<T> make (T a,
                         T b)
    { Vector<T> v (2);
      v [0] = a;
      v [1] = b;
      return v; 
    }*/
	
	
  bool find (const T &value,
             size_t &index) const
	  // Output: index: valid if (bool)Return
	  { for (index = 0; index < P::size (); index++)
	      if (P::at (index) == value)
	        return true;
	    return false;
	  }
  typename Vector<T>::const_iterator constFind (const T &value) const
	  { typename Vector<T>::const_iterator it = P::begin (); 
	    while (it != P::end ()) 
	      if (*it == value)
	        break;
	      else
	        it++;
	    return it;
	  }
  typename Vector<T>::iterator find (const T &value)
	  { typename Vector<T>::iterator it = P::begin (); 
	    while (it != P::end ()) 
	      if (*it == value)
	        break;
	      else
	        it++;
	    return it;
	  }
  bool contains (const T &value) const
    { return constFind (value) != P::end (); }
  size_t count (const T &value) const
    { size_t n = 0;
      for (const T& t : *this)
        if (value == t)
          n++;
      return n;
    }
  size_t binSearch (const T &value,
                    bool exact = true) const
    // Return: if exact then NO_INDEX or vec[Return] = value else min {i : vec[i] >= value}
    // Requires: *this is sorted ascending
    { if (P::empty ())
    	  return NO_INDEX;
    	size_t lo = 0;  // vec.at(lo) <= value
    	size_t hi = P::size () - 1;  
    	// lo <= hi
    	if (value < P::at (lo))
    	  return exact ? NO_INDEX : lo;
    	if (P::at (hi) < value)
    	  return NO_INDEX;
    	// at(lo) <= value <= at(hi)
    	for (;;)
    	{
	    	const size_t m = (lo + hi) / 2;
	    	if (   P::at (m) == value
	    		  || P::at (m) <  value
	    		 )
	    		if (lo == m)  // hi in {lo, lo + 1}
	    			break;
	    		else
	    		  lo = m;
	      else
	      	hi = m;
	    }
	    if (P::at (lo) == value)
	    	return lo;
	    if (! exact || P::at (hi) == value)
	    	return hi;
	    return NO_INDEX;
    }

  Vector<T>& operator<< (const T &value)
    { P::push_back (value);
    	return *this;
    }
  Vector<T>& operator<< (const vector<T> &other)
    { P::insert (P::end (), other. begin (), other. end ());
    	return *this;
    }
  void setAll (const T &value)
    { ITER (typename P, it, *this)
    	  *it = value;
    }
  void eraseAt (size_t index)
    { eraseMany (index, index + 1); }
  void eraseMany (size_t from,
                  size_t to)
    { P::erase ( P::begin () + (ptrdiff_t) from
    	         , P::begin () + (ptrdiff_t) to
    	         ); 
    }
  void reverse ()
    { for (size_t i = 0; i < P::size (); i++)
    	{ const size_t j = P::size () - 1 - i;
    		if (i >= j)
    			break;
    	  swap (P::at (i), P::at (j));
    	}
    }
  void randomOrder (ulong seed)
		{ Rand rand (seed);
			ITER (typename P, it, *this)
	      swap (*it, P::at ((size_t) rand. get ((ulong) P::size ())));
		}
  void sortBubble ()  
    // Fast if *this is pre-sorted
    { FOR_START (size_t, i, 1, P::size ())
		    FOR_REV (size_t, j, i)
		      if (P::at (j + 1) > P::at (j))
        	  swap (P::at (j), P::at (j + 1));
		      else
		      	break;
    }
  T pop (size_t n = 1)
    { T t = T ();
      while (n)
      { t = P::at (P::size () - 1);
    	  P::pop_back ();
        n--;
      }
    	return t;
    }
};



template <typename T /* : Root */>
struct VectorPtr : Vector <const T*>
{
private:
	typedef  Vector <const T*>  P;
public:


  VectorPtr ()
	  {}	  
	explicit VectorPtr (size_t n, 
	                    const T* value = nullptr)
	  : P (n, value)
	  {}
	VectorPtr (const VectorPtr<T> &x)
	  : P ()
	  { *this = x; }
	template <typename U>
  	VectorPtr (const VectorPtr<U> &x)
  	  : P ()
  	  { P::reserve (x. size ());
  	    insertAll (*this, x);
  	  }	  
  static VectorPtr<T> make (const T* a)
    { VectorPtr<T> v (1);
      v [0] = a;
      return v; 
    }
  static VectorPtr<T> make (const T* a,
                            const T* b)
    { VectorPtr<T> v (2);
      v [0] = a;
      v [1] = b;
      return v; 
    }


  VectorPtr<T>& operator<< (const T* value)
    { return static_cast <VectorPtr<T>&> (P::operator<< (value)); }
  VectorPtr<T>& operator<< (const vector<const T*> &other)
    { return static_cast <VectorPtr<T>&> (P::operator<< (other)); }
	void deleteData ()
	  {	CONST_ITER (typename P, it, *this)
			  delete *it;
			P::clear ();  
	  }
  void erasePtr (size_t index)
    { delete P::at (index);
      P::eraseAt (index);
    }
  void sortBubble ()
    { FOR_START (size_t, i, 1, P::size ())
		    FOR_REV (size_t, j, i)
		      if (* P::at (j + 1) > * P::at (j))
        	  swap (P::at (j), P::at (j + 1));
		      else
		      	break;
    }
};



template <typename T /* : Root */>
struct VectorOwn : VectorPtr <T>
{
private:
	typedef  VectorPtr <T>  P;
public:

  VectorOwn ()
	  {}	  
	VectorOwn (const VectorOwn<T> &x)
	  : P ()
	  { *this = x; }
	VectorOwn (const VectorPtr<T> &x)
	  : P ()
	  { P::operator= (x); }
	VectorOwn<T>& operator= (const VectorOwn<T> &x)
	  { P::deleteData ();
	  	P::reserve (x. size ());
	  	CONST_ITER (typename P, it, x)
	  	  P::push_back (static_cast <const T*> ((*it)->copy ()));
	  	return *this;
	  }
 ~VectorOwn ()
    { P::deleteData (); }
};



struct DisjointCluster
// Cormen, Leiserson, Rivest, Introduction to Algorithms, p. 449
{
protected:
  DisjointCluster* parentDC;
    // !nullptr
    // Tree
    // = this <=> root
  size_t rankDC;
    // Upper bound on the height of *this
    // (Height = max. # arcs between *this and a leaf)
public:

protected:
	DisjointCluster ()
	  { init (); }
public:
	
	void init ()
    { rankDC = 0;
			parentDC = this;
		}
	void merge (DisjointCluster &other);
  DisjointCluster* getDisjointCluster ();
};



template <typename T>
struct Heap : Root, Nocopy
// Priority queue
// Heap property: comp(&arr[parent(index)],&arr[index]) >= 0
// More operations than in std::priority_queue
{
private:
  Vector<T> arr;
    // Elements are not owned by arr
  const CompareInt comp;
    // !nullptr
  typedef void (*SetHeapIndex) (T &item, size_t index);
    // Example: item.heapIndex = index
  const SetHeapIndex setHeapIndex;
    // Needed to invoke increaseKey()
public:


  struct Error : runtime_error
    { explicit Error (const string &str) : runtime_error (("Heap: " + str)) {} };


  explicit Heap (const CompareInt &comp_arg,
					       const SetHeapIndex &setHeapIndex_arg = nullptr,
					       size_t toReserve = 0)
    : comp (comp_arg)
    , setHeapIndex (setHeapIndex_arg)
    { arr. reserve (toReserve); }


  bool empty () const
    { return arr. empty (); }
  Heap& operator<< (T item)
    { arr << item;
      increaseKey (arr. size () - 1);
      return *this;
    }
  T increaseKey (size_t index)
    { T item = arr [index];
      size_t p;
      while (index && comp (& arr [p = parent (index)], & item) < 0)
      { assign (arr [p], index);
        index = p;
      }
      assign (item, index);
      return item;
    }
  T decreaseKey (size_t index)
    { T item = arr [index];
      heapify (index, arr. size ());
      return item;
    }
  T getMaximum () const
    { if (arr. empty ()) 
    	  throw Error ("getMaximum");
      return arr [0];
    }
  void deleteMaximum ()
    // Time: O(1) amortized
    { if (arr. empty ()) 
    	  throw Error ("deleteMaximum");
      T item = arr. back ();
      arr. pop_back ();
      if (arr. empty ())
        return;
      assign (item, 0);
      reinsertMaximum ();
    }
  void reinsertMaximum ()
    // Do this if the key of getMaximum() is changed
    { heapify (0, arr. size ()); }
  void sort ()
    { if (arr. empty ())
        return;
      for (size_t i = arr. size () - 1; i > 0; i--)
      { swap (0, i);
        heapify (0, i);
      }
    }
private:
  size_t parent (size_t index) const
    { if (! index) throw Error ("parent");
      return (index + 1) / 2 - 1;
    }
  size_t left (size_t index) const
    { return 2 * (index + 1) - 1; }
  size_t right (size_t index) const
    { return left (index) + 1; }
  void assign (T item,
               size_t index)
    { arr [index] = item;
      if (setHeapIndex)
        setHeapIndex (item, index);
    }
  void swap (size_t i,
             size_t j)
    { T item = arr [i];
      assign (arr [j], i);
      assign (item, j);
    }
  void heapify (size_t index,
                size_t maxIndex)
    // Requires: Heap property holds for all index1 < maxIndex except parent(index1) == index
    { if (maxIndex > arr. size ()) throw Error ("heapify: maxIndex");
      if (index >= maxIndex)       throw Error ("heapify: index");
      for (;;)
      { size_t extr = index;
        const size_t l = left (index);
        if (   l < maxIndex
            && comp (& arr [extr], & arr [l]) < 0)
          extr = l;
        const size_t r = right (index);
        if (   r < maxIndex
            && comp (& arr [extr], & arr [r]) < 0)
          extr = r;
        if (extr == index)
          break;
        swap (index, extr);
        index = extr;
      }
    }
public:

  // Test
  static void testStr ()
    { Heap <string> heap (strComp);
      heap << "Moscow" << "San Diego" << "Los Angeles" << "Paris";
      while (! heap. empty ())  
      { cout << heap. getMaximum () << endl;
        heap. deleteMaximum ();
      }
    }
private:
  static int strComp (const void* s1,
                      const void* s2)
    { const string& s1_ = * static_cast <const string*> (s1);
      const string& s2_ = * static_cast <const string*> (s2);
      if (s1_ < s2_) return -1;
      if (s1_ > s2_) return  1;
                     return  0;
    }
};



template <typename T>
struct Set : set<T>
{
private:
	typedef  set<T>  P;
public:
	bool universal;
	  // true => empty()
	

	Set ()
	  : universal (false)
	  {}
	explicit Set (bool universal_arg)
	  : universal (universal_arg)
	  {}
	Set (const Set<T> &other)
	  : P ()
	  { *this = other; }
	template <typename U>
	  Set (const map<T,U> &other)
	    : universal (false)
	    { for (typename map<T,U>::const_iterator it = other. begin (); it != other. end (); it++)
	        P::insert (it->first);
	    }
	Set<T>& operator= (const Set<T> &other)
	  { universal = other. universal;
	  	return static_cast <Set<T>&> (P::operator= (other)); 
	  }
  bool operator== (const Set<T> &other) const
    { return universal
               ? other. universal ? true : false
               : other. universal
                 ? false
                 :    P::size () == other. size ()
                   && contains (other);
    }


  bool contains (const T& el) const
    { return universal || P::find (el) != P::end (); }
  T front () const
    { return * P::begin (); }
  T back () const
    { return * std::prev (P::end ()); }
  template <typename U /* : T */>
    bool contains (const Set<U> &other) const
      { if (universal)
    	    return true;
    	  if (other. universal)
    		  return false;
    	  if (other. size () > P::size ())
    	    return false;
        for (const U& u : other)
          if (! contains (u))
            return false;
        return true;
      }
  bool intersects (const Set<T> &other) const
     { if (universal && other. universal)
     	   return true;
     	 if (universal)
     	   return ! other. empty ();
     	 if (other. universal)
     	   return ! P::empty ();
     	 if (P::size () < other. size ())
         return intersects_ (other);
       return other. intersects_ (*this);
     }
private:
  bool intersects_ (const Set<T> &other) const
    { for (const T& t : *this)
        if (other. contains (t))
          return true;
      return false;
    }
public:

  Set<T>& operator<< (const T &el)
    { if (! universal)
    	  P::insert (el);
    	return *this;
    }
  Set<T>& operator<< (const Set<T> &other)
    { if (! universal)
    	{ if (other. universal)
    	  { P::clear ();
    	  	universal = true;
    	  }
    	  else
    	    P::insert (other. begin (), other. end ());
    	}
    	return *this;
    }
  template <typename From>
    void insertAll (const From &from)
      { P::insert (from. begin (), from. end ()); }
  Set<T>& checkUnique (const T& el)
    { ASSERT (! contains (el));
      return operator<< (el);
    }
	void intersect (const Set<T> &other) 
		{ if (other. universal)
			  return;
			if (universal)
			{ operator= (other);
				return;
			}
      for (Iter <Set<T> > iter (*this); iter. next (); )
				if (! other. contains (*iter))
					iter. erase ();
		}
	size_t intersectSize (const Set<T> &other) const
	  // Return: universal <=> SIZE_MAX
		{ if (other. universal)
			  return universal ? SIZE_MAX : P::size ();
			if (universal)
				return other. size ();
		  size_t n = 0;
		  CONST_ITER (typename Set<T>, it, *this)
				if (other. contains (*it))
					n++;
			return n;
		}
  size_t setMinus (const Set<T> &other)
    { ASSERT (! universal);
    	size_t n = 0;
    	if (other. universal)
    	{ n = P::size ();
    		P::clear ();
    	}
    	else
	    //CONST_ITER (typename Set<T>, it, other)
	    	for (const T t : other)
	        n += P::erase (t);
      return n;
    }
};


template <typename T>
  void setMove (Set<T>* from,
	              Set<T>* to,
	              T el)
    { if (from == to)
    	  return;
    	IMPLY (from, ! from->universal);
    	IMPLY (to,   ! to  ->universal);
    	if (from)
    	  { EXEC_ASSERT (from->erase (el) == 1); }
    	if (to)
    	  { EXEC_ASSERT (to->insert (el). second); }
    }



struct DiGraph : Root
// Directed graph
// Bi-directional access
// n = number of nodes
// m = number of arcs
{
  struct Arc; 
  

  struct Node : Root, protected DisjointCluster
  {
    friend struct DiGraph;
    
    const DiGraph* graph;
  private:
    List<Node*>::iterator graphIt;
      // In graph->nodes
  public:
    List<Arc*> arcs [2 /*bool: out*/];

    Node* scc;  
      // Strongly-connected component
      // = root of SCC-subtree of DFS tree
      // scc->orderDfs <= orderDfs 
    size_t orderDfs;
      // Order by depth-first search
      // 0 <=> not visited by DFS
  private:
    bool inStack;  // auxiliary
      // => orderDfs
  public:

    explicit Node (DiGraph &graph_arg)
			: graph (nullptr)
			, scc (nullptr) 
			, orderDfs (0)
			, inStack (false)
			{ attach (graph_arg); }
    Node (const Node &other)
      : Root (other)
      , DisjointCluster (other)
      , graph (nullptr)
      , scc      (other. scc)
      , orderDfs (other. orderDfs)
      , inStack  (other. inStack)
      {}
      // To be followed by: attach()
    Node* copy () const
      { return new Node (*this); }
   ~Node ();
      // Remove this from graph 
      // Time: O(m) for all nodes
    void qc () const;
    void saveText (ostream& os) const;
      // Invokes: getName(), saveContent()
  protected:
    virtual void saveContent (ostream &/*os*/) const 
      {}
  public:

  	void attach (DiGraph &graph_arg);
      // Requires: !graph; no Arc's
      // Invokes: graph_arg.nodes.push_back(this)
      // Time: O(1)
    virtual string getName () const
      // Return: !empty()
	    { ostringstream oss;
		  	oss << this;
			  return oss. str ();
		  }
    bool isIncident (const Node* n,
                     bool out) const;
      // Return: n is among arcs[out]->node[out]
    VectorPtr<Node> getNeighborhood (bool out) const;
    VectorPtr<Node> getNeighborhood () const
      { return getNeighborhood (false) << getNeighborhood (true); }
      // Return: !contains(this)
    VectorPtr<DiGraph::Node> getChildren () const
      { return getNeighborhood (false); }
    void deleteNeighborhood (bool out);
  private:
    Node* setScc (size_t &visitedNum,
                  stack<Node*, vector<Node*> > &sccStack);
      // If node n is reachable from this then
      //   the SCC of n is reachable from this and the SCC is a subtree of DFS tree
      // Output: scc, orderDfs in nodes reachable from this
      // Tarjan's alogorithm
      // Return: n s.t. n->inStack and there is a path from this to n
      // Requires: n->inStack <=> there is a path from n to this
      // Time: O(n + m) for all nodes

    virtual void contractContent (const Node* /*from*/) {}
      // Required time: O(1)
  public:      
  	Node* getConnectedComponent ()
  	  { return static_cast <Node*> (getDisjointCluster ()); }
    void contract (Node* from);
      // Update: this: No parallel arcs
      // Invokes: contractContent(from), Arc::contractContent(), delete from
      // Requires: from != this
      //           No parallel arcs
      // Time: O(n + m log n) for all nodes
    void remove ();
      // Output: graph = nullptr
      // Requires: No Arc's
      // Invokes: list::erase()
  };


  struct Arc : Root
  {
    friend struct DiGraph;
    friend struct Node;
   
    Node* node [2 /*bool: out*/];
      // !nullptr
  private:
    List<Arc*>::iterator arcsIt [2 /*bool: out*/];
      // in node [b] -> arcs [! b]
  public:

    Arc (Node* start,
         Node* end)
      { node [false] = nullptr;
      	node [true]  = nullptr;
      	attach (start, end);
      }
  private:
  	void attach (Node* start,
                 Node* end);
      // Adds *this to the graph
      //   node[i]->arcs[!i].push_back(this)
      // Requires: !node[i]
      // Time: O(1)
  public:
  	Arc (const Arc& other)
  	  : Root (other)
  	  { node [false] = nullptr;
      	node [true]  = nullptr;
      }
      // To be followed by: attach()
    Arc* copy () const
      { return new Arc (*this); }
   ~Arc ();
      // Remove this from node->graph
      // Time: O(1)

    virtual void saveContent (ostream &/*os*/) const 
      {}
  private:
    virtual void contractContent (const Arc* /*from*/) 
      {}
      // Required time: O(1)
  };


  List<Node*> nodes;
    // size() == n


  DiGraph ()
    {}
  typedef  map <const Node* /*old*/, Node* /*new*/>  Old2new;
  DiGraph (const DiGraph &other)
    { Old2new old2new;
    	init (other, old2new);
    }
  DiGraph (const DiGraph &other,
           Old2new &old2new)
    { init (other, old2new); }
private:
  void init (const DiGraph &other,
             Old2new &old2new);
    // Output: old2new
public:
  DiGraph* copy () const
    { return new DiGraph (*this); }
 ~DiGraph ();
    // Invokes: Node::delete
    // Time: O(n + m)
  void qc () const;
  void saveText (ostream &os) const;   


  void connectedComponents ();
    // Output: Node::getConnectedComponent()  
    // Invokes: DisjointCluster::init()
  void scc (); 
    // Output: Node::{scc,orderDfs}
    // Invokes: Node::setScc()
    // Time: O(n + m)
  void contractScc ();
    // Output: DAG
    // Requires: After scc()
    // Invokes: Node::contract()
    // Time: O(n + m log n)
  Set<const Node*> getEnds (bool out) const;
    // Input: out: false - roots
    //             true  - leaves
  const Node* getRoot (bool out) const
		{ const Set<const DiGraph::Node*> ends (getEnds (out));
			if (ends. size () == 1)
			  return * ends. begin ();
			return nullptr;
		}
  typedef  map <const Node*, const Node*>  Node2Node;  
    // !nullptr
  static Node2Node reverse (const Node2Node& old2new);
  void borrowArcs (const Node2Node &node2node,
                   bool parallelAllowed);
    // Input: node2node: other node to *this node
    // Time: |node2node| log|node2node| outdegree_max outdegree'_max,
    //          where outdegree_max = max(outdegree(node2node.keys()),
    //                outdegree'_max = max(outdegree(node2node.values())
};



struct Tree : DiGraph
{
	struct Node : DiGraph::Node
	{
		Node (Tree &tree,
		      Node* parent_arg)
			: DiGraph::Node (tree)
			{ setParent (parent_arg); }
		  // Input: parent_arg: may be nullptr
   	void saveText (ostream &os) const;
   	  // Invokes: getName(), saveContent(), getSaveSubtreeP()

    virtual bool getSaveSubtreeP () const 
      { return true; }
	  virtual double getParentDistance () const
	    { return -1; }
	    // Return: -1 || >= 0
	  virtual string getNewickName () const
	    { return getName (); }
		static string name2newick (const string &s);
	private:
	  void printNewick_ (ostream &os,
	                     bool internalNames) const;
	    // Input: os.setprecision
	    // Invokes: getParentDistance(), getNewickName(), name2newick()
  public:
    void printNewick (ostream &os,
                      bool internalNames) const
      { printNewick_ (os, internalNames);
      	os << ';';
      }
    // Input: internalNames <=> print name at each internal node
	  const Tree& getTree () const
  	  { return * static_cast <const Tree*> (graph); }
		bool isLeaf () const
		  { return arcs [false]. empty (); }
		const Node* getParent () const
			{ return arcs [true]. empty () ? nullptr : static_cast <Node*> (arcs [true]. front () -> node [true]); }
		  // Return: nullptr <=> root
		const Node* getSuperParent (size_t height) const;
		  // Return: may be nullptr
		  // getSuperParent(0) = this
		void setParent (Node* newParent);
		  // Update: *newParent
		  //         getTree()->root if !newParent
    struct TipName : Root
    { string name; 
      size_t depth; 
      TipName ()
        : depth (0)
        {}
      TipName (const string &name_arg,
               size_t depth_arg)
        : name (name_arg)
        , depth (depth_arg)
        {}
      void saveText (ostream &os) const
        { os << name << '\t' << depth; }
    };
    TipName getTipName () const;
      // Return: identification of *this by a tip
		size_t getDepth () const
		  { if (const Node* parent_ = getParent ())
		  		return parent_->getDepth () + 1;
		  	return 0;
		  }
		size_t getHeight () const;
		  // Return: 0 <=> isLeaf()
		bool descendentOf (const Node* ancestor) const
		  { if (! ancestor)
		  	  return true;
		  	if (this == ancestor)
		  		return true;
		  	if (const Node* parent_ = getParent ())
		  		return parent_->descendentOf (ancestor);
		  	return false;
		  }
		const Node* getPrevAncestor (const Node* ancestor) const
		  { if (! ancestor)
		      return nullptr;
		    const Node* parent_ = getParent ();
		    if (parent_ == ancestor)
		      return this;
		    if (! parent_)
		      return nullptr;
		    return parent_->getPrevAncestor (ancestor);
		  }
		size_t getSubtreeSize () const;
		  // Does not count *this
		  // Return: 0 <=> isLeaf()
		size_t getLeavesSize () const;
    void getLeaves (VectorPtr<Node> &leaves) const;
      // Update: leaves
		const Node* getClosestLeaf (size_t &leafDepth) const;
		  // Return: !nullptr
		  // Output: leafDepth; nullptr <=> Return = this
    const Node* getOtherChild (const Node* child) const;
      // Return: May be nullptr; != child
      // Requires: getChildren().size() <= 2
    const Node* leftmostDescendent () const;
    const Node* rightmostDescendent () const;
	  void childrenUp ();
	    // Children->setParent(getParent())
	    // Post-condition: arcs[false].empty()
	  Node* isTransient () const
	    { return arcs [false]. size () == 1 
	    	        ? static_cast <Node*> (arcs [false]. front () -> node [false]) 
	    	        : nullptr; 
	    }
	    // Return: Single child of *this
	  void isolateChildrenUp ();
	    // Invokes: childrenUp(), remove()
	  Tree::Node* isolateTransient ()
			{ Tree::Node* transient = isTransient ();
				if (transient)
					isolateChildrenUp ();
			  return transient;
			}
	  bool deleteTransient ()
	    { if (! isolateTransient ())
		 			return false;
		    delete this;
        return true;
	    }
	  void deleteSubtree ();
	    // Postcondition: isLeaf()
	  const Node* makeRoot ();
	    // Redirect Arc's so that this = getTree()->root
	    // Return: old getTree()->root, !nullptr
    void getArea (uint distance,
                  VectorPtr<Tree::Node> &area,
                  VectorPtr<Tree::Node> &boundary) const
      { getArea_ (distance, nullptr, area, boundary); }
      // Output: area: connected Node's with one root, distinct
      //         boundary: distinct
      //         area.contains(boundary)
  private:
    void getArea_ (uint distance,
                   const Tree::Node* prev,
                   VectorPtr<Tree::Node> &area,
                   VectorPtr<Tree::Node> &boundary) const;
      // Update: area, bounday
      //         area.contains(boundary)
  public:
    template <typename Compare>
    	void sort (const Compare &compare)
  			{ VectorPtr<DiGraph::Node> children (getChildren ());
  				Common_sp::sort (children, compare);
  				CONST_ITER (VectorPtr<DiGraph::Node>, it, children)
  				{	Node* s = const_static_cast <Node*> (*it);
  					s->setParent (const_cast <Node*> (s->getParent ()));  // To reorder arcs[false]
  				  s->sort (compare);
  				}
  			}
	};
	const Node* root;
	  // nullptr <=> nodes.empty()


  Tree ()
    : root (nullptr)
    {}
  void qc () const;
	void saveText (ostream &os) const
	  { if (root)
	  	  root->saveText (os);
      os << endl;
    }


  void printAsn (ostream &os) const;	
    // http://www.ncbi.nlm.nih.gov/tools/treeviewer/biotreecontainer/
  static const Node* getLowestCommonAncestor (const Node* n1,
                                              const Node* n2);
    // Return: nullptr <=> !n1 || !n2
  static const Node* getLowestCommonAncestor (const VectorPtr<Node> &nodeVec);
    // Return: nullptr <= nodeVec.empty()
    // Input: nodeVec: may be nullptr
  static Set<const Node*> getParents (const VectorPtr<Node> &nodeVec);
    // Return: !nullptr, !contains(getLowestCommonAncestor(nodeVec)), contains(nodeVec)
    // Invokes: getLowestCommonAncestor(nodeVec)
  static Set<const Node*> getPath (const Node* n1,
                                   const Node* n2)
    { return move (getParents (VectorPtr<Node>::make (n1, n2))); }

  void setRoot ();
    // Output: root
  size_t deleteTransients ();
    // Return: # Node's delete'd

  template <typename Compare>
    void sort (const Compare &compare)
      { if (root)
      	  const_cast <Node*> (root) -> sort (compare); 
      }
private:
	static bool compare_std (const DiGraph::Node* a,
	                         const DiGraph::Node* b);

public:
  void sort ()
    { sort (compare_std); }
};



struct Progress : Nocopy
{
private:
	static size_t beingUsed;
	  // Singleton
public: 
	uint n_max;
	  // 0 <=> unknown
	bool active;
	uint n;
	string step;
	uint displayPeriod;
	
	explicit Progress (uint n_max_arg = 0,
	                   uint displayPeriod_arg = 1)
	  : n_max (n_max_arg)
	  , active (! beingUsed && displayPeriod_arg)
	  , n (0)
	  , displayPeriod (displayPeriod_arg)
	  { if (active) 
	  	  beingUsed++; 
	  }
 ~Progress ()
    { if (active)
    	{ if (! uncaught_exception ())
    	  { report ();
    	    cerr << endl;
    	  }
    	  beingUsed--;
    	}
    }
    
  void operator() (const string& step_arg = string ())
    { n++;
    	step = step_arg;
    	if (   active 
    		  && n % displayPeriod == 0
    		 )
    	  report ();
    }
private:
	void report () const
	  { cerr << '\r';
    #ifndef _MSC_VER
      cerr << "\33[2K";
    #endif
	    cerr << n; 
	  	if (n_max)
	  		cerr << " / " << n_max;
	  	if (! step. empty ())
	  		cerr << ' ' << step;
	  	cerr << ' ';
	  }
public:
	static void disable ()
	  { beingUsed++; }
	  
	  
	struct Start : Nocopy  
	{
	private:
		AutoPtr <Progress>& prog;
	public:
		explicit Start (AutoPtr<Progress> &prog_arg,
							      uint n_max = 0,
						        uint displayPeriod = 1);
	~Start ()
	   { prog. reset (); }
	};
};



struct Input : Root, Nocopy
{
protected:
	AutoPtr <char> buf;
	ifstream ifs;
public:
	bool eof;
	  // Init: false
	uint lineNum;
	  // # lines read
protected:
	Progress prog;
public:


protected:	
	Input (const string &fName,
      	 size_t bufSize,
         uint displayPeriod);
};
	


struct LineInput : Input
{
	string line;
	  // Current line

	
	explicit LineInput (const string &fName,
          	          size_t bufSize = 100 * 1024,
          	          uint displayPeriod = 0)
    : Input (fName, bufSize, displayPeriod)
    {}


	bool nextLine ();
  	// Output: eof, line
  	// Update: lineNum
	string getString ()
	  { string s; 
	  	while (nextLine ())
	  	{ if (! s. empty ())
	  			s += "\n";
	  	  s += line;
	  	}
	  	return s;
	  }
	void setVector (Vector<string> &vec)
	  { while (nextLine () && ! line. empty ())
	  	  vec << line;
	  }
};
	


struct ObjectInput : Input
{
	explicit ObjectInput (const string &fName,
          	            size_t bufSize = 100 * 1024,
          	            uint displayPeriod = 0)
    : Input (fName, bufSize, displayPeriod)
    {}


	bool next (Root &row);
	  // Output: row
  	// Update: lineNum
};
	


struct CharInput : Input
{
  uint charNum;
private:
  bool ungot;
  bool eol;
public:

	
	explicit CharInput (const string &fName,
              	      size_t bufSize = 100 * 1024,
              	      uint displayPeriod = 0)
    : Input (fName, bufSize, displayPeriod)
    , charNum ((uint) -1)
    , ungot (false)
    , eol (false)
    {}
  

	char get ();
	  // Output: eof
	  // Update: lineNum, charNum
	void unget ();
	  // Requires: To be followed by get()
	  

  struct Error : runtime_error
  { explicit Error (const CharInput &in,
		                const string &what_arg = string ()) 
			: runtime_error ("Error at line " + toString (in. lineNum + 1) 
		                   + ", pos. " + toString (in. charNum + 1)
		                   + (what_arg. empty () ? string () : (": " + what_arg + " is expected"))
		                  )
	    {}
	};
};
	


struct Token : Root
{
	enum Type {eSystem, eText, eNumber};
	Type type;
	string name;
	  // May be: !goodName(name)
	uint num;
	  // Valid if type = eNumber

	  
	Token ()
	  { clear (); }
	Token (const string& name_arg,
	       Type type_arg)
	  : type (type_arg)
	  , name (name_arg)
	  , num (0)
	  {}
	Token (CharInput &in,
	       Type expected)
	  { readInput (in, expected); }
	explicit Token (CharInput &in)
	  { readInput (in); }
	void saveText (ostream &os) const
	  { switch (type)
	  	{ case eSystem: os << name;                 break;
	  		case eText:   os << '\"' << name << '\"'; break;
	  		case eNumber: os << num;                  break;
	  	}
	  }

	  
	bool empty () const
	  { return type == eSystem && name. empty (); }
	void clear ()
	  { name. clear ();
	  	type = eSystem;
	    num = (uint) -1;
	  }
	  // empty()
	void readInput (CharInput &in);
	  // Postcondition: !empty()
	  // type != eNumber => !name.empty()
	  // type = eNumber => num != UINT_MAX
	void readInput (CharInput &in,
	               Type expected);
	  // Invokes: readInput(in)
};



struct OFStream : ofstream
{
	OFStream ()
	  {}
	OFStream (const string &dirName,
	          const string &pathName,
	          const string &extension)
	  { open (dirName, pathName, extension); }
	void open (const string &dirName,
	           const string &pathName,
	           const string &extension);
	  // Input: !pathName.empty()
};


template <typename T>
  OFStream& operator<< (OFStream &ofs,
                        const List<T> &ts)
    { for (auto& t : ts)
        ofs << t << endl;
      return ofs;
    }


template <typename T>
  OFStream& operator<< (OFStream &ofs,
                        const Vector<T> &ts)
    { for (auto& t : ts)
        ofs << t << endl;
      return ofs;
    }




struct Csv : Root
// Line of Excel .csv-file
{
private:
  const string &s;
  size_t pos;
public:

  
  Csv (const string &s_arg)
  : s (s_arg)
  , pos (0)
  {}
  
  
  bool goodPos () const
    { return pos < s. size (); }
  string getWord ();
    // Return: Next word
    // Requires: goodPos()
private:
  void findChar (char c)
    { while (goodPos () && s [pos] != c)
        pos++;
    }
};

  
  
void csvLine2vec (const string &line,
                  Vector<string> &words);
  // Output: words
  // Invokes: Csv



struct TabDel
// Usage: {<<field;}* str();
{
private:
  ostringstream tabDel;
public:
  
  TabDel ()
    {}
    
  template <typename T>
    TabDel& operator<< (const T &field)
      { if (tabDel. tellp ())  
          tabDel << '\t'; 
        tabDel << field; 
        return *this; 
      }    
  string str () const
    { return tabDel. str (); }
};




// Json


struct JsonNull;
struct JsonInt;
struct JsonDouble;
struct JsonString;
struct JsonArray;
struct JsonMap;
  

struct Json : Root, Nocopy  // Heaponly
{
protected:
  Json (JsonContainer* parent,
        const string& name);
  Json ()
    {};
public:  
  virtual void print (ostream& os) const = 0;
  
  virtual const JsonNull* asJsonNull () const
    { return nullptr; }  
  virtual const JsonInt* asJsonInt () const
    { return nullptr; }  
  virtual const JsonDouble* asJsonDouble () const
    { return nullptr; }  
  virtual const JsonString* asJsonString () const
    { return nullptr; }  
  virtual const JsonArray* asJsonArray () const
    { return nullptr; }  
  virtual const JsonMap* asJsonMap () const
    { return nullptr; }  

protected:
  static string toStr (const string& s)
    { return "'" + to_c (s) + "'"; }    
  static Token readToken (istream &is);
  static void parse (istream& is,
                     const Token& firstToken,
                     JsonContainer* parent,
                     const string& name);
public:
    
  int getInt () const;
    // Requires: JsonInt
  double getDouble () const;
    // Requires: JsonDouble
  string getString () const;
    // Requires: JsonString
  const Json* at (const string& name_arg) const;
    // Requires: JsonMap
  const Json* at (size_t index) const;
    // Requires: JsonArray
  size_t getSize () const;
    // Requires: JsonArray
};


struct JsonNull : Json
{
  explicit JsonNull (JsonContainer* parent,
                     const string& name = noString)
    : Json (parent, name)
    {}    
  void print (ostream& os) const
    { os << "null"; }

  const JsonNull* asJsonNull () const
    { return this; }  
};


struct JsonInt : Json
{
  int n;
  
  JsonInt (int n_arg,
           JsonContainer* parent,
           const string& name = noString)
    : Json (parent, name)
    , n (n_arg)
    {}
  void print (ostream& os) const
    { os << n; }

  const JsonInt* asJsonInt () const
    { return this; }  
};


struct JsonDouble : Json
{
  double n;
  uint decimals;

  JsonDouble (double n_arg,
              uint decimals_arg,
              JsonContainer* parent,
              const string& name = noString)
    : Json (parent, name)
    , n (n_arg)
    , decimals (decimals_arg)
    {}
  void print (ostream& os) const
    { os << fixed; os. precision ((streamsize) decimals); os << n; }

  const JsonDouble* asJsonDouble () const
    { return this; }  
};


struct JsonString : Json
{
  string s;

  JsonString (const string& s_arg,
              JsonContainer* parent,
              const string& name = noString)
    : Json (parent, name)
    , s (s_arg)
    {}
  void print (ostream& os) const
    { os << toStr (s); }

  const JsonString* asJsonString () const
    { return this; }  
};


struct JsonContainer : Json
{
protected:
  JsonContainer (JsonContainer* parent,
                 const string& name)
    : Json (parent, name)
    {}
  JsonContainer ()
    {}
public:  
};


struct JsonArray : JsonContainer
{
  friend struct Json;
private:
  VectorOwn<Json> data;
public:

  explicit JsonArray (JsonContainer* parent,
                      const string& name = noString)
    : JsonContainer (parent, name)
    {}
private:
  JsonArray (istream& is,
             JsonContainer* parent,
             const string& name);
public:
  void print (ostream& os) const;

  const JsonArray* asJsonArray () const
    { return this; }
};


struct JsonMap : JsonContainer
{
  friend struct Json;
private:
  typedef  map <string, const Json*>  Map;
  Map data;
public:
  
  explicit JsonMap (JsonContainer* parent,
                    const string& name = noString)
    : JsonContainer (parent, name)
    {}
  JsonMap ();
    // Output: jRoot = this
  explicit JsonMap (const string &fName);
private:
  JsonMap (istream& is,
           JsonContainer* parent,
           const string& name)
    : JsonContainer (parent, name)
    { parse (is); }
  void parse (istream& is);
public:
 ~JsonMap ();
  void print (ostream& os) const;

  const JsonMap* asJsonMap () const
    { return this; }
};



extern JsonMap* jRoot;



//

struct Chronometer : Nocopy
{
  const string name;
  ostream& os;
  clock_t time;

  Chronometer (const string &name_arg,
               ostream& os_arg)
    : name (name_arg)
    , os (os_arg)
    , time (0)
    {}
 ~Chronometer ()
    { os << name << ": Duration: ";       
      os << fixed; os. precision (2); os << (double) time / CLOCKS_PER_SEC << " sec." << endl; 
    }
    
  struct Measure 
  {
  private:
    Chronometer& chr;
    const clock_t start;
  public:
    explicit Measure (Chronometer &chr_arg)
      : chr (chr_arg)
      , start (clock ())
      {}
   ~Measure ()
      { chr. time += clock () - start; }
  };
};



struct Offset
{
private:
	static size_t size;
public:
	static const size_t delta;

	Offset ()
  	{ size += delta; }
 ~Offset ()
  	{ size -= delta; }

  static void newLn (ostream &os) 
    { os << endl << string (size, ' '); }
};



class ONumber
{
	ostream &o;
  const streamsize prec_old;
	const ios_base::fmtflags flags_old;
public:
	ONumber (ostream &o_arg,
	         streamsize precision,
	         bool scientific_arg)
	  : o (o_arg)
	  , prec_old (o. precision ())
	  , flags_old (o. flags ())
	  { if (scientific_arg)
	  	  o << scientific;
	  	else
	  		o << fixed;
      o. precision (precision);
	  }
 ~ONumber ()
    { o. flags (flags_old); 
      o. precision (prec_old); 
    }
};




void exec (const string &cmd);




struct Application : Singleton<Application>, Root
// Usage: int main (argc, argv) { Application app; return app. run (argc, argv); }
{  
private:
  // arg's are unique
  Vector<string> positionalArgs;
  Vector<string> positionalValues;
    // size() = positionalArgs.size()
    // Default: string()
  map<string/*arg*/,string/*value*/> keyArgs;
  map<string/*arg*/,bool/*value*/> flagArgs;
    // Default: false
public:
  

protected:
  Application ()
    { addKey ("log");
      addKey ("verbose", "0");
      addKey ("json");
      addFlag ("noprogress");
    }
    // Invoke: addKey(), addFlag(), addPositional()
  void addKey (const string &key, 
               const string &defValue = string ());
  void addFlag (const string &flag);
  void addPositional (const string &positional);
public:


private:
  void setKey (const string &key, 
               const string &value);
    // Update: keyArgs[key] = value
  void setFlag (const string &flag);
    // Update: flagArgs[flag] = true
  void setPositional (const string &value);
    // Update: positionalValues << value
protected:
  string getArg (const string &arg) const;
    // Input: keyArgs, positionalArgs, positionalValues
  bool getFlag (const string &flag) const;
    // Input: flagArgs
public:
  void instruction () const;
    // Invokes: ERROR_MSG
  int run (int argc, 
           const char* argv []);
    // Update: keyArgs, flagArgs, positionalArgs
    // Invokes: body()
private:
  virtual void body () const = 0;
};



}



#endif

