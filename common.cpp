// common.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"

#include <sstream>
#include <cstring>
#ifndef _MSC_VER
  #include <dirent.h>
#endif
#include <signal.h>  



namespace Common_sp
{
 


bool Chronometer::enabled = false;
bool qc_on = false;



namespace 
{

void segmFaultHandler (int /*sig_num*/)
{
  signal (SIGSEGV, SIG_DFL); 
  errorExit ("Segmentation fault", true); 
}

}



const vector<bool> Bool {false, true};



bool initCommon ()
{
  MODULE_INIT

  signal (SIGSEGV, segmFaultHandler);  
        
#ifdef _MSC_VER
  #pragma warning (disable : 4127)
#endif
  ASSERT (numeric_limits<char>::is_signed);
  ASSERT (sizeof (long) >= 4);
#ifdef _MSC_VER
  #pragma warning (default : 4127)
#endif

  ASSERT (SIZE_MAX == std::numeric_limits<size_t>::max ());

  return true;
}


namespace
{
  const bool init_ = initCommon ();
}




vector<string> programArgs;
string programName;
ostream* logPtr = nullptr;



void errorExit (const char* msg,
                bool segmFault)
// alloc() may not work
{ 
	ostream* os = logPtr ? logPtr : & cout; 
	
	// time ??
#ifndef _MSC_VER
	const char* hostname = getenv ("HOSTNAME");
	const char* pwd = getenv ("PWD");
#endif
	*os << endl
      << "ERROR" << endl
      << msg << endl
    #ifndef _MSC_VER
	    << "HOSTNAME: " << (hostname ? hostname : "?") << endl
	    << "PWD: " << (pwd ? pwd : "?") << endl
    #endif
	    << "Progam name: " << programName << endl
	    << "Command line:";
	 FOR (size_t, i, programArgs. size ())
	   *os << " " << programArgs [i];
	 *os << endl;
//system (("env >> " + logFName). c_str ());

  os->flush ();

  if (segmFault)
    abort ();
  exit (1);
}



namespace
{
	
uint powInt_ (uint a,
              uint b)
// Input: a: !0, != 1
{
	if (! b)
		return 1;
	if (b == 1)
		return a;
	const uint half = b / 2;
	const uint res = powInt_ (a, half);
	return res * res * (divisible (b, 2) ? 1 : a);
}
	
}


uint powInt (uint a,
             uint b)
{
	if (a)
		if (a == 1)
			return 1;
		else
			return powInt_ (a, b);
	else
		if (b)
			return 0;
		else
			throw runtime_error ("powInt: 0^0");
}




// string

const string noString;



bool isRight (const string &s,
              const string &right)
{
  if (s. size () < right. size ())
    return false;
  return s. substr (s. size () - right. size ()) == right;
}



bool trimPrefix (string &s,
                 const string &prefix)
{
  if (! isLeft (s, prefix))
    return false;
  s. erase (0, prefix. size ());
  return true;
}



bool trimSuffix (string &s,
                 const string &suffix)
{
  if (! isRight (s, suffix))
    return false;
  s. erase (s. size () - suffix. size ());
  return true;
}



void trimSuffixNonAlphaNum (string &s)
{
  for (;;)
  {
    trimTrailing (s);
    if (s. empty ())
      break;
    const size_t pos = s. size () - 1;
    if (   isAlpha (s [pos])
        || isDigit (s [pos])
       )
      break;
    s. erase (pos);
  }
}



bool trimTailAt (string &s,
                 const string &tailStart)
{
  const size_t pos = s. find (tailStart);
  const bool trimmed = pos != string::npos;
  if (trimmed)
    s. erase (pos);
  return trimmed;
}



bool goodName (const string &name)
{
  if (name. empty ())
    return false;
  if (name [0] == ' ')
    return false;
  if (*(name. end () - 1) == ' ')
    return false;

  for (const char c : name)
    if (! printable (c))
      return false;
      
  return true;
}



bool strBlank (const string &s)
{
  for (const char c : s)
    if (! isSpace (c))
      return false;
  return true;
}



void strUpper (string &s)
{
  for (char& c : s)
    c = toUpper (c);
}



void strLower (string &s)
{
  for (char& c : s)
    c = toLower (c);
}



bool isUpper (const string &s)
{
  for (const char c : s)
    if (! isUpper (c))
    	return false;
  return true;
}



bool isLower (const string &s)
{
  for (const char c : s)
    if (! isLower (c))
    	return false;
  return true;
}



string::const_iterator stringInSet (const string &s,
                    	           	  const string &charSet)
{
  CONST_ITER (string, it, s)
    if (! charInSet (*it, charSet))
      return it;

  return s. end ();
}



size_t strCountSet (const string &s,
 		                const string &charSet)
{
  size_t n = 0;
  for (const char c : s)
	  if (charInSet (c, charSet))		
	    n++;
  return n;
}



void strDeleteSet (string &s,
		               const string &charSet)
{
  FOR_REV (size_t, i, s. size ())
    if (charInSet (s [i], charSet))
    	s. erase (i, 1);
}



void trimLeading (string &s)
{
  size_t i = 0;
  for (; i < s. size () && isSpace (s [i]); i++)
    ;
  s. erase (0, i);
}



void trimTrailing (string &s)
{
	size_t i = s. size ();
	while (i)
		if (isSpace (s [i - 1]))
			i--;
		else
			break;
  s. erase (i);
}



size_t containsWord (const string& hay,
                     const string& needle)
{
  size_t pos = string::npos;
  for (;;)
  {
    pos = hay. find (needle, pos == string::npos ? 0 : pos + 1);
    if (pos == string::npos)
      break;
    const size_t pos1 = pos + needle. size ();
    if (   (pos  == 0            || ! isLetter (hay [pos - 1]))
        && (pos1 == hay. size () || ! isLetter (hay [pos1]))
       )
      break;
  }
  return pos;  
}



void replace (string &s,
              char from,
              char to)
{ 
  for (char& c : s)
	  if (c == from)
	  	c = to;
}



void replace (string &s,
              const string &fromChars,
              char to)
{
  for (char& c : s)
	  if (charInSet (c, fromChars))
	  	c = to;
}



void replaceStr (string &s,
                 const string &from,
                 const string &to)
{
  if (from == to)
    return;
    
  for (;;)
  {
    const size_t pos = s. find (from);
    if (pos == string::npos)
      break;
  #if 0
    s = s. substr (0, pos) + to + s. substr (pos + from. size ());
  #else
    s. replace (pos, from. size (), to);
  #endif
  }
}
  
  
  
string to_c (const string &s)
{
  string r;
  for (const char c : s)
    if (c == '\n')
      r += "\\n";
    else
    {
      if (charInSet (c, "\"\'\\"))
        r += '\\';
      r += c;
    }   
  return r;
}



void collapseSpace (string &s)
{
  string s1;
  do 
  {
    s1 = s;
    replaceStr (s, "  ", " ");
  }
  while (s1 != s);
}
               


  
string str2streamWord (const string &s,
                       size_t wordNum)
{
  istringstream iss (s);
  string word;
  FOR (size_t, i, wordNum + 1)
    if (iss. eof ())
    	return string ();
    else
      iss >> word;
  return word;
}
                


string str2sql (const string &s)
{
	string r = "'";
  for (const char c : s)
	{
	  if (c == '\'')
	  	r += "'";
	  r += c;
	}
	r += "'";
	
	return r;
}



string sql2escaped (const string &s)
{
  string s1;
  for (const char c : s)
  {
    if (charInSet (c, "[]*_%\\"))
      s1 += '\\';
    s1 += c;
  }
  
  return s1;
}



string findSplit (string &s,
                  char c)
{
	const size_t pos = s. find (c);
	if (pos == string::npos)
	{
		const string s1 (s);
		s. clear ();
		return s1;
	}
	const string before (s. substr (0, pos));
	s. erase (0, pos + 1);
	return before;
}



string rfindSplit (string &s,
                   char c)
{
	const size_t pos = s. rfind (c);
	if (pos == string::npos)
	{
		const string s1 (s);
		s. clear ();
		return s1;
	}
	const string after (s. substr (pos + 1));
	s. erase (pos);
	return after;
}



List<string> str2list (const string &s,
                       char c) 
{
	List<string> res;
	string s1 (s);
	while (! s1. empty ())
	  res << findSplit (s1, c);
	return res;
}



string list2str (const List<string> &strList,
                 const string &sep) 
{
	string s;
	for (const string& e : strList)
	{
		if (! s. empty ())
			s += sep;
	  s += e;
	}
	return s;
}



bool fileExists (const string &fName)
{
  ifstream f (fName. c_str ());
  bool ok = f. good ();
#ifndef _MSC_VER
  if (ok)
    ok = ! directoryExists (fName);
#endif

  return ok;
}



#ifndef _MSC_VER
bool directoryExists (const string &dirName)
{
  DIR* dir = opendir (dirName. c_str ());
  const bool yes = (bool) (dir);
  if (yes)
  {
    if (closedir (dir))
      ERROR;
  }
  return yes;
}
#endif



//

size_t strMonth2num (const string& month)
{
  if (isDigit (month [0]))
  {
    const int m = str2<int> (month);
    ASSERT (m >= 1);
    ASSERT (m <= 12);
    return (size_t) (m - 1);
  }
  
	size_t i = NO_INDEX;
	     if (month == "Jan")
		i = 0;
  else if (month == "Feb")
		i = 1;
  else if (month == "Mar")
		i = 2;
  else if (month == "Apr")
		i = 3;
  else if (month == "May")
		i = 4;
  else if (month == "Jun")
		i = 5;
  else if (month == "Jul")
		i = 6;
  else if (month == "Aug")
		i = 7;
  else if (month == "Sep")
		i = 8;
  else if (month == "Oct")
		i = 9;
  else if (month == "Nov")
		i = 10;
  else if (month == "Dec")
		i = 11;
  else 
  	ERROR;
  	
  return i;
}




// istream

bool getChar (istream &is,
              char &c)
{
  ASSERT (is. good ());

  const int i = is. get ();
  c = EOF;
  if (is. eof ())
  {
    ASSERT (i == c);
    return false;
  }
  ASSERT (i >= 0 && i <= 255);
  c = static_cast<char> (i);

  return true;
}



void skipLine (istream &is)
{
#if 1
  char c;
  while (getChar (is, c) && c != '\n')  // UNIX
    ;
#else
	char c = '\0';
	while (! is. eof () && c != '\n')  // UNIX
  {	
  	ASSERT (is. good ());
	  c = (char) is. get ();
	}
#endif
}



void readLine (istream &is,
               string &s)
{
  s. clear ();
  for (;;)
  {	
    char c;
    if (! getChar (is, c))
      break;
	  if (c == '\n')  // UNIX
	  	break;
	  s += c;
	}
}



string getToken (istream &is,
                 const string &skip,
                 const string &delimeters)
{
  // Skipping skip
  for (;;)
  {
    IMPLY (! is. eof (), is. good ());
    const char c = (char) is. get ();
    if (is. eof ())
      return string ();
    ASSERT (c);
    if (charInSet (c, skip))
      continue;
    if (charInSet (c, delimeters))
      return string (1, c);
    is. unget ();
    break;
  }  
  
  string token;
  for (;;)
  {
    IMPLY (! is. eof (), is. good ());
    const char c = (char) is. get ();
    if (   is. eof () 
        || charInSet (c, skip)
        || charInSet (c, delimeters)
       )
    {
      is. unget ();
      break;
    }
    ASSERT (c);
    token += c;
  }
  ASSERT (! token. empty ());

  return token;
}

 


// verbose

namespace
{
	int verbose_ = 0;
}


bool verbose (int inc)
{ 
	return verbose_ + inc > 0;
}



Verbose::Verbose (int verbose_arg)
: verbose_old (verbose_)
{
	verbose_ = verbose_arg;
}


Verbose::Verbose ()
: verbose_old (verbose_)
{
	verbose_++;
}


Verbose::~Verbose ()
{
	verbose_ = verbose_old;
}



Unverbose::Unverbose ()
{
  verbose_--;
}


Unverbose::~Unverbose ()
{
  verbose_++;
}




// Root

void Root::saveFile (const string &fName) const
{
	if (fName. empty ())
		return;
  
  OFStream f ("", fName, "");
  saveText (f);
}




// Rand

const long Rand::max_ = 2147483647;



ulong Rand::get (ulong max)
{
  ASSERT (max > 0);
  ASSERT (max <= (ulong) max_);
  
  const long up = max_ - (max_ % (long) max);  // Multiple of max
  do run ();
    while (seed >= up);  // P(stay in cycle) < 0.5
  return (ulong) (seed % (long) max);
}



void Rand::qc () const
{
  if (! qc_on)
    return;
  ASSERT (seed > 0);
  ASSERT (seed < max_);
}



void Rand::run ()
{
	qc ();
  
  const long ia = 16807;
  const long iq = 127773;
  const long ir = 2836;

  const long k = seed / iq;
  seed = ia * (seed - k * iq) - ir * k;
  if (seed < 0)
  	seed += max_;
}




// Named

Named::Named (const string& name_arg)
: name (name_arg) 
{
#ifndef NDEBUG
  if (! goodName (name))
    ERROR_MSG ("Bad name: \"" + name_arg + "\"");
#endif
}



void Named::qc () const
{
  if (! qc_on)
    return;
  Root::qc ();
    
  ASSERT (goodName (name));
}



// DisjointCluster

void DisjointCluster::merge (DisjointCluster &other)
// Trying to keep rank's small
{ 
  DisjointCluster* root1 =        getDisjointCluster ();
  DisjointCluster* root2 = other. getDisjointCluster ();
	if (root1->rankDC > root2->rankDC)
		root2->parentDC = root1;
	else
	{
		root1->parentDC = root2;
		if (root1->rankDC == root2->rankDC)
			root2->rankDC ++;
	}
}



DisjointCluster* DisjointCluster::getDisjointCluster ()
{ 
  ASSERT (parentDC);
	if (this == parentDC)
	  return this;
	parentDC = parentDC->getDisjointCluster ();
	return parentDC;
}




//////////////////////////////// DiGraph /////////////////////////////////////

// DiGraph::Node

void DiGraph::Node::attach (DiGraph &graph_arg) 
{
  ASSERT (! graph);
#ifndef NDEBUG
	for (const bool b : Bool)
	  ASSERT (arcs [b]. empty ());
#endif
  
  graph = & graph_arg;
  graph_arg. nodes << this;
  graphIt = graph_arg. nodes. end (); 
  graphIt--;
}

 

DiGraph::Node::~Node ()
{
	for (const bool b : Bool)
    deleteNeighborhood (b);
  if (graph)
    const_cast <DiGraph*> (graph) -> nodes. erase (graphIt);
}

 

void DiGraph::Node::qc () const
{
  if (! qc_on)
    return;
  Root::qc ();
    
  ASSERT (*graphIt == this);
}



void DiGraph::Node::saveText (ostream &os) const
{
  os << getName ();
  saveContent (os);
  
  if (orderDfs)
    os << "  DFS_order = " << orderDfs;
  if (scc)
  {
    os << "  SCC: " << scc->getName ();
  }

  os << endl;  
	for (const bool b : Bool)
  {
    os << "  " << (b ? "Out" : "In") << ":" << endl;
    for (const Arc* arc : arcs [b])
    {
      os << "    " << arc->node [b] -> getName ();
      os << ": ";
      arc->saveContent (os);
      os << endl;
    }
  }
}

 

bool DiGraph::Node::isIncident (const DiGraph::Node* n,
                                bool out) const
{
	ASSERT (n);
	ASSERT (n->graph == graph);

  for (const Arc* arc : arcs [out])
    if (arc->node [out] == n)
    	return true;
  return false;
}



VectorPtr<DiGraph::Node> DiGraph::Node::getNeighborhood (bool out) const
{
	VectorPtr<Node> s;  s. reserve (arcs [out]. size ());
  for (const Arc* arc : arcs [out])
    s << arc->node [out];
	return s;
}

 

void DiGraph::Node::deleteNeighborhood (bool out)
{
  for (Iter <List<Arc*> > iter (arcs [out]); iter. next (); )
    delete *iter;  // changes arcs[i]
}



DiGraph::Node* DiGraph::Node::setScc (size_t &visitedNum,
                                      stack<DiGraph::Node*, vector<DiGraph::Node*> > &sccStack)
// Tarjan's alogorithm:
/*
  Input: Graph G = (V, E), Start node v0
  
  index = 0                       // DFS node number counter 
  S = empty                       // An empty stack of nodes
  tarjan(v0)                      // Start a DFS at the start node
  
  procedure tarjan(v)
    v.index = index               // Set the depth index for v
    v.lowlink = index
    index = index + 1
    S.push(v)                     // Push v on the stack
    forall (v, v') in E do        // Consider successors of v 
      if (v'.index is undefined)  // Was successor v' visited? 
        tarjan(v')                // Recurse
        v.lowlink = min(v.lowlink, v'.lowlink)
      elseif (v' in S)            // Is v' on the stack?
        v.lowlink = min(v.lowlink, v'.lowlink)
    if (v.lowlink == v.index)     // Is v the root of an SCC?
      print "SCC:"
      repeat
        v' = S.pop
        print v'
      until (v' == v)
*/
{
  if (orderDfs)
  {
    if (inStack)
      return scc;
    else
      return nullptr;
  }
  ASSERT (! inStack);
  
  visitedNum++;
  orderDfs = visitedNum;

  sccStack. push (this);
  inStack = true;

  scc = this;

  for (Arc* arc : arcs [true])
    if (Node* lowNode = arc->node [true] -> setScc (visitedNum, sccStack))
      if (lowNode->orderDfs < scc->orderDfs)
        scc = lowNode;

  ASSERT (scc->orderDfs <= orderDfs);  
  if (scc->orderDfs < orderDfs)  
  {
    ASSERT (scc != this);
    return scc;
  }
  
  // this is the root of a SCC
  ASSERT (scc->orderDfs == orderDfs);  
  ASSERT (scc == this);
  for (;;)
  {
    Node* n = sccStack. top ();
    sccStack. pop ();
    n->inStack = false;
    if (n == this)
      break;
    n->scc = this;
  }

  return nullptr;
}

 

void DiGraph::Node::contract (Node* from)
{
  ASSERT (from);
  ASSERT (this != from);
  
  contractContent (from);

	for (const bool b : Bool)
  {  
    map <Node*, Arc*> m;
    for (Arc* arc : arcs [b])
      m [arc->node [b]] = arc;

  #if 1      
    for (Iter <List<Arc*> > iter (from->arcs [b]); iter. next (); )
    {
      Arc* arc = *iter;
      if (m. find (arc->node [b]) == m. end ())
      {
        iter. erase ();
        arc->node [! b] = this;
        arcs [b]. push_back (arc);
        arc->arcsIt [! b] = arcs [b]. end ();
        arc->arcsIt [! b] --;
      }
      else
        m [arc->node [b]] -> contractContent (arc);
    }
  #else
    List<Arc*>::iterator it = from->arcs [b]. begin ();
    while (it != from->arcs [b]. end ())
    {
      List<Arc*>::iterator next = it;
      next++;
      
      Arc* arc = *it;
      if (m. find (arc->node [b]) == m. end ())
      {
        from->arcs [b]. erase (it);
        arc->node [! b] = this;
        arcs [b]. push_back (arc);
        arc->arcsIt [! b] = arcs [b]. end ();
        arc->arcsIt [! b] --;
      }
      else
        m [arc->node [b]] -> contractContent (arc);
        
      it = next;
    }
  #endif
  }
  
  delete from;
}



void DiGraph::Node::remove ()
{
#ifndef NDEBUG
	for (const bool b : Bool)
	  ASSERT (arcs [b]. empty ());
#endif
  const_cast <DiGraph*> (graph) -> nodes. erase (graphIt);  
  graph = nullptr;
}


 

// DiGraph::Arc

void DiGraph::Arc::attach (Node* start,
                           Node* end) 
{      
  ASSERT (start);
  ASSERT (end);
  ASSERT (start->graph == end->graph);
  ASSERT (! node [false]);
  ASSERT (! node [true]);

  node [false] = start;
  node [true]  = end;

	for (const bool b : Bool)
  {
    node [b] -> arcs [! b]. push_back (this);
    arcsIt [b] = node [b] -> arcs [! b]. end ();
    arcsIt [b] --;
  }
}

 

DiGraph::Arc::~Arc ()
{
	for (const bool b : Bool)
    node [b] -> arcs [! b]. erase (arcsIt [b]);
}

 


// DiGraph

void DiGraph::init (const DiGraph &other,
	                  Old2new &old2new)
{
  ASSERT (old2new. empty ());

  for (const Node* node : other. nodes)
  {
  	Node* n = node->copy ();
  	n->attach (*this);
  	old2new [node] = n;
  } 
  ASSERT (old2new. size () == other. nodes. size ());


  for (const Node* node : other. nodes)
  {
  	Node* n = old2new [node];
  	ASSERT (n);
  	
  	EXEC_ASSERT (n->parentDC = old2new [static_cast <Node*> (n->parentDC)]);
  	
  	if (n->scc)
  	{
	  	n->scc = old2new [n->scc];
	  	ASSERT (n->scc);
	  }
  	
	  for (const Arc* arc : node->arcs [true])
	  {
	  	Arc* a = arc->copy ();
	  	a->attach (n, old2new [arc->node [true]]);
	  }
  }   
}



DiGraph::~DiGraph ()
{
  while (! nodes. empty ())
  {
    Node* n = nodes. front ();
    if (! n)
      errorExit ("DiGraph::Node is nullptr");
    delete n;
  }
}



void DiGraph::qc () const
{
  if (! qc_on)
    return;

#ifndef NDEBUG
  Set<const Node*> nodes_;
  Set<const Arc*> arcs_ [2];
  Set<string> names;
  size_t arcs = 0;
  for (const Node* node : nodes)
  {
    ASSERT (node);
    ASSERT (node->graph == this);
    nodes_ << node;
    node->qc ();
  	for (const bool b : Bool)
      for (const Arc* arc : node->arcs [b])
      {
        ASSERT (arc);
        arcs_ [b] << arc;
        arcs++;
        if (b)
        {
    	  	Unverbose unv;
          arc->qc ();
        }
      }
    if (names. contains (node->getName ()))
    {
      cout << "Duplicate name: " << node->getName () << endl;
      ERROR;  
    }
    names << node->getName ();
  }
  ASSERT (nodes. size () == nodes_. size ());
	for (const bool b : Bool)
    ASSERT (2 * arcs_ [b]. size () == arcs);
  ASSERT (arcs_ [false] == arcs_ [true]);
#endif
}

 

void DiGraph::saveText (ostream &os) const
{
  for (const Node* node : nodes)
  {
    node->saveText (os);
    os << endl;
  }
}

 

void DiGraph::connectedComponents ()
{
  for (Node* node : nodes)
	  node->DisjointCluster::init ();
  for (Node* node : nodes)
  	for (const bool b : Bool)
      for (Arc* arc : node->arcs [b])
	    	arc->node [b] -> DisjointCluster::merge (* arc->node [! b]);
}



void DiGraph::scc ()
{
  size_t visitedNum = 0;
  vector<Node*> vec (nodes. size ());
  vec. clear ();
  stack<Node*, vector<Node*> > sccStack (vec);
  for (Node* node : nodes)
  {
  //os << node->getName ();  
    if (! node->orderDfs)
      node->setScc (visitedNum, sccStack);
    ASSERT (sccStack. empty ());
  }
}

 

void DiGraph::contractScc ()
{
#if 1
  for (Iter <List<Node*> > iter (nodes); iter. next (); )
  {
    Node* n = *iter;
    ASSERT (n->scc);
    if (n->scc != n)
      n->scc->contract (n);
  }
#else
  List<Node*>::iterator it = nodes. begin ();
  while (it != nodes. end ())
  {
    List<Node*>::iterator next = it;
    next++;
    
    Node* n = *it;
    ASSERT (n->scc);
    if (n->scc != n)
      n->scc->contract (n);
      
    it = next;
  }
#endif
}



VectorPtr<DiGraph::Node> DiGraph::getEnds (bool out) const
{
  VectorPtr<Node> s;
  for (const Node* node : nodes)
    if (node->arcs [out]. empty ())
      s << node;
  return s;
}



DiGraph::Node2Node DiGraph::reverse (const Node2Node& old2new)
{
  Node2Node new2old;
  for (const auto it : old2new)
    new2old [it. second] = it. first;    
  return new2old;
}



void DiGraph::borrowArcs (const Node2Node &node2node,
                          bool parallelAllowed)
{
#ifndef NDEBUG
  const DiGraph* otherGraph = nullptr;
#endif
  for (const auto it : node2node)
  {
    const Node* other = it. first;
    const Node* from  = it. second;
    ASSERT (other);
    ASSERT (from);
  #ifndef NDEBUG
    if (otherGraph)
      { ASSERT (otherGraph == other->graph); }
    else
      otherGraph = other->graph;
  #endif
    ASSERT (this == from->graph);
    ASSERT (otherGraph != this);
    const VectorPtr<Node> otherNeighborhood (other->getNeighborhood (true));
    for (const Node* otherNeighbor : otherNeighborhood)
      if (const Node* to = findPtr (node2node, otherNeighbor))
        if (parallelAllowed || ! from->isIncident (to, true))
          new Arc ( const_cast <Node*> (from)
                  , const_cast <Node*> (to)
                  );
  }
}




// Tree::TreeNode

void Tree::TreeNode::qc () const
{
  if (! qc_on)
    return;
  DiGraph::Node::qc ();
    
  ASSERT (! (isLeafType () && isInteriorType ()));
  IMPLY (isLeaf (), isLeafType ());
  IMPLY (isInteriorType () && getParent (), getParent () -> isInteriorType ());
}
  


void Tree::TreeNode::saveText (ostream &os) const
{
  os << getName () << ": ";
  saveContent (os);

  
  bool saveSubtreeP = false;
	for (const Arc* arc : arcs [false])
	  if (static_cast <const TreeNode*> (arc->node [false]) -> getSaveSubtreeP ())
	  {
	    saveSubtreeP = true;
	    break;
	  }

  if (saveSubtreeP)
  {
  	Offset ofs;
  	for (const Arc* arc : arcs [false])
  	{ 
  		Offset::newLn (os);
  	  static_cast <TreeNode*> (arc->node [false]) -> saveText (os);
  	}
  }
  else
    if (! isLeaf ())
      os << " (" << getLeavesSize () << ")";
}



void Tree::TreeNode::printNewick_ (ostream &os,
	                                 bool internalNames,
	                                 bool minimalLeafName) const
{
  // Cf. saveText() ??

	if (isLeaf ())
		os << name2newick (getNewickName (minimalLeafName));
	else
	{
		os << "(";
		bool first = true;
		for (const Arc* arc : arcs [false])
		{
			const TreeNode* n = static_cast <TreeNode*> (arc->node [false]);
			if (! first)
			  os << ",";
			n->printNewick_ (os, internalNames, minimalLeafName);
			first = false;
		}
	  
		os << ")";
		if (internalNames)
			os << name2newick (getNewickName (minimalLeafName));
	}
	
	const double dist = getParentDistance ();
	if (dist == dist && dist != -1)
		os << ":" << fixed << dist;
}



string Tree::TreeNode::name2newick (const string &s) 
{
  string s1 (s);
  replace (s1, "\"\' ():;,[]<>=", '_');
  return s1. substr (0, 100 /*50*/ /*PAR*/);
    // Newick allowes the name length <= 50 characters, but http://www.trex.uqam.ca/ breaks with 50
}



const Tree::TreeNode* Tree::TreeNode::getSuperParent (size_t height) const
{
  const TreeNode* n = this;
  FOR (size_t, i, height)
  {
  	ASSERT (n);
    n = n->getParent ();
  }
  return n;
}



void Tree::TreeNode::setParent (TreeNode* newParent)
{ 
	ASSERT (newParent != this);
	
	if (! arcs [true]. empty ())
	{
		Arc* a = arcs [true]. front ();
		delete a;
  }
  ASSERT (arcs [true]. empty ());

	if (newParent)
	{
	  new Arc (this, newParent);
	  ASSERT (getParent () == newParent);
	}
	else
		const_cast <Tree&> (getTree ()). root = this;
}



Tree::TreeNode::TipName Tree::TreeNode::getTipName () const
{ 
  if (isLeaf ()) 
    return TipName (getName (), 0);
    
  TipName tn_best;
	for (const DiGraph::Arc* arc : arcs [false])
	{
	  const TipName tn = static_cast <Tree::TreeNode*> (arc->node [false]) -> getTipName ();
	  if (tn_best. name. empty () || tn_best. name > tn. name)
	    tn_best = tn;
	}
	tn_best.depth ++;

	return tn_best;
}



size_t Tree::TreeNode::getHeight () const
{
  size_t n = 0;
	for (const Arc* arc : arcs [false])
	  maximize (n, 1 + static_cast <const TreeNode*> (arc->node [false]) -> getHeight ());
	return n;
}



size_t Tree::TreeNode::getInteriorHeight () const
{
  ASSERT (isInteriorType ());
  
  size_t n = 0;
	for (const Arc* arc : arcs [false])
	{
	  const TreeNode* node = static_cast <const TreeNode*> (arc->node [false]);
	  if (node->isInteriorType ())
	    maximize (n, 1 + node->getInteriorHeight ());
	}
	return n;
}



void Tree::TreeNode::getBifurcatingInteriorBranching (size_t &bifurcatingInteriorNodes,
                                                      size_t &branches) const
{
  ASSERT (isInteriorType ());
  ASSERT (/*a*/ arcs [false]. size () >= 2);  // transient nodes are prohibited

  // Make *this bifurcating
  branches                 += arcs [false]. size () - 2;  // b
  bifurcatingInteriorNodes += arcs [false]. size () - 2;  // n

  size_t interiorArcs = 0;  // i
	for (const Arc* arc : arcs [false])
	{
	  const TreeNode* node = static_cast <const TreeNode*> (arc->node [false]);
	  if (node->isInteriorType ())
	  {
	    interiorArcs++;
	    node->getBifurcatingInteriorBranching (bifurcatingInteriorNodes, branches);
	  }
	}  
	
	branches += interiorArcs;
	if (interiorArcs)
	  bifurcatingInteriorNodes++;
	ASSERT (branches >= bifurcatingInteriorNodes);

  ASSERT (branches <= 2 * bifurcatingInteriorNodes);  
/*Proof:
  a >= 2.
  i <= a.
  For each parent node:
  b = a - 2 + i.
  n = a - 2 + (bool) i.
  x := 2 n - b = a - 2 + 2 (bool) i - i = (a - i) + 2 ((bool) i - 1).
  if i > 0 then x = a - i >= 0.
  if i = 0 then x = a - 2 >= 0.
*/
}



size_t Tree::TreeNode::getSubtreeSize (bool countLeaves) const
{
	size_t n = 0;
	for (const Arc* arc : arcs [false])
	{
	  const TreeNode* child = static_cast <TreeNode*> (arc->node [false]);
    if (! countLeaves && child->arcs [false]. empty ())
      continue;
	  n += 1 + child->getSubtreeSize (countLeaves);
	}
	return n;
}



double Tree::TreeNode::getSubtreeLength () const
{
	double len = 0;
	for (const Arc* arc : arcs [false])
	{
	  const TreeNode* node = static_cast <TreeNode*> (arc->node [false]);
	  len += node->getParentDistance () + node->getSubtreeLength ();
	}
	return len;
}



void Tree::TreeNode::setLeaves () 
{
	leaves = 0;
	for (const Arc* arc : arcs [false])
	{
	  TreeNode* child = static_cast <TreeNode*> (arc->node [false]);
	  child->setLeaves ();
	  ASSERT (child->leaves);
	  leaves += child->leaves;
	}
	if (! leaves)
	  leaves = 1;
}



size_t Tree::TreeNode::getLeavesSize () const
{
	size_t n = 0;
	for (const Arc* arc : arcs [false])
	  n += static_cast <TreeNode*> (arc->node [false]) -> getLeavesSize ();
	return max<size_t> (n, 1);
}



namespace
{
  struct ChildFreq
  {
    Tree::TreeNode* node;
    size_t freq;
    ChildFreq (Tree::TreeNode* node_arg,
               size_t freq_arg)
      : node (node_arg)
      , freq (freq_arg)
      { ASSERT (node);
        ASSERT (freq);
      }
    static bool compare (const ChildFreq& a,
                         const ChildFreq& b)
      { LESS_PART (b, a, freq);
        LESS_PART (a, b, node);
        return false;
      }
  };
}



void Tree::TreeNode::children2frequentChild (double rareProb)
{
  ASSERT (rareProb >= 0);
  ASSERT (rareProb < 0.5);  // => in a bifurcating node at least one child is frequentChild
  ASSERT (frequentChild);  
  
  Vector<ChildFreq> childFreqs;
  {  
    const VectorPtr<DiGraph::Node> children (getChildren ());
    childFreqs. reserve (children. size ());
    for (const DiGraph::Node* child : children)
    {
      const TreeNode* node = static_cast <const TreeNode*> (child);
      childFreqs << ChildFreq (const_cast <TreeNode*> (node), node->getLeavesSize ());  // Time = O(|nodes|^2) ??
    }
  }
  Common_sp::sort (childFreqs, ChildFreq::compare);  // Bifurcatization
  
  size_t sum = 0;
#ifndef NDEBUG
  size_t freq_prev = numeric_limits<size_t>::max ();
#endif
  for (ChildFreq& cf : childFreqs)
  {
  #ifndef NDEBUG
    ASSERT (freq_prev >= cf. freq);
    freq_prev = cf. freq;
  #endif
    sum += cf. freq;
    ASSERT (sum);
    if ((double) cf. freq / (double) sum < rareProb)
      continue;
    cf. node->frequentChild = true;      
    cf. node->children2frequentChild (rareProb);    
  }
}
  


void Tree::TreeNode::getLeaves (VectorPtr<TreeNode> &leafVec) const
{
  if (arcs [false]. empty ())
    leafVec << this;
  else
  	for (const Arc* arc : arcs [false])
  	  static_cast <TreeNode*> (arc->node [false]) -> getLeaves (leafVec);
}



const Tree::TreeNode* Tree::TreeNode::getClosestLeaf (size_t &leafDepth) const
{
	const TreeNode* leaf = nullptr;
  leafDepth = SIZE_MAX;
	for (const Arc* arc : arcs [false])
	{
		size_t depth1;
		const TreeNode* leaf1 = static_cast <TreeNode*> (arc->node [false]) -> getClosestLeaf (depth1);
		if (minimize (leafDepth, depth1))
			leaf = leaf1;
	}
	
	if (leaf)
		leafDepth++;
	else
	{
		leaf = this;
		leafDepth = 0;
	}

	ASSERT (! leafDepth == (leaf == this));
		
	return leaf;
}



const Tree::TreeNode* Tree::TreeNode::getOtherChild (const TreeNode* child) const
{
  ASSERT (child);
  ASSERT (child->getParent () == this);

  const TreeNode* otherChild = nullptr;
	for (const Arc* arc : arcs [false])
	{
		const TreeNode* n = static_cast <TreeNode*> (arc->node [false]);
	  if (n != child)
	  {
	  	ASSERT (! otherChild);
	  	otherChild = n;
	  }
	}
  return otherChild;
}



const Tree::TreeNode* Tree::TreeNode::getLeftmostDescendent () const
{
  const TreeNode* n = this;
  while (! n->isLeaf ())
    n = static_cast <TreeNode*> (n->arcs [false]. front () -> node [false]);
  return n;
}



const Tree::TreeNode* Tree::TreeNode::getRightmostDescendent () const
{
  const TreeNode* n = this;
  while (! n->isLeaf ())
    n = static_cast <TreeNode*> (n->arcs [false]. back () -> node [false]);
  return n;
}



void Tree::TreeNode::childrenUp ()
{
  const VectorPtr<DiGraph::Node> children (getChildren ());
	for (const DiGraph::Node* node : children)
	{	
		TreeNode* n = const_static_cast <TreeNode*> (node);
		n->setParent (const_cast <TreeNode*> (getParent ()));  
	}
	ASSERT (arcs [false]. empty ());
}



void Tree::TreeNode::isolateChildrenUp ()
{ 
  childrenUp ();
	if (! arcs [true]. empty ())
	{ 
		Arc* a = arcs [true]. front ();
		delete a;
  }
  remove ();
}



void Tree::TreeNode::deleteSubtree ()
{
	for (const DiGraph::Arc* arc : arcs [false])
		static_cast <TreeNode*> (arc->node [false]) -> deleteSubtree ();
	while (! arcs [false]. empty ())
	{
		TreeNode* n = static_cast <TreeNode*> (arcs [false]. front() -> node [false]);
		delete n;
	}
}



const Tree::TreeNode* Tree::TreeNode::makeRoot ()
{
	const TreeNode* root_old = getTree (). root;
	ASSERT (root_old);
	
	TreeNode* parent_new = nullptr;
	TreeNode* node = this;
	while (node)
	{
		TreeNode* parent_old = const_cast <TreeNode*> (node->getParent ());
	  node->setParent (parent_new);
	  parent_new = node;
	  node = parent_old;
	}
	ASSERT (arcs [true]. empty ());
	ASSERT (this == getTree (). root);
	ASSERT (root_old != getTree (). root);
	
	return root_old;
}



void Tree::TreeNode::getArea_ (uint radius,
                               const Tree::TreeNode* prev,
                               VectorPtr<Tree::TreeNode> &area,
                               VectorPtr<Tree::TreeNode> &boundary) const
{
  area << this;

  size_t degree = (size_t) (prev ? 1 : 0);  // Degree of *this in area
  if (radius)
  {
    const TreeNode* parent_ = getParent ();
    if (parent_ && parent_ != prev)
    {
      parent_->getArea_ (radius - 1, this, area, boundary);
      degree++;
    }
    for (const Arc* arc : arcs [false])
    {
      const TreeNode* child = static_cast <TreeNode*> (arc->node [false]);
      if (child != prev)
      {
        child->getArea_ (radius - 1, this, area, boundary);
        degree++;
      }
    }
  }
  
  if (degree <= 1)
    boundary << this;
}




// Tree

void Tree::qc () const
{
  if (! qc_on)
    return;
	DiGraph::qc ();

#ifndef NDEBUG		
  Set<string> names;
  bool transient = false;
	for (const DiGraph::Node* node : nodes)
	{
	  if (node->arcs [true]. size () > 1)
	  {
	    cout << "Multiple parents of " << node->getName () << endl;
	    const VectorPtr<DiGraph::Node> neighbors (node->getNeighborhood (true));
	    ASSERT (neighbors. size () >= 2);
	    for (const DiGraph::Node* neighbor : neighbors)
	      cout << " " << neighbor->getName ();
	    cout << endl << endl;
	    if (verbose ())
	      print (cout);
	    ERROR;
	  }
	  const TreeNode* n = static_cast <const TreeNode*> (node);
	  if (n->isLeaf ())
	  {
	    const string newickName (TreeNode::name2newick (n->getNewickName (true)));
      if (names. contains (newickName))
      {
        cout << "Duplicate name: " << newickName << endl;
        ERROR;
      }
      names << newickName;
    }
    if (n->isTransient ())
      transient = true;
	}
#endif
	ASSERT (! root == nodes. empty ());
	IMPLY (root, getRoot (true) == root);
	IMPLY (root, (nodes. size () > 1) == ! root->isLeaf ());
	
#ifndef NDEBUG		
	IMPLY (! transient, nodes. size () <= 2 * root->getLeavesSize () - 1);
#endif
}



namespace 
{
  typedef  map<size_t, string> AsnFeatures;
  void printAsnFeatures (ostream &os,
                         const AsnFeatures &features)
  {
    if (features. empty ())
      return;
    os << ",\n\
      features {\n";
    bool first = true;
    for (const auto it : features)
    {
      if (! first)
        os << ',';
      os << "\
        {\n\
          featureid " << it. first << ",\n\
          value \"" << it. second << "\"\n\
        }";
      first = false;
    }
    os << "\n\
      }";
  }
}



void Tree::printAsn (ostream &os) const
{
  map<const DiGraph::Node*, size_t/*index*/> node2index;
  size_t index = 0;
  for (const DiGraph::Node* n : nodes)
  {
    node2index [n] = index;
    index++;
  }    
  ASSERT (node2index. size () == nodes. size ());

  os << "BioTreeContainer ::= {fdict {\n\
    {\n\
      id 0,\n\
      name \"label\"\n\
    },\n\
    {\n\
      id 1,\n\
      name \"dist\"\n\
    }\n\
  },\n\
  nodes {";

  bool first = true;
  for (const DiGraph::Node* n : nodes)
  {
    const TreeNode* node = static_cast <const TreeNode*> (n);
    const TreeNode* parent = node->getParent ();
    if (! first)
      os << ',';
    os << "\n\
    {\n\
      id " << node2index [n];
    if (parent)
      os << ",\n\
      parent " << node2index [parent];
    AsnFeatures features;
    if (node->isLeaf ())
      features [0] = n->getName ();
    if (parent)
    {
      ostringstream oss;
      oss << node->getParentDistance ();
      features [1] = oss. str ();
    }
    printAsnFeatures (os, features);
    os << "\n\
    }";
    first = false;
  }
  os << "\n\
  }\n\
}\n\
";
}



void Tree::printArcLengths (ostream &os) const
{
  for (const DiGraph::Node* n : nodes)
  {
    const TreeNode* node = static_cast <const TreeNode*> (n);
    if (n == root)
      continue;
  	const double dist = node->getParentDistance ();
  	if (dist == dist && dist != -1 && dist > 0)
  	{
  		os << node->getLcaName () << " " << dist << " " << node->getRootDistance ();
  		double distPar = 0;
      if (const TreeNode* parent = node->getParent ())
        if (parent != root)
        {
        	distPar = parent->getParentDistance ();
        	if (distPar == distPar && distPar != -1)
        	{
        	  ASSERT (distPar > 0);
        	  os << " " << log (distPar) - log (dist);
        	}
        }
      if (! distPar)
     	  os << " ?";
      os << endl;
    }
  }  
}



double Tree::getAveArcLength () const
{
  double len = 0;
  size_t n = 0;
  for (const DiGraph::Node* node : nodes)
  {
    if (node == root)
      continue;
    const double arcLen = static_cast <const TreeNode*> (node) -> getParentDistance ();
    ASSERT (arcLen >= 0);
    len += arcLen;
    n++;
  }
  return len / (double) n;
}



Tree::Patristic::Patristic (const TreeNode* leaf1_arg, 
                            const TreeNode* leaf2_arg,
                            double distance_arg)
: leaf1 (leaf1_arg)
, leaf2 (leaf2_arg)
, distance (distance_arg)
{
  ASSERT (leaf1);
  ASSERT (leaf2);
  ASSERT (distance == distance);  // != NAN
  ASSERT (leaf1->graph == leaf2->graph);
  ASSERT (leaf1->getName () != leaf2->getName ());
  if (leaf1->getName () > leaf2->getName ())
    swap (leaf1, leaf2);
}



namespace 
{

typedef  map <const Tree::TreeNode*, double>  Leaf2dist;
 


Vector<Tree::Patristic> node2leafDistances (const Tree::TreeNode* node,
                                            Leaf2dist &leaf2dist) 
// Output: leaf2dist
{
  ASSERT (node);
  ASSERT (leaf2dist. empty ());
  
  Vector<Tree::Patristic> res;
  if (node->isLeaf ())
    leaf2dist [node] = 0;
  else
		for (const Tree::Arc* arc : node->arcs [false])
		{
			const Tree::TreeNode* n = static_cast <Tree::TreeNode*> (arc->node [false]);
			Leaf2dist nodeLeaf2dist;
			res << node2leafDistances (n, nodeLeaf2dist);
			ASSERT (! nodeLeaf2dist. empty ());
			const double dist = n->getParentDistance ();
			ASSERT (dist == dist);  // != NAN
			for (auto it : nodeLeaf2dist)
			  nodeLeaf2dist [it. first] += dist;
		  for (const auto it1 : leaf2dist)
			  for (const auto it2 : nodeLeaf2dist)
			    res << Tree::Patristic (it1. first, it2. first, it1. second + it2. second);
			for (const auto it : nodeLeaf2dist)
			  leaf2dist [it. first] = it. second;
	  }
    
  return res;
}

}



Vector<Tree::Patristic> Tree::getLeafDistances () const
{
  Leaf2dist leaf2dist;
  return node2leafDistances (root, leaf2dist);
}



size_t Tree::countInteriorNodes () const
{ 
  size_t n = 0;
  for (const DiGraph::Node* node : nodes)
    if (static_cast <const TreeNode*> (node) -> isInteriorType ())
      n++;
  return n;
} 



double Tree::getBifurcatingInteriorBranching () const
{ 
  if (! (root && root->isInteriorType ()))
    return -1;
    
  size_t bifurcatingInteriorNodes = 0;
  size_t branches = 0;
  root->getBifurcatingInteriorBranching (bifurcatingInteriorNodes, branches);
  const double branching = (double) branches / (double) bifurcatingInteriorNodes;
  ASSERT (branching >= 1);
  ASSERT (branching <= 2);
  return branching;
}



size_t Tree::countInteriorUndirectedArcs () const
{ 
  size_t n = countInteriorNodes ();
  if (root->isInteriorType ())
    n--;

  const VectorPtr<DiGraph::Node> children (root->getChildren ());
  switch (children. size ())
  {
    case 1: if (static_cast <const TreeNode*> (children [0]) -> isInteriorType ())
              n--;
            break;
    case 2: // root is transient in the undirected tree
      {
        size_t m = 0;
        for (const DiGraph::Node* child : children)
          if (static_cast <const TreeNode*> (child) -> isInteriorType ())
            m++;
        if (m)
          n--;
      }
  }

  return  n;
}



namespace 
{

bool getParentsOrTarget (const Tree::TreeNode* from,
			                   const Tree::TreeNode* target,
			                   VectorPtr<Tree::TreeNode> &parents) 
// Return: true <=> from->descendentOf(target)
{
	ASSERT (from)
	ASSERT (target)
  ASSERT (from->graph == target->graph);
  ASSERT (parents. empty ());
  
	while (from)
	{
		if (from == target)
			return true;
		parents << from;
		from = from->getParent ();
	}
	
	return false;
}

}



const Tree::TreeNode* Tree::getLowestCommonAncestor (const TreeNode* n1,
	                                                   const TreeNode* n2) 
{
  IMPLY (n1 && n2, n1->graph == n2->graph);
  
  if (   ! n1 
  	  || ! n2
  	 )
  	return nullptr;
  	
  if (n1 == n2)
    return n1;
	
	static VectorPtr<TreeNode> vec1;  
	vec1. clear ();
	if (getParentsOrTarget (n1, n2, vec1))
		return n2;
	ASSERT (! vec1. empty ());

	static VectorPtr<TreeNode> vec2; 
	vec2. clear ();
	if (getParentsOrTarget (n2, n1, vec2))
		return n1;
	ASSERT (! vec2. empty ());
	
	const TreeNode* m = nullptr;
	size_t i1 = vec1. size () - 1;
	size_t i2 = vec2. size () - 1;
	ASSERT (vec1 [i1] == vec2 [i2]);
	while (vec1 [i1] == vec2 [i2])
	{
		m = vec1 [i1];
		ASSERT (i1);
		ASSERT (i2);
		i1--;
		i2--;
	}
	ASSERT (m);
	return m;
}



const Tree::TreeNode* Tree::getLowestCommonAncestor (const VectorPtr<TreeNode> &nodeVec) 
{
  if (nodeVec. empty ())
    return nullptr;
    
	const TreeNode* n = nodeVec [0];
	FOR_START (size_t, i, 1, nodeVec. size ())
    n = getLowestCommonAncestor (n, nodeVec [i]);
	return n;
}



Set<const Tree::TreeNode*> Tree::getParents (const VectorPtr<TreeNode> &nodeVec) 
{
	Set<const TreeNode*> s;  
  const TreeNode* lca = getLowestCommonAncestor (nodeVec);
  for (const TreeNode* n : nodeVec)
	  while (n != lca)
	  {
	  	ASSERT (n);
		  s << n;
		  n = n->getParent ();
		}

	return s;
}



VectorPtr<Tree::TreeNode> Tree::getPath (const TreeNode* n1,
                                         const TreeNode* n2)
{ 
  ASSERT (n1);
  ASSERT (n2);
  ASSERT (n1->graph == n2->graph);
  
  if (n1 == n2)
    return VectorPtr<Tree::TreeNode> ();
  
	static VectorPtr<TreeNode> vec1;  
	vec1. clear ();
	if (getParentsOrTarget (n1, n2, vec1))
		return vec1;
	ASSERT (! vec1. empty ());

	static VectorPtr<TreeNode> vec2; 
	vec2. clear ();
	if (getParentsOrTarget (n2, n1, vec2))
		return vec2;
	ASSERT (! vec2. empty ());

	size_t i1 = vec1. size () - 1;
	size_t i2 = vec2. size () - 1;
	ASSERT (vec1 [i1] == vec2 [i2]);
	while (vec1 [i1] == vec2 [i2])
	{
    vec1. pop_back ();
    vec2. pop_back ();
		ASSERT (i1);
		ASSERT (i2);
		i1--;
		i2--;
	}

  VectorPtr<TreeNode> res;
  res. reserve (vec1. size () + vec2. size ());
  res << vec1 << vec2;
  
  return res;
}



void Tree::setFrequentChild (double rareProb)
{ 
  if (! root)
    return;
    
  for (DiGraph::Node* node : nodes)
    static_cast <TreeNode*> (node) -> frequentChild = false;
  const_cast <TreeNode*> (root) -> frequentChild = true;
  const_cast <TreeNode*> (root) -> children2frequentChild (rareProb);
  
  VectorPtr<TreeNode> transients;  transients. reserve (nodes. size ());
  for (const DiGraph::Node* node : nodes)
  {
    const TreeNode* tn = static_cast <const TreeNode*> (node);
    if (! tn->frequentChild)
      continue;
    const VectorPtr<DiGraph::Node> children (tn->getChildren ());
    size_t freqChildren = 0;
    for (const DiGraph::Node* child : children)
      if (static_cast <const TreeNode*> (child) -> frequentChild)
        freqChildren++;
    if (freqChildren == 1)  // transient, not leaf
      transients << tn;
  }  
  for (const TreeNode* tn : transients)
    const_cast <TreeNode*> (tn) -> frequentChild = false;
}



void Tree::setFrequentDegree (double rareProb)
{ 
  ASSERT (rareProb >= 0);
  ASSERT (rareProb < 0.3);
  
  if (! root)
    return;
    
  const size_t allLeaves = root->getLeavesSize ();
#if 0
  Binomial bin;
  bin. setParam ((int) allLeaves, rareProb);
  const double pValue = 0.01;  // PAR
#endif
  for (DiGraph::Node* node : nodes)
  {
    size_t degree = 0;
    {
      size_t sum = 0;
      const VectorPtr<DiGraph::Node> children (node->getChildren ());
      for (const DiGraph::Node* child : children)
      {
        const TreeNode* tn = static_cast <const TreeNode*> (child);
        const size_t leaves = tn->getLeavesSize ();  // Time = O(|nodes|^2) ??
        ASSERT (leaves);
        if ((double) leaves / (double) allLeaves >= rareProb)
      //if (1 - bin. cdf ((int) leaves - 1) <= pValue)
          degree++;
        sum += leaves;
      }
      if (node != root)
      {
        ASSERT (allLeaves > sum);
        const size_t leaves = allLeaves - sum;
        if ((double) leaves / (double) allLeaves >= rareProb)
      //if (1 - bin. cdf ((int) leaves - 1) <= pValue)
          degree++;
      }
    }
    static_cast <TreeNode*> (node) -> frequentDegree = degree;
  }

  // Leaves
  for (const DiGraph::Node* node : nodes)
  {
    const TreeNode* tn = static_cast <const TreeNode*> (node);
    if (! tn->isLeafType ())
      continue;
    ASSERT (tn->frequentDegree == 1);
    const TreeNode* parent = tn;
    while (parent && parent->frequentDegree == 1)
      parent = parent->getParent ();
    if (parent && parent->frequentDegree < 3)  // 0 or 2
      const_cast <TreeNode*> (tn) -> frequentDegree = 0;
  }  
}



void Tree::setRoot ()
{
  root = nullptr;
  for (const DiGraph::Node* node : nodes)
    if (! static_cast <const TreeNode*> (node) -> getParent ())
    {
      ASSERT (! root);
      const_static_cast <TreeNode*> (node) -> setParent (nullptr);
    }
  IMPLY (! nodes. empty (), root);
}



size_t Tree::deleteTransients ()
{
	size_t n = 0;
 	for (List<DiGraph::Node*>::const_iterator it = nodes. begin (); 
 		   it != nodes. end ();
 		  )
 	{
 		TreeNode* node = static_cast <TreeNode*> (*it);
 		it++;
 		if (node->deleteTransient ())
      n++;
  }  
  return n;
}



size_t Tree::restrictLeaves (const Set<string> &leafNames,
                             bool deleteTransientAncestor)
{
  size_t n = 0;
 	Vector<DiGraph::Node*> nodeVec;  nodeVec. reserve (nodes. size ());
 	insertAll (nodeVec, nodes);
 	for (DiGraph::Node* node_ : nodeVec)  
 	{
 	  const TreeNode* node = static_cast <const TreeNode*> (node_);
    if (node->isLeafType ())
      if (! leafNames. contains (node->getName ()))
      {
        deleteLeaf (const_cast <TreeNode*> (node), deleteTransientAncestor);
        n++;
      }
  }
  return n;
}



bool Tree::compare_std (const DiGraph::Node* a,
	                      const DiGraph::Node* b)
{
	ASSERT (a);
	ASSERT (b);
	ASSERT (a->graph == b->graph);
	const TreeNode* a_ = static_cast <const TreeNode*> (a);
	const TreeNode* b_ = static_cast <const TreeNode*> (b);

	LESS_PART (*b_, *a_, getLeavesSize ());
  LESS_PART (*a_, *b_, getLeftmostDescendent  () -> getName ());
  LESS_PART (*a_, *b_, getRightmostDescendent () -> getName ());

  return false;
}




// Progress

size_t Progress::beingUsed = 0;



Progress::Start::Start (AutoPtr<Progress> &prog_arg,
									      uint n_max,
								        uint displayPeriod)
: prog (prog_arg)
{
	ASSERT (! prog. get ());
	prog. reset (new Progress (n_max, displayPeriod));
}




// Input

Input::Input (const string &fName,
	            size_t bufSize,
	            uint displayPeriod)
: buf (new char [bufSize])
, ifs (fName. c_str ())
, eof (false)
, lineNum (0)
, prog (0, displayPeriod)  
{ 
  if (! ifs. good ())
    ERROR_MSG ("Bad file: " + fName);
  EXEC_ASSERT (ifs. rdbuf () -> pubsetbuf (buf. get (), (long) bufSize));   
}
 



// LineInput

bool LineInput::nextLine ()
{ 
	if (eof)  
	{ 
		line. clear ();
	  return false;
	}
	
	readLine (ifs, line);
  trimTrailing (line); 
	eof = ifs. eof ();
	lineNum++;

	const bool end = line. empty () && eof;
	if (! end)
		prog ();
		
	return ! end;
}




// ObjectInput

bool ObjectInput::next (Root &row)
{ 
	row. clear ();

	if (eof)
	  return false;

	row. read (ifs);
	lineNum++;

 	eof = ifs. eof ();
  if (eof)
  {
  	ASSERT (row. empty ());
  	return false;
  }

	prog ();
	
  ASSERT (ifs. peek () == '\n');
	row. qc ();

  skipLine (ifs);

	return true;
}



// CharInput

char CharInput::get ()
{ 
  ungot = false;
  
	const char c = (char) ifs. get ();

	eof = ifs. eof ();
	ASSERT (eof == (c == (char) EOF));

	if (eol)
	{ 
	  lineNum++;
		charNum = 0;
		prog ();
  }
  else
	  charNum++;

	eol = (eof || c == '\n');

	return /*eof ? (char) EOF :*/ c;
}



void CharInput::unget ()
{ 
	ASSERT (! ungot);
	
	ifs. unget (); 
	charNum--;
	ungot = true;
}



string CharInput::getLine ()
{
  string s;
  while (! eof)
  {
    const char c = get ();
    if (eol)
      break;
    s += c;
  }
  return s;
}




// Token

void Token::readInput (CharInput &in)
{
	clear ();

  // Skip spaces
	char c = '\0';
	do { c = in. get (); }
	  while (! in. eof && isSpace (c));
	if (in. eof)
	{
	  ASSERT (empty ());
		return;  
  }
		
  charNum = in. charNum;
	if (c == quote)
	{
		type = eText;
		for (;;)
		{ 
			c = in. get (); 
			if (in. eol)
				throw CharInput::Error (in, "ending quote");
			if (c == quote)
				break;
			name += c;
		}
	}
	else if (isDigit (c))
	{
		type = eNumber;
		while (! in. eof && isDigit (c))
		{ 
			name += c;
			c = in. get (); 
		}
		if (! in. eof)
			in. unget ();
		num = str2<uint> (name);
		ASSERT (num != numeric_limits<uint>::max ());
	}
	else if (isLetter (c))
	{
		type = eName;
		while (! in. eof && isLetter (c))
		{ 
			name += c;
			c = in. get (); 
		}
		if (! in. eof)
			in. unget ();
	}
	else 
	{
	  ASSERT (type == eDelimiter);
		name = c;
	}	
	ASSERT (! empty ());
}



void Token::qc () const
{
  if (! qc_on)
    return;
  if (! empty ())
  {
  	ASSERT (! name. empty ());
    ASSERT (! contains (name, quote));
  	IMPLY (type != eText, ! contains (name, ' '));
  	IMPLY (type != eNumber, num == noNum);
  	const char c = name [0];
  	switch (type)
  	{ 
  	  case eText:      break;
  		case eNumber:    ASSERT (isDigit (c)); 
  		                 ASSERT (num != noNum);
  		                 break;
  		case eName:      ASSERT (isLetter (c) && ! isDigit (c)); 
  		                 break;
  		case eDelimiter: ASSERT (name. size () == 1); 
  		                 break;
  		default: throw runtime_error ("Unknown type");
  	}
  }
}




// OFStream

void OFStream::open (const string &dirName,
					           const string &pathName,
					           const string &extension)
{ 
	ASSERT (! is_open ());	
	ASSERT (! pathName. empty ());
	
	string name;
	if (! dirName. empty ())
	  name = dirName + "/";
	name += pathName;
	if (! extension. empty ())
		name += "." + extension;
	
	ofstream::open (name. c_str ());
/*	
	if (eof ())
		cout << "eof" << endl;
	if (fail ())
		cout << "fail" << endl;
	if (bad ())
		cout << "bad" << endl;
*/	
	ASSERT (good ());	
}




// Csv
 
string Csv::getWord ()
{ 
  ASSERT (goodPos ());
  
  size_t start = pos;
  size_t stop = NO_INDEX;
  if (s [pos] == '\"')
  {
    pos++;
    start = pos;
    findChar ('\"');
    ASSERT (s [pos] == '\"');
    stop = pos;
  }
  
  findChar (',');
  if (stop == NO_INDEX)
    stop = pos;
  pos++;
  
  return s. substr (start, stop - start);
}


  
  
void csvLine2vec (const string &line,
                  StringVector &words)
{
  words. clear ();
  Csv csv (line);
  string s;
  while (csv. goodPos ())
  {
    s = csv. getWord ();
    trim (s);
    words << s;
  }
}




// Json

Json::Json (JsonContainer* parent,
            const string& name)
{ 
  ASSERT (parent);
  
  if (const JsonArray* jArray = parent->asJsonArray ())
  {
    ASSERT (name. empty ());  
    const_cast <JsonArray*> (jArray) -> data << this;
  }
  else if (const JsonMap* jMap = parent->asJsonMap ())
  {
    ASSERT (! name. empty ());  
    ASSERT (! contains (jMap->data, name));
    const_cast <JsonMap*> (jMap) -> data [name] = this;
  }
}



Token Json::readToken (istream &is)
{
  const string delim ("[]{},:");
  
  string s;
  bool spaces = true;
  bool isC = false;
  while (is)
  {
    char c;
    is. get (c);
    const bool isDelim = charInSet (c, delim);
    if (spaces)
      if (isSpace (c))
        continue;
      else
      {
        if (isDelim)
          return Token (string (1, c), Token::eDelimiter);
        spaces = false;
        isC = c == '\'';
      }
    else
    {
      // Opposite to toStr()
      if (isC)
      {
        if (c == '\'')
          break;
        else
        if (c == '\\')
        {
          is. get (c);
          if (c == 'n')
            c = '\n';
        }
      }
      else
        if (isSpace (c) || isDelim)
        {
          is. unget ();
          break;
        }
    }
      
    s += c;
  }
  
  if (isC)
    s. erase (0, 1);
  else
    ASSERT (! s. empty ());
  
  return Token (s, isC ? Token::eText : Token::eNumber);
}



void Json::parse (istream& is,
                  const Token& firstToken,
                  JsonContainer* parent,
                  const string& name)
{
  string s (firstToken. name);
  strLower (s);
  
  if (firstToken. type == Token::eDelimiter && s == "{")
    new JsonMap (is, parent, name);
  else if (firstToken. type == Token::eDelimiter && s == "[")
    new JsonArray (is, parent, name);
  else if (firstToken. type == Token::eText && s == "null")
    new JsonNull (parent, name);
  else if (firstToken. type == Token::eNumber)
  {
    if (   contains (s, '.')
        || contains (s, 'e')
       )
    {
      uint decimals = 0;
      size_t pos = s. find ('.');
      if (pos != string::npos)
      {
        pos++;
        while (   pos < s. size () 
               && isDigit (s [pos])
              )
        {
          pos++;
          decimals++;
        }
      }
      new JsonDouble (str2<double> (s), decimals, parent, name);
    }
    else
      new JsonInt (str2<int> (s), parent, name); 
  }
  else
    new JsonString (firstToken. name, parent, name);
}



int Json::getInt () const
{ 
  if (! this || asJsonNull ())
    throw runtime_error ("undefined");
  if (const JsonInt* j = asJsonInt ())
    return j->n;
  throw runtime_error ("Not a JsonInt");
}



double Json::getDouble () const
{ 
  if (! this)
    throw runtime_error ("undefined");
  if (asJsonNull ())
    return NAN;
  if (const JsonDouble* j = asJsonDouble ())
    return j->n;
  throw runtime_error ("Not a JsonDouble");
}



string Json::getString () const
{ 
  if (! this || asJsonNull ())
    throw runtime_error ("undefined");
  if (const JsonString* j = asJsonString ())
    return j->s;
  throw runtime_error ("Not a JsonString");
}



const Json* Json::at (const string& name_arg) const
{ 
  if (! this)
    throw runtime_error ("undefined");
  if (asJsonNull ())
    return nullptr;
  if (const JsonMap* j = asJsonMap ())
    return findPtr (j->data, name_arg);
  throw runtime_error ("Not a JsonMap");
}



const Json* Json::at (size_t index) const
{ 
  if (! this)
    throw runtime_error ("undefined");
  if (asJsonNull ())
    return nullptr;
  if (const JsonArray* j = asJsonArray ())
  {
    if (index >= j->data. size ())
      throw runtime_error ("Index out of range");
    else
      return j->data [index];
  }
  throw runtime_error ("Not a JsonArray");
}        



size_t Json::getSize () const
{ 
  if (! this)
    throw runtime_error ("undefined");
  if (asJsonNull ())
    return 0;
  if (const JsonArray* j = asJsonArray ())
    return j->data. size ();
  throw runtime_error ("Not a JsonArray");
}        



// JsonArray


JsonArray::JsonArray (istream& is,
                      JsonContainer* parent,
                      const string& name)
: JsonContainer (parent, name)
{
  bool first = true;
  for (;;)
  {
    Token token (readToken (is));
    if (token. isDelimiter (']'))
      break;
    if (! first)
    {
      ASSERT (token. isDelimiter (','));
      token = readToken (is);
    }
    parse (is, token, this, string());
    first = false;
  }
}



void JsonArray::print (ostream& os) const
{ 
  os << "[";
  bool first = true;
  for (const Json* j : data)
  {
    ASSERT (j);
    if (! first)
      os << ",";
    j->print (os);
    first = false;
  }
  os << "]";
}



// JsonMap

JsonMap::JsonMap ()
{
  ASSERT (! jRoot);
  jRoot = this;
}



JsonMap::JsonMap (const string &fName)
{
  ifstream ifs (fName. c_str ());
  if (! ifs. good ())
    ERROR_MSG ("Bad file: " + fName);
  const Token token (readToken (ifs));
  ASSERT (token. isDelimiter ('{'));
  parse (ifs);
}



void JsonMap::parse (istream& is)
{
  ASSERT (data. empty ());
  
  bool first = true;
  for (;;)
  {
    Token token (readToken (is));
    if (token. isDelimiter ('}'))
      break;
    if (! first)
    {
      ASSERT (token. isDelimiter (','));
      token = readToken (is);
    }
    ASSERT (! token. name. empty ());
    const Token colon (readToken (is));
    ASSERT (colon. isDelimiter (':'));
    Json::parse (is, readToken (is), this, token. name);
    first = false;
  }
}



JsonMap::~JsonMap ()
{ 
  for (auto it : data)
    delete it. second;
}



void JsonMap::print (ostream& os) const
{ 
  os << "{";
  bool first = true;
  for (const auto it : data)
  {
    if (! first)
      os << ",";
    os << toStr (it. first) << ":";
    const Json* j = it. second;
    ASSERT (j);
    j->print (os);
    first = false;
  }
  os << "}";
}




JsonMap* jRoot = nullptr;




// Offset

size_t Offset::size = 0;
const size_t Offset::delta = 2;




//

void exec (const string &cmd)
{
  ASSERT (! cmd. empty ());
  EXEC_ASSERT (system (cmd. c_str ()) == 0);
}




// Application

Application::Arg::Arg (const string &name_arg,
                       const string &description_arg)
: Named (name_arg)
, description (description_arg)
{
  trim (name);
  ASSERT (! name. empty ());
  ASSERT (! contains (name, ' '));
  
  ASSERT (! description. empty ());
}



void Application::Key::saveText (ostream &os) const
{
  os << "[-" << name; 
  if (! flag)
  {
    os << " ";
    const bool quoted = value. empty () || contains (value, ' ');
    if (quoted)
      os << "\"";
    os << value;
    if (quoted)
      os << "\"";
  }
  os << "]";
}



void Application::addKey (const string &name, 
                          const string &argDescription,
                          const string &defaultValue)
{
  ASSERT (! contains (args, name));
  keys << Key (name, argDescription, defaultValue);
  args [name] = & keys. back ();
}



void Application::addFlag (const string &name,
                           const string &argDescription)
{
  ASSERT (! contains (args, name));
  keys << Key (name, argDescription);
  args [name] = & keys. back ();
}



void Application::addPositional (const string &name,
                                 const string &argDescription)
{
  ASSERT (! contains (args, name));
  positionals << Positional (name, argDescription);
  args [name] = & positionals. back ();
}



string Application::getArg (const string &name) const
{
  if (contains (args, name))  
    return args. at (name) -> value;
  throw runtime_error ("Parameter \"" + name + "\" is not found");
}



bool Application::getFlag (const string &name) const
{
  const string value (getArg (name));
  const Key* key = args. at (name) -> asKey ();
  if (! key || ! key->flag)
    throw runtime_error ("Parameter \"" + name + "\" is not a flag");
  return value == "true";
}



string Application::getInstruction () const
{
  string instr (description);
  instr += "\nUsage: " + programName;
  for (const Positional& p : positionals)
    instr += " " + p. str ();
  for (const Key& key : keys)
    instr += " " + key. str ();
  
  instr += string ("\n") + "Help:  " + programName + " -help";

  return instr;
}



string Application::getHelp () const
{
  string instr (getInstruction ());
  instr += "\nParameters:";
  const string par ("\n  ");
  for (const Positional& p : positionals)
    instr += par + p. str () + ": " + p. description;
  for (const Key& key : keys)
    instr += par + key. str () + ": " + key. description;;
  
  return instr;
}



int Application::run (int argc, 
                      const char* argv []) 
{
	try
  { 
    for (int i = 0; i < argc; i++)  
      programArgs. push_back (argv [i]);
    ASSERT (! programArgs. empty ());
  
      
    // positionals, keys
    bool first = true;
    posIt = positionals. begin ();
    Key* key = nullptr;
    Set<string> keysRead;
    for (string s : programArgs)
    {
      if (first)
      {
        programName = rfindSplit (s, fileSlash);
        ASSERT (! programName. empty ());
      }
      else
      {
        if (! s. empty () && s [0] == '-')
        {
          if (key)
            errorExitStr ("Key with no value: " + key->name + "\n" + getInstruction ());
          const string name (s. substr (1));
          if (name == "help")
          {
            cout << getHelp () << endl;
            exit (1);
          }
          if (! contains (args, name))
            errorExitStr ("Unknown key: " + name + "\n" + getInstruction ());
          key = const_cast <Key*> (args [name] -> asKey ());
          if (! key)
            errorExitStr (name + " is not a key\n" + getInstruction ());
          if (keysRead. contains (name))
            errorExitStr ("Parameter \"" + name + "\" is used more than once");
          else
            keysRead << name;
          if (key->flag)
          {
            key->value = "true";
            key = nullptr;
          }
        }
        else
          if (key)
          {
            ASSERT (! key->flag);
            key->value = s;
            key = nullptr;
          }
          else
          {
            if (posIt == positionals. end ())
              errorExitStr ("Too many positional arguments\n" + getInstruction ());
            (*posIt). value = s;
            posIt++;
          }
      }
      first = false;
    }
  
  
    const string logFName = getArg ("log");
  	ASSERT (! logPtr);
    if (! logFName. empty ())
  		logPtr = new OFStream ("", logFName, "");
  
  	if (getFlag ("qc"))
  		qc_on = true;

  	Verbose vrb (str2<int> (getArg ("verbose")));
  	
  	if (getFlag ("noprogress"))
  		Progress::disable ();
  	if (getFlag ("profile"))
  		Chronometer::enabled = true;
  
  	const string jsonFName = getArg ("json");
  	ASSERT (! jRoot);
  	if (! jsonFName. empty ())
  	{
  		new JsonMap ();
  	  ASSERT (jRoot);
  	}
  
  
    if (programArgs. size () == 1)
    {
      cout << getInstruction () << endl;
      exit (1);
    }
    
    if (posIt != positionals. end ())
      errorExitStr ("Too few positional arguments\n" + getInstruction ());


  	body ();

  
  	if (! jsonFName. empty ())
  	{
  	  ASSERT (jRoot);
  		OFStream f ("", jsonFName, "");
      jRoot->print (f);
      delete jRoot;
      jRoot = nullptr;
    }
  
  	if (! logFName. empty ())
  	{
  	  delete logPtr;
  	  logPtr = nullptr;
  	  if (remove (logFName. c_str ()))
  	  {
  	    cout << "Cannot remove log file \"" << logFName << "\"" << endl;
  	    abort ();
  	  }
    }
	}
	catch (const std::exception &e) 
	{ 
	  errorExit (e. what ());
  }


  return 0;
}




}


