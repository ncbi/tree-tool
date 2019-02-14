// setMinus.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



#define SORT 0



namespace 
{


struct Item
{
	bool num {false};
	enum Position {SMALLEST, NORMAL, BIGGEST};
	Position pos {SMALLEST};
	// Valid if pos == NORMAL
  string s;  // Valid if !num
	size_t i {0};     // Valid if  num
	

	explicit Item (bool num_arg)
	  : num (num_arg)
	  {}
	Item& operator= (const Item &item)
		{
			ASSERT (num == item. num);
			pos = item. pos;
			if (pos == NORMAL)
			{
				if (num)
					i = item. i;
				else
				{
					s = item. s;
					ASSERT (! s. empty ());
				}
			}
			return *this;
		}


	string name () const
	  { ASSERT (pos == NORMAL);
	  	return num ? toString (i) : s; 
	  }
	bool operator== (const Item &item) const
	  { ASSERT (num == item. num);
	  	if (         pos != NORMAL
	  		  || item. pos != NORMAL
	  		 )
	  		return false;
	  	if (num)
	  		return i == item. i;
	  	return s == item. s;	  	
	  }
	bool operator< (const Item &item) const
	  { ASSERT (num == item. num);
	  	switch (pos)
	  	{
	  		case SMALLEST: return true;
	  		case BIGGEST:  return false;
	  		case NORMAL:
			  	switch (item. pos)
			  	{
			  		case SMALLEST: return false;
			  		case BIGGEST:  return true;
			  		case NORMAL: 
					  	if (num)
					  		return i < item. i;
					  	return s < item. s;
			  	}
			}
	  	ERROR;
	  	return false;
	  }
	void read (const string &fName,
	           istream &is,
	           Item &itOld)
	  // Update: itOld
	  { ASSERT (pos != BIGGEST);
	  	for (;;)
	  	{
	  		if (is. eof ())
	  		{
	  			pos = BIGGEST;
	  			break;
	  		}
		  	readLine (is, s);
		    trim (s);  
		  	if (s. empty ())
		  		continue;
	  		pos = NORMAL;
	  		break;
		  }
		  ASSERT (pos != SMALLEST);
		  if (pos == NORMAL)
		  {
  		  if (num)
		  	  i = str2<size_t> (s);
		  	if (*this == itOld)
		  		throw runtime_error ("File: " + fName + ": Not unique: " + name ());
		  	if (*this < itOld)
		  		throw runtime_error ("File: " + fName + ": Not ordered: " + itOld. name () + " and " + name ());
      }
	  	itOld = *this;
	  }
  int compare (const Item &item) const
	  { if (*this == item)
	  	  return 0;
	  	if (*this < item)
	  		return -1;
	  	return 1;
	  }
};



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print the set-theoretic minus: <list1> \\ <list2>\n\
Time: O(L1 " /*"ln(L1) "*/ "+ L2 " /*"ln(L2) "*/ ") where L1=|list1|, L2=|list2|")
  	{
  	  addPositional ("list1", "File with ordered set of items");
  	  addPositional ("list2", "File with ordered set of items to be removed from list1");
  	  addFlag ("number", "Items are numbers");
  	}


#if SORT
	bool getMatch (const vector<string> &vec,
			         	 uint &index,
			         	 const string &target)
	// Update: index
	// Requires: List is sorted ascending
	{
	  ASSERT (! target. empty ());
	
	  for (; index < vec. size (); index++)
	  {
	    const int cmp = target. compare (vec. at (index));
	    if (cmp < 0)
	      return false;
	    if (cmp == 0)
	      return true;
	  }
		
	  return false;
	}
#endif



	void body () const final
	{
		const string list1Name = getArg ("list1");
		const string list2Name = getArg ("list2");
		const bool num = getFlag ("number");


  #if SORT
		vector <string> vec1 (readList (list1Name));
		sortAll (vec1);
		
		vector <string> vec2 (readList (list2Name));
		sortAll (vec2);
		
		uint j = 0;
		CONST_ITER (vector <string>, it, vec1)
		  if (! getMatch (vec2, j, *it))
		    cout << *it << endl;          
  #else
	  ifstream f1 (list1Name);
	  ifstream f2 (list2Name);
	  Item it1 (num);
	  Item it2 (num);
	  Item it1old (num);
	  Item it2old (num);
	  while (! f1. eof ())
	  {
	  	it1. read (list1Name, f1, it1old);

      int cmp;
      for (;;)
		  {
		  	cmp = it1. compare (it2);
		    if (cmp <= 0)
		    	break;
        if (f2. eof ())
       	  break;
  	  	it2. read (list2Name, f2, it2old);
		  }

	    if (cmp < 0)
		    cout << it1. name () << endl;
	  }
	#endif
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}
