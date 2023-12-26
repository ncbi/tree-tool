// setMinus.cpp

/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* Author: Vyacheslav Brover
*
* File Description:
*   Set-theoretic minus
*
*/


#undef NDEBUG

#include "common.hpp"
using namespace Common_sp;
#include "version.inc"

#include "common.inc"



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
      version = VERSION;
  	  addPositional ("list1", "File with ordered set of items");
  	  addPositional ("list2", "File with ordered set of items to be removed from list1");
  	  addFlag ("number", "Items are numbers");
  	  addFlag ("subset", "Break if <list1> is not a subset of <list2>");
  	}


	void body () const final
	{
		const string list1Name = getArg ("list1");
		const string list2Name = getArg ("list2");
		const bool num         = getFlag ("number");
		const bool subset      = getFlag ("subset");


	  IFStream f1 (list1Name);
	  IFStream f2 (list2Name);
	  constexpr size_t buf_size = 16 * 1024 * 1024;  // PAR
	  char* buf1 = new char [buf_size];  // Not delete'd
	  char* buf2 = new char [buf_size];  // Not delete'd
    QC_ASSERT (f1. rdbuf () -> pubsetbuf (buf1, buf_size)); 
    QC_ASSERT (f2. rdbuf () -> pubsetbuf (buf2, buf_size)); 
	  Item it1 (num);
	  Item it2 (num);
	  Item it1old (num);
	  Item it2old (num);
	  while (! f1. eof ())
	  {
	  	it1. read (list1Name, f1, it1old);

      int cmp = 0;
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
	    {
	      if (subset)
	        throw runtime_error (it1. name ());
        cout << it1. name () << endl;
		  }
	  }
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}
