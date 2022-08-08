// fixed2tsv.cpp

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
*   Convert a fixed-field text file to a tsv-table
*
*/

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



namespace
{
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert a fixed-field text file to a tsv-table")
  	{
      version = VERSION;
  	  addPositional ("in", "Fixed-field text file");
  	  addPositional ("field_start", "1-based start positions of fieldas");
  	  addFlag ("skip_comment", "Skip comment starting with '#'");
  	}
  	
  	
 
	void body () const final
	{
		const string fName       = getArg ("in");
		const string fieldFName  = getArg ("field_start");
		const bool   skipComment = getFlag ("skip_comment");
		
		
		Vector<size_t> starts;
		{
		  ifstream f (fieldFName);
		  for (;;)
		  {
		    size_t n = no_index;
		    f >> n;
		    if (n == no_index)
		      break;
		    if (starts. empty ())
		    {
		      if (n != 1)
		        throw runtime_error ("First field must start with 1");
		    }
		    else
		      if (n <= starts. back ())
		        throw runtime_error ("Field start " + to_string (n) + " is before " + to_string (starts. back ()));
		    starts << n;
		  }
		}
		if (starts. empty ())
		  throw runtime_error ("No field start positions");
		// Make 0-based
		for (size_t& i : starts)
		  i--;
		
		
		
		{
  		LineInput f (fName, 1);  // PAR
  		while (f. nextLine ())
  		{
  		  if (skipComment && ! f. line. empty () && f. line [0] == '#')
  		    continue;
  		  if (contains (f. line, '\t'))
  		    throw runtime_error ("Tab character in line:\n" + f. line);
  		  if (f. line. size () <= starts. back ())
  		    throw runtime_error ("Two short line:\n" + f. line);
  		  StringVector s;
  		  FOR_REV (size_t, i, starts. size ())
  		  {
  		    s << f. line. substr (starts [i]);
  		    f. line. erase (starts [i]);
  		  }
  		  ASSERT (s. size () == starts. size ());
  		  FOR_REV (size_t, i, s. size ())
  		  {
  		    trim (s [i]);
  		    cout << s [i];
  		    if (i)
  		      cout << '\t';
  		  }
  		  cout << endl;
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



