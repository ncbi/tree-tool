// splitList.cpp

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
*   Split a list into subsets
*
*/


#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace 
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Partition <in> into parts of size <size> sequentially.\n\
Parts are named <out_dir>/<i>")
  {
	  addPositional ("in", "Text file");
	  addPositional ("size", "# lines in one part");
	  addPositional ("out_dir", "Output directory");
	  addKey ("start_part", "Start number of a part", "1");	  
	}


	void body () const final
	{
		const string in       = getArg ("in");
		const streamsize size = str2<streamsize> (getArg ("size"));
		const string out_dir  = getArg ("out_dir");
		const uint start_part = str2<uint> (getArg ("start_part"));
		ASSERT (size > 0);
		ASSERT (start_part >= 1);

    LineInput inF (in, 1024 * 1024);  // PAR
    OFStream outF;
    uint part = 0;
    streamsize n = size;  // # lines in outF
    bool writing = false;
    Progress prog;
    while (inF. nextLine ())
    {
    	if (n == size)
    	{
        prog ();
        if (outF. is_open ())
    		  outF. close ();
    		part++;
    		if (part >= start_part)
    		  writing = true;
    		if (writing)
    		  outF. open (out_dir, toString (part), "");
    		n = 0;
    	}
    	if (writing)
   	    outF << inF. line << endl;
   	  ASSERT (outF. good ());
    	n++;
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
