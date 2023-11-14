// efffectiveSize.cpp

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
*   Print effective number of different classes
*
*/


#undef NDEBUG

#include "common.hpp"
using namespace Common_sp;
#include "version.inc"

#include "common.inc"



namespace 
{
	
	
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print effective number of different classes")
    {
      version = VERSION;
  	  addPositional ("in", "Text file");
  	  addKey ("num", "Approximate number of different classes", "1000");
    }



	void body () const final
	{
	  const string inFName =          getArg ("in");
	  const size_t num     = (size_t) arg2uint ("num");
	  QC_ASSERT (num > 0);
	  
	  
    unordered_map<string,size_t> name2num;
    {
      LineInput f (inFName);
      name2num. rehash (num * 2);
      while (f. nextLine ())
      {
        trim (f. line);
        name2num [f. line] ++;
      }
    }
    
    size_t n = 0;
    for (const auto& it : name2num)
      n += it. second;

    if (verbose ())
    {
      cout << "# Items: " << n << endl;
      cout << "# Different classes: " << name2num. size () << endl;
    }

    double s = 0.0;    
    for (const auto& it : name2num)
    {
      const double p = (double) it. second / (double) n;
      s += p * p;
    }
    
    cout << 1.0 / s << endl;
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



