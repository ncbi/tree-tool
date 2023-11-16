// count.cpp

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
*   Print statistics of a list of numbers
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
using namespace DM_sp;
#include "../version.inc"

#include "../common.inc"




namespace 
{
	
	
void report (const string &attr,
             Real value)
{	
  unique_ptr<const ONumber> on;
  if (isInteger (value))
    on. reset (new ONumber (cout, 0, false));
  cout << attr << '\t' << value << endl; 
}

	


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Statistics of a sequence of numbers from cin", false)
    {
      version = VERSION;
    }



	void body () const final
	{
		MeanVar mv;
		string s;
		while (cin >> s)
		{
			const Real x = str2real (s);
			if (! isNan (x))
				mv << x;
		}
		
		report ("count",  mv. n);
		report ("mean",   mv. getMean ());
		report ("var",    mv. getVar ());
		report ("SD",     mv. getSD ());
	  report ("sum",    mv. s);
		report ("meanSD", mv. getSD () / sqrt (mv. n));
	  report ("min",    mv. v_min);
	  report ("max",    mv. v_max);
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



