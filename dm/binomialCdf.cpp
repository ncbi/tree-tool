// binomialCdf.cpp

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
*   C.d.f. of a binomial distribution
*
*/


#undef NDEBUG

#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Binomial distribution CDF")
  	{
  	  addPositional ("N", "Number of coin tosses");	
  	  addPositional ("P", "Probability of heads");	
  	  addPositional ("M", "Number of heads");
  	  addFlag("reverse", "1 - CDF");
  	  addFlag("add_pmf", "Add p.m.f.");
  	}



	void body () const final
	{
		const uint n       = str2<uint> (getArg ("N"));
		const Prob p       = str2real (getArg ("P"));
		const uint m       = str2<uint> (getArg ("M"));
		const bool reverse = getFlag ("reverse");
		const bool add_pmf = getFlag ("add_pmf");


    Binomial bin;
    bin. setParam ((int) n, p);
    
    Prob res = bin. cdf (m);
    if (reverse)
      res = 1 - res;
    cout << res << endl;    

    if (add_pmf)
      cout << "pmf = " << bin. pdf (m) << endl;
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


