// numeric_test.cpp

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
*   Test of numeric.{hpp,cpp}
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
using namespace DM_sp;
#include "version.inc"



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Test numeric.cpp")
  	{
  	  version = VERSION;
  	  addPositional ("go", "Go");
  	}



	void body () const final
	{  
  #if 0
    cout << "sizeof(float) = " << sizeof (float) << endl;
    cout << "sizeof(double) = " << sizeof (double) << endl;
  #endif

    ASSERT_EQ (zeta (1.5),    2.612, 1e-3);
    ASSERT_EQ (zeta (1.5, 1), 2.612, 1e-3);
    ASSERT_EQ (zeta (2),   1.645, 1e-3);
    ASSERT_EQ (zeta (3),   1.202, 1e-3);
    ASSERT_EQ (zeta (4),   1.0823, 1e-4);
    ASSERT_EQ (zeta (2, 2), 0.645, 1e-3);
    
    {
	    Sum sum;
	    const uint n = 1000000;
	    FOR (uint, i, n)
	      sum. add (1.0 / n);
	    ASSERT_EQ (sum. get (), 1, 1e-6);
	  }

    {
	    SumLn sum;
	    const uint n = 1000000;
	    FOR (uint, i, n)
	      sum. addExp (- log (n));
	    ASSERT_EQ (exp (sum. getLn ()), 1, 1e-6);
 	  }

  #if 0
    #define REAL2STR(x)  cout << #x << ": " << real2str (x,6) << endl;
 	  REAL2STR (0);
 	  REAL2STR (1);
 	  REAL2STR (0.0001);
 	  REAL2STR (1e5);
 	  REAL2STR (0.000123456);
 	  #undef REAL2STR
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


