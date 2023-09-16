// conv_comb.cpp

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
*   Convex combinination of two dissimilarities
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace
{



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convex combinination of two dissimilarities")
    {
      version = VERSION;
  	  addPositional ("in1", "File with lines: <obj1> <obj2> <dissimilarity>");
  	  addPositional ("in2", "File with lines: <obj1> <obj2> <dissimilarity>. Object pairs must match <in1>");
  	  addPositional ("coeff", "Coefficient for dissimilarity in <in1>, 0..1");
  	}



	void body () const final
	{
		const string fName1 = getArg ("in1");
		const string fName2 = getArg ("in2");
		const Real coeff    = str2real (getArg ("coeff"));
		QC_ASSERT (coeff >= 0.0);
		QC_ASSERT (coeff <= 1.0);
		
		
    LineInput in1 (fName1);  
    LineInput in2 (fName2);  
    Istringstream iss1;
    Istringstream iss2;
    string name1_1, name1_2, dissimS1,
           name2_1, name2_2, dissimS2;
    while (in1. nextLine ())
    {
      if (! in2. nextLine ())
        throw runtime_error ("File <in2> has fewer lines than file <in1>");

      iss1. reset (in1. line);
      name1_1. clear ();
      name1_2. clear ();
      dissimS1. clear ();
      iss1 >> name1_1 >> name1_2 >> dissimS1;

      iss2. reset (in2. line);
      name2_1. clear ();
      name2_2. clear ();
      dissimS2. clear ();
      iss2 >> name2_1 >> name2_2 >> dissimS2;

      try
      {
        if (name1_1 != name2_1)
          throw runtime_error ("First name " + strQuote (name1_1) + " in " + fName1 + " is different from first name " + strQuote (name2_1) + " in " + fName2);
        if (name1_2 != name2_2)
          throw runtime_error ("Second name in " + strQuote (name1_2) + " in " + fName1 + " is different from second name " + strQuote (name2_2) + " in " + fName2);
          
        cout         << name1_1 
             << '\t' << name1_2
             << '\t' << coeff * str2real (dissimS1) + (1.0 - coeff) * str2real (dissimS2)
             << endl;
      }
      catch (const exception &e)
      {
        throw runtime_error ("File " + fName1 + ": " + in1. lineStr () + "\n" 
                             "File " + fName2 + ": " + in2. lineStr () + "\n" 
                             + e. what ()
                            );
      }
    }
    if (in2. nextLine ())
      throw runtime_error ("File <in2> has more lines than file <in1>");
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);

}



