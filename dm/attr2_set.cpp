// attr2_set.cpp

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
*   Set the values of a two-way attribute for specified pairs of objects
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "matrix.hpp"
#include "dataset.hpp"
using namespace DM_sp;
#include "../version.inc"

#include "../common.inc"



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print a " + dmSuff + "-file setting the values of a symmetric two-way attribute for specified pairs of objects")
    {
      version = VERSION;
      addPositional ("file", dmSuff + "-file");
      addPositional ("attr", "2-way attribute name");
      addPositional ("pair_value", "File with pairs of objects with new values");
    }
	
	
	
	void body () const final
	{
		const string inFName   = getArg ("file");
		const string attrName  = getArg ("attr");
		const string pairFName = getArg ("pair_value");


    Dataset ds (inFName);    
	  ds. qc ();
    const Attr* attr_ = ds. name2attr (attrName);
    QC_ASSERT (attr_);
    const RealAttr2* attr = attr_->asRealAttr2 ();
    QC_ASSERT (attr);
    
    LineInput f (pairFName);
    string name1, name2;
    Istringstream iss;    
    while (f. nextLine ())
    {
      iss. reset (f. line);
      Real value = NaN;
      iss >> name1 >> name2 >> value;
      const size_t i = ds. getName2objNum (name1);
      if (i == no_index)
        throw runtime_error ("Object " + strQuote (name1) + " is not in " + inFName);
      const size_t j = ds. getName2objNum (name2);
      if (j == no_index)
        throw runtime_error ("Object " + strQuote (name2) + " is not in " + inFName);
      var_cast (attr->matr). putSymmetric (i, j, value);
    }
	  ds. qc ();
    
    ds. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


