// attr2_2pairs.cpp

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
*   Convert a two-way attribute to pairs
*
*/


#undef NDEBUG
#include "../../common.inc"

#include "../../common.hpp"
using namespace Common_sp;
#include "../numeric.hpp"
#include "../matrix.hpp"
#include "../dataset.hpp"
using namespace DM_sp;
#include "../../version.inc"



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print pairs for an Attr2: <objName1>\t<objName2>\t<value>")
    {
      version = VERSION;
  	  addPositional ("file", dmSuff + "-file");
  	  addPositional ("attr2Name", "Attribute name of an object-object table in the " + dmSuff + "-file");
  	  addFlag ("symmetric", "Attribute is symmetric, print only lines where <objName1> <= <objName2>");
  	  addFlag ("diagonal", "Print the diagonal elements");
  	}



	void body () const final
	{
		const string fName     = getArg ("file");
		const string attr2Name = getArg ("attr2Name");
		const bool symmetric   = getFlag ("symmetric");
		const bool diagonal    = getFlag ("diagonal");
		
		
    Dataset ds (fName);
    
    const RealAttr2* dist = nullptr;
    {
      const Attr* attr = ds. name2attr (attr2Name);
      ASSERT (attr);
      dist = attr->asRealAttr2 ();
    }
    if (! dist)
    	throw runtime_error ("Two-way attribute " + attr2Name + " is not found");

    const Matrix& matr = const_cast <RealAttr2*> (dist) -> matr;
    
    FOR (size_t, col, ds. objs. size ())
    FOR (size_t, row, ds. objs. size ())
      if (! symmetric || ds. objs [row] -> name <= ds. objs [col] -> name)
        if (row != col || diagonal)
          cout         << ds. objs [row] -> name 
               << '\t' << ds. objs [col] -> name 
               << '\t' << matr. get (false, row, col)
               << endl;
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



