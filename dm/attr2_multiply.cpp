// attr2_multiply.cpp

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
*   Multiply a two-way attribute by a coefficient
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "matrix.hpp"
#include "dataset.hpp"
using namespace DM_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{


const string clusterAttrName = "Cluster";



struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Multiply a two-way attribute and print the resulting data set")
    {
      version = VERSION;
  	  addPositional ("file", dmSuff + "-file");
  	  addPositional ("attr2Name", "Attribute name of an object-object table in the " + dmSuff + "-file");
  	  addPositional ("coefficient", "Coeffciient to multiply <attr2Name> by");
  	}



	void body () const final
	{
		const string fName     = getArg ("file");
		const string attr2Name = getArg ("attr2Name");
		const Real coefficient = str2real (getArg ("coefficient"));
		ASSERT (! isNan (coefficient));
		
		
    Dataset ds (fName);
    
    // dist
    RealAttr2* dist = nullptr;
    {
      const Attr* attr = ds. name2attr (attr2Name);
      ASSERT (attr);
      dist = const_cast <RealAttr2*> (attr->asRealAttr2 ());
    }
    ASSERT (dist);
    
    dist->decimals += (uint) max<long> (0, - DM_sp::round (log10 (coefficient)));
    
    Matrix& matr = dist->matr;
    FOR (size_t, row, ds. objs. size ())
    FOR (size_t, col, ds. objs. size ())
      matr. putProd (false, row, col, coefficient);
      
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



