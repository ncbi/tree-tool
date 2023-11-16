// sim2dist.cpp

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
*   Convert a similarity to distance
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


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert similarity to distance")
    {
      version = VERSION;
  	  addPositional ("file", dmSuff + "-file without the extension");
  	  addPositional ("attrName", "Attribute name of real-valued object-object table in the " + dmSuff + "-file");
  	  addFlag ("sqr", "Make squared distances");
  	}



	void body () const
	{
		const string fName    = getArg ("file");
		const string attrName = getArg ("attrName");
		const bool makeSqr    = getFlag ("sqr");
		
		
    Dataset ds (fName);
    
    // sim
    const RealAttr2* sim = nullptr;
    {
      const Attr* attr = ds. name2attr (attrName);
      ASSERT (attr);
      sim = attr->asRealAttr2 ();
    }
    ASSERT (sim);    
    ASSERT (! sim->asPositiveAttr2 ());
    
	
	  // Check sim->matr  
	  {
      Real maxCorrection;
      size_t row_bad, col_bad;
      const_cast <RealAttr2*> (sim) -> matr. symmetrize (maxCorrection, row_bad, col_bad);
      if (maxCorrection > 2 * pow (10, - (Real) sim->decimals))
        ds. comments << "maxCorrection = " + toString (maxCorrection) + " at " + ds. objs [row_bad] -> name + ", " + ds. objs [col_bad] -> name;
    }

	  {
  	  size_t row;
  	  size_t col;
  	  if (sim->matr. existsMissing (false, row, col))
      {
        cout << attrName << '[' << ds. objs [row] -> name << "," << ds. objs [col] -> name << "] is missing" << endl;
        exit (1);
      }  
    }
	

    auto dist = new PositiveAttr2 (sim->name + "_dist", ds, sim->decimals); 
    dist->matr = sim->matr;
    dist->matr. similarity2sqrDistance ();
    if (! makeSqr)
    {
      dist->matr. sqrtAll ();
      dist->decimals += 4;  // PAR
    }
    
    delete sim;
    
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



