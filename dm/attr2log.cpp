// attr2log.cpp

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
*   Logarithmize an attribute
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
    : Application ("Print a " + dmSuff + "-file adding log(attr) named \"log_\" + attr")
    {
      addPositional ("file", dmSuff + "-file");
      addPositional ("attr", "Attribute name");
      addFlag ("inf2nan", "Convert -INF to NaN");
    }
	
	
	
	void body () const final
	{
		const string inFName  = getArg ("file");
		const string attrName = getArg ("attr");
		const bool inf2nan    = getFlag ("inf2nan");


    Dataset ds (inFName);    
    const Attr* attr_ = ds. name2attr (attrName);
    ASSERT (attr_);
    const RealAttr1* attr = attr_->asRealAttr1 ();
    ASSERT (attr);
    
    auto logAttr = new RealAttr1 ("log_" + attrName, ds, attr->decimals + 1);
    FOR (size_t, objNum, ds. objs. size ())
    {
      const Real r = (*attr) [objNum];
      if (isNan (r))
        continue;
      ASSERT (r >= 0);
      Real res = log (r);
      if (inf2nan && res == -INF)
      	res = NaN;
      (*logAttr) [objNum] = res;
    }
    
    ds. print (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


