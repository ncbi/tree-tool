// attr2_2paup.cpp

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
*   Convert a two-way attribute to the PAUP* format
*
*/


#undef NDEBUG

#include "../../common.hpp"
using namespace Common_sp;
#include "../matrix.hpp"
#include "../dataset.hpp"
using namespace DM_sp;
#include "../../version.inc"

#include "../../common.inc"



namespace 
{


struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Print an Attr2 in PAUP* format")
    {
      version = VERSION;
      
  	  addPositional ("file", dmSuff + "-file");
  	  addPositional ("attr2Name", "Attribute name of an object-object table in the " + dmSuff + "-file");
  	  addPositional ("map", "Output map file for nw_rename");
  	}



	void body () const final
	{
		const string fName     = getArg ("file");
		const string attr2Name = getArg ("attr2Name");
		const string mapFName  = getArg ("map");		
		
		
    Dataset ds (fName);
    
    // dist
    const RealAttr2* dist = nullptr;
    {
      const Attr* attr = ds. name2attr (attr2Name);
      ASSERT (attr);
      dist = attr->asRealAttr2 ();
    }
    ASSERT (dist);

    Matrix& matr = const_cast <RealAttr2*> (dist) -> matr;
    
    OFStream f (mapFName);
    
    cout << "#NEXUS" << endl;
    cout << "begin distances;" << endl;
    cout << "dimensions ntax = " << ds. objs. size () << ";" << endl;
    cout << "format triangle = both;" << endl;
    cout << "matrix";
    ONumber on (cout, 5, false);
    FOR (size_t, row, ds. objs. size ())
    {
      const string name (("X" + toString (row + 1) + "          "). substr (0, 10));
      cout << endl << name;
      f << name << ' ' << ds. objs [row] -> name << endl;
      FOR (size_t, col, ds. objs. size ())
        cout << ' ' << matr. get (false, row, col);
    }
    cout << ';' << endl << "end;" << endl;
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



