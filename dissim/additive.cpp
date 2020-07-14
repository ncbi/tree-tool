// additive.cpp

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
*   Check additivity of a dissimilarity
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace 
{


Real getAdditivity (const Matrix &matr,
                    size_t x,
                    size_t y,
                    size_t u,
                    size_t v)
// Return: <= 1.0 <=> additive; may be NaN
{
  // QC
  Vector<size_t> vec (4);
  vec [0] = x;
  vec [1] = y;
  vec [2] = u;
  vec [3] = v;
  vec. sort ();
  ASSERT (vec. isUniq ());
  
  const Real xy = matr. get (false, x, y);
  const Real xu = matr. get (false, x, u);
  const Real xv = matr. get (false, x, v);
  const Real yu = matr. get (false, y, u);
  const Real yv = matr. get (false, y, v);
  const Real uv = matr. get (false, u, v);
  
  const Real s = xy + xu + xv + yu + yv + uv;
  if (isNan (s) || s == INF)
    return NaN;
    
  ASSERT (xy >= 0.0);
  ASSERT (xu >= 0.0);
  ASSERT (xv >= 0.0);
  ASSERT (yu >= 0.0);
  ASSERT (yv >= 0.0);
  ASSERT (uv >= 0.0);
    
  return (xy + uv) / max (xu + yv, yu + xv);
}



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Check additivity of a dissimilarity")
  	{
  	  version = VERSION;
  		// Input
  	  addPositional ("file", dmSuff + "-file without the extension");
  	  addPositional ("attrName", "Attribute name of a distance in the <file>");
  	}



	void body () const final
	{
		const string fName    = getArg ("file");
		const string attrName = getArg ("attrName");
		
		
    Dataset ds (fName);
    
    // dist
    const PositiveAttr2* dist = nullptr;
    {
      const Attr* attr = ds. name2attr (attrName);
      QC_ASSERT (attr);
      dist = attr->asPositiveAttr2 ();
    }
    QC_ASSERT (dist);
    
    
    Vector<size_t> objs (ds. objs. size (), no_index);
    FFOR (size_t, i, objs. size ())
      objs [i] = i;
    objs. randomOrder ();
    size_t quartets = 0;
    size_t additive = 0;
    for (size_t i = 0; i + 4 < objs. size (); i += 4)
    {
      const Real additivity = getAdditivity ( dist->matr
                                            , objs [i + 0]
                                            , objs [i + 1]
                                            , objs [i + 2]
                                            , objs [i + 3]
                                            );
      if (isNan (additivity))
        continue;
      // Print distribution of additivity ??
      quartets++;
      if (additivity <= 1.0)
        additive++;
    }
    ASSERT (additive <= quartets);
    cout << additive << " / " << quartets << " = " << Real (additive) / Real (quartets) << endl;
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



