// bool2dissim.cpp

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
*   Add a dissimilarity attribute and print the resulting Data Master file
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
  

struct WeightAttr
{
  const ExtBoolAttr1* attr;
  Real weight;
};


  
Vector<WeightAttr> weightAttrs;



struct ObjPair 
{
  size_t row;
  size_t col;
  Real diff;
};



void computeObjPair (size_t from, 
                     size_t to, 
                     Notype /*&res*/,
                     Vector<ObjPair> &objPairs)
{
  Progress prog (to - from);  
  FOR_START (size_t, i, from, to)
  {
    prog ();
    ObjPair& objPair = objPairs [i];
    objPair. diff = 0.0;
    for (const WeightAttr& wa : weightAttrs)
      if (   (* wa. attr) [objPair. row] != UBOOL
          && (* wa. attr) [objPair. col] != UBOOL
          && (* wa. attr) [objPair. row] != (* wa. attr) [objPair. col]
         )
        objPair. diff += wa. weight;
  }
}
  
	
	
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Add a dissimilarity attribute and print the resulting Data Master file")
    {
  	  version = VERSION;
  	  addPositional ("data", dmSuff + "-file with Boolean attributes"); 
  	  addPositional ("weight", "Attribute weight file where each line has format: <attr_name> <weight, >= 0>");
  	  addPositional ("dissim", "Name of two-way dissimilarity attribute to add (Hamming distance)"); 
    }


	
	void body () const final
  {
	  const string dsFName        = getArg ("data");
    const string weightFName    = getArg ("weight");
    const string dissimAttrName = getArg ("dissim"); 


    Dataset ds (dsFName);
    ds. qc ();

    {
      LineInput f (weightFName);
      string name;
      Real weight;
      Istringstream iss;
      while (f. nextLine ())
      {
        iss. reset (f. line);
        weight = NaN;
        iss >> name >> weight;
        QC_ASSERT (weight >= 0);
        const Attr* attr = ds. name2attr (name);
        QC_ASSERT (attr);
        if (const ExtBoolAttr1* boolAttr = attr->asExtBoolAttr1 ())
          weightAttrs << WeightAttr {boolAttr, weight};
      }
    }

	  auto dist = new PositiveAttr2 (dissimAttrName, ds, 6);  // PAR
	  Vector<ObjPair> objPairs;  objPairs. reserve ((ds. objs. size () * (ds. objs. size () - 1)) / 2);
	  FOR (size_t, row, ds. objs. size ())
	  {
	  	dist->matr. putDiag (row, 0);
	    FOR_START (size_t, col, row + 1, ds. objs. size ())
	      objPairs << ObjPair {row, col, 0};
	  }
    vector<Notype> notypes;
	  arrayThreads (computeObjPair, objPairs. size (), notypes, ref (objPairs));
    for (const ObjPair& objPair : objPairs)
	    dist->matr. putSymmetric (objPair. row, objPair. col, objPair. diff);
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



