// snp2dissim.cpp

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
*   Add a SNP dissimilarity attribute and print the resulting Data Master file
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/dataset.hpp"
#include "../dm/optim.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace 
{
  

struct MutAttr
{
  const ExtBoolAttr1* attr;
  uint mutations;
};


  
Vector<MutAttr> mutAttrs;



struct ObjPair 
{
  size_t row;
  size_t col;
  Real dissim;
};



struct Func : Func1
{
  size_t row;
  size_t col;
  
  Func (size_t row_arg,
        size_t col_arg)
    : row (row_arg)
    , col (col_arg)
    {}
  
  Real f (Real t) final
    {
      ASSERT (t >= 0.0);
      if (t == 0.0)
        return -INF;
      Real s = 0.0;
      for (const MutAttr& ma : mutAttrs)
        if (   ma. mutations
            && (* ma. attr) [row] != UBOOL
            && (* ma. attr) [col] != UBOOL
           )
        {
          const Real rate = (Real) ma. mutations;
          const Real a = exp (- 2.0 * t * rate);
          ASSERT (a >= 0.0);
          ASSERT (a < 1.0);
          const Prob p = (* ma. attr) [row] == (* ma. attr) [col]
                           ? a * rate / (1.0 + a)
                           : - rate / (1.0 / a - 1.0);  // = - a * rate (1.0 - a)
          s += p;
        }
      return s;
    }  
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
    Func f (objPair. row, objPair. col);    
    objPair. dissim = f. findZero (0.0, 1.0, 1e-5);  // PAR  // 1.0 = tree length
      // f. findZeroPositive (0.05, 1e-6);  // PAR
  }
}
  
	
	
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Add a SNP dissimilarity attribute and print the resulting Data Master file")
    {
  	  version = VERSION;
  	  addPositional ("data", dmSuff + "-file with Boolean attributes"); 
  	  addPositional ("mutations", "Attribute file where each line has format: <attr_name> <# mutations>");
  	  addPositional ("dissim", "Name of two-way dissimilarity attribute to add (Hamming distance)"); 
    }


	
	void body () const final
  {
	  const string dsFName        = getArg ("data");
    const string mutationsFName = getArg ("mutations");
    const string dissimAttrName = getArg ("dissim"); 


    Dataset ds (dsFName);
    ds. qc ();

    {
      LineInput f (mutationsFName);
      string name;
      uint mutations;
      Istringstream iss;
      while (f. nextLine ())
      {
        iss. reset (f. line);
        mutations = (uint) -1;
        iss >> name >> mutations;
        QC_ASSERT (mutations != (uint) -1);
        const Attr* attr = ds. name2attr (name);
        QC_ASSERT (attr);
        if (const ExtBoolAttr1* boolAttr = attr->asExtBoolAttr1 ())
          mutAttrs << MutAttr {boolAttr, mutations};
      }
    }

	  auto dist = new PositiveAttr2 (dissimAttrName, ds, 6);  // PAR
	  Vector<ObjPair> objPairs;  objPairs. reserve ((ds. objs. size () * (ds. objs. size () - 1)) / 2);
	  FOR (size_t, row, ds. objs. size ())
	  {
	  	dist->matr. putDiag (row, 0.0);
	    FOR_START (size_t, col, row + 1, ds. objs. size ())
	      objPairs << ObjPair {row, col, 0.0};
	  }
    vector<Notype> notypes;
	  arrayThreads (computeObjPair, objPairs. size (), notypes, ref (objPairs));
    for (const ObjPair& objPair : objPairs)
	    dist->matr. putSymmetric (objPair. row, objPair. col, objPair. dissim);
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



