// printDistTree.cpp

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
*   Print a distance tree
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "distTree.hpp"
using namespace DistTree_sp;
#include "../version.inc"



namespace 
{



struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Print a tree made by makeDistTree")
	{
	  version = VERSION;
    // Input
	  addPositional ("input_tree", "Tree file");
	  addKey ("data", dmSuff + "-file without " + strQuote (dmSuff) + " to read object comments");
	  addKey ("dissim_attr", "Dissimilarity attribute name in the <data> file");
	  addKey ("variance", "Dissimilarity variance: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]); 
	  addKey ("variance_power", "Power for -variance pow; > 0", "NaN");
	  addKey ("name_match", "File with lines: <name_old> <tab> <name_new>, to replace leaf names");
	  addKey ("decimals", "Number of decimals in arc lengths", toString (dissimDecimals));
	  addFlag ("order", "Order subtrees by the number of leaves descending");
    // Output
	  addKey ("format", "newick|dm (Data Master format, makeDistTree output)|ASNT (textual ASN.1)", "newick");
	  addFlag ("ext_name", "Extended leaf names for newick");
	}



	void body () const final
  {
	  const string input_tree     = getArg ("input_tree");
	  const string dataFName      = getArg ("data");
	  const string dissimAttrName = getArg ("dissim_attr");
	               varianceType   = str2varianceType (getArg ("variance"));  // Global    
	               variancePower  = str2real (getArg ("variance_power"));    // Global
	  const string name_match     = getArg ("name_match");
	  const size_t decimals       = str2<size_t> (getArg ("decimals"));
  	const bool order            = getFlag ("order");
	  const string format         = getArg ("format");
  	const bool ext_name         = getFlag ("ext_name");

    if (input_tree. empty ())
      throw runtime_error ("-input_tree must be present");
  //if (dataFName. empty () != dissimAttrName. empty ())
    //throw runtime_error ("The both data file and the dissimilarity attribute must be present or absent");
		if (! isNan (variancePower) && varianceType != varianceType_pow)
		  throw runtime_error ("-variance_power requires -variance pow");
		if (isNan (variancePower) && varianceType == varianceType_pow)
		  throw runtime_error ("-variance_power is needed by -variance pow");
		if (variancePower <= 0.0)
		  throw runtime_error ("-variance_power must be positive");
		      

    DistTree tree (input_tree, dataFName, dissimAttrName, string());
    tree. qc ();    
    if (order)
      tree. sort ();
    if (! dataFName. empty ())
      tree. setLeafNormCriterion ();
    tree. qc ();    
    
    if (! name_match. empty ())
    {
      LineInput f (name_match, 10 * 1024, 1000);  // PAR
      string name_old, name_new;
      while (f. nextLine ())
      {
        name_new = f. line;
        name_old = findSplit (name_new, '\t');
        QC_ASSERT (! name_old. empty ());
        QC_ASSERT (! name_new. empty ());
        if (const Leaf* leaf = findPtr (tree. name2leaf, name_old))
          var_cast (leaf) -> name = name_new;
      #if 0
        else
          throw runtime_error ("Object '" + name_old + "' does not exist");
      #endif
      }
    }

   	cout << fixed << setprecision ((int) decimals);  
    if (format == "newick")
      tree. printNewick (cout, false, ! ext_name);
    else if (format == "dm")
      tree. saveText (cout);
    else if (format == "ASNT")
      tree. printAsn (cout);
    else
      throw runtime_error ("Unknown format " + format);
	}
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


