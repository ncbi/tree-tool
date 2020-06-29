// distTree_new.cpp

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
*   Add new objects to a distance tree
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
		: Application ("Find location of new objects in a distance tree.\n\
Update: <incremental distance tree directory>/search/")
		{
		  version = VERSION;
		  
		  // Input
		  addPositional ("data", "Directory with data");
		  addFlag ("init", "Initialize search");
		  
  	//addKey ("dissim_power", "Power to raise dissimilarity in", "1");

  	  addKey ("variance", "Dissimilarity variance function: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]);
  	  addKey ("variance_power", "Power for -variance pow; >= 0", "NaN");
  	  addFlag ("variance_dissim", "Variance is computed off dissimilarities");
  	  addKey ("variance_min", "Min. dissimilarity variance; to be added to the computed variance", "0");
  	    	  
		  addKey ("name", "Name of the object");
		  addKey ("dissim", "File of the format: <obj1> <obj2> <dissimilarity>");
		  
		  // Output
		  addKey ("request", "Output file of the format: <obj1> <obj2>");
		  addKey ("leaf", "Output file of the format: <obj_new> <obj1>-<obj2> <leaf_len> <arc_len>");
		}
	
	
	
	void body () const final
  {
	  const string dataDir       = getArg ("data");
	  const bool   init          = getFlag ("init");

	             //dissim_power        = str2real (getArg ("dissim_power"));      // Global

	               varianceType    = str2varianceType (getArg ("variance"));  // Global
	               variancePower   = str2real (getArg ("variance_power"));    // Global
	               variance_min    = str2real (getArg ("variance_min"));      // Global
	//const bool   variance_dissim =           getFlag ("variance_dissim");  

	  const string name          = getArg ("name");
	  const string dissimFName   = getArg ("dissim");
	  const string requestFName  = getArg ("request");
	  const string leafFName     = getArg ("leaf");
	   
    if (! isRight (dataDir, "/"))
      throw runtime_error (strQuote (dataDir) + " must end with '/'");

		if (! isNan (variancePower) && varianceType != varianceType_pow)
		  throw runtime_error ("-variance_power requires -variance pow");
		if (isNan (variancePower) && varianceType == varianceType_pow)
		  throw runtime_error ("-variance_power is needed by -variance pow");
		if (variancePower <= 0.0)
		  throw runtime_error ("-variance_power must be positive");
      
    ASSERT (name. empty () == dissimFName.  empty ());
    ASSERT (name. empty () == requestFName. empty ());
    ASSERT (name. empty () == leafFName.    empty ());


    if (verbose ())
    {
      DistTree::printParam (cout);
      cout << endl;
    }

    DistTree tree (dataDir, string (), false, false, false);
    tree. qc ();     

    if (verbose ())
    {
      tree. printInput (cout);
      cout << endl;
    }
    
    if (name. empty ())
    {
      const string newDir (dataDir + "search/");
      FileItemGenerator fig (1, true, newDir);  // PAR
  	  string item;
  	  while (fig. next (item))
      {
        const NewLeaf nl (tree, newDir, item, init);
        nl. qc ();
      }
    }
    else
    {
      const NewLeaf nl (tree, name, dissimFName, leafFName, requestFName, init);
      nl. qc ();
    }
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


