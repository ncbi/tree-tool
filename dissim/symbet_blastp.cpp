// symbet_blastp.cpp  

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
*   Dissimilarity by BLASTP symmetric best hits
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "evolution.hpp"
using namespace DM_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{
  
  
struct Hit : Named
{
  size_t positive {0};
  size_t nident {0};
  
  Hit (const string &name_arg,
       size_t positive_arg,
       size_t nident_arg)
    : Named (name_arg)
    , positive (positive_arg)
    , nident (nident_arg)
    {
      QC_ASSERT (! name. empty ());
      QC_ASSERT (! contains (name, ' '));
      QC_ASSERT (nident);
      QC_ASSERT (positive >= nident);
    }
  Hit () = default;
  Hit (Hit&& other) = default;
  Hit& operator= (Hit&& other) = default;
  
  
  bool operator< (const Hit &other) const
    { LESS_PART (*this, other, positive);
      LESS_PART (*this, other, nident);
      LESS_PART (*this, other, name);  // Tie resolution
      return false;
    }
};



void read_blastp (const string &fName,
                  unordered_map<string/*qseqid*/,Hit> &name2hit)
{
  ASSERT (name2hit. empty ());
  name2hit. rehash (100000);  // PAR
  LineInput f (fName);
  Istringstream iss;
  while (f. nextLine ())
  {
    iss. reset (f. line);
    string qseqid, sseqid;
    size_t positive = 0;
    size_t nident = 0;
    iss >> qseqid >> sseqid >> positive >> nident;
    QC_ASSERT (nident);
    Hit hit (sseqid, positive, nident);
    Hit& old = name2hit [qseqid];
    if (old < hit)
      old = std::move (hit);
  }
}



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print dissimilarity by BLASTP symmetric best hits")
    {
      version = VERSION;
  	  addPositional ("blastp1", "File 1 with BLASTP hits: qseqid sseqid positive nident");
  	  addPositional ("blastp2", "File 2 with BLASTP hits: qseqid sseqid positive nident");
  	  addPositional ("n1", "Number of proteins in genome 1");
  	  addPositional ("n2", "Number of proteins in genome 2");
  	}



	void body () const final
	{
		const string fName1 = getArg ("blastp1");
		const string fName2 = getArg ("blastp2");
		const uint n1       = arg2uint ("n1");
		const uint n2       = arg2uint ("n2");
		
		
	  unordered_map<string/*qseqid*/,Hit>	name2hit_1;
	  unordered_map<string/*qseqid*/,Hit>	name2hit_2;

	  read_blastp (fName1, name2hit_1);
	  read_blastp (fName2, name2hit_2);
	  
	  if ((size_t) n1 < name2hit_1. size ())
	    throw runtime_error ("Number of proteins in genome 1 (" + to_string (n1) + ") is less than the number of qseqids in " + fName1 + " (" + to_string (name2hit_1. size ()) + ")");
	  if ((size_t) n2 < name2hit_2. size ())
	    throw runtime_error ("Number of proteins in genome 2 (" + to_string (n2) + ") is less than the number of qseqids in " + fName2 + " (" + to_string (name2hit_2. size ()) + ")");
	  	  
	  size_t n = 0;
	  for (const auto& it : name2hit_1)
	    if (const Hit* hit = findPtr (name2hit_2, it. second. name))
	      if (it. first == hit->name)
	        n++;
 	    // n = sum of Hit::positive or Hit::nident ??
	        
	  if (verbose ())
	  {
  	  PRINT (name2hit_1. size ());
  	  PRINT (name2hit_2. size ());
  	  PRINT (n);
  	}
	  cout << intersection2dissim ( (Real) n1
	                              , (Real) n2
	                              , (Real) n
	                              , 50.0   // PAR
	                              , 0.5    // PAR
	                              , true)  
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



