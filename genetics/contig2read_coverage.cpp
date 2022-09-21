// contig2read_coverage.cpp

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
*   Print average read coverage per contig and flag noise
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
#include "../tsv/tsv.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace 
{
	
	
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print a tsv-table with columns: Contig\tLength\tCoverage_avg\tNoise")
	  {
      version = VERSION;
		  addPositional ("in", "Tab-delimited file where each line corresponds to one nucleotide and has format: <contig>\t<coverage>");
		  addPositional ("noise", "noise threshold (0..1)");
	  }

  
  
	void body () const final
  {
	  const string inFName        = getArg ("in");
	  const Prob   noiseThreshold = str2real (getArg ("noise"));
	  QC_ASSERT (isProb (noiseThreshold));
	  

    MeanVar totalMv;    
    map<string/*contig*/,MeanVar> mvs;
    {
      LineInput f (inFName);
      Istringstream iss;
      string contig;
      while (f. nextLine ())
      {
        iss. reset (f. line);
        size_t coverage;
        iss >> contig >> coverage;
        mvs [contig] << (Real) coverage;
        totalMv << (Real) coverage;
      }
    }
    
    TextTable tt;
    {
      tt. header << TextTable::Header ("Contig") 
                 << TextTable::Header ("Length") 
                 << TextTable::Header ("Coverage_avg")
                 << TextTable::Header ("Noise");
      const Real totalMean = totalMv. getMean ();
      StringVector row;
      row << string () << to_string (totalMv. n) << to_string (totalMean) << "0";
      tt. rows << move (row);
      for (const auto& it : mvs)
      {
        const MeanVar& mv = it. second;
        const Real mean = mv. getMean ();
        row. clear ();
        row << it. first << to_string (mv. n) << to_string (mean) << (mean < noiseThreshold * totalMean ? "1" : "0");
        tt. rows << move (row);
      }
    }
    tt. qc ();
    tt. saveText (cout);
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



