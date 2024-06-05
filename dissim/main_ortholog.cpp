// main_ortholog.cpp

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
*   Find main common ortholog for a set of protein FASTA files
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "../genetics/seq.hpp"
using namespace Seq_sp;
#include "align.hpp"
#include "../version.inc"

#include "../common.inc"



namespace 
{


struct Genome : Named
{ 
  VectorOwn<Peptide> peptides;
  const Peptide* main {nullptr};
    // In peptides
    
    
  explicit Genome (const string &fName)
    : Named (getFileName (fName))
    {
      Multifasta faIn (fName, true, 0);
      size_t len = 0;
		  while (faIn. next ())
		  {
		    auto pep = new Peptide (faIn, 1000/*PAR*/, false);
		    pep->ambig2X ();
		    pep->qc ();
		    peptides << pep;
		    if (maximize (len, pep->seq. size ()))
		      main = pep;
		  }
		  if (! main)
		    throw runtime_error ("Empty FASTA file: " + fName);
		}
};

	
	
bool blosum62 = false;



Real getDissim (const Peptide* p1,
                const Peptide* p2)
{
  ASSERT (p1);
  ASSERT (p2);
  Align_sp::Align al (*p1, *p2, false/*global*/, 0/*match_len_min*/, blosum62);
	al. setAlignment (p1->seq, p2->seq); 
	al. qc ();
  return al. getDistance (Align_sp::Align::dist_rel_min_edit);
}
	
	

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find main common ortholog for a set of protein FASTA files") 
    {
      version = VERSION;
  	  // Input
  	  addPositional ("dir", "Directory with protein FASTA files");
  	  addFlag ("blosum62", "Use BLOSUM62, otherwise PAM30");
  	  // Output
  	  addPositional ("out", "Output FASTA fule with main orthologous proteins");
    }


	
	void body () const final
  {
	  const string dirName  = getArg ("dir");
	               blosum62 =               getFlag ("blosum62");
	  const string outFName = getArg ("out");
                                
                                
    Vector<Genome> genomes;
    {
      DirItemGenerator gen (1, dirName, false);
   	  string item;
   	  while (gen. next (item))
   	  {            
   	    Genome genome (dirName + "/" + item);
   	    genomes << std::move (genome);
   	  }
   	}
   	if (genomes. size () <= 1)
   	  throw runtime_error ("Too few genomes");
    
    // Init: Genome::main
    {
      Progress prog (genomes. size ());
      const Genome* prev = nullptr;
      for (const Genome& g : genomes)
      {
        prog (g. name);
        if (prev)
        {
          ASSERT (prev != & g);
          Real dissim_min = inf;
          for (const Peptide* p1 : g. peptides)
            for (const Peptide* p2 : prev->peptides)
              if (minimize (dissim_min, getDissim (p1, p2)))
              {
                var_cast (g).      main = p1;
                var_cast (prev) -> main = p2;
              }
          prev = nullptr;
        }
        else
          prev = & g;
      }
    }

    FOR (size_t, i, 10)  // PAR
    {
      cerr << endl << "Iteration " << i + 1 << endl;
      size_t changes = 0;
      Real dissim_total = 0.0;
      {
        Progress prog (genomes. size ());
        for (const Genome& g1 : genomes)
        {
          prog (g1. name);
          Real dissim_sum_min = inf;
          const Peptide* main = nullptr;
          for (const Peptide* p : g1. peptides)
          {
            Real dissim_sum = 0.0;          
            for (const Genome& g2 : genomes)
              if (& g2 != & g1)
              {
                const Real dissim = getDissim (p, g2. main);
                dissim_sum += dissim;
              }
            if (minimize (dissim_sum_min, dissim_sum))
              main = p;
          }
          if (! main)
            throw runtime_error ("Cannot align " + g1. name);
          dissim_total += dissim_sum_min;
          if (g1. main != main)
          {
            var_cast (g1). main = main;
            changes++;
          }
        }
      }
      cerr << "Changes: " << changes << endl;
      cerr << "Average dissimilarity: " << dissim_total / (Real) genomes. size () << endl;
      if (! changes)
        break;
    }
    
    
    OFStream f (outFName);
    for (const Genome& g : genomes)
    {
      const Peptide* p = g. main;
      ASSERT (p);
      var_cast (p) -> name = g. name + "-" + p->getId ();
      p->saveText (f);
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



