// blastmat.cpp

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
*   Convert a BLAST similarity matrix to a distance matrix, analyze
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;




namespace 
{
  
  
  
struct ThisApplication : Application
{
	ThisApplication ()
	  : Application ("Convert a BLAST similarity matrix to a distance matrix, analyze")
	  {
		  addPositional ("mat", "File with a BLAST score matrix"); 
		  addPositional ("gap", "Gap cost, >= 0");
		  addPositional ("gap_open", "Gap open cost, >= 0");
	  }


	
	void body () const final
  {
	  const string matName  = getArg ("mat");
	  const double gap_sim      = - str2<double> (getArg ("gap"));
	  const double gap_open_sim = - str2<double> (getArg ("gap_open"));
	  
	  
    const SubstMat mat (matName);
    {
      Unverbose unv;
      if (verbose ())
        mat. saveText (cout);
    }
    mat. qc ();
    
    if (verbose ())
      mat. saveText (cout);

    double dist [SubstMat::sim_size + 1] [SubstMat::sim_size + 1];  // (chars + gap) * (chars + gap)
    
    double dist_ave = 0.0;
    FOR (size_t, i, SubstMat::sim_size)
      if (mat. goodChar (i))
        FOR (size_t, j, SubstMat::sim_size)
          if (mat. goodChar (j))
          {
            dist [i] [j] = mat. getSubstitutionDist (i, j);  
            if (i != j)
              dist_ave += dist [i] [j];
          }
    dist_ave /= (double) (SubstMat::sim_size * (SubstMat::sim_size - 1));
    if (verbose ())
      cout << "dist_ave (not on diagonal) = " << dist_ave << endl;
    
    FOR (size_t, i, SubstMat::sim_size)
    {
      dist [i] [SubstMat::sim_size] = mat. getDeletionDist (i, gap_sim);
      dist [SubstMat::sim_size] [i] = dist [i] [SubstMat::sim_size];
    }
    dist [SubstMat::sim_size] [SubstMat::sim_size] = 0.0;

    // QC: dist[][]
    FOR (size_t, i, SubstMat::sim_size + 1)
    {
      if (! charInSet (char (i), peptideAlphabet))
        continue;
      // dist = 0
      QC_ASSERT (! dist [i] [i]);
      FOR (size_t, j, i)
      {
        if (! charInSet (char (j), peptideAlphabet))
          continue;
        // Symmetricity
        QC_ASSERT (dist [i] [j] == dist [j] [i]);
        // dist = 0
        if (dist [i] [j] <= 0.0)
          cout << "dist(" << char(i) << ',' << char(j) << ") <= 0" << endl;
        // Triangle inequality
        FOR (size_t, k, SubstMat::sim_size + 1)
          if (k == SubstMat::sim_size || charInSet (char (k), peptideAlphabet))
            if (dist [i] [j] > dist [i] [k] + dist [k] [j] - (k == SubstMat::sim_size ? gap_open_sim : 0.0))  // ??
            {
              const char ck = (k == SubstMat::sim_size ? '-' : char (k));
              cout << "dist(" << char(i) << ',' << char(j) << ") > dist(" << char(i) << ',' << ck << ") + dist(" << ck << ',' << char(j) << ")" << endl;
            }
      }
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



