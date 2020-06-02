// prot2dissim.cpp

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
*   Align and compute dissimilarity for a pair of proteins
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../genetics/seq.hpp"
using namespace Seq_sp;



namespace 
{
	

Peptide readPeptide (const string &fName)
{
  LineInput in (fName);
  EXEC_ASSERT (in. nextLine ());
  const Peptide pep (in, 1024 * 1024, false);  // PAR
  pep. qc ();
  
  return pep;
}



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Align and compute dissimilarity for a pair of proteins")
    {
  	  addPositional ("FASTA1", "FASTA file 1 with a protein sequence");
  	  addPositional ("FASTA2", "FASTA file 1 with a protein sequence");
  	  addPositional ("mat", "Amino acid substution matrix");
  	  addPositional ("gap_open_cost", "Gap open cost, >= 0");
  	  addPositional ("gap_cost", "Gap cost, >= 0");
  	  addKey ("semiglobal_match_len_min", "Min. match length. Valid for semiglobal alignment; 0 - global alignment", "0");
    }


	
	void body () const final
  {
	  const string inFName1    = getArg ("FASTA1");
	  const string inFName2    = getArg ("FASTA2");
	  const string matFName    = getArg ("mat");
	  const double gapOpenCost = str2<double> (getArg ("gap_open_cost"));
	  const double gapCost     = str2<double> (getArg ("gap_cost"));
	  const size_t semiglobalMatchLen_min = str2<size_t> (getArg ("semiglobal_match_len_min"));
    
    if (gapOpenCost < 0.0)
      throw runtime_error ("-gap_open_cost must be non-negative");
    if (gapCost < 0.0)
      throw runtime_error ("-gap_cost must be non-negative");
    

    const Peptide pep1 (readPeptide (inFName1));
    const Peptide pep2 (readPeptide (inFName2));
    const SubstMat substMat (matFName);


		const Align align (pep1, pep2, substMat, gapOpenCost, gapCost, semiglobalMatchLen_min);
		align. qc ();
	  cout << "Score = " << align. score << endl;  
	  cout << "Self_score1 = " << align. self_score1 << endl;  
	  cout << "Self_score2 = " << align. self_score2 << endl;  

		cout << "Identity = " << align. matches << '/' << align. transformations. size () 
		     << " (" << ((double) align. matches / (double) align. transformations. size () * 100) << "%)" << endl;
		cout << "Min. edit distance = " << align. getMinEditDistance () << endl;
		cout << endl;
		align. printAlignment (pep1. seq, pep2. seq, 60);  // PAR
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



