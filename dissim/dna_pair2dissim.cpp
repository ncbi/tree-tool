// dna_pair2dissim.cpp

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
*   Align and compute dissimilarity for a pair of DNA sequences
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "../genetics/seq.hpp"
using namespace Seq_sp;
#include "align.hpp"
using namespace Align_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{
	
	
const Dna* readDna (const string &fName)
{
  LineInput in (fName);
  EXEC_ASSERT (in. nextLine ());
  const auto* dna = new Dna (in, 1024 * 1024, false);  // PAR
  dna->qc ();
  
  return dna;
}

	

struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Align and compute dissimilarity for pairs of DNA sequences")
    {
      version = VERSION;
  	  addPositional ("pairs", "File with pairs of DNA file names");
  	  addPositional ("dnaDir", "Directory with DNA sequences");
  	  addPositional ("align_len_min", "Min. sequence length");
  	  addPositional ("out", "Output file with lines: <DNA1> <DNA2> <dissimilarity>. <DNA1> < <DNA2>. Lines match the file <pairs>");
  	  addFlag       ("large", "DNA sequences are in subdirectoirs of <dnaDir> which are hash-codes modulo " + to_string (small_hash_class_max) + " of the DNA file names");
  	  addFlag       ("global", "Global alignment, otherwise semiglobal");
  	  addFlag       ("relative", "Dissimilarity relative to sequence length");
  	  addFlag       ("diff", "Dissimilarity is the number of different nucleotides");
  	  addKey        ("band", "Alignment band; 0 - no band", "0");
  	  addKey        ("coeff", "Coefficient to multiply the dissimilarity by", "1.0");  // ??
  	  addKey        ("power", "Power to raise the dissimilarity in",  "1.0");
  	  addKey        ("name_new", "Name of a sequence which is not in <dnaDir>");
  	  addKey        ("file_new", "File with the sequence <name_new>");
  	  addFlag       ("print", "Print alignment");
    }


	
	void body () const final
  {
	  const string pairsFName    =               getArg ("pairs");
	  const string dnaDir        =               getArg ("dnaDir");
	  const string outFName      =               getArg ("out");
	  const bool   large         =               getFlag ("large");
	  const bool   global        =               getFlag ("global");
	  const size_t align_len_min = str2<size_t> (getArg ("align_len_min"));
	  const bool   relative      =               getFlag ("relative");
	  const bool   diff          =               getFlag ("diff");
	  const size_t band          = str2<size_t> (getArg ("band"));
	  const Real   coeff         = str2real     (getArg ("coeff"));
	  const Real   power         = str2real     (getArg ("power"));
	  const string name_new      =               getArg ("name_new");
	  const string file_new      =               getArg ("file_new");
	  const bool   printP        =               getFlag ("print");
	  QC_ASSERT (coeff > 0.0);
	  QC_ASSERT (power > 0.0);
	  QC_ASSERT (! (relative & diff));
	  QC_ASSERT (name_new. empty () == file_new. empty ());	  
	  QC_IMPLY (band, diff);
    
    
    OFStream output (outFName);
    ONumber on (output, 6, true);  // PAR
    map<string/*fName*/,const Dna*> name2dna;  // not delete'd
    PairFile in (pairsFName, false, false, 1);  // PAR
    while (in. next ())
    {
      const string h1 (large ? to_string (str2hash_class (in. name1, false)) + "/" : "");
      const string h2 (large ? to_string (str2hash_class (in. name2, false)) + "/" : "");
      const string fName1 (in. name1 == name_new ? file_new : (dnaDir + "/" + h1 + in. name1));
      const string fName2 (in. name2 == name_new ? file_new : (dnaDir + "/" + h2 + in. name2));
      if (! contains (name2dna, in. name1))  name2dna [in. name1] = readDna (fName1);
      if (! contains (name2dna, in. name2))  name2dna [in. name2] = readDna (fName2);
      const Dna& dna1 = * name2dna [in. name1];
      const Dna& dna2 = * name2dna [in. name2];

      unique_ptr<Align_sp::Align> align (new Align_sp::Align (dna1, dna2, ! global, global ? 0 : align_len_min, band));
			if (diff)
			{
			  align->setAlignment (dna1. seq, dna2. seq);
			  if (band && align->getDiff () >= band)
			  {
			    align. reset (new Align_sp::Align (dna1, dna2, ! global, global ? 0 : align_len_min, 0));
  			  align->setAlignment (dna1. seq, dna2. seq);
  			}
			}
			
			if (printP)
			{
			  cout << in. name1 << ' ' << in. name2 << endl;
				align->printDistances (cout);
				cout << endl;
				align->printAlignment (cout, 60);  // PAR
			}
			
					
			const Real d = relative 
			                 ? align->getDissim () 
			                 : diff
			                   ? (Real) align->getDiff ()
			                   : align->getMinEditDistance ();
		  const Real dissim = coeff * pow (d, power);  
      output << in. name1 << '\t' << in. name2 << '\t' << dissim << endl;
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



