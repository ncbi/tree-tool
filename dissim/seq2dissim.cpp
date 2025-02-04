// seq2dissim.cpp

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
*   Align and compute dissimilarity for a pair of sequences.
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "../genetics/seq.hpp"
using namespace Seq_sp;
#include "align.hpp"
#include "../version.inc"

#include "../common.inc"



namespace 
{
  
  
Seq* readSeq (const string &fName,
              bool aa)
{
  LineInput in (fName);
  EXEC_ASSERT (in. nextLine ());
  Seq* seq = nullptr;
  if (aa)
  {
    auto pep = new Peptide (in, 10 * 1024, true);  
    pep->pseudo = true;
    seq = pep;
  }
  else
    seq = new Dna (in, 1024 * 1024, true);  
  QC_ASSERT (seq);
  seq->unSparse ();  
  seq->qc ();
  
  return seq;
}



bool ambigClear (const string &seq,  // sparse
                 size_t mutStart, 
                 size_t mutSize,
                 size_t noambigWindow,
                 bool aa)
{
  ASSERT (noambigWindow != no_index);

  for (size_t i = mutStart, n = 0; i-- > 0 && n < noambigWindow; )
  {
    if (isAmbig (seq [i], aa))
      return false;
    if (seq [i] != '-')
      n++;
  }

  for (size_t i = mutStart + mutSize, n = 0; seq [i] && n < noambigWindow; i++)
  {
    if (isAmbig (seq [i], aa))
      return false;
    if (seq [i] != '-')
      n++;
  }
            
  return true;
}



struct ThisApplication final : Application
{
  
  ThisApplication ()
    : Application ("Align and compute dissimilarity for a pair of sequences.\nPrint: ref_match: <target start>-<target stop> (human coordinates)")
    {
      version = VERSION;
  	  addPositional ("target", "FASTA file 1 with a sequence (target)");
  	  addPositional ("reference", "FASTA file 1 with a sequence (reference)");
  	  addKey ("prot_name", "Protein name. Empty string means that sequences are DNA");
  	  addFlag ("global", "Global alignment, otherwise semiglobal");
  	  addKey ("match_len_min", "Min. match length. Valid for semiglobal alignment", "60");
  	  addKey ("noambig", "Min. window around a mutation with no ambiguity. -1: allow ambiguities", "-1");
  	  addKey ("mutation", "File for mutations");
	  	addFlag ("ambig_start", "Replace deletions at the start of <target> by the ambiguity character");
           	  //ambig_end ??
      addFlag ("blosum62", "Use substitution matrix BLASUM62, otherwise PAM30. Valid for protein sequences");
      addFlag ("alignment", "Print alignment");
    }


	
	void body () const final
  {
	  const string targetFName    = getArg ("target");
	  const string refFName       = getArg ("reference");
	  const string protName       = getArg ("prot_name");
	  const bool global           = getFlag ("global");
	        size_t match_len_min  = str2<size_t> (getArg ("match_len_min"));
	  const size_t noambig        = getArg ("noambig") == "-1" ? no_index : str2<size_t> (getArg ("noambig"));
	  const string mutFName       = getArg ("mutation");
    const bool   ambig_start    = getFlag ("ambig_start");
    const bool   blosum62       = getFlag ("blosum62");
    const bool   printAlignment = getFlag ("alignment");
    

    if (global)
    	match_len_min = 0;
    else
    	if (match_len_min == 0)
    		throw runtime_error ("match_len_min cannot be 0 for a semiglobal alignment");

    const bool aa = ! protName. empty ();
    
    if (! aa && blosum62)
      throw runtime_error ("-blasum62 requires protein sequences");
      

    const unique_ptr<const Seq> targetSeq (readSeq (targetFName, aa));
    const unique_ptr<const Seq> refSeq    (readSeq (refFName,    aa));
    ASSERT (targetSeq);
    ASSERT (refSeq);

    unique_ptr<Align_sp::Align> align;
    if (aa)
      align. reset (new Align_sp::Align (* targetSeq->asPeptide (), * refSeq->asPeptide (), ! global, match_len_min, blosum62));
    else
      align. reset (new Align_sp::Align (* targetSeq->asDna (), * refSeq->asDna (), ! global, match_len_min, 0 /*band*/));
    ASSERT (align);
		
		if (verbose ())
		  align->saveText (cout);
				
		align->setAlignment (targetSeq->seq, refSeq->seq);

 		size_t targetStart = 0;
		size_t targetStop = targetSeq->seq. size ();
    {		
  		while (align->sparse2 [targetStart] == '-')
  	  {
  	    ASSERT (align->sparse1 [targetStart] != '-');
  		  targetStart++;
  		}
  		
  		size_t i = align->sparse2. size ();
  		while (align->sparse2 [i - 1] == '-')
  	  {
  	    ASSERT (align->sparse1 [i] != '-');
  	    ASSERT (i);
  	    i--;
  	    ASSERT (targetStop);
  		  targetStop--;
  		}
    }
		if (targetStart >= targetStop)
		  throw runtime_error ("No alignment");
		cout << "ref_match: " << targetStart + 1 << ' ' << targetStop << endl;
  		

		if (! mutFName. empty ())
		{
		  OFStream f (mutFName);
		  size_t refPos = 0;
  		size_t mismatchStart = no_index;  
  		size_t refStart = no_index;
  		FFOR (size_t, i, align->sparse1. size ())
  	  {
  		  if (align->sparse1 [i] == align->sparse2 [i])
  		  {
  		    if (mismatchStart != no_index)
  		    {
  		      ASSERT (mismatchStart < i);
  		      const size_t len = i - mismatchStart;
  		      string ref    (align->sparse2. substr (mismatchStart, len));
  		      string allele (align->sparse1. substr (mismatchStart, len));
  		      replaceStr (ref,    "-", "");
  		      replaceStr (allele, "-", "");
  		      ASSERT (ref != allele);
  		      ASSERT (refStart <= refSeq->seq. size ());
  		      const bool outside =     ref. empty ()
                  		           && (   refStart == 0 
                  		               || refStart == refSeq->seq. size ()
                  		              );
            if (! outside)
            {
              if (ambig_start && refStart == 0)
              {
                ASSERT (! ref. empty ());
                allele = string (ref. size (), aa ? 'X' : 'n');
              }
    		      const Mutation mut (protName, refStart, ref, allele);
    		      mut. qc ();
    		      if (   noambig == no_index 
    		          || (! mut. ambig && ambigClear (align->sparse1, mismatchStart, len, noambig, aa))
    		         )
    		      f << mut << endl;
    		    }
  		      mismatchStart = no_index;
  		      refStart = no_index;
  		    }
  		  }
  		  else
  		    if (mismatchStart == no_index)
  		    {
  		      mismatchStart = i;
  		      refStart = refPos;
  		    }
  		  if (align->sparse2 [i] != '-')
  		    refPos++;
  		}
  		ASSERT (refPos == refSeq->seq. size ());
    }


    cout << endl;
    align->printStats (cout);
    align->printDistances (cout);
    if (printAlignment)
		  align->printAlignment (cout, 60);  // PAR
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



