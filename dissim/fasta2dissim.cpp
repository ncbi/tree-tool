// fasta2dissim.cpp

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
*   Create a dataset with a dissimilarity attribute of a collection of DNA or proteins or print: <obj1> <obj2> <dissim> <score> <self-score1> <self-score2>
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/matrix.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "../genetics/seq.hpp"
using namespace Seq_sp;
#include "align.hpp"



namespace 
{
	

const string attrName ("dissim");



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Create a dataset with a dissimilarity attribute \"" + attrName + "\" of a collection of DNA or proteins or print: <obj1> <obj2> <dissim> <score> <self-score1> <self-score2>")
    {
  	  // Input
  	  addPositional ("Multi_FASTA", "Multi-FASTA file with DNA or protein sequences");
  	  addFlag ("aa", "Amino-acid (protein) sequences, otherwise DNA");
  	  addFlag ("global", "Global alignment, otherwise semiglobal");
  	  addKey ("match_len_min", "Min. match length. Valid for semiglobal alignment", "60");
  	  addFlag ("blosum62", "Use BLOSUM62, otherwise PAM30");
  	  addFlag ("mismatch_frac", "Distance is mismatch fraction");
  	  addKey ("power", "Raise raw dissimilarity to this power", "1");
  	  addKey ("coeff", "Multiply raw dissimilarity by this coefficient", "1");
  	  addKey("class", "File with classes of the sequencess (in the same order as Multi_FASTA)");
  	  // Output
  	  addKey ("dataset", "Output dataset file with dissimilarity matrix without " + dmSuff);  
    }


	
	void body () const final
  {
	  const string inFName       =               getArg ("Multi_FASTA");
	  const string dsFName       =               getArg ("dataset");
	  const bool aa              =               getFlag ("aa");
	  const bool global          =               getFlag ("global");
	        size_t match_len_min = str2<size_t> (getArg ("match_len_min"));
	  const bool blosum62        =               getFlag ("blosum62");
	  const bool mismatch_frac   =               getFlag ("mismatch_frac");
	  const Real power           = str2real     (getArg ("power"));
	  const Real coeff           = str2real     (getArg ("coeff"));
	  const string classF        =               getArg ("class");
    QC_ASSERT (power > 0);
    QC_ASSERT (coeff > 0);
    
    if (global)
    	match_len_min = 0;
    else
    	if (match_len_min == 0)
    		throw runtime_error ("match_len_min cannot be 0 for a semiglobal alignment");    		

    QC_ASSERT (! (blosum62 && mismatch_frac));


    Dataset ds;  

    VectorOwn<Seq> seqs;  seqs. reserve (1000);  // PAR
    {
		  Multifasta faIn (inFName, aa);
		  while (faIn. next ())
		  {
		  	const Seq* seq = nullptr;
		  	if (aa)
		  	{
			    auto pep = new Peptide (faIn, 1000/*PAR*/, false);
			    pep->ambig2X ();
			    pep->pseudo = true;
			    seq = pep;
		  	}
		  	else
			    seq = new Dna (faIn, 1000/*PAR*/, false);
			  ASSERT (seq);
			  seq->qc ();
		    seqs << seq;
		    string s (seq->name);
		    const size_t objNum = ds. appendObj (findSplit (s));
		    trim (s);
		    const_cast <Obj*> (ds. objs [objNum]) -> comment = s;
		  }
		}
    ASSERT (seqs. size () == ds. objs. size ());
    if (ds. objs. empty ())
    	throw runtime_error ("No DNA sequences");
    if (ds. objs. size () < 2)
    	throw runtime_error ("Too few DNA sequences");
    	
    auto lenAttr = new IntAttr1 ("Length", ds);
    FFOR (size_t, objNum, seqs. size ())
    {
      const string& seq = seqs [objNum] -> seq;
      size_t size = seq. size ();
      ASSERT (size);
      if (aa && isRight (seq, "*"))
        size--;
      (*lenAttr) [objNum] = (int) size;
    }
			
    if (! classF. empty ())
    {
      StringVector seqClass (classF, (size_t) 1000);  // PAR
      if (seqClass. size () != seqs. size ())
        throw runtime_error ("# Sequences does not match # classes");
      for (string& s : seqClass)
        trim (s);
      //
      auto classAttr = new NominAttr1 ("Class", ds);
      FFOR (size_t, objNum, seqClass. size ())
        (*classAttr) [objNum] = classAttr->category2index (seqClass [objNum]);
    }

	  auto dissimAttr = new PositiveAttr2 (attrName, ds, 6);  // PAR
	  {
		  Progress prog;
		  FFOR (size_t, row, seqs. size ())
		  {
		    const Seq* seq1 = seqs [row];
		    ASSERT (seq1);
		  	prog (seq1->getId ());
		    IMPLY (! aa, seq1->asDna ());
		    IMPLY (aa, seq1->asPeptide ());
		  	dissimAttr->matr. put (false, row, row, 0);
		    FFOR_START (size_t, col, row + 1, seqs. size ())
		    {
			    const Seq* seq2 = seqs [col];
			    ASSERT (seq2);
			    IMPLY (! aa, seq2->asDna ());
			    IMPLY (aa, seq2->asPeptide ());
			    unique_ptr<const Align_sp::Align> align;
			    if (const Dna* dna1 = seq1->asDna ())
						align. reset (new Align_sp::Align (*dna1, * seq2->asDna (),     ! global, match_len_min, 0));
					else if (const Peptide* pep1 = seq1->asPeptide ())
						align. reset (new Align_sp::Align (*pep1, * seq2->asPeptide (), ! global, match_len_min, blosum62));
					else
						ERROR;
					ASSERT (align. get ());
					const Real dissim_raw = mismatch_frac
					                          ? align->getMismatchFrac ()
					                          : /*aa 
					                             ? align->getDissim ()   // Can be NaN
					                             :*/ align->getMinEditDistance ();
					const Real dissim = coeff * pow (dissim_raw, power);
					if (dsFName. empty ())
						cout         << seq1->getId ()
						     << '\t' << seq2->getId ()
						     << '\t' << dissim
						     << '\t' << align->score
						     << '\t' << align->self_score1
						     << '\t' << align->self_score2
						     << endl;
					else
  				  dissimAttr->matr. putSymmetric (row, col, dissim); 
		    }
		  }
		}
		
		if (! dsFName. empty ())
		{
	    OFStream os ("", dsFName, dmExt); 
	    ds. saveText (os);
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



