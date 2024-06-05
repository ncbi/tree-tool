// align.cpp

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
*   Alignment utilities
*
*/


#undef NDEBUG

#include "align.hpp"

#include "nw/nw_aligner.hpp"
//#include "../nw/mm_aligner.hpp"
#include "nw/nw_band_aligner.hpp"

#include "evolution.hpp"
using namespace DM_sp;

#include "../common.inc"



#define WU_BLASTN 



namespace Align_sp
{


namespace 
{

string peptide2stnd (const Peptide &pep)
{
  string s (pep. seq);
  replace (s, 'J', 'X');
  replace (s, 'U', 'X');
  replace (s, 'O', 'X');
  ASSERT (s. size () == pep. seq. size ());
  return s;
}



int pep2selfScore (const SNCBIFullScoreMatrix &mat,
	                 const Peptide &pep,
	                 size_t start,
	                 size_t stop)
{
  ASSERT (start <= stop);
  ASSERT (stop <= pep. seq. size ());
  
  const string seq (peptide2stnd (pep));
  
	int score = 0;
#ifndef NDEBUG
	bool ambig = false;
#endif
	FOR_START (size_t, i, start, stop)
	{
	  const char c = seq [i];
	  if (c == '-')
	    continue;
  #ifndef NDEBUG
	  if (pep. isAmbiguous (c))
	    ambig = true;
	#endif
		const int num = (int) c;
	  score += mat. s [num] [num];
	}
	IMPLY (! ambig, score >= 0);	

	return score;
}



#if 0
int pep2score (const SNCBIFullScoreMatrix &mat,
               int gap_open,
	             const Peptide &pep1,
	             const Peptide &pep2)
{
  ASSERT (gap_open <= 0);
  ASSERT (pep1. sparse);
  ASSERT (pep2. sparse);
  ASSERT (pep1. seq. size () == pep2. seq. size ());
  
	int score = 0;
	FOR_START (size_t, i, start, stop)
	{
    const char c1 = pep1. seq [i];
    const char c2 = pep2. seq [i];
    if (   c1 == '-' 
        && c2 == '-'
       )
      ;
    else if (   c1 == '-' 
             || c2 == '-'
            )
      score += gap_open;
    else
    {
  		const int num1 = (int) c1;
  		const int num2 = (int) c2;
  	  score += mat. s [num1] [num2];
    }

	}
	ASSERT (score >= 0);	

	return score;
}
#endif

}




// Align

Align::Align (const Peptide &pep1,
	            const Peptide &pep2,
	            bool semiglobal_arg,
	            size_t match_len_min,
	            bool blosum62)
: prot (true)
, semiglobal (semiglobal_arg)
{
  ASSERT (! pep1. sparse);
  ASSERT (! pep2. sparse);
  
  /*           ------ default -----
    matrix     gap_open  gap_extent   
    --------   --------  ----------   
		BLOSUM90         10           1
	  BLOSUM80         10           1
	  BLOSUM62         11           1 
	  BLOSUM50         13           2
		BLOSUM45         14           2
	  PAM250           14           2
	  PAM70            10           1
	  PAM30             9           1
	  IDENTITY         15           2
  */
  
  // PAR
	int gap_open   = 1;  
	int gap_extent = 1;
	if (blosum62)
	{
	  gap_open   = -11;  
		gap_extent =  -2; 
	}
	else
	{
	  gap_open   = -8;  
		gap_extent = -2;  
	}
  ASSERT (gap_open <= 0);
  ASSERT (gap_extent < 0);

#if 0
	SNCBIFullScoreMatrix mat;
	NCBISM_Unpack (blosum62 ? & SNCBIPackedScoreMatrix NCBISM_Blosum62 : & SNCBIPackedScoreMatrix NCBISM_Pam30, & mat);
#endif

	CNWAligner al ( peptide2stnd (pep1)
	              , peptide2stnd (pep2)
	              , blosum62 ? & NCBISM_Blosum62 : & NCBISM_Pam30
	              );
  al. SetWg (gap_open);
  al. SetWs (gap_extent);
	al. SetEndSpaceFree (semiglobal, semiglobal, semiglobal, semiglobal);
	score = al. Run ();
  tr = al. GetTranscriptString ();  
  
	finish (pep1, pep2, match_len_min, Peptide::stdMinComplexity);
		
  const SNCBIFullScoreMatrix& mat = al. GetScoreMatrix ();
	self_score1 = pep2selfScore (mat, pep1, start1, stop1);
	self_score2 = pep2selfScore (mat, pep2, start2, stop2);
}



Align::Align (const Dna &dna1,
	            const Dna &dna2,
	            bool semiglobal_arg,
	            size_t match_len_min,
	          //bool fast,
	            size_t band)
: prot (false)
, semiglobal (semiglobal_arg)
{
  ASSERT (! dna1. sparse);
  ASSERT (! dna2. sparse);
//IMPLY (band, ! fast);
  
  unique_ptr<CNWAligner> al;
/*if (fast)
  {
    // P(non-optimal) ~= 1e-4
    al. reset (new CMMAligner ());
  	al->EnableMultipleThreads ();
  }
  else*/ if (band)
  {
    auto al_ = new CBandAligner ();
    al. reset (al_);
    al_->SetBand (band);
  }
  else 
    al. reset (new CNWAligner ());
  ASSERT (al. get ());

#ifdef WU_BLASTN
	constexpr int match_score    =  5;
	constexpr int mismatch_score = -4;
	constexpr int gap_open       =  -1;  // was: 0 ??
	constexpr int gap_extent     = -10;
#else	
  // NCBI BLASTN
	constexpr int match_score    =  2;
	constexpr int mismatch_score = -3;
	constexpr int gap_open       = -5;
	constexpr int gap_extent     = -2;
#endif
  static_assert (match_score > 0, "match_score");
  static_assert (mismatch_score < 0, "mismatch_score");
  static_assert (gap_open <= 0, "gap_open");
  static_assert (gap_extent < 0, "gap_extent");	
  
  al->SetWm (match_score);
  al->SetWms (mismatch_score);
  al->SetWg (gap_open);
  al->SetWs (gap_extent);
  al->SetScoreMatrix (nullptr);
	
	{
		string seq1 (dna1. seq);
		string seq2 (dna2. seq);	
		strUpper (seq1);
		strUpper (seq2);
		al->SetSequences (seq1, seq2);
	}

	al->SetEndSpaceFree (semiglobal, semiglobal, semiglobal, semiglobal);
	score = al->Run ();		
  tr = al->GetTranscriptString();  	
	finish (dna1, dna2, match_len_min, Dna::stdMinComplexity);
	
	self_score1 = match_score * (int) (stop1 - start1);
	self_score2 = match_score * (int) (stop2 - start2);
	ASSERT (self_score1 >= 0);
	ASSERT (self_score2 >= 0);
}



void Align::finish (const Seq &seq1,
	                  const Seq &seq2,
	                  size_t match_len_min,
	                  double stdMinComplexity)
{
  ASSERT ((bool) match_len_min == semiglobal);
  ASSERT (! semiglobal_ends);

  const string& s1 = seq1. seq;
  const string& s2 = seq2. seq;

  ASSERT (! s1. empty ());
  ASSERT (! s2. empty ());
  ASSERT (! tr. empty ());

  matches       = strCountSet (tr, "M");
  substitutions = strCountSet (tr, "R");
  insertions    = strCountSet (tr, "I");
  deletions     = strCountSet (tr, "D");
  ASSERT (matches + substitutions + insertions + deletions == tr. size ());
  // Global alignment
  ASSERT (s1. size () == size1 ());
  ASSERT (s2. size () == size2 ());
  
	stop1 = s1. size ();
	stop2 = s2. size ();
  if (semiglobal)
  {
  	for (const char* c = & tr [0];                        *c == 'I'; c++)  { semiglobal_ends++; start2++; }
  	for (const char* c = & tr [0];                        *c == 'D'; c++)  { semiglobal_ends++; start1++; }
  	for (const char* c = & tr [tr. size () - 1]; stop2 && *c == 'I'; c--)  { semiglobal_ends++; stop2--; }
  	for (const char* c = & tr [tr. size () - 1]; stop1 && *c == 'D'; c--)  { semiglobal_ends++; stop1--; }
  }
	ASSERT (start1 <= stop1);
	ASSERT (start2 <= stop2);
	
	lowComplexity1 = seq1. getComplexityInt (start1, stop1) < stdMinComplexity;
	lowComplexity2 = seq2. getComplexityInt (start2, stop2) < stdMinComplexity;

  if (semiglobal)
  	if (   stop1 - start1 < match_len_min
	      || stop2 - start2 < match_len_min
	      || lowComplexity1 
	      || lowComplexity2
	     )
	  badMatch = true;
}



void Align::printStats (ostream &os) const
{ 
  os << "alignment_length:\t" << getAlignmentSize () << endl;
  
  #define P(x)  { os << #x << ":\t" << (x) << endl; }   
  P (score);
	P (matches);
  P (substitutions);
  P (insertions);
  P (deletions);
  os << "start1:\t" << start1 + 1 << endl;
  P (stop1);
  os << "start2:\t" << start2 + 1 << endl;
  P (stop2);
  P (self_score1);
  P (self_score2);
  P (lowComplexity1);
  P (lowComplexity2);
  P (badMatch);
  #undef P
  
  os << endl;
}



void Align::setAlignment (const string &seq1,
	                        const string &seq2) 
{
  ASSERT (! contains (seq1, '-'));
  ASSERT (! contains (seq2, '-'));

  sparse1.   reserve (tr. size ());
  sparse2.   reserve (tr. size ());
  consensus. reserve (tr. size ());
  size_t i1 = 0;
  size_t i2 = 0;
  for (const char c : tr)
  {
  	char c1 = '-';
  	char c2 = '-';
  	char cons = ' ';
  	switch (c)
  	{
  		case 'M': cons = '|';      // Match
  		          c1 = seq1 [i1];  
  			        c2 = seq2 [i2];
  			        break;
  		case 'R': c1 = seq1 [i1];  // Replacement
  			        c2 = seq2 [i2];
  			        break;
  		case 'I': c2 = seq2 [i2];  // Insertion
  			        break;
  		case 'D': c1 = seq1 [i1];  // Deletion
  			        break;
  		default : ERROR_MSG ("Unknown alignment trace: " + to_string (c));
  	}
    sparse1 += c1; 
  	sparse2 += c2; 
  	consensus += cons;
  	if (c1 != '-')
  		i1++;
  	if (c2 != '-')
  		i2++;
  }
  ASSERT (i1 == seq1. size ());
  ASSERT (i2 == seq2. size ());
  ASSERT (sparse1. size () == tr. size ());
  ASSERT (sparse2. size () == tr. size ());
}



void Align::printAlignment (ostream &os,
                            size_t line_len) const
{
  ASSERT (line_len);
  
  size_t qEnd = 0;
  size_t sEnd = 0;
  for (size_t i = 0; i < tr. size (); i += line_len)
  {
    const string qStr (sparse1.substr (i, line_len));
    const string sStr (sparse2.substr (i, line_len));
    qEnd += (qStr. size () - strCountSet (qStr, "-"));
    sEnd += (sStr. size () - strCountSet (sStr, "-"));
    os << qStr << ' ' << qEnd << endl;
    os << consensus. substr (i, line_len) << endl;
    os << sStr << ' ' << sEnd << endl;
    os << endl;
  }
}



Real Align::getDissim () const
{	
  if (self_score1 <= 0)
    return NaN;
  if (self_score2 <= 0)
    return NaN;
  
  if (badMatch)
  	return NaN;
	if (score <= 0)
		return inf;
		
  // Heuristic
	return intersection2dissim (self_score1, self_score2, score, 0, 0.5, true);  // PAR
}



Real Align::getMinEditDistance () const
{	
  ASSERT (self_score1 >= 0);
  ASSERT (self_score2 >= 0);
  
  if (badMatch)
  	return inf;
  const Real dist = self_score1 + self_score2 - 2 * score;
  ASSERT (dist >= 0.0);
	return dist;
}



size_t Align::getDiff () const
{
  ASSERT (! sparse1. empty ());
  ASSERT (sparse1. size () == sparse2. size ());

  size_t start = 0;
  size_t stop  = sparse1. size ();
  if (semiglobal)
  {
    while (   sparse1 [start] == '-'
           || sparse2 [start] == '-'
          )
      start++;
    while (   sparse1 [stop - 1] == '-'
           || sparse2 [stop - 1] == '-'
          )
    {
      ASSERT (stop);
      stop--;
    }
  }
  ASSERT (start < stop);
  
  size_t n = 0;
  FOR_START (size_t, i, start, stop)
    if (prot)
    {
      if (! aaMatch ( sparse1 [i]
                    , sparse2 [i]
                    )
         )
        n++;
    }
    else
    {
      if (! nucleotideMatch ( sparse1 [i]
                            , sparse2 [i]
                            )
         )
        n++;
    }
  
  return n;
}


const StringVector Align::distanceNames {"dissim", "min_edit", "rel_min_edit", "diff", "mismatch_frac"};
  


void Align::printDistances (ostream &os) const
{ 
  for (const Distance d : {dist_dissim, dist_min_edit, dist_rel_min_edit, dist_diff, dist_mismatch_frac})  
    os << distanceNames [d] << ":\t" << getDistance (d) << endl;
	os << "identity:\t" << getIdentity () << endl;

  os << endl;
}



}
