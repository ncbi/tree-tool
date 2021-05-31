// align.hpp

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



#ifndef ALIGN_HPP
#define ALIGN_HPP


#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "../genetics/seq.hpp"
using namespace Seq_sp;



namespace Align_sp
{
	
	
struct Align : Root
// Global alignment
// Needleman-Wunsch algorithm 
{
  const bool semiglobal;
  
	// Output:
	int score {0};
	string tr;
	  // "Transcript"
	size_t matches {0};
	  // = identities
	size_t substitutions {0};
	size_t insertions {0};
		// In seq1
	size_t deletions {0};
		// = insertions in seq2
	
	// Positions in seq1/seq2
	size_t start1 {0};
	size_t start2 {0};
  size_t stop1 {0};
	size_t stop2 {0};
	
	int self_score1 {0};
	int self_score2 {0};
	bool lowComplexity1 {false};
	bool lowComplexity2 {false};
	
	bool badMatch {false};
	
	// Alignment
  string sparse1;    
  string sparse2;    
  string consensus;  

	
  // Input: match_len_min: valid if semiglobal, otherwise 0
  // !semiglobal = global
	Align (const Peptide &pep1,
	       const Peptide &pep2,
	       bool semiglobal_arg,
	       size_t match_len_min,
	       bool blosum62);
	  // Input: blosum62: false <=> PAM30
	Align (const Dna &dna1,
	       const Dna &dna2,
	       bool semiglobal_arg,
	       size_t match_len_min,
	       size_t band);
	  // Input: band => use banded aligner
	void saveText (ostream &os) const override
	  { os << tr << endl
	  	   <<        score 
	  	   << ' ' << matches 
	  	   << ' ' << substitutions 
	  	   << ' ' << insertions 
	  	   << ' ' << deletions 
	  	   << ' ' << start1 << '-' << stop1 
	  	   << ' ' << start2 << '-' << stop2
	  	   << ' ' << self_score1
	  	   << ' ' << self_score2
	  	   << ' ' << lowComplexity1
	  	   << ' ' << lowComplexity2
	  	   << ' ' << badMatch
	  	   << endl;
	  }
private:
	void finish (const Seq &seq1,
	             const Seq &seq2,
	             size_t match_len_min,
	             double stdMinComplexity);
public:
	
	
	// Sizes of the input sequences
  size_t size1 () const
    { return tr. size () - insertions; }
  size_t size2 () const
    { return tr. size () - deletions; }
  // Return: !isNan() => >= 0
	Real getDissim () const;
	Real getMinEditDistance () const;
	Real getMismatchFrac () const
	  { return (Real) (tr. size () - matches) / (Real) tr. size (); }
  //
  void setAlignment (const string &seq1,
	                   const string &seq2);
    // Output: Alignment
  // Input: Alignment
  void printAlignment (size_t line_len) const;
	size_t getDiff () const;
};
	


}



#endif
