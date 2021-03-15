// blast_merge.cpp

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
*   Merge BLASTP HSPs
*
*/
   
   
#undef NDEBUG 
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
#include "../dm/matrix.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace
{
  
  
void findBoundary (const Vector<bool> &match,
                   size_t start,
                   size_t stop,
                   size_t &boundary,
                   Prob &leftIdent,
                   Prob &rightIdent)
// Output: boundary: between start and stop
//                   start of right sequence
{
  ASSERT (start <= stop);
  ASSERT (stop <= match. size ());
  
  const size_t len = stop - start;

  size_t ident = 0;
  FOR_START (size_t, i, start, stop)
    if (match [i])
      ident++;
  ASSERT (ident <= len);
  
  boundary = start;
  leftIdent = NaN;
  rightIdent = (Real) ident / (Real) len;

  Real chi2_max = 0.0;
  Matrix mat (2, 0.0);
  mat. put (false, 1, 0, (Real) ident);
  mat. put (false, 1, 1, (Real) (len - ident));
  FOR_START (size_t, i, start + 1, stop)
  {
    const Real chi2 = mat. getChi2 ();
    if (maximize (chi2_max, chi2))
    {
      boundary = i;
      leftIdent  = mat. get (false, 0, 0) / mat. sumRow (false, 0);
      rightIdent = mat. get (false, 1, 0) / mat. sumRow (false, 1);
    }
    const size_t col = ! match [i];
    mat. putInc (false, 0, col,  1.0);
    mat. putInc (false, 1, col, -1.0);
  }
}
  
  

struct Blastp 
{
  string targetName;
  string refName; 
  size_t length {0}, nident {0}  // aa
       , targetStart {0}, targetStop {0}, targetLen {0}
       ,    refStart {0},    refStop {0},    refLen {0};
    // Positions are 0-based
    // start < stop
  string targetSeq;
  string refSeq;
  
  bool targetStrand {false}; 
  Prob ident_frac {NaN};
  

  explicit Blastp (const string &line)
    {
      static Istringstream iss;
	    iss. reset (line);
	    iss >> targetName >> refName >> length >> nident >> targetStart >> targetStop >> targetLen >> refStart >> refStop >> refLen >> targetSeq >> refSeq;
	  // format:  qseqid      sseqid    length    nident         qstart         qend         qlen      sstart      send      slen         sseq    qseq
    // blastp:  ...         ...          663       169              2          600          639           9       665       693          ...
	    ASSERT (! targetSeq. empty ());	
	    ASSERT (targetSeq. size () == refSeq. size ());
	    ASSERT (length == targetSeq. size ());

	    ASSERT (refStart < refStop);  

	    ASSERT (targetStart != targetStop);
	    targetStrand = targetStart < targetStop;  
	    if (! targetStrand)
	      swap (targetStart, targetStop);
	      
	    ASSERT (refStart >= 1);
	    ASSERT (targetStart >= 1);
	    ASSERT (refStart < refStop);
	    ASSERT (targetStart < targetStop);
	    refStart--;
	    targetStart--;

      ident_frac = (Real) nident / (Real) length;
    }    
//Blastp (Blastp &&) = default;
//Blastp& operator= (Blastp &&) = default;
  void saveText (ostream &os) const
    { os         << targetName 
         << '\t' << refName 
         << '\t' << length 
         << '\t' << nident 
         << '\t' << (targetStrand ? targetStart + 1 : targetStop)
         << '\t' << (targetStrand ? targetStop : targetStart + 1)
         << '\t' << targetLen 
         << '\t' << refStart + 1
         << '\t' << refStop 
         << '\t' << refLen 
         << '\t' << targetSeq 
         << '\t' << refSeq
         << endl;
    }
  void saveSummary (ostream &os) const
    { const ONumber on (os, 2, false);
      os         << targetName 
         << '\t' << refName 
         << '\t' << ident_frac * 100.0      
         << '\t' << pTargetCoverage ()
         << '\t' << pRefCoverage ()
         << endl;
    }
    
    
  size_t targetCoverage () const
    { return targetStop - targetStart; }
  size_t refCoverage () const
    { return refStop - refStart; }
  Prob pTargetCoverage () const
    { return 100.0 * (Real) targetCoverage () / (Real) targetLen; }
  Prob pRefCoverage () const
    { return 100.0 * (Real) refCoverage () / (Real) refLen; }
  bool operator< (const Blastp &other) const
    { 
      LESS_PART (other, *this, refCoverage ());
      LESS_PART (*this, other, refStart); 
      LESS_PART (*this, other, refName); 
      return false;
    }
    

  void trim (Prob ident_max)
    {
      ASSERT (isProb (ident_max));
      
      Vector<bool> match (targetSeq. size ());
	    FFOR (size_t, i, targetSeq. size ())
	      match [i] = (targetSeq [i] == refSeq [i]);
	      
	    size_t start = 0;
	    size_t stop = match. size ();
	    for (;;)
	    {
	      while (start < stop && ! match [start])
	        start++;
	      while (start < stop && ! match [stop - 1])
	        stop--;
  	    size_t boundary = start;
  	    Prob leftIdent = NaN;
  	    Prob rightIdent = NaN;
  	    findBoundary (match, start, stop, boundary, leftIdent, rightIdent);
  	    if (   isNan (leftIdent)
  	        || isNan (rightIdent)
  	       )
  	      break;
  	    ASSERT (isProb (leftIdent));
  	    ASSERT (isProb (rightIdent));
  	    if (leftIdent < ident_max)
  	    {
  	      ASSERT (start < boundary)
	        start = boundary;
  	    }
  	    else if (rightIdent < ident_max)
  	    {
	        ASSERT (stop > boundary);
          stop = boundary;
  	    }
  	    else
  	      break;
  	  }
  	  ASSERT (start <= stop);
  	  ASSERT (stop <= targetSeq. size ());
  	  
  	  length = stop - start;
  	  
  	  nident = 0;
  	  FOR_START (size_t, i, start, stop)
  	    if (match [i])
  	      nident++;
  	      
      ident_frac = (Real) nident / (Real) length;

  	  FOR (size_t, i, start)
  	  {
  	    if (targetSeq [i] != '-')
  	      targetStart++;
  	    if (refSeq [i] != '-')
  	      refStart++;
  	  }
  	  
  	  FOR_REV_END (size_t, i, stop, targetSeq. size ())
  	  {
  	    if (targetSeq [i] != '-')
  	      targetStop--;
  	    if (refSeq [i] != '-')
  	      refStop--;
  	  }
  	  
  	  targetSeq. erase (stop);
  	  refSeq.    erase (stop);
  	  targetSeq. erase (0, start);
  	  refSeq.    erase (0, start);
	    ASSERT (targetSeq. size () == refSeq. size ());
	    ASSERT (targetSeq. size () == length);
    }
    
    
  bool refMerge (const Blastp &other) 
    {
      ASSERT (& other != this);
      ASSERT (targetName == other. targetName);
      ASSERT (refName    == other. refName);
      ASSERT (targetLen  == other. targetLen);
      ASSERT (refLen     == other. refLen);
      ASSERT (refCoverage () >= other. refCoverage ());
      
      if (targetStrand != other. targetStrand)
        return false;
      if (refStop < other. refStart)
        return false;
      if (other. refStop < refStart)
        return false;
      if (targetStop < other. targetStart)
        return false;
      if (other. targetStop < targetStart)
        return false;
        
      minimize (refStart, other. refStart);
      maximize (refStop,  other. refStop);  
      minimize (targetStart, other. targetStart);
      maximize (targetStop,  other. targetStop);  
      
      minimize (ident_frac, other. ident_frac);
      
      length = 0;
      nident = 0;
      targetSeq. clear ();
      refSeq.    clear ();
        
      return true;
    }
};



void process (Vector<Blastp> &als,
              Prob ident_max,
              size_t coverage_min)
// Update: als
{
  if (als. size () > 1)
  {
    // QC
    for (Blastp& blastp : als)
    {
      ASSERT (blastp. targetName == als. front (). targetName);
      ASSERT (blastp. refName    == als. front (). refName);
      ASSERT (blastp. targetLen  == als. front (). targetLen);
      ASSERT (blastp. refLen     == als. front (). refLen);
    }
      
    for (Blastp& blastp : als)
      blastp. trim (ident_max);

    for (Iter<Vector<Blastp>> it (als); it. next ();)
      if (it->refCoverage () < coverage_min)
        it. erase ();          
    
    als. sort ();
    
    FOR (size_t, i, als. size ())
      FOR_REV_END (size_t, j, i + 1, als. size ())
        if (als [i]. refMerge (als [j]))
          als. eraseAt (j);
  }  

	for (const Blastp& blastp : als)
	  blastp. saveSummary (cout);
}



} // namespace
 



// ThisApplication

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("For multiple BLASTP HSPs with the same queries and subjects: trim bad ends and merge the HSPs\nPrint: qseqid sseqid ident_frac p_target_coverage p_ref_coverage")
    {
      version = VERSION;
      // Input
      const string blastFormat ("qseqid sseqid length nident qstart qend qlen sstart send slen sseq qseq");
      addPositional ("blastp", "blastp output in the format: " + blastFormat + "\nSorted by qseqid, sseqid"); 
      addPositional ("ident_max", "Max. identity of bad ends of an HSP");
      addPositional ("coverage_min", "Min. coverage of a trimmed HSP");
    }



  void body () const final
  {
    const string blastpFName  =               getArg ("blastp");
    const Prob   ident_max    = str2real     (getArg ("ident_max"));
    const size_t coverage_min = str2<size_t> (getArg ("coverage_min"));
    
    
    QC_ASSERT (isProb (ident_max));
    

    LineInput f (blastpFName);
    Vector<Blastp> als;
	  while (f. nextLine ())
      try
  	  {
  	    Blastp al (f. line);
  	    if (   ! als. empty () 
  	        && ! (   als. front (). targetName == al. targetName
  	              && als. front (). refName    == al. refName
  	             )
  	       )
  	    {
    	    QC_ASSERT (als. front (). targetName <= al. targetName);
    	    QC_IMPLY (als. front (). targetName == al. targetName, als. front (). refName < al. refName);
  	      process (als, ident_max, coverage_min);
  	      als. clear ();
  	    }
  	    als << move (al);
  	  }
		  catch (...)
		  {
		  	cout << f. line << endl;
		  	throw;
		  }
		process (als, ident_max, coverage_min);		
  }
};




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



