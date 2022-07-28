// blast_best_hits.cpp

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
*   Find best BLAST HSPs
*
*/
   
   
#undef NDEBUG 
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



namespace
{
  
  
struct Blast : Root
{
  string name;
  string targetName;
  string refName; 
  size_t length {0}, nident {0}  // aa
       , targetStart {0}, targetStop {0}
       , refStart {0}, refStop {0}, refLen {0};
    // Positions are 0-based
    // start < stop
  string refTitle;
  bool targetStrand {false}; 
  bool refStrand {false}; 
  

  Blast () = default;
  Blast (const Blast&) = default;
  Blast (Blast&&) = default;
  Blast& operator= (Blast &&) = default;
  Blast (const string &line,
         bool nameP)
    {
      static Istringstream iss;
	    iss. reset (line);
	    if (nameP)
	    {
	      iss >> name;
	      QC_ASSERT (! name. empty ());
	    }
	    iss >> targetName >> refName >> length >> nident >> targetStart >> targetStop >> refStart >> refStop >> refLen; 
	  // format:  qseqid      sseqid    length    nident         qstart         qend     sstart      send       slen

      processStartStop (targetStart, targetStop, targetStrand);
      processStartStop (refStart, refStop, refStrand);
      QC_ASSERT (refStop <= refLen);
	    
	    refTitle = line. substr ((size_t) iss. tellg () + 1);  // stitle
	    trim (refTitle);	    
	    trimPrefix (refTitle, refName + " ");
	    trim (refTitle);
    }    
private:
  static void processStartStop (size_t &start,
                                size_t &stop,
                                bool &strand)
    { QC_ASSERT (start != stop);
	    strand = start < stop;  
	    if (! strand)
	      swap (start, stop);	      
	    QC_ASSERT (start >= 1);
	    ASSERT (start < stop);
	    start--;
    }
public:
  void qc () const final
    { if (! qc_on)
        return;
      QC_ASSERT (! targetName. empty ());
      QC_ASSERT (! refName. empty ());
      QC_ASSERT (nident <= length);
      QC_ASSERT (nident);
      QC_ASSERT (targetStart < targetStop);
    }
  void saveText (ostream &os) const final
    { if (! name. empty ())
        os << name << '\t';
      os         << targetName 
         << '\t' << refName 
         << '\t' << length 
         << '\t' << nident 
         << '\t' << (targetStrand ? targetStart + 1 : targetStop)
         << '\t' << (targetStrand ? targetStop : targetStart + 1)
         << '\t' << (refStrand ? refStart + 1 : refStop)
         << '\t' << (refStrand ? refStop : refStart + 1)
         << '\t' << refLen
         << '\t' << refTitle
         << endl;
    }
  bool empty () const final
    { return targetName. empty (); }
    

  double ident_frac () const
    { return (double) nident / (double) length; }
  size_t refcoverage () const
    { return refStop - refStart; }
  double refcoverage_frac () const
    { return (double) refcoverage () / (double) refLen; }
  void saveSummary (ostream &os) const
    { const ONumber on (os, 2, false);
      if (! name. empty ())
        os << name << '\t';
      os         << targetName 
         << '\t' << refName 
         << '\t' << ident_frac () * 100.0      
         << '\t' << targetStart + 1
         << '\t' << targetStop
         << '\t' << targetStrand  
         << '\t' << refStart + 1
         << '\t' << refStop
         << '\t' << refStrand  
         << '\t' << refLen  
         << '\t' << refcoverage_frac () * 100.0      
         << '\t' << refTitle
         << endl;
    }
  bool operator< (const Blast &other) const
    { LESS_PART (*this, other, targetStart);
      LESS_PART (other, *this, targetStop);
      return false;
    }
};



void process (Vector<Blast> &als)
{
  als. sort ();

  const Blast* prev = nullptr;
  for (const Blast& al : als)
  {
    if (! prev || prev->targetStop < al. targetStop)
	    al. saveSummary (cout);
	  if (prev)
	  {
	    ASSERT (prev->name       == al. name);
	    ASSERT (prev->targetName == al. targetName);
	  }
	  prev = & al;
	}

  als. clear ();
}



} // namespace
 



// ThisApplication

static const string header ("qseqid\tsseqid\tpident\tqstart\tqend\tqstrand\tsstart\tsend\tsstrand\tslen\tpscoverage\tstitle");



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find BLAST HSPs with largest coverage\nPrint tsv-file with columns: " + header)
    {
      version = VERSION;
      // Input
      const string blastFormat ("qseqid sseqid length nident qstart qend sstart send slen stitle");
      addPositional ("blast", "BLAST output in the format: " + blastFormat + "\nGrouped by name, qseqid"); 
      addPositional ("pident_min", "Min. identity percent of HSPs");
      addPositional ("coverage_min", "Min. coverage of HSP reference to report");
      addPositional ("pcoverage_min", "Min. coverage percent of HSP reference to report");
      addFlag ("name", "First column of <blast> is the query name");
    }



  void body () const final
  {
    const string blastFName        =               getArg ("blast");
    const double ident_frac_min    = str2<double> (getArg ("pident_min")) / 100.0;
    const size_t coverage_min      = str2<size_t> (getArg ("coverage_min"));    
    const double coverage_frac_min = str2<double> (getArg ("pcoverage_min")) / 100.0;  
    const bool nameP               =               getFlag ("name");  
    
    QC_ASSERT (ident_frac_min >= 0.0);
    QC_ASSERT (ident_frac_min <= 1.0);
    
    QC_ASSERT (coverage_frac_min >= 0.0);
    QC_ASSERT (coverage_frac_min <= 1.0);


    LineInput f (blastFName);
    try
    {
      cout << '#';
      if (nameP)
        cout << "name\t";
      cout << header << endl;
      Vector<Blast> als;
  	  while (f. nextLine ())
  	  {
  	    Blast al (f. line, nameP);
  	    al. qc ();
	    //QC_IMPLY (! als. empty (), als. back (). targetName <= al. targetName);
  	    if (al. ident_frac () <= ident_frac_min)
  	      continue;
  	    if (al. refcoverage () < coverage_min)
  	      continue;
  	    if (al. refcoverage_frac () < coverage_frac_min)
  	      continue;
  	    if (   ! als. empty () 
  	        && (   als. back (). name != al. name 
  	            || als. back (). targetName != al. targetName)
  	       )
  	      process (als);
	      als << move (al);
  	    ASSERT (al. empty ());
  	  }
      process (als);
  	}
	  catch (...)
	  {
	    cout << "Line # " << f. lineNum << endl;
	  	cout << f. line << endl;
	  	throw;
	  }
  }
};




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



