// blastp2exons.cpp

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
*   Find exons
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{
  
  
struct Align
{
  string seq;
    // aa
  Vector<size_t> matches;
  Vector<bool> intron;
  // seq.size() = matches.size() = intron.size()
  Set<string> seqids;
  
  
  Align () = default;
  explicit Align (const size_t len)
    : seq (string (len, ' '))
    , matches (len, 0)
    , intron (len, false)
    { QC_ASSERT (len); }
    
    
  bool empty () const
    { return seq. empty (); }
  void setAa (size_t pos, 
              char aa)
    { if (seq [pos ] == ' ')
        seq [pos] = aa;
      QC_ASSERT (seq [pos] == aa);
    }
  void printStat () const
    { FFOR (size_t, i, seq. size ())
	      cout << i + 1 << '\t' << seq [i] << '\t' << matches [i] << '\n';
    }  
  void setIntron (size_t window,
                  double threshold)
    { FFOR (size_t, i, seq. size ())
        intron [i] = (center2matches (i, window) < threshold * (double) seqids. size ());
    }
  void printExons (ostream &os) const
    { FFOR (size_t, i, seq. size ())
        if (   seq [i] != ' ' 
            && ! intron [i]
           )
          os << seq [i];
      os << '\n';
    }
private:
  double interval2matches (int start,
                           size_t stop) const
    { ASSERT (start < (int) stop);
      double s = 0.0;
      double n = 0;
      FOR_START (size_t, i, (size_t) max (0, start), min (matches. size (), stop))
      { s += (double) matches [i];
        n += 1.0;
      }
      return s / n;
    }
  double center2matches (size_t center,
                         size_t window) const
    { const int start = (int) center - (int) window / 2;
      return interval2matches (start, size_t (start + (int) window));
    }
};

  
  
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Find exons")
    {
      version = VERSION;
  	  addPositional ("blastp", "blastp output in format: qseqid sseqid qstart qend qlen qseq sseq");
  	  addKey ("stat", "Output statistcis file");
    }



	void body () const final
  {
	  const string blastpFName = getArg ("blastp");
	  const string statFName   = getArg ("stat");
  
  
    map<string/*qseqid*/, Align> seqid2align;
	  {
  	  LineInput in (blastpFName, 10000);  // PAR
  	  Istringstream iss;
  	  while (in. nextLine ())
  	  {
  	    string qseqid, sseqid, qseq, sseq;
  	    size_t qstart, qend, qlen;
  	    iss. reset (in. line);
  	    iss >> qseqid >> sseqid >> qstart >> qend >> qlen >> qseq >> sseq;
  	    QC_ASSERT (! sseq. empty ());  // in.line is not truncated
  	    QC_ASSERT (qlen);
  	    QC_ASSERT (qstart);
  	    qstart--;
  	    QC_ASSERT (qstart < qend);
  	    QC_ASSERT (qend <= qlen);
  	    QC_ASSERT (qseq. size () == sseq. size ());
  	    QC_ASSERT (qseq. size () >= qend - qstart);  // Can contain '-'
  	    QC_ASSERT (qseq. front () != '-');
  	    QC_ASSERT (qseq. back  () != '-');
  	    if (qseqid == sseqid)
  	      continue;
  	    Align& align = seqid2align [qseqid];
  	    if (align. empty ())
  	      align = std::move (Align (qlen));
  	    align. seqids << sseqid;
  	    QC_ASSERT (align. seq. size () == qlen);
  	    size_t i = qstart;
  	    FFOR (size_t, j, qseq. size ())
  	      if (qseq [j] != '-')
    	    { 
    	      if (qseq [j] == sseq [j])
    	        align. matches [i] ++;
    	      align. setAa (i, qseq [j]);
    	      i++;
    	    }
  	  }
  	}
  	
  	
  	
  	unique_ptr<OFStream> statF;
    if (! statFName. empty ())
    	statF. reset (new OFStream (statFName));
  	for (auto& it : seqid2align)
	  {
	    Align& align = it. second;
    	if (statF)
    	{
  	    *statF << it. first << '\n'; 
  	    align. printStat ();
  	    *statF << '\n';
  	  }
  	  cout << '>' << it. first << '\n';   	  
  	  align. setIntron (30, 0.2);  // PAR
  	  align. printExons (cout);
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



