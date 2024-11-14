// tblastn2orf.cpp

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
*   Find best annotation and print proteins
*
*/


#undef NDEBUG

#include "../common.hpp"
#include "../graph.hpp"
using namespace Common_sp;
#include "seq.hpp"  
using namespace Seq_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{
  
  
SubstMat sm;
unique_ptr<OFStream> gapF;


struct Exon;
map <Pair<string>/*ref,contig*/, VectorPtr<Exon>> contig2exons;



struct Intron;
  

  
// Merge with other struct Alignment's ??
struct Exon final : DiGraph::Node
{
  // Input
	string qseqid;  // reference protein
	string sseqid;  // contig accession
	size_t qstart {0}, qend {0};  // aa
	size_t sstart {0}, send {0};  // nt
	  // sstart < send
	string qseq;
	string sseq;  
	  // Can contain '*'
	Strand strand {0};
	AlignScore score {0};
	// Output
	AlignScore totalScore {0};
	const Intron* bestIntron {nullptr};


  Exon (DiGraph &graph_arg,
        const string &qseqid_arg,
      	const string &sseqid_arg,
      	size_t qstart_arg,
      	size_t qend_arg,
      	size_t sstart_arg,
      	size_t send_arg,
      	const string &qseq_arg,
      	const string &sseq_arg,
      	Strand strand_arg)
    : DiGraph::Node (graph_arg)
    , qseqid (qseqid_arg)
    , sseqid (sseqid_arg)
    , qstart (qstart_arg)
    , qend (qend_arg)
    , sstart (sstart_arg)
    , send (send_arg)  
    , qseq (qseq_arg)
    , sseq (sseq_arg)
    , strand (strand_arg)
    { 
      finish (); 
    }	
    
    
  Exon (DiGraph &graph_arg,
        const string &line)
    : DiGraph::Node (graph_arg)
    {
	  	QC_ASSERT (! line. empty ());
	  	{
  	  	istringstream iss (line);
  		  iss >> qseqid >> sseqid >> qstart >> qend >> sstart >> send >> qseq >> sseq;
  		}
	  	QC_ASSERT (! sseq. empty ());  // line is not truncated

		  QC_ASSERT (qstart);
		  qstart--;

		  QC_ASSERT (sstart);
		  QC_ASSERT (send);
		  QC_ASSERT (sstart != send);
		  strand = 1;
		  if (sstart < send)
		  	sstart--;
		  else
		  {
		  	send--;
		  	strand = -1;
		  	swap (sstart, send);
		  }  		  

      finish ();
    }
    
    
private:
  void finish ();
    // Output: bestIntron
public:
  void saveText (ostream &os) const final;
  
  
  void qc () const final
    {
      if (! qc_on)
        return;
        
      DiGraph::Node::qc ();
        
	    QC_ASSERT (! qseqid. empty ());	    
	    QC_ASSERT (! sseqid. empty ());	    

		  QC_ASSERT (qstart < qend);
		  QC_ASSERT (sstart < send);
	    QC_ASSERT (divisible ((uint) (send - sstart), 3));	  

		  QC_ASSERT (qseq. size () == sseq. size ());		  
	  	QC_ASSERT (! qseq. empty ()); 
	    QC_ASSERT (qLen () <= qseq. size ());
	    QC_ASSERT (sLen () <= sseq. size ());
	    
	    // Needed for Intron::Intron 
	    QC_ASSERT (qseq. front () != '-');
	    QC_ASSERT (qseq. back  () != '-');
	    QC_ASSERT (sseq. front () != '-');
	    QC_ASSERT (sseq. back  () != '-');
	    
	    QC_ASSERT (score >= 0);
	    QC_ASSERT (totalScore >= 0);
    }
    
    
  size_t qLen () const
    { return qend - qstart; }
  size_t sLen () const
    { return (send - sstart) / 3; }
    // aa
  size_t qCenter () const  // --> center of match density ??
    { return (qstart + qend) / 2; }
  bool arcable (const Exon &next) const
    {
      ASSERT (qseqid == next. qseqid);
      ASSERT (sseqid == next. sseqid);
      
      // PAR
      constexpr size_t intron_max = 10000;  // nt 
      
      if (strand != next. strand)
        return false;
    //if (same frame and overlap) return false;  // ??
      if (qCenter () >= next. qCenter ())
        return false;  // => DAG
      if (qend + 20 < next. qstart)
        return false;
      if (   strand == 1 
          && send + intron_max < next. sstart
         )
        return false;
      if (   strand == -1 
          && next. sstart + intron_max < send
         )
        return false;
         
      return true;
    }
    
    
  void setBestIntron ();
    // Update: totalScore
    // Output: bestIntron
  string getSeq (size_t start) const;
};



struct Intron final : DiGraph::Arc
{
  AlignScore score {0};
    // Partial intron score
    // Minimized
  // In Exon::sseq
  size_t prev_end {no_index};
  size_t next_start {0};
  
  
  Intron (Exon* prev,
          Exon* next)
    : DiGraph::Arc (prev, next)
    , prev_end (prev->qseq. size ())
    {
      ASSERT (prev);
      ASSERT (next);
      
      const size_t start = max (next->qstart, prev->qCenter ());
      const size_t end   = min (prev->qend,   next->qCenter ());
      if (start >= end)
        return;
      
      // [start, end) = overlap of prev->qseq and next->qseq  
      const size_t len = end - start;  
      ASSERT (len);
      ASSERT (len <= prev->qLen ());  
      ASSERT (len <= next->qLen ());  
            
      Vector<AlignScore> prevScores;  prevScores. reserve (len + 1);
      {
        size_t pos = prev->qstart;
        FFOR (size_t, i, prev->qseq. size ())
          if (prev->qseq [i] != '-')
          {
            if (between (pos, start, end))
              prevScores << sm. char2score ( prev->qseq [i]
                                           , prev->sseq [i]
                                           );
            pos++;
          }
      }
      prevScores << 0;
      ASSERT (prevScores. size () == len + 1);
      
      Vector<AlignScore> nextScores;  nextScores. reserve (len + 1);
      nextScores << 0;
      {
        size_t pos = next->qstart;
        FFOR (size_t, i, next->qseq. size ())
          if (next->qseq [i] != '-')
          {
            if (between (pos, start, end))
              nextScores << sm. char2score ( next->qseq [i]
                                           , next->sseq [i]
                                           );
            pos++;
          }
      }
      ASSERT (nextScores. size () == len + 1);
      
      FOR_REV (size_t, i, len)
        prevScores [i] += prevScores [i + 1];
      FOR_START (size_t, i, 1, len + 1)
        nextScores [i] += nextScores [i - 1];
        
      size_t bestSuff = no_index;
      score = score_inf; 
      FOR (size_t, i, len + 1)
        if (minimize (score, prevScores [i] + nextScores [i]))
          bestSuff = i;
      ASSERT (score != score_inf);
      ASSERT (bestSuff != no_index);
      ASSERT (bestSuff <= len);
      
      const size_t split = start + bestSuff;

      // prev_end
      {
        size_t pos = prev->qstart;
        FFOR (size_t, i, prev->qseq. size ())
          if (prev->qseq [i] != '-')
          {
            if (pos == split)
            {
              prev_end = i;
              break;
            }
            pos++;
          }
      }
      ASSERT (prev_end <= prev->qseq. size ());
      
      // next_start
      {
        size_t pos = next->qstart;
        FFOR (size_t, i, next->qseq. size ())
          if (next->qseq [i] != '-')
          {
            if (pos == split)
            {
              next_start = i;
              break;
            }
            pos++;
          }
      }
      ASSERT (next_start <= next->qseq. size ());
    }
    
    
  void saveText (ostream &os) const final
    {
      os << '\t' << prev_end 
         << '\t' << next_start
         << '\t' << score 
         << '\t';
      ASSERT (node [true]);
      static_cast <const Exon*> (node [true]) -> saveText (os);
    }
    
    
  AlignScore getTotalScore ()
    {
      const Exon* next = static_cast <const Exon*> (node [true]);
      ASSERT (next);
      var_cast (next) -> setBestIntron ();  // DAG => no loop 
      return next->score - score;
    }
};



void Exon::finish ()
{
  ASSERT (graph);
  ASSERT (score == 0);
  ASSERT (! bestIntron);
  
  size_t gap = 0;
  size_t qstart_new = qstart;
  size_t sstart_new = sstart;
  size_t send_new   = send;
  FFOR (size_t, i, qseq. size ())
  {
    if (qseq [i] == '-')
    {
      gap++;
      QC_ASSERT (sseq [i] != '-');
    }
    else
    {
      QC_ASSERT (qseq [i] != '*');
      // sseq: prohibit '*' in a non-intron gap ??
      if (gap >= 15)  // PAR  
      {
        auto exon = new Exon ( * var_cast (graph)
                             , qseqid
                             , sseqid
                             , qstart_new
                             , qend
                             , strand == 1 ? sstart_new : sstart
                             , strand == 1 ? send       : send_new
                             , qseq. substr (i)
                             , sseq. substr (i)
                             , strand
                             );
        qend = qstart_new;
        const size_t intron_len = gap * 3;
        if (strand == 1)
        {
          QC_ASSERT (sstart_new > intron_len);
          send = sstart_new - intron_len;
        }
        else
          sstart = send_new + intron_len;
        ASSERT (i >= gap);
        qseq. erase (i - gap);
        sseq. erase (i - gap);
        bestIntron = new Intron (this, exon);  
        break;
      }
      qstart_new++;
      gap = 0;
    }

    // sstart_new, send_new
    if (sseq [i] != '-')
    {
      if (strand == 1)
        sstart_new += 3;
      else
      {
        QC_ASSERT (send_new >= 3);
        send_new -= 3;
      }
    }
  }
  
  FFOR (size_t, i, qseq. size ())
    score += sm. char2score (qseq [i], sseq [i]);
  QC_ASSERT (score > 0);

  const Pair<string> p (qseqid, sseqid);
  contig2exons [p] << this;
  
  qc ();
}



void Exon::saveText (ostream &os) const 
{
  os         << qstart + 1 << '(' << sstart + 1 << ')'
     << '\t' << qend   + 1 << '(' << send   + 1 << ')'
     << '\t' << (int) strand
     << '\t' << score 
     << '\t' << totalScore;
  if (bestIntron)
    bestIntron->saveText (os);
}



void Exon::setBestIntron ()
{
  ASSERT (score > 0);
  ASSERT (totalScore >= 0);
  
  if (totalScore > 0)
    return;

  if (arcs [true]. empty ())
    totalScore = score;
  else
  {
    if (bestIntron)
      totalScore = score + var_cast (bestIntron) -> getTotalScore ();
    else
      for (const DiGraph::Arc* arc : arcs [true])
      {
        ASSERT (arc);
        ASSERT (arc->node [false] == this);
        const Intron* intron = static_cast <const Intron*> (arc);
        if (maximize (totalScore, score + var_cast (intron) -> getTotalScore ()))
          bestIntron = intron;
      }
    ASSERT (bestIntron);
  }
  QC_ASSERT (totalScore > 0);
}



string Exon::getSeq (size_t start) const
{
  if (gapF)
  {
    size_t gap = 0;
    for (const char c : qseq)
      if (c == '-')
        gap++;
      else
      {
        if (gap > 3)  // PAR
          *gapF << gap << '\n';
        gap = 0;
      }
  }
  
  // s
  const size_t end = bestIntron ? bestIntron->prev_end : sseq. size (); 
  ASSERT (start <= end);
  ASSERT (end <= sseq. size ());
  string s = sseq. substr (start, end - start);
  replaceStr (s, "-", noString);
  
  if (bestIntron)
  {
    const Exon* next = static_cast <const Exon*> (bestIntron->node [true]);
    ASSERT (next);
    s += next->getSeq (bestIntron->next_start);
  }
  
  return s;
}




struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Find best annotation and print proteins")
    {
      version = VERSION;
  	  addPositional ("tblastn", "tblastn output in format: qseqid sseqid qstart qend sstart send qseq sseq. qseqid = <variant>-<gene>, where <gene> has no dashes");
  	  addKey ("gap_stats", "Output file with gap lengths");
    }



	void body () const final
  {
	  const string tblastnFName = getArg ("tblastn");
	  const string gap_stats    = getArg ("gap_stats");
  
  
    sm = SubstMat (execDir + "/matrix/BLOSUM62");  // PAR
    sm. qc ();
  
  
	  DiGraph graph;
	  {
  	  LineInput in (tblastnFName);
  	  while (in. nextLine ())
  	    try { new Exon (graph, in. line); }
  	      catch (const exception &e)
  	        { throw runtime_error (string (e. what ()) + "\n" + in. lineStr ()); }
  	}
  	graph. qc ();
  	
  	
  	map<string/*gene*/,const Exon*> gene2exon;
  	for (const auto& it : contig2exons)
  	{
  	  ASSERT (! it. second. empty ());
  	  
  	  // new Intron
  	  for (const Exon* next : it. second)
  	  {
  	    ASSERT (next);
    	  for (const Exon* prev : it. second)
    	    if (prev == next)
    	      break;
    	    else
    	      if (prev->arcable (*next))
    	        new Intron (var_cast (prev), var_cast (next));
    	}
    	
    	AlignScore totalScore_max = 0;
    	const Exon* bestExon = nullptr;
  	  for (const Exon* exon : it. second)
  	  {
  	    var_cast (exon) -> setBestIntron ();
  	    ASSERT (exon->totalScore > 0);
  	    if (maximize (totalScore_max, exon->totalScore))
  	      bestExon = exon;
  	  }
  	  ASSERT (bestExon);  	  
  	  
  	  const string& ref  = it. first. first;
  	  const string& subj = it. first. second;
  	  
  	  if (verbose ())
  	  {
    	  cerr << ref << '\t' << subj << '\t';
    	  bestExon->saveText (cerr);
    	  cerr << endl;
    	}
    	
    	const size_t dash = ref. rfind ('-');
    	QC_ASSERT (dash != string::npos);
    	auto& exonIt = gene2exon [ref. substr (dash + 1)];
    	if (   ! exonIt 
    	    || exonIt->totalScore < bestExon->totalScore
    	   )
    	  exonIt = bestExon;
    }
  	graph. qc ();


    if (! gap_stats. empty ())
      gapF. reset (new OFStream (gap_stats));
          
    
    for (const auto& it : gene2exon)
    {
      const string s (it. second->getSeq (0));
      ASSERT (! s. empty ());
      if (verbose ())
      {
        cerr << it. first << '\t';
        it. second->saveText (cerr);
        cerr << '\t' << s << endl;
      }
      const Dna dna (it. first, s, false);
      dna. qc ();
      dna. saveText (cout);
    }
  }
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



