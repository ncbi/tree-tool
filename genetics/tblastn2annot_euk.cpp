// tblastn2annot_euk.cpp

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
*   Find best gene matches and print proteins
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
  bool nextExonFixed {false};
  const bool prevExonFixed {false};
	// Output
	bool bestIntronSet {false};
	const Intron* bestIntron {nullptr};
	AlignScore totalScore {- score_inf};


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
    , prevExonFixed (true)
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
    // Output: score, bestIntron
  void trimHangingDashes ();
public:
  void saveText (ostream &os) const final;
  void qc () const final;
    
    
  size_t qLen () const
    { return qend - qstart; }
  size_t sLen () const
    { return (send - sstart) / 3; }
    // aa
  // --> center of match density ??
  size_t qCenter () const  
    { return (qstart + qend) / 2; }
  size_t sCenter () const  
    { return (sstart + send) / 2; }
  //
  bool arcable (const Exon &next) const
    {
      ASSERT (qseqid == next. qseqid);
      ASSERT (sseqid == next. sseqid);
      ASSERT (! nextExonFixed);
      ASSERT (! next. prevExonFixed);
      
      // PAR
      constexpr size_t intron_max = 30000;  // nt 
      
      if (strand != next. strand)
        return false;
    //if (same frame and overlap) return false;  // ??
    
      // => DAG
      if (qCenter () >= next. qCenter ())
        return false;  
      if (strand == 1)
      {
        if (sCenter () >= next. sCenter ())
          return false;  
      }
      else
        if (next. sCenter () >= sCenter ())
          return false;  

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
    // Update: bestIntronSet, totalScore
    // Output: bestIntron
  string getSeq (size_t start) const;
  size_t pos2q (size_t pos) const
    {
      ASSERT (pos <= qseq. size ());
      size_t j = qstart;
      FFOR (size_t, i, qseq. size () + 1)
      {
        ASSERT (j >= qstart);
        ASSERT (j <= qend);
        if (i == pos)
          return j;
        if (qseq [i] != '-')
          j++;
      }
      ERROR;
    }
  size_t pos2s (size_t pos) const
    {
      ASSERT (pos <= sseq. size ());
      size_t j = (strand == 1 ? sstart : send);
      FFOR (size_t, i, sseq. size () + 1)
      {
        ASSERT (j >= sstart);
        ASSERT (j <= send);
        if (i == pos)
          return j;
        if (sseq [i] != '-')
        {
          if (strand == 1)
            j += 3;
          else
          {
            ASSERT (j >= 3);
            j -= 3;
          }
        }
      }
      ERROR;
    }
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
    , prev_end (prev->sseq. size ())
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


  void qc () const final
    {
      if (! qc_on)
        return;
      DiGraph::Arc::qc ();

      const Exon* prev = static_cast <const Exon*> (node [false]);
      const Exon* next = static_cast <const Exon*> (node [true]);
      ASSERT (prev);
      ASSERT (next);
      QC_ASSERT (prev->qseqid == next->qseqid);
      QC_ASSERT (prev->sseqid == next->sseqid);
      QC_ASSERT (prev->strand == next->strand);

      QC_ASSERT (prev_end   <= prev->qseq. size ());
      QC_ASSERT (next_start <= next->qseq. size ());
    }
    
    
  void saveText (ostream &os) const final
    {
      os << "  prev_end: " << prev_end 
         << "  next_start: " << next_start
         << "  score: " << score;
      ASSERT (node [true]);

      Offset ofs;
      Offset::newLn (os);
      os << "Exon:";
      static_cast <const Exon*> (node [true]) -> saveText (os);
    }
    
    
  AlignScore getTotalScore ()
    {
      const Exon* next = static_cast <const Exon*> (node [true]);
      ASSERT (next);
      var_cast (next) -> setBestIntron ();  // DAG => no loop 
      return next->totalScore - score;
    }
};



void Exon::finish ()
{
  ASSERT (graph);
  ASSERT (score == 0);
  ASSERT (! bestIntron);
  ASSERT (! nextExonFixed);
  
  trimHangingDashes ();

  size_t gap = 0;
  bool pseudo = false;
  size_t qstart_new = qstart;
  size_t sstart_new = sstart;
  size_t send_new   = send;
  FFOR (size_t, i, qseq. size ())
  {
    if (qseq [i] == '-')
    {
      gap++;
      QC_ASSERT (sseq [i] != '-');
      if (sseq [i] == '*')
        pseudo = true;
    }
    else
    {
      QC_ASSERT (qseq [i] != '*');
      if (pseudo || gap >= 15)  // PAR  
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
        nextExonFixed = true;
        trimHangingDashes ();
        break;
      }
      qstart_new++;
      gap = 0;
      pseudo = false;
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
//QC_ASSERT (score > 0);

  const Pair<string> p (qseqid, sseqid);
  contig2exons [p] << this;
  
  qc ();
}



void Exon::trimHangingDashes ()
{
  while (   ! qseq. empty ()
         && (   qseq. front () == '-'
             || sseq. front () == '-'
            )
        )
  {
    qseq. erase (0, 1);
    sseq. erase (0, 1);
    qstart++;
    if (strand == 1)
      sstart += 3;
    else
    {
      ASSERT (send >= 3);
      send -= 3;
    }
  }
  
  while (   ! qseq. empty ()
         && (   qseq. back () == '-'
             || sseq. back () == '-'
            )
        )
  {
    qseq. erase (qseq. size () - 1);
    sseq. erase (sseq. size () - 1);
    qend--;
    if (strand == 1)
    {
      ASSERT (send >= 3);
      send -= 3;
    }
    else
      sstart += 3;
  }
}



void Exon::qc () const 
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
        
  QC_ASSERT (pos2q (0)             == qstart);
  QC_ASSERT (pos2q (qseq. size ()) == qend);
  QC_ASSERT (pos2s (0)             == (strand == 1 ? sstart : send));
  QC_ASSERT (pos2s (sseq. size ()) == (strand == 1 ? send   : sstart));
    
  QC_IMPLY (nextExonFixed, bestIntron);
  QC_IMPLY (bestIntron, arcs [true]. find (var_cast (bestIntron)) != no_index);

//QC_ASSERT (score >= 0);
//QC_ASSERT (totalScore >= 0);
}



void Exon::saveText (ostream &os) const 
{
  os        << qstart + 1 << '(' << sstart + 1 << ')'
     << "-" << qend   + 1 << '(' << send   + 1 << ')'
     << "  strand: " << (int) strand
     << "  score: " << score 
     << "  totalScore: " << totalScore;
  if (bestIntron)
  {
    Offset ofs;
    Offset::newLn (os);
    os << "Intron:";
    bestIntron->saveText (os);
  }
}



void Exon::setBestIntron ()
{
  if (bestIntronSet)
    return;    
  bestIntronSet = true;

  ASSERT (nextExonFixed == (bool) bestIntron);

  if (arcs [true]. empty ())
  {
    ASSERT (! bestIntron);
    totalScore = score;
  }
  else
  {
    if (bestIntron)
      totalScore = score + var_cast (bestIntron) -> getTotalScore ();
    else
    {
      totalScore = - score_inf;
      for (const DiGraph::Arc* arc : arcs [true])
      {
        ASSERT (arc);
        ASSERT (arc->node [false] == this);
        const Intron* intron = static_cast <const Intron*> (arc);
        if (maximize (totalScore, score + var_cast (intron) -> getTotalScore ()))
          bestIntron = intron;
      }
    }
    ASSERT (bestIntron);
  }
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
    : Application ("Find best gene matches and print proteins")
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
  	  LineInput in (tblastnFName, 1000);  // PAR
  	  while (in. nextLine ())
  	    try { new Exon (graph, in. line); }
  	      catch (const exception &e)
  	        { throw runtime_error (string (e. what ()) + "\n" + in. lineStr ()); }
  	}
  	graph. qc ();
  	
  	
  	map<string/*gene*/,const Exon*> gene2exon;
  	{
  	  Progress prog (contig2exons. size ());
    	for (const auto& it : contig2exons)
    	{
    	  prog ();
    	  ASSERT (! it. second. empty ());
    	  
    	  // new Intron
    	  for (const Exon* next : it. second)
    	  {
    	    ASSERT (next);
    	    if (next->prevExonFixed)
    	      continue;
      	  for (const Exon* prev : it. second)
      	    if (prev == next)
      	      break;
      	    else
      	      if (   ! prev->nextExonFixed 
      	          && prev->arcable (*next)
      	         )
      	        new Intron (var_cast (prev), var_cast (next));
      	}
      	
      	AlignScore totalScore_max = - score_inf;
      	const Exon* bestExon = nullptr;
    	  for (const Exon* exon : it. second)
    	  {
    	    var_cast (exon) -> setBestIntron ();
    	    if (exon->totalScore <= 0)  // PAR ??
    	      continue;
    	    if (maximize (totalScore_max, exon->totalScore))
    	      bestExon = exon;
    	  }
    	  if (! bestExon)
    	    continue;
    	  
    	  const string& ref  = it. first. first;
    	  const string& subj = it. first. second;
    	  
    	#if 0
    	  if (   ref  == "IDAU-4471"
    	      && subj == "1890392768"
    	     )
      	  for (const Exon* exon : it. second)
      	  {
      	    exon->saveText (cout);
    	      cout << endl;
    	    }
    	#endif
    	  
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
    }
  	graph. qc ();


    if (! gap_stats. empty ())
      gapF. reset (new OFStream (gap_stats));
          
    
  	{
  	  Progress prog (gene2exon. size ());
      for (const auto& it : gene2exon)
      {
        prog ();
        const Exon* exon = it. second;
        ASSERT (exon);
        const string s (exon->getSeq (0));
        ASSERT (! s. empty ());
        if (verbose ())
        {
          cerr << it. first << '\t';
          exon->saveText (cerr);
          cerr << '\t' << s << endl;
        }

        string ref        (" ref="        + exon->qseqid + ":");
        string contig     (" contig="     + exon->sseqid + ":");
        string ref_hsp    (" ref_hsp="    + exon->qseqid + ":");
        string contig_hsp (" contig_hsp=" + exon->sseqid + ":");
        const AlignScore totalScore = exon->totalScore;
        size_t start = 0;
        for (;;)
        {
          ASSERT (exon);
          const size_t end = exon->bestIntron ? exon->bestIntron->prev_end : exon->sseq. size ();
          ref        += to_string (exon->pos2q (start) + 1) + "-" + to_string (exon->pos2q (end));  
          contig     += to_string (exon->pos2s (start) + 1) + "-" + to_string (exon->pos2s (end));  
          ref_hsp    += to_string (exon->qstart + 1) + "-" + to_string (exon->qend);  
          contig_hsp += to_string (exon->sstart + 1) + "-" + to_string (exon->send);
          if (! exon->bestIntron)
            break;
          start = exon->bestIntron->next_start;
          exon = static_cast <const Exon*> (exon->bestIntron->node [true]);
          ref        += ",";
          contig     += ",";
          ref_hsp    += ",";
          contig_hsp += ",";
        }
        const string name (it. first + ref + contig + ref_hsp + contig_hsp + " strand=" + to_string (exon->strand == 1 ? 1 : 0) + " score=" + to_string (int (totalScore)));
        Peptide pep (name, s, false);
        pep. pseudo = true;
        pep. qc ();
        pep. saveText (cout);
      }
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



