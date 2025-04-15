// tblastn2marker_euk.cpp

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
  
  
//unique_ptr<OFStream> gapF;
char delim {'\0'};

// PAR
constexpr double complexity_min = 3.0;  
static_assert (complexity_min >= Peptide::stdMinComplexity);

map <Pair<string>/*ref,contig*/, VectorPtr<Exon>> contig2exons;



void process (DiGraph &graph,
              const string &gene,
              double threshold)
{
  graph. qc ();
  if (graph. empty ())
    return;

  // contig2exons
  for (const DiGraph::Node* n : graph. nodes)
  {
    const Exon* exon = static_cast <const Exon*> (n);
    ASSERT (exon);
    if (contains (exon->sseq, '*'))
      continue;  
    const Peptide pep ("x", exon->sseq, false);
    if (pep. getComplexity () < complexity_min)  
      continue;
    const Pair<string> p (exon->qseqid, exon->sseqid);
    contig2exons [p] << exon;
  }
    
  ASSERT (! gene. empty ());
	const string prefix (gene + ifS (delim, string (1, delim)));
	
	const Exon* bestExon_gene = nullptr;
	for (const auto& it : contig2exons)
	{
	  ASSERT (! it. second. empty ());
	  const string& ref  = it. first. first;
	  const string& subj = it. first. second;
	  ASSERT (isLeft (ref, prefix));
	  
  	const Exon* bestExon = exons2bestInitial (it. second);
	  if (! bestExon)
	    continue;
	  
	#if 0  
	  if (   ref  == "VUSY-4471"
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
  	
  	if (   ! bestExon_gene 
  	    || bestExon_gene->totalScore < bestExon->totalScore
  	   )
  	  bestExon_gene = bestExon;
  }
	graph. qc ();
  if (! bestExon_gene)
    return;


  const string s (bestExon_gene->getSeq (0));
  ASSERT (! s. empty ());
  if (verbose ())
  {
    cerr << gene << '\t';
    bestExon_gene->saveText (cerr);
    cerr << '\t' << s << endl;
  }

  const AlignScore totalScore = bestExon_gene->totalScore;
  if (totalScore / (double) s. size () >= threshold)
  {
    string ref        (" ref="        + bestExon_gene->qseqid + ":");
    string contig     (" contig="     + bestExon_gene->sseqid + ":");
    string ref_hsp    (" ref_hsp="    + bestExon_gene->qseqid + ":");
    string contig_hsp (" contig_hsp=" + bestExon_gene->sseqid + ":");
    size_t start = 0;
    const Exon* exon = bestExon_gene;
    for (;;)
    {
      ASSERT (exon);
      const size_t end = exon->bestIntron ? exon->bestIntron->prev_end : exon->sseq. size ();
      ref        += to_string (exon->pos2q (start) + 1) + "-" + to_string (exon->pos2q (end));  
      contig     += to_string (exon->pos2s (start) + 1) + "-" + to_string (exon->pos2s (end));  
      ref_hsp    += to_string (exon->qstart        + 1) + "-" + to_string (exon->qend);  
      contig_hsp += to_string (exon->sstart        + 1) + "-" + to_string (exon->send);
      if (! exon->bestIntron)
        break;
      start = exon->bestIntron->next_start;
      exon = static_cast <const Exon*> (exon->bestIntron->node [true]);
      ref        += ",";
      contig     += ",";
      ref_hsp    += ",";
      contig_hsp += ",";
    }
    const string name (gene + ref + contig + ref_hsp + contig_hsp + " strand=" + to_string (exon->strand == 1 ? 1 : 0) + " score=" + to_string (int (totalScore)));
    Peptide pep (name, s, false);
    pep. pseudo = true;
    pep. qc ();
    if (pep. getComplexity () >= complexity_min)  
      pep. saveText (cout);
  }

  
  graph. clear ();
  contig2exons. clear ();
}



struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Find best protein matches and print proteins")
    {
      version = VERSION;
  	  addPositional ("tblastn", "tblastn output in format: qseqid sseqid qstart qend sstart send qseq sseq. Ordered by qseqid. If -delimeter then qseqid = <gene><delimiter><variant>, where <gene> has no <delimiter>.");
  	  addKey ("delimiter", "Delimiter chracater separating gene names from variant names in <tblastn>");
  	  addKey ("threshold", "Min. protein score to length ratio to be reported", "4");
  	  addKey ("matrix", "Protein matrix", "BLOSUM62");
  	//addKey ("gap_stats", "Output file with gap lengths");
    }



	void body () const final
  {
	  const string tblastnFName = getArg ("tblastn");
	  const string delimiterS   = getArg ("delimiter"); 
	  const double threshold    = arg2double ("threshold");
	  const string matrix       = getArg ("matrix");
	//const string gap_stats    = getArg ("gap_stats");
  
  
    if (! delimiterS. empty ())
    {
      QC_ASSERT (delimiterS. size () == 1);
      delim = delimiterS [0];
    }
  
    const SubstMat sm (execDir + "/matrix/" + matrix);  
    sm. qc ();
  
  #if 0
	  if (! gap_stats. empty ())
      gapF. reset (new OFStream (gap_stats));
  #endif
  
    DiGraph graph;  // of Exon*
	  {
  	  LineInput in (tblastnFName, 10000);  // PAR
  	  string gene_prev;
  	  while (in. nextLine ())
  	  {
      	const size_t pos = in. line. find (delim ? delim : '\t');
      	QC_ASSERT (pos != string::npos);
      	const string gene (in. line. substr (0, pos));
      	if (gene. empty ())
  	      throw runtime_error ("No gene\n" + tblastnFName + ": " + in. lineStr ()); 
      	if (gene < gene_prev)
  	      throw runtime_error ("Alphabetically disordered gene\n" + tblastnFName + ": " + in. lineStr ()); 
      	if (gene_prev != gene)
      	{
      	  process (graph, gene_prev, threshold);
      	  gene_prev = gene;
      	}
  	    try { new Exon (graph, sm, in. line); }
  	      catch (const exception &e)
  	        { throw runtime_error (string (e. what ()) + "\n" + tblastnFName + ": " + in. lineStr ()); }
  	  }
   	  process (graph, gene_prev, threshold);
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



