// tblastn2frameshift.cpp

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
*   Find frame shuifts
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
  
  
void process (DiGraph &graph,
              VectorPtr<Exon> &exons)
{
	graph. qc ();
  if (exons. empty ())
    return;
	
	const Exon* initExon = Exon::exons2bestInitial (exons);
	ASSERT (initExon);
  const Exon* exon = initExon;
  for (;;)
  {
    ASSERT (exon);
    if (verbose ())
    {
      exon->saveText (cout);
      cout << endl;
    }
    const size_t end = exon->bestIntron ? exon->bestIntron->prev_end : exon->qseq. size ();
    if (! exon->bestIntron)
      break;
    const Exon* next = static_cast <const Exon*> (exon->bestIntron->node [true]);
    ASSERT (next);
    if ((exon->sstart % 3) != (next->sstart % 3))
      cout         << exon->qseqid
           << '\t' << exon->sseqid
           << '\t' << "fs_" + to_string (exon->pos2q (end)) 
                      + "_" + to_string (exon->pos2s (end)) 
                      + "_" + to_string (exon->strand == 1 ? 1 : 0)
         << '\n';
    exon = next;
  }

  graph. clear ();
  exons. clear ();
}
  
  
}



// ThisApplication

struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Find frame shifts using tblastn output.\nOutput: qseqid sseqid fs_<qpos>_<spos>_<sstrand(1/0)>, where <qpos> and <spos> are 0-based")
    {
      version = VERSION;
      // Input
      addPositional ("tblastn", "tblastn output in the format: qseqid sseqid qstart qend sstart send qseq sseq. Ordered by qseqid, sseqid"); 
  	  addKey ("matrix", "Protein matrix", "BLOSUM62");
    }



  void body () const final
  {
    const string blastFName = getArg ("tblastn");
	  const string matrix     = getArg ("matrix");


    const SubstMat sm (execDir + "/matrix/" + matrix);  
    sm. qc ();
  
    DiGraph graph;  // of Exon*
    VectorPtr<Exon> exons;
	  {
  	  LineInput in (blastFName, 10000);  // PAR
  	  string qseqid_prev;
  	  string sseqid_prev;
  	  while (in. nextLine ())
  	    try 
  	    { 
  	      auto exon = new Exon (graph, sm, in. line); 
  	      if (! (   qseqid_prev == exon->qseqid
  	             && sseqid_prev == exon->sseqid
  	            )
  	         )
  	      {
  	        process (graph, exons);
            if (qseqid_prev > exon->qseqid)
              throw runtime_error ("Non-sorted qseqid");
            if (   qseqid_prev == exon->qseqid 
                && sseqid_prev >  exon->sseqid
               )
              throw runtime_error ("Non-sorted sseqid");
  	      }
  	      exons << exon;
  	      qseqid_prev = exon->qseqid;
  	      sseqid_prev = exon->sseqid;
  	    }
	      catch (const exception &e)
	        { throw runtime_error (string (e. what ()) + "\n" + blastFName + ": " + in. lineStr ()); }
  	}
  	process (graph, exons);
  }
};




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



