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
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "../version.inc"



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find best annotation and print proteins")
    {
      version = VERSION;
  	  addPositional ("dna", "DNA multi-FASTA file used fro TBLASTN");
  	  addPositional ("tblastn", "tblastn output in format: qseqid sseqid length positive qstart qend sstart send slen sseq");
  	  addPositional ("gencode", "Genetic code");
    }



	void body () const final
  {
    const string dnaFName     = getArg ("dna");
	  const string tblastnFName = getArg ("tblastn");
	  const Gencode gencode     = (Gencode) str2<int> (getArg ("gencode"));
  
  
	  map<string/*contig*/,Vector<Cds>> contig2cdss;
	  {
  	  LineInput in (tblastnFName);
  	  while (in. nextLine ())
  	  {
  	  	string qseqid;  // reference protein
  	  	string sseqid;  // contig accession
  	  	string sseq;
  	  	size_t length, positive, qstart, qend /* with '*' */;  // aa
  	  	size_t sstart, send, slen;  // na
  	  	{
    	  	QC_ASSERT (! in. line. empty ());
    	  	istringstream iss (in. line);
    		  iss >> qseqid >> sseqid >> length >> positive >> qstart >> qend >> sstart >> send >> slen >> sseq;
    		}
  	  	replaceStr (sseq, "-", "");
  	  	QC_ASSERT (! sseq. empty ());
  		  QC_ASSERT (length);
  		  QC_ASSERT (positive);
  		  QC_ASSERT (positive <= length);
  		  QC_ASSERT (qstart);
  		  QC_ASSERT (sstart);
  		  QC_ASSERT (send);
  		  QC_ASSERT (sstart != send);
  		  QC_ASSERT (sstart <= slen);
  		  QC_ASSERT (send   <= slen);		  
  		  qstart--;
  		  QC_ASSERT (qend >= qstart + positive);
  		  Strand strand = 1;
  		  if (sstart < send)
  		  	sstart--;
  		  else
  		  {
  		  	send--;
  		  	strand = -1;
  		  }
        {
    	  	const size_t sMatchLen = sstart < send ? send - sstart : sstart - send;
    	  	QC_ASSERT (sMatchLen >= 3 * positive);
    	    QC_ASSERT (divisible ((uint) sMatchLen, 3));	  	
    	  }  	  		  	
  	  	const double positivesFrac = (double) positive / (double) length;
  	  	if (positivesFrac < 0.85)  // PAR
  	  		continue;

  	  	PeptideOrf bestOrf;
  	    {
  		  	const Peptide pep ("pep", sseq, false);
  		  	const Vector<PeptideOrf> orfs (pep. getPeptideOrfs (sstart, strand, true, false, 20));  // PAR
  		  	size_t goodSize_max = 0;
  		  	for (const PeptideOrf& orf : orfs)
  		  	{
  		  	  orf. qc ();
      		  ASSERT (orf. stop <= length);
  		  	  if (   orf. good (20/*PAR*/) 
  		  	  	  && maximize (goodSize_max, orf. size ())
  		  	  	 )
  		  	  	bestOrf = orf;
  		  	}
  		  	// Non-PeptideOrf::good() PeptideOrf's: must be significantly better than the good() ones
  		  	size_t size_max = 0;
  		  	for (const PeptideOrf& orf : orfs)
  		  	  if (   orf. size () >= max ((size_t) ((double) goodSize_max * 1.2/*PAR*/), goodSize_max + 20/*PAR*/)
  		  	  	  && maximize (size_max, orf. size ())
  		  	  	 )
  		  	  	bestOrf = orf;
  		  }		  
  		  if (bestOrf. empty ())
  		    continue;
  	  		  		  
  		  contig2cdss [sseqid] << move (Cds (bestOrf. cdsStart (), bestOrf. cdsStop (), qseqid, positivesFrac));
  		  // add all ORFs of size >= 150 aa ??
  	  }
  	}
  	
  	
  	for (auto& it : contig2cdss)
  	{
      Vector<Cds>& cdss = it. second;
      DnaAnnot da;    
      {
        cdss. sort ();
      	const Cds* prev = nullptr;
  	    for (const Cds& cds : cdss)
  				if (! prev || ! cds. worse (*prev))
  				{
  			    da. cdss << cds;
  			    prev = & cds;
  			  }
	    }
	    cdss. clear ();
      const Cds* cds = da. run ();
      while (cds)
      {
        cds->qc ();
        cdss << *cds;
      	cds = cds->bestPrev;
      }
  	}
  	
  	
	  Multifasta faIn (dnaFName, false);
	  while (faIn. next ())
	  {
	    Dna dna (faIn, 1e6/*PAR*/, false);
    	strLower (dna. seq);
		  dna. qc ();
		  if (const Vector<Cds>* cdss = findPtr (contig2cdss, dna. getId ()))
		    for (const Cds& cds : *cdss)
		    {
		      QC_ASSERT (cds. right () <= dna. seq. size ());
		      string seq (dna. seq. substr (cds. left (), cds. size ()));
		      if (! cds. strand ())
		        reverseDna (seq);
		      const Dna orf (dna. getId () + ":" + to_string (cds. start_human ()) + ".." + to_string (cds. stop_human ()), seq, false);  // ??
		      size_t translationStart = 0;
	        const Peptide pep (orf. makePeptide (1, gencode, true, true, translationStart));
	        if (! translationStart)
	          pep. saveText (cout);
		    }
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



