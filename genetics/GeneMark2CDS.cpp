// GeneMark2CDS.cpp

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
*   Save CDSs predicted by GeneMark
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "../version.inc"

#include "../common.inc"




namespace 
{
	
	
struct Exon : Root
{
  size_t start {0};
  size_t stop {0};
  size_t codonStart {0};
    // Incomplete first codon length
  
  
  Exon (size_t start_arg,
        size_t stop_arg,
        size_t codonStart_arg)
    : start (start_arg)
    , stop (stop_arg)
    , codonStart (codonStart_arg)
    {
    	QC_ASSERT (start < stop);
    }
  void qc () const override
    {
      if (! qc_on)
        return;
      QC_ASSERT (start < stop);
      QC_ASSERT (codonStart < 3);
      QC_ASSERT (getLen () > codonStart);
    }
  void saveText (ostream &os) const override
    { os << start << ".." << stop << '(' << codonStart << ')'; }


  size_t getLen () const
    { return stop - start; }
	string getSeq (const string &contigDna) const
	  { return contigDna. substr (start, getLen ()); }
};

	
	
struct EuCds : Root
// Eukaryotic CDS
{
	const uint gene;
	const string contig;
	  // Name
	const bool strand;
	Vector<Exon> exons;
	size_t start {no_index};
	size_t stop {no_index};
private:
	bool cut {false};
public:
	
	
	EuCds (uint gene_arg,
	       const string &contig_arg,
	       const string &strand_arg)
	  : gene (gene_arg)
	  , contig (contig_arg)
	  , strand (strand_arg == "+")
	  { 
	  	QC_ASSERT (strand_arg == "+" || strand_arg == "-");
	  }	     
	void qc () const override
	  { 
	    if (! qc_on)
	      return;
	    QC_ASSERT (gene);
	  	QC_ASSERT (! contig. empty ());
	  	QC_ASSERT (! contains (contig, ' '));
	  	size_t prev_stop = 0;
	  	for (const Exon& exon : exons)
	  	{
	  	  exon. qc ();
	  	  QC_ASSERT (exon. start >= prev_stop);
	  	  prev_stop = exon. stop;
	  	}
	  	QC_IMPLY (start != no_index && stop != no_index, start < stop);
	  	QC_ASSERT (cut);
  	  const Exon* prev = nullptr;
  	  if (strand)
  	  	for (const Exon& exon : exons)
  	  	{
  	  	  QC_IMPLY (! prev, ! exon. codonStart);
  	  	  QC_IMPLY (prev, (prev->getLen () - prev->codonStart + exon. codonStart) % 3 == 0);
  	  	  prev = & exon;
  	  	}
  	  else
  	  	FOR_REV (size_t, i, exons. size ())
  	  	{
  	  	  const Exon& exon = exons [i];
  	  	  QC_IMPLY (! prev, ! exon. codonStart);
  	  	  QC_IMPLY (prev, (prev->getLen () - prev->codonStart + exon. codonStart) % 3 == 0);
  	  	  prev = & exon;
  	  	}
	  }
	void saveText (ostream &os) const override
	  { os         << gene 
	  	   << '\t' << contig
	  	   << '\t' << strand
	  	   << '\t' << start
	  	   << '\t' << stop
	  	   << '\t';
	  	save (os, exons, '\t');
	  	os << endl;
	  }
	 
	  
	char strand2char () const
	  { return strand ? '+' : '-'; }
	bool trunc5 () const
	  { return (strand ? start : stop) == no_index; }
	bool trunc3 () const
	  { return (strand ? stop : start) == no_index; }
	void cutExons ()
	  { ASSERT (! cut);
	    for (Iter<Vector<Exon>> it (exons); it. next (); )
	    {
	      if (start != no_index)
	      {
	        if (start >= it->stop)
	        {
	          QC_ASSERT (! it->codonStart);
	          it. erase ();
	          continue;
	        }
	        else if (start >= it->start)
	        {
	          QC_IMPLY (strand, ! it->codonStart);
	          it->start = start;
	        }
	      }
	      if (stop != no_index)
	      {
	        if (stop <= it->start)
	        {
	          QC_ASSERT (! it->codonStart);
	          it. erase ();
	        }
	        else if (stop < it->stop)
	        {
	          QC_IMPLY (! strand, ! it->codonStart);
	          it->stop = stop;
	        }
	      }
	    }
	    
	    if (! exons. empty ())
	    {
	      if (start == no_index && strand)
	      {
	        Exon& exon = exons. front ();
	        exon. start += exon. codonStart;
	        exon. codonStart = 0;
	      }
	      if (stop == no_index && ! strand)
	      {
	        Exon& exon = exons. back ();
	        exon. stop -= exon. codonStart;
	        exon. codonStart = 0;
	      }
	    }
	    cut = true;
	  }
	Dna getDna (const Dna &contigDna) const
	  { ASSERT (cut);
	    string s;
	  	for (const Exon& exon : exons)
	  		s += exon. getSeq (contigDna. seq);
	  	Dna dna (contig + "." + to_string (gene) + " trunc5:" + to_string ((int) trunc5 ()) + " trunc3:" + to_string ((int) trunc3 ()), s, false);
	  	if (! strand)
	  		dna. reverse ();
	  	return dna;
	  }
};
  



struct ThisApplication final : Application
{
  ThisApplication ()
  : Application ("Save CDSs predicted by GeneMark; FASTA description contains: trunc5:(0/1) trunc3:(0/1) ")
  { 
    version = VERSION;
	  addPositional ("fasta", "FASTA DNA file");
	  addPositional ("annot", "GeneMark file with annotations");
	  addFlag ("gtf", "<annot> is a GeneMark .gtf-file, otherwise a GeneMark .lst-file");
	  addKey ("cds", "Output FASTA DNA file of CDSs");
	  addKey ("prot", "Output FASTA DNA file of proteins");
	  addKey ("prot_len_min", "Min. protein length", "20");
	  addKey ("prot_len_max", "Min. protein length", "50000");
	  addKey ("gencode", "NCBI genetic code", "0");
	  addFlag ("complete", "Save only non-truncated");
	  addKey ("ambig", "Min. number of ambiguities to discard", "1");
	//addFlag ("noerror", "Do not abort on errors");
	}



	void body () const final
	{
    const string fastaFName   = getArg ("fasta");  
    const string annotFName   = getArg ("annot"); 
    const bool gtf            = getFlag ("gtf"); 
    const string cdsFName     = getArg ("cds"); 
    const string protFName    = getArg ("prot"); 
    const size_t prot_len_min = str2<size_t> (getArg ("prot_len_min"));
    const size_t prot_len_max = str2<size_t> (getArg ("prot_len_max"));
    const Gencode gencode     = (Gencode) str2<int> (getArg ("gencode"));
    const bool complete       = getFlag ("complete");
    const size_t ambig        = str2<size_t> (getArg ("ambig"));
  //const bool noerror        = getFlag ("noerror");
    
    
    if (! gencode)
    	throw runtime_error ("Bad genocde");
        
    
    VectorOwn<EuCds> cdss;  cdss. reserve (10000);  // PAR
    if (gtf)
    {
    	uint gene_prev = 0;
    	LineInput f (annotFName, 1000);  // PAR
 	    EuCds* cds = nullptr;
    	while (f. nextLine ())
    	  try
	    	{
	    	        string contig      = findSplit (f. line, '\t');
	    	  const string method      = findSplit (f. line, '\t');
	    	  const string type        = findSplit (f. line, '\t');
	    	        size_t start       = str2<size_t> (findSplit (f. line, '\t'));
	    	  const size_t stop        = str2<size_t> (findSplit (f. line, '\t'));
	    	  const string dot         = findSplit (f. line, '\t');
	    	  const string strand      = findSplit (f. line, '\t');
	    	  const string frame       = findSplit (f. line, '\t');
	    	  const string gene_id_txt = findSplit (f. line, ' ');
	    	        string gene_idS    = std::move (f. line);
	  	    
	  	    const size_t spacePos = contig. find (' ');
	  	    if (spacePos != string::npos)
	  	      contig. erase (spacePos);

	  	    QC_ASSERT (isLeft (method, "GeneMark.hmm"));  // PAR
	  	    QC_ASSERT (start <= stop);
	  	    QC_ASSERT (start >= 1);
	  	    QC_ASSERT (gene_id_txt == "gene_id");
	  	    if (type == "exon" /*"CDS"*/)
	  	    	continue;	
          if (type == "gene")
            continue;
          if (type == "mRNA")
            continue;
          if (type == "intron")
            continue;
	
	  	    start--;
	  	    EXEC_ASSERT (trimPrefix (gene_idS, "\""));
	  	    const uint gene = str2<uint> (findSplit (gene_idS, '_'));
	  	  //ASSERT (isLeft (descr, "g\"; transcript_id \""));
	  	    QC_ASSERT (gene);
	  	    QC_ASSERT (gene == gene_prev || gene == gene_prev + 1);
	  	    if (gene != gene_prev)
	  	    {
	  	      if (cds)
	  	      {
	  	        cds->cutExons ();
	  	        cds->qc ();
	  	      }
	  	      cds = new EuCds (gene, contig, strand);
	  	    	cdss << cds;
	  	    }
	  	    ASSERT (cds);
	  	    QC_ASSERT (cds->contig == contig);
	  	    QC_ASSERT (strand == string (1, cds->strand2char ()));
	  	    if (type == "CDS" /*"exon"*/)
	  	    {
	  	    	QC_ASSERT (dot == "." /*"0"*/);
	  	      QC_ASSERT (frame != ".");
	  	      const size_t codonStart = str2<size_t> (frame);
	  	      QC_ASSERT (codonStart < 3);
	  	    	cds->exons << Exon (start, stop, codonStart);
	  	    }
	  	    else if (   type == "start_codon"
	  	    	       || type == "stop_codon"
	  	    	      )
	  	    {
	  	    	QC_ASSERT (dot == ".");
	  	    	QC_ASSERT (frame == "0");
	  	    	QC_ASSERT (stop - start == 3);
	  	    	if ((type == "start_codon") == cds->strand)
	  	    	  cds->start = start;
	  	    	else
	  	    	  cds->stop = stop;
	  	    }
	  	    else
	  	    	throw runtime_error ("Unknown type " + strQuote (type));
	
	  	    gene_prev = gene;
	      }
	      catch (const exception &e)
	      {
	      	throw runtime_error (  f. lineStr () + "\n" 
	      	                     + (cds ? cds->contig + "." + to_string (cds->gene) + "\n" : "")
	      	                     + e. what ());
	      }
      if (cds)
      {
        cds->cutExons ();
        cds->qc ();
      }
    }
    else
    {
    	LineInput f (annotFName, 1000);  // PAR
    	while (f. nextLine () && ! f. line. empty ())
    	  ;
    	uint n = 0;
    	EuCds* cds = nullptr;
    	for (;;)
    	  try
	    	{
	    		if (! f. expectPrefix ("FASTA definition line: ", true))
	    			break;
	    		const string contig (f. line);
	    		f. expectPrefix ("Predicted genes", false);
	    		f. expectPrefix ("   Gene    Strand    LeftEnd    RightEnd       Gene     Class", false);
	    		f. expectPrefix ("    #                                         Length", false);
	  	    Istringstream iss ;
	    		while (f. nextLine () && ! f. line. empty ())
	    		{
		  	    iss. reset (f. line);
		  	    uint gene;
		        size_t len, geneClass /*??*/;
		        string strand, leftEndS, rightEndS;
		  	    iss >> gene >> strand >> leftEndS >> rightEndS >> len >> geneClass;
		  	    n++;
		  	    QC_ASSERT (n == gene);
		  	    cds = new EuCds (gene, contig, strand);
		  	    Exon exon (str2<size_t> (leftEndS) - 1, str2<size_t> (rightEndS), 0 /*??*/);
		  	    if (trimPrefix (leftEndS,  "<"))
		  	      cds->start = exon. start;
		  	    if (trimPrefix (rightEndS, ">"))
		  	      cds->stop  = exon. stop;
		  	    cds->exons << exon;
		  	    QC_ASSERT (exon. getLen () == len);
		  	    cds->cutExons ();
		  	    cds->qc ();
		  	    cdss << cds;
	    		}
	    	}
	      catch (const exception &e)
	      {
	      	throw runtime_error (  f. lineStr () + "\n" 
	      	                     + (cds ? cds->contig + "." + to_string (cds->gene) + "\n" : "")
	      	                     + e. what ());
	      }
    }
    cout << "# CDSs: " << cdss. size () << endl;
    
    
    map <string/*contig*/, VectorPtr<EuCds>> contig2cdss;
    size_t good = 0;
    for (const EuCds* cds : cdss)
    {
      if (cds->exons. empty ())
        continue;
    	if (   complete
    	    && (   cds->trunc5 ()
    		      || cds->trunc3 ()
    		     )
    		 )
    		continue;
  		good++;
  		contig2cdss [cds->contig] << cds;
  	}
    cout << "# Good CDSs: " << good << endl;


	  size_t longs = 0;
    {
    	unique_ptr<OFStream> cdsF;
    	if (! cdsFName. empty ())
    		cdsF. reset (new OFStream (cdsFName));
    		
    	unique_ptr<OFStream> protF;
    	if (! protFName. empty ())
    		protF. reset (new OFStream (protFName));
    		
		  Multifasta f (fastaFName, false, 100);  // PAR
		  while (f. next ())
		  {
		    const Dna contigDna (f, 10000, false);  // PAR
		    contigDna. qc ();
        const VectorPtr<EuCds>& contigCdss = contig2cdss [contigDna. getId ()];
        for (const EuCds* cds : contigCdss)
        {
          ASSERT (cds);
        	const Dna cdsDna (cds->getDna (contigDna));
        	ASSERT (! cdsDna. seq. empty ());
        	cdsDna. qc ();
        	try 
        	{
	        	Peptide pep (cdsDna. cds2prot (gencode, cds->trunc5 (), cds->trunc3 (), true, true /*'tar' etc.*/));
	        	pep. name += " " + cdsDna. getDescription (false);
	        	try { pep. qc (); }
    	        catch (const exception &e)
    	        {
    	          throw runtime_error (e. what () + string ("\n") + pep. str ());
    	        }	        	
	        	if (pep. seq. size () < prot_len_min)
	        		continue;
	        	if (pep. seq. size () > prot_len_max)
	        		continue;
	        	if (pep. getXs () >= ambig)
	        		continue;
	        	if (cdsF)
	        	  cdsDna. saveText (*cdsF);
	        	if (protF)
	        		pep. saveText (*protF);
	        	longs++;
	        }
	        catch (const exception &e)
	        {
 	          throw runtime_error (e. what () + string ("\n") + cds->str () + "\n" + cdsDna. str ());
	        }
	      }
      }
    }
    cout << "# Long, clean CDSs: " << longs << endl;
  }  
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



