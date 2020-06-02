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
*   Save full-length CDSs predicted by GeneMark
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;




namespace 
{
	
	
struct Exon : Root
{
  size_t start;
  size_t stop;
  
  Exon (size_t start_arg,
        size_t stop_arg)
    : start (start_arg)
    , stop (stop_arg)
    {
    	ASSERT (start < stop);
    }
  void saveText (ostream &os) const override
    { os << start << ".." << stop; }

  size_t getLen () const
    { return stop - start; }
	string getSeq (const string &contigDna) const
	  { return contigDna. substr (start, getLen ()); }
};

	
	
struct EuCds
{
	uint gene;
	string contig;
	bool strand;
	Vector<Exon> exons;
	bool has_start {false};
	bool has_stop {false};
	
	EuCds (uint gene_arg,
	       const string &contig_arg,
	       const string &strand_arg)
	  : gene (gene_arg)
	  , contig (contig_arg)
	  , strand (strand_arg == "+")
	  { ASSERT (gene);
	  	ASSERT (! contig. empty ());
	  	ASSERT (! contains (contig, ' '));
	  //ASSERT (! contains (contig, '.'));
	  	ASSERT (strand_arg == "+" || strand_arg == "-");
	  }
	     
	void print () const
	  { cout         << gene 
	  	     << '\t' << contig
	  	     << '\t' << strand
	  	     << '\t' << exons
	  	     << endl;
	  }
	char strand2char () const
	  { return strand ? '+' : '-'; }
	Dna getDna (const Dna &contigDna) const
	  { string s;
	  	for (const Exon& exon : exons)
	  		s += exon. getSeq (contigDna. seq);
	  	Dna dna (contig + "." + toString (gene), s, false);
	  	if (! strand)
	  		dna. reverse ();
	  	return dna;
	  }
};
  



struct ThisApplication : Application
{
  ThisApplication ()
  : Application ("Save full-length CDSs predicted by GeneMark")
  { 
	  addPositional ("fasta", "FASTA DNA file");
	  addPositional ("annot", "GeneMark file with annotations");
	  addFlag ("gtf", "<annot> is a GeneMark .gtf-file, otherwise a GeneMark .lst-file");
	  addKey ("cds", "Output FASTA DNA file of CDSs");
	  addKey ("prot", "Output FASTA DNA file of proteins");
	  addKey ("min_prot_len", "Min. protein length", "20");
	  addKey ("gencode", "NCBI genetic code", "0");
	  addFlag ("noerror", "Do not abort on errors");
	}



	void body () const final
	{
    const string fastaFName = getArg ("fasta");  
    const string annotFName = getArg ("annot"); 
    const bool gtf          = getFlag ("gtf"); 
    const string cdsFName   = getArg ("cds"); 
    const string protFName  = getArg ("prot"); 
    const size_t min_prot_len = str2<size_t> (getArg ("min_prot_len"));
    const Gencode gencode   = (Gencode) str2<int> (getArg ("gencode"));
    const bool noerror      = getFlag ("noerror");
    
    
    if (! gencode)
    	throw runtime_error ("Bad genocde");
        
    
    Vector<EuCds> cdss;  cdss. reserve (10000);  // PAR
    if (gtf)
    {
    	uint gene_prev = 0;
    	LineInput f (annotFName, 1000 * 1024, 1000);  // PAR
 	    Istringstream iss;
    	while (f. nextLine ())
    	  try
	    	{
	  	    iss. reset (f. line);
	        string contig, method, type, dot, strand, frame, gene_id_txt, gene_idS;
	        size_t start, stop;
	  	    iss >> contig >> method >> type >> start >> stop >> dot >> strand >> frame >> gene_id_txt >> gene_idS;
	  	    ASSERT (method == "GeneMark.hmm");  // PAR
	  	    ASSERT (start <= stop);
	  	    ASSERT (start >= 1);
	  	    ASSERT (gene_id_txt == "gene_id");
	  	    if (type == "CDS")
	  	    	continue;	
	
	  	    start--;
	  	    EXEC_ASSERT (trimPrefix (gene_idS, "\""));
	  	    const uint gene = str2<uint> (findSplit (gene_idS, '_'));
	  	  //ASSERT (isLeft (descr, "g\"; transcript_id \""));
	  	    ASSERT (gene);
	  	    ASSERT (gene == gene_prev || gene == gene_prev + 1);
	  	    if (gene != gene_prev)
	  	    	cdss << EuCds (gene, contig, strand);
	  	    EuCds* cds = & cdss. back ();
	  	    ASSERT (cds->contig == contig);
	  	    ASSERT (strand == string (1, cds->strand2char ()));
	  	    if (type == "exon")
	  	    {
	  	    	ASSERT (dot == "0");
	  	    	ASSERT (frame == ".");
	  	    	cds->exons << Exon (start, stop);
	  	    }
	  	    else if (   type == "start_codon"
	  	    	       || type == "stop_codon"
	  	    	      )
	  	    {
	  	    	ASSERT (dot == ".");
	  	    	ASSERT (frame == "0");
	  	    	ASSERT (stop - start == 3);
	  	    	if (type == "start_codon")
	  	    	  cds->has_start = true;
	  	    	else
	  	    	  cds->has_stop = true;
	  	    }
	  	    else
	  	    	throw runtime_error ("Unknown type " + type);
	
	  	    gene_prev = gene;
	      }
	      catch (...)
	      {
	      	cout << "Line # " << f. lineNum << endl;
	      	throw;
	      }
    }
    else
    {
    	LineInput f (annotFName, 1000 * 1024, 1000);  // PAR
    	while (f. nextLine () && ! f. line. empty ())
    	  ;
    	uint n = 0;
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
		  	    ASSERT (n == gene);
		  	    EuCds cds (gene, contig, strand);
		  	    cds. has_start = ! trimPrefix (leftEndS,  "<");
		  	    cds. has_stop  = ! trimPrefix (rightEndS, ">");
		  	    const Exon exon (str2<size_t> (leftEndS) - 1, str2<size_t> (rightEndS));
		  	    cds. exons << exon;
		  	    ASSERT (exon. getLen () == len);
		  	    cdss << cds;
	    		}
	    	}
	    	catch (...)
	    	{
	      	cout << "Line # " << f. lineNum << endl;
	    		throw;
	    	}
    }
    cout << "# CDSs: " << cdss. size () << endl;
    
    
    map <string/*contig*/, VectorPtr<EuCds>> contig2cdss;
    size_t good = 0;
    for (const EuCds& cds : cdss)
    	if (   cds. has_start
    		  && cds. has_stop
    		  && ! cds. exons. empty ()
    		 )
    	{
    		good++;
    		contig2cdss [cds. contig] << & cds;
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
    		
		  Multifasta f (fastaFName, false);
		  f. prog. active = false;
		  while (f. next ())
		  {
		    const Dna contigDna (f, 10000, false);  // PAR
		    contigDna. qc ();
        const VectorPtr<EuCds>& contigCdss = contig2cdss [contigDna. name];
        for (const EuCds* cds : contigCdss)
        {
        	const Dna cdsDna (cds->getDna (contigDna));
        	ASSERT (! cdsDna. seq. empty ());
        	cdsDna. qc ();
        	try 
        	{
	        	const Peptide pep (cdsDna. cds2prot (gencode, true));
	        	pep. qc ();
	        	if (pep. seq. size () < min_prot_len)
	        		continue;
	        	if (pep. getXs ())
	        		continue;
	        	if (cdsF)
	        	  cdsDna. saveText (*cdsF);
	        	if (protF)
	        		pep. saveText (*protF);
	        	longs++;
	        }
	        catch (...)
	        {
	        	cds->print ();
	        	cdsDna. print (cout);
	        	if (! noerror)
	        	  throw;
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



