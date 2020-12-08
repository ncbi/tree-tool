// fasta2hash.cpp

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
*   Print hash codes for full-length CDSs or proteins
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
  
  

const StringVector gene_finders {"prodigal", "GeneMark"};


  
void reducePep (Peptide &pep)
{
  static const string s_Blast_Alphabet10 ("ABCBBFGHIIKIIBOPKKAAUIWXFB");  // By I.Tolstoy
  for (char& c : pep. seq)
    c = s_Blast_Alphabet10 [(size_t) (c - 'A')];
}



hash<string> str_hash;

unique_ptr<OFStream> targetSeqF;



void addHash (const string &s,
              const string &seqName,
              Vector<size_t> &hashes,
              Vector<size_t> &targetHashes)
{
	ASSERT (! s. empty ());
	const size_t h = str_hash (s);
  hashes << h;
  if (   targetSeqF. get ()
      && targetHashes. containsFast (h)
     )
  	*targetSeqF << seqName << endl;
}
  

  
struct ThisApplication : Application
{
  ThisApplication ()
  : Application ("Print hash codes for full-length CDSs or proteins")
  { 
  	// Input
	  addPositional ("in", "Multi-FASTA file");
	  addKey ("gene_finder", "Gene finder used to find CDSs: " + gene_finders. toString (", "));
	  addFlag ("cds", "Are input sequences CDSs? If not then proteins");
	  addFlag ("translate", "Translate CDSs to proteins");
	  addFlag ("reduce", "Reduce 20 amino acids to 13 characters");
	  addKey ("kmer", "Cut the sequences into k-mers of this size. 0 - do not cut", "0");
	  addKey ("min_prot_len", "Min. protein length of input sequences", "0");
	  addKey ("max_hashes", "Max. number of hashes to print. 0 - print all", "0");
	  addKey ("target_hashes", "Target hashes");
	  // Output
	  addPositional ("out", "Output file of sorted unique hashes. Sorting is numeric.");
	  addKey ("out_prot", "File to save proteins with hashes");
	  addKey ("target_seq", "File to save the sequence names of the hashes in -target_hashes");
	}



	void body () const final
	{
    const string in            =               getArg ("in");  
    const string gene_finder   =               getArg ("gene_finder");
    const bool cds             =               getFlag ("cds");
    const bool translate       =               getFlag ("translate");
    const bool reduce          =               getFlag ("reduce");
    const size_t kmer          = str2<size_t> (getArg ("kmer"));
    const size_t prot_len_min  = str2<size_t> (getArg ("min_prot_len"));
    const size_t max_hashes    = str2<size_t> (getArg ("max_hashes"));
    const string target_hashes =               getArg ("target_hashes");
    const string out           =               getArg ("out"); 
    const string out_prot      =               getArg ("out_prot");
    const string target_seq    =               getArg ("target_seq");
    QC_ASSERT (! out. empty ()); 
    QC_IMPLY (translate, cds);
    QC_IMPLY (reduce, ! cds || translate);
    QC_IMPLY (! out_prot. empty (), ! cds || translate);
    QC_IMPLY (! target_seq. empty (), ! target_hashes. empty ());
  //if (! gene_finders. contains (gene_finder))
    //throw runtime_error ("Uknown gene_finder: " + gene_finder);
    

    Vector<size_t> targetHashes;
    if (! target_hashes. empty ())
    {
      LineInput f (target_hashes);  
  	  while (f. nextLine ())
    	  targetHashes << str2<size_t> (f. line);
    	cout << "# Target hashes: " << targetHashes. size () << endl;
    }
    targetHashes. sort ();
    QC_ASSERT (targetHashes. isUniq ());
    
    if (! target_seq. empty ())
      targetSeqF. reset (new OFStream (target_seq));
    

    const size_t seq_len_min = cds && ! translate 
                                 ? prot_len_min * 3 
                                 : prot_len_min;
        
    
    static_assert (sizeof (size_t) == 8, "Size of size_t must be 8 bytes");
    
    Vector<size_t> hashes;  hashes. reserve (10000);   // PAR
    size_t sequences = 0;
    unique_ptr<OFStream> protF (out_prot. empty () ? nullptr : new OFStream (out_prot));
    {
		  Multifasta f (in, ! cds);
		  f. prog. active = false;
		  while (f. next ())
		  {
		    string s;
		    string seqName;
		    if (cds)
		    {
  		    const Dna dna (f, 10000, false);  // PAR
  		    seqName = dna. name;
		    	try { dna. qc (); }
	  		    catch (...)
		  		  {
	    	    	if (gene_finder == "GeneMark")
		    	      continue;
	   	    		throw;
		  		  }
  		    if (gene_finder == "prodigal" && ! contains (dna. name, ";partial=00"))  
  		      continue;
  		    if (dna. getXs ())
  		      continue;
  		    if (translate)
  		    	try 
	  		    {
	    		    Peptide pep (dna. cds2prot (11, true));  // gencode - PAR ??
	    		    pep. qc ();
						  ASSERT (! pep. getXs ());  
	    	      strUpper (pep. seq);  // Start codons are lowercased
		  		    if (protF. get ())
		  		    	pep. saveText (*protF);
	    	      if (reduce)
	    	        reducePep (pep);
	    	      s = pep. seq;
	    	    }
	    	    catch (...)
	    	    {	    	    	
	    	    	if (gene_finder == "prodigal")
	    	    		throw;
	    	      continue;
	    	    }
  		    else
  		      s = dna. seq;
	  		}
    		else
    		{
  		    Peptide pep (f, 1000, false);  // PAR
  		    seqName = pep. name;
  		    pep. qc ();
  		    if (gene_finder == "prodigal" && ! contains (pep. name, ";partial=00"))  
  		      continue;
  		    pep. trimStop ();
   	      strUpper (pep. seq);  // Start codons are lowercased
  		    if (pep. getXs ())
  		      continue;
  		    if (contains (pep. seq, '*'))
  		      throw runtime_error ("Protein " + pep. name + " contains a stop codon");
  		    if (protF. get ())
  		    	pep. saveText (*protF);
  	      if (reduce)
   	        reducePep (pep);
		      s = pep. seq;
    		}
	      if (s. size () < seq_len_min)
	      	continue;

	      sequences++;
        if (kmer)
          FOR (size_t, i, s. size ())
          {
            const string sub (s. substr (i, kmer));
            if (sub. size () < kmer)
              break;
            addHash (sub, seqName, hashes, targetHashes);
          }
        else
          addHash (s, seqName, hashes, targetHashes);
	    }
    }
    cout << "Good sequences: " << sequences << endl;
    cout << "All hashes: " << hashes. size () << endl;

    hashes. sort ();
    hashes. uniq ();
    cout << "Unique hashes: " << hashes. size () << endl;
    
    {
      OFStream f (out);
      if (max_hashes)
      {
        const size_t size = min (hashes. size (), max_hashes);
        FOR (size_t, i, size)
          f << hashes [i] << endl;
      }
      else
        for (const auto& h : hashes)
          f << h << endl;
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



