// prot_collection2dissim.cpp

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
*   Compute dissimilarity between collections of proteins
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/dataset.hpp"
#include "align.hpp"
#include "evolution.hpp"
using namespace DM_sp;
#include "../genetics/seq.hpp"
using namespace Seq_sp;
#include "../version.inc"




namespace 
{
	
	
struct Fasta
{
  typedef  unordered_map<string/*Peptede::getId()*/,Peptide>  Peptides; 
  Peptides peptides;
  

  Fasta () = default;


  void readPeptides (const string &fName,
                     bool sparse,
                     size_t num
                  //,const StringVector& pamNames
                    )
    {
      peptides. rehash (num);
      Unverbose unv;
      Multifasta faIn (fName, true);
      while (faIn. next ())
      {
        Peptide pep (faIn, 1000/*PAR*/, sparse);
      //pep. ambig2X ();
        pep. pseudo = true;
        pep. qc ();
        const string id (pep. getId ());
      #if 0  // Some HMMs may be discarded by positiveAverage: fix positiveAverage.cpp to keep them!
        if (! pamNames. contains (id))
          throw runtime_error ("Unknown component name: " + strQuote (id));
      #endif
        peptides [id] = move (pep);
      }
    }
};
	

	
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Compute dissimilarity between collections of proteins")
    {
      version = VERSION;
  	  // Input
  	  addPositional ("prot_info", string ("Format: ") + PositiveAverageModel::format + ", where <attrname>'s are the protein identifiers in <in>");
  	  addPositional ("in", "File with pairs of protein FASTA file names. Protein sequences may contain stop codons");
  	  addFlag ("blosum62", "Use BLOSUM62, otherwise PAM30");
  	  // Output
  	  addPositional("out", "Output list of dissimilarities for each identifier in <prot_list> if -separate, or one dissimilarity otherwise");  
  	  addFlag ("separate", "Print scores for each target protein separately");
  	  addKey ("raw_power", "Raise raw (not averaged) dissimilarities to this power", "1.0");
  	  addKey ("coeff", "Coefficient to multiply the averaged dissimilarity by", "1.0");
    }


	
	void body () const final
  {
	  const string prot_infoFName = getArg ("prot_info");
	  const string inFName        = getArg ("in");
	  const bool blosum62         = getFlag ("blosum62");

	  const string outFName       = getArg ("out");
	  const bool separate         = getFlag ("separate");
	  const Real raw_power        = str2real (getArg ("raw_power"));
	  const Real coeff            = str2real (getArg ("coeff"));

    if (raw_power <= 0.0)
      throw runtime_error ("-raw_power must be positive");    
    if (separate && raw_power != 1.0)
    	throw runtime_error ("-raw_power must be 1 if -separate");
	  if (coeff <= 0.0)
      throw runtime_error ("-coeff must be positive");    
    if (separate && coeff != 1.0)
    	throw runtime_error ("-coeff must be 1 if -separate");


    PositiveAverageModel pam (prot_infoFName, ! separate);  
    cout << "# Targets: " << pam. components. size () << endl;
    pam. qc ();
    
    StringVector pamNames;
    for (PositiveAverageModel::Component& comp : pam. components)
    	pamNames << comp. name;
    pamNames. sort ();
    QC_ASSERT (pamNames. isUniq ());    


    Vector<Pair<string>> pairs;  pairs. reserve (10000);  // PAR
    map<string,Fasta> file2fasta;
    {
      PairFile pf (inFName, false, false);
      while (pf. next ())
      {
        if (! contains (file2fasta, pf. name1))
          file2fasta [pf. name1] = move (Fasta ());
        if (! contains (file2fasta, pf. name2))
          file2fasta [pf. name2] = move (Fasta ());
        pairs << move (Pair<string> (move (pf. name1), move (pf. name2)));
      }
    }

    for (auto& it : file2fasta)
      it. second. readPeptides (it. first, false, pam. components. size () /*, pamNames*/);

        
    Progress prog (pairs. size ());
    OFStream out (outFName);
    const ONumber on (out, 6, true);  // PAR
    for (const auto& it : pairs)
    {
      prog ();
      out         << it. first 
          << '\t' << it. second 
          << '\t';
      ASSERT (contains (file2fasta, it. first));
      ASSERT (contains (file2fasta, it. second));
      const Fasta::Peptides& peptides1 = file2fasta [it. first].  peptides; 
      const Fasta::Peptides& peptides2 = file2fasta [it. second]. peptides; 
      ASSERT (peptides1. bucket_count () >= pam. components. size ());
      ASSERT (peptides2. bucket_count () >= pam. components. size ());
     	TabDel td (6, false);  // PAR
     	pam. clearValues ();
      for (PositiveAverageModel::Component& comp : pam. components)
      	if (   contains (peptides1, comp. name)
      		  && contains (peptides2, comp. name)
      		 )
      	{
      		const Peptide* pep1 = findPtr (peptides1, comp. name);
      		const Peptide* pep2 = findPtr (peptides2, comp. name);
      		ASSERT (pep1);
      		ASSERT (pep2);
      		
      		Real dissim = NaN;
      		try
      		{
    				Align_sp::Align al (*pep1, *pep2, false, 0, blosum62);  // PAR
    	      dissim = al. getDissim ();  
    	      if (verbose ())
    	      {
    	        al. setAlignment (pep1->seq, pep2->seq); 
    	        al. printAlignment (60);  // PAR
    	        cout << "dissim = " << dissim << endl;
    	      }
    	    }
    	    catch (const exception &e)
    	    {
    	      throw runtime_error (pep1->str () + "\n" + pep2->str () + "\n" + e. what ());
    	    }
  	      IMPLY (! isNan (dissim), dissim >= 0.0);
  	      
  	    	if (separate)
     	  		td << dissim;
     	  	else
  	    	  comp. setValue (pow (dissim, raw_power));
      	}
      	else
      	{
  	    	if (separate)
     	  		td << "nan";
      	}
      pam. qc ();


    	if (separate)
      	out << td. str ();
      else
      {
    		out << coeff * pam. get ();
    		if (verbose ())
    		  pam. saveText (cout);
      }
    		
      out << endl;
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



