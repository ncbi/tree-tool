// mutation_tab.cpp

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
*   Print mutations table in a tab-delimited fromat
*
*/


#undef NDEBUG 
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



namespace 
{
	
	
void removeMutualGaps (string &a,
                       string &b)
// Due to muscle
{
	ASSERT (a. size () == b. size ());
	
  size_t i = 0;
  while (i < a. size ())
  {
  	if (   a [i] == '-'
  		  && b [i] == '-'
  		 )
  	{
  		a. erase (i, 1);
  		b. erase (i, 1);
  	}
  	else
  		i++;
  }
}

	
	
void cleanAlignment (string &a,
                     const string &b)
{
	ASSERT (! a. empty ());
	ASSERT (a. size () == b. size ());
	ASSERT (a. front () == '*');
	ASSERT (a. back ()  == '*');
	ASSERT (b. front () == '*');
	ASSERT (b. back ()  == '*');
	
	size_t start = 0;
	while (a [start] == '*')
	  start++;
	
	size_t stop = a. size ();
	while (stop && a [stop - 1] == '*')
	  stop--;
	
	ASSERT (start);
  ASSERT (stop < a. size ());
	ASSERT (start < stop);
	
	start--;
	while (   a [start] == '*'
	       && b [start] != '*'
	      )
	{
		a [start] = '-';
		start--;
	}

	while (   a [stop] == '*'
	       && b [stop] != '*'
	      )
	{
		a [stop] = '-';
		stop++;
	}

	size_t j = 0;
	FFOR_START (size_t, i, 1, a. size ())
	{
		if (   a [i]     == '-'
			  && a [i - 1] != '-'
			 )
			j = i;
	  if (   a [i]     != '-'
	  	  && a [i - 1] == '-'
	  	  && b [j] == a [i]
	  	 )
	  {
	  	ASSERT (j < i);
	  	ASSERT (a [j] == '-');
	    swap (a [j], a [i]);	  	
	    j++;
	  }
	}	  	
}

	

string refSeqWhole;



struct BlastAlignment 
{
  // BLASTP global alignment
  string accession;   
  string metadata;
  Vector<string> mutations;
    // size() = refSeq.size()
  

  explicit BlastAlignment (const string &line)
    {
      string targetSeq, refSeq;
      {
		    istringstream iss (line);
		    iss >> metadata >> targetSeq >> refSeq;
		  }
			ASSERT (targetSeq. size () == refSeq. size ());

      accession = findSplit (metadata, '|');
      ASSERT (! accession. empty ());
      ASSERT (! metadata. empty ());
      replace (metadata, '|', '\t');
      replaceStr (metadata, "NULL", "");
      trimSuffix (metadata, "\t");
      
      removeMutualGaps (targetSeq, refSeq);            
      cleanAlignment (targetSeq, refSeq);
      cleanAlignment (refSeq, targetSeq);
      removeMutualGaps (targetSeq, refSeq);            
      
      const bool refSeqWhole_empty = refSeqWhole. empty ();
      
      size_t len = 0;
      for (const char c : refSeq)
      	if (   c != '*'
      		  && c != '-'
      		 )
      	{
      		len++;
      		if (refSeqWhole_empty)
      			refSeqWhole += c;
      	}
      ASSERT (refSeqWhole. size () == len);
      		
      mutations. resize (len);
      
      size_t pos = 0;
      string prefix;
      FFOR (size_t, i, targetSeq. size ())
      {
      	if (   targetSeq [i] == '-'
      	    && refSeq    [i] == '-'
      	   )
      	{
      		cout << targetSeq << endl;
      		cout << refSeq << endl;
      		ERROR;
      	}
        if (targetSeq [i] == refSeq [i])
        	if (prefix. empty ())
        		;
        	else
        	{
        		ASSERT (! pos);
        		ASSERT (mutations [pos]. empty ());
        		mutations [pos] = prefix + targetSeq [i];
        		prefix. clear ();
        	}        	
        else
        {
        	ASSERT (targetSeq [i] != '*');
        	ASSERT (refSeq    [i] != '*');
        	if (refSeq [i] == '-')
        		if (pos)
	        	{
	        		ASSERT (i);
	        		if (mutations [pos - 1]. empty ())
		        	{
		        		ASSERT (targetSeq [i - 1] != '-');
		        		ASSERT (targetSeq [i - 1] != '*');
		        		mutations [pos - 1] = targetSeq [i - 1];	        	
		        	}
	        	  mutations [pos - 1] += targetSeq [i];
		        }
		        else
		        	prefix += targetSeq [i];
        	else
        	  mutations [pos] += targetSeq [i];
        }
        if (   refSeq [i] != '*'
        	  && refSeq [i] != '-'
        	 )
        	pos++;
      }
      ASSERT (pos == len);
    }
};




// ThisApplication

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print mutations table in a tab-delimited fromat")
    {
      version = VERSION;
      addPositional ("blastp", "blastp output in the format: qseqid qseq sseq. qseqid has format: accession|.... Sequences should have flanking stop codons");
      addPositional ("ref_accession", "Accession of the reference protein");
      addPositional ("metadata_headers", "Column names for metadata, separated by '|'");
    }



  void body () const final
  {
    const string blastpFName      = getArg ("blastp");
    const string refAccession     = getArg ("ref_accession");  
          string metadata_headers = getArg ("metadata_headers");
          
    replace (metadata_headers, '|', '\t');
    
    
	  Vector<BlastAlignment> blastAls;
    {
      LineInput f (blastpFName, 1);  // PAR
  	  while (f. nextLine ())
  	    try { blastAls << BlastAlignment (f. line); }
  	      catch (...)
  	      {
  	      	cout << f. line << endl;
  	      	throw;
  	      }
  	}
  	ASSERT (! blastAls. empty ());
  	
  	const size_t len = blastAls [0]. mutations. size ();
  	ASSERT (len == refSeqWhole. size ());
  	
  	Vector<size_t> positions;  positions. reserve (len);
  	{
	  	Vector<bool> hasMutation (len, false);
	  	for (const BlastAlignment& al : blastAls)
	  	{
	  		ASSERT (al. mutations. size () == len);
				FOR (size_t, i, len)
				  if (! al. mutations [i]. empty ())
				  {
				  	hasMutation [i] = true;
			  	  ASSERT (al. accession != refAccession);  
				  }
	  	}	  	
			FOR (size_t, i, len)
			  if (hasMutation [i])
	  	    positions << i;
	  }  	    
  	if (positions. empty ())
  		return;
  		
  	cout << "Accession\t" << metadata_headers;
    for (size_t pos : positions)
    	cout << '\t' << pos + 1;
    cout << endl;
    
    bool refFound = false;
  	for (const BlastAlignment& al : blastAls)
  		if (al. accession == refAccession)
	  	{
	  		cout << al. accession << '\t' << al. metadata;
		    for (size_t pos : positions)
		    	cout << '\t' << refSeqWhole [pos];
		    cout << endl;	  		
	  		refFound = true;
	  		break;
	  	}
	  ASSERT (refFound);
  	
  	for (const BlastAlignment& al : blastAls)
  		if (al. accession != refAccession)
	  	{
	  		cout << al. accession << '\t' << al. metadata;
	  		bool hasMutation = false;
		    for (size_t pos : positions)
		    {
		    	const string& mut = al. mutations [pos];
		    	IMPLY (contains (mut, '-'), mut. size () == 1);
		    	cout << '\t' << (mut == "-" ? "Del" : mut);
		    	if (! mut. empty ())
		    		hasMutation = true;
		    }
		    ASSERT (hasMutation);
		    cout << endl;	  		
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



