// genbank_text.cpp

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
*   GenBank format parsing
*
*/


#undef NDEBUG
#include "../common.inc"

#include <cmath>

#include "genbank_text.hpp"
using namespace Common_sp;

 

namespace Asn_sp
{


// GenbankText

GenbankText::Record::Record (const string &name_arg,
  	                         bool isFeature_arg)
: Named (name_arg)
, isFeature (isFeature_arg)
{ 
	if (verbose ())
		cout << name << endl;
	ASSERT (! name. empty ()); 
#if 0
	if (isFeature)
	{ 
		ASSERT (isLower (name)); 
	}
#if 0
	else
		ASSERT (isUpper (name));
#endif
#endif
}



void GenbankText::Record::addValue (const string &s)
{ 
	if (s. empty ())
	  return;
	if (! value. empty ())
		value += " ";
	value += s;
}



string GenbankText::Record::getValue () const
{ 
	if (! isFeature)
	  return value;
	ASSERT (isQuoted (value));
	return unQuote (value);
}



GenbankText::GenbankText (LineInput &f)
: genomeCoverage (NAN)
{
  records. reserve (1024);  // PAR
  
	const size_t nameLen = 12;
	Record* record = nullptr;
	bool featuresP = false;
	while (   f. nextLine () 
	       && f. line != "//"
	       && f. line != "ORIGIN"
	      )
	{
		string name (f. line. substr (0, nameLen));
		trim (name);
		string value;
		if (f. line. size () >= nameLen)
		  value = f. line. substr (nameLen);
		trim (value);

  	if (featuresP)
  	{
  		if (name == "source")
  		{
  		//ASSERT (value == "1..50");
  			record = nullptr;
  		}
  		else if (name. empty ())
  		{
  			if (value [0] == '/')
  			{
  				const size_t pos = value. find ('=');
  				if (pos == string::npos)
  				{
	  				name = value. substr (1);
  					value = "\"1\"";
  				}
  				else
  				{
	  				name = value. substr (1, pos - 1);
	  				value = value. substr (pos + 1);
	  			}
  				records << std::move (Record (name, true));
  				record = & records. back ();
  			}
  		}
  		else
  			featuresP = false;
  	}
  		
  	if (! featuresP)
  	{
	  	if (name == "FEATURES")
	  	{
  			featuresP = true;
  			ASSERT (value == "Location/Qualifiers");
  			record = nullptr;
  		}
	    else if (! name. empty ())
	    {
  			records << std::move (Record (name, false));
 				record = & records. back ();
 		  }
  	}
  		
  	if (record)
  	{
		  record->addValue (value);
		  if (record->name == "COMMENT")
		  {
  		  const string assemblerS      ("Assembly Method       :: ");
  		  const string sequencerS      ("Sequencing Technology :: ");
  		  const string genomeCoverageS ("Genome Coverage       :: ");
		         if (isLeft (value, assemblerS))
		      assembler = value. substr (assemblerS. size ());
		    else if (isLeft (value, sequencerS))
		      sequencer = value. substr (sequencerS. size ());
		    else if (isLeft (value, genomeCoverageS))
		    {
		      string s (value. substr (genomeCoverageS. size ()));
		      trim (s);
		      strLower (s);
		      string prefix (findSplit (s));
		      if (   prefix != "unknown"
		          && prefix != "low"
		          && prefix != "not"
		         )
		      {
  		      if (prefix [prefix. size () - 1] == 'x')
  		        prefix = prefix. substr (0, prefix. size () - 1);
  		      if (   prefix [0] == '<'
  		          || prefix [0] == '>'
  		         )
  		        prefix = prefix. substr (1);
  		      if (isDigit (prefix [0]))
  		        genomeCoverage = atof (prefix. c_str ());
  		    }
		    }
		  }
		}
	}
}


void GenbankText::saveText (ostream &os) const
{
	for (const Record& rec : records)
	  rec. saveText (os);
	os << "Sequencer: " << sequencer << endl;
	os << "Assembler: " << assembler << endl;
	os << "Coverage:  " << genomeCoverage << endl;
}



string GenbankText::name2value (const string &name,
                                bool isUniq) const
{
	const Record* r = nullptr;
	for (const Record& rec : records)
	  if (rec. name == name)
	  {
	  	if (r)
	  	{
	  		cout << "Duplicate name " << name << endl;
	  		ERROR;
	  	}
	  	r = & rec;
	  	if (! isUniq)
	  		break;
	  }
  return r ? r->getValue () : string ();
}



string GenbankText::keyword2name (const string &keyword,
                                  bool caseSensitive,
                                  string &content) const
{
  string keyword1 (keyword);
  if (! caseSensitive)
    strUpper (keyword1);
  for (const Record& rec : records)
  {
    content = rec. getValue ();
    if (! caseSensitive)
      strUpper (content);
    if (contains (content, keyword1))
      return rec. name;
  }
  
  return string ();
}



}
