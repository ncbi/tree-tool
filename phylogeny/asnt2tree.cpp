// asnt2tree.cpp

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
*   Convert a tree in textual ASN.1 format to internal format
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "distTree.hpp"
//#include "asn.hpp"
using namespace DistTree_sp;
#include "../version.inc"



namespace 
{


OFStream* featureF = nullptr;
bool featureHeader = false;



namespace 
{

struct Asn : Root
{
	CharInput in;
	  // ASN.1 text
  string title;    
  StringVector fDict;
    // Index: feature id

  // Current node
  constexpr static const uint no_id {numeric_limits<uint>::max ()};
  uint id {0};
  uint parent {no_id};
  StringVector features;
    // size() = fDict.size()

private:
	Token last;
public:
    

protected:
  explicit Asn (const string &fName)
		: in (fName, 100000)  // PAR
		{}
public:
		
		
  void asnRead ();
    // Invokes: processNode()
  virtual void processNode () = 0;
  bool expectTry (const string& text);
  void expect (const string& text)
    { if (! expectTry (text))
      	in. error (strQuote (text));
    }
//void printFeatures (ostream &os) const;
};




struct Brace : Nocopy
{
private:
	Asn &asn;
public:
	explicit Brace (Asn &asn_arg)
	  : asn (asn_arg)
	  { asn. expect ("{"); }
 ~Brace ()
	  { if (! asn. expectTry ("}"))
      	errorExit (asn. in. errorText (strQuote ("}")). c_str ());
	  }
};




struct AsnList : Root
{
	Asn &asn;
	

protected:
	explicit AsnList (Asn &asn_arg)
	  : asn (asn_arg)
	  {}
public:

	
	void asnRead ()
		{ expect ("{");
			bool first = true;
			while (! expectTry ("}"))
			{
				if (first)
					first = false;
				else
			  	expect (",");
		  	const Brace br (asn);
		  	readElement ();
			}
		}
protected:
	virtual void readElement () = 0;
  bool expectTry (const string& text)
    { return asn. expectTry (text); }
  void expect (const string& text)
    { asn. expect (text); }
};




struct FeatureDictionary : AsnList
{
	explicit FeatureDictionary (Asn &asn_arg)
	  : AsnList (asn_arg)
	  {}

private:
	void readElement () final
	  { expect ("id");
	  	const Token id (asn. in, Token::eInteger, false, false);
	  	if (id. n < 0)
	  	  asn. in. error ("Non-negative integer");
	  	const size_t n = (size_t) id. n;
	    expect (",");
	    expect ("name");
	    const Token text (asn. in, Token::eText, false, false);
	    asn. fDict. resize (n + 1);
	    ASSERT (asn. fDict [n]. empty ());
	    asn. fDict [n] = text. name;
	  }
};




struct Features : AsnList
{
	explicit Features (Asn &asn_arg)
	  : AsnList (asn_arg)
	  {}

private:
	void readElement () final
	  { expect ("featureid");
	  	const Token id (asn. in, Token::eInteger, false, false);
	  	if (id. n < 0)
	  	  asn. in. error ("Non-negative integer");
	  	const size_t n = (size_t) id. n;
	    expect (",");
	    expect ("value");
	  	const Token value (asn. in, Token::eText, false, true);
	  	ASSERT (asn. features [n]. empty ());
	  	asn. features [n] = value. name;
	  }
};




struct Nodes : AsnList
{
	Progress prog;
	
	
	explicit Nodes (Asn &asn_arg)
	  : AsnList (asn_arg)
	  , prog (0, 1000)  // PAR
	  {}


private:
	void readElement () final
	  { 
	    expect ("id");
	  	{	
	  	  const Token id (asn. in, Token::eInteger, false, false);
  	  	if (id. n < 0)
  	  	  asn. in. error ("Non-negative integer");
	  	  asn. id = (uint) id. n;
	  	  prog (to_string (asn. id));  
	  	}
      asn. parent = Asn::no_id;
      asn. features. clear ();
      asn. features. resize (asn. fDict. size ());
	    while (expectTry (","))
	      if (expectTry ("parent"))
		    { 
		      const Token parent (asn. in, Token::eInteger, false, false);
    	  	if (parent. n < 0)
    	  	  asn. in. error ("Non-negative integer");
	  	  	asn. parent = (uint) parent. n;
	  	  }
		    else if (expectTry ("features"))
		    {
			    Features features (asn);
			    features. asnRead ();
			  }
	    asn. processNode ();
	  }
};




// Asn

void Asn::asnRead ()
{
	const Token titleToken (in, Token::eName, false, false);
	title = titleToken. name;
	expect (":");
	expect (":");
	expect ("=");

  expect ("{");
	
	expect ("fdict");
  FeatureDictionary fDictList (*this);
  fDictList. asnRead ();

	expect (",");    	
	expect ("nodes");
	Nodes nodes (*this);
	nodes. asnRead ();
	
	// "label", "user"
}



bool Asn::expectTry (const string& text)
{ 
	ASSERT (! text. empty ());
	
	if (last. empty ()) 
		last = Token (in, true, true);
	ASSERT (! last. empty ());
	if (verbose ())
	  PRINT (last);
  const bool yes =    (   last. type == Token::eName
                       || last. type == Token::eDelimiter
                      )
  	               && last. name == text;
  if (yes)
  	last. clear ();
  return yes;
}



#if 0
void Asn::printFeatures (ostream &os) const
{
	ASSERT (fDict. size () == features. size ());
	FOR (size_t, i, fDict. size ())
	  if (! features [i]. empty ())
	  	os << fDict [i] << ": " << features [i] << endl;
}
#endif


}




struct TreeAsn : Asn
{
  DistTree& tree;
  map<size_t/*id*/,Steiner*> id2steiner;
  

  TreeAsn (const string &fName,
           DistTree &tree_arg)
    : Asn (fName)
    , tree (tree_arg)
    {
      ASSERT (tree. nodes. empty ());
    }

    
  void processNode () final
    { 
      constexpr size_t labelNum   =  0;
    	constexpr size_t distNum    =  1;
    	QC_ASSERT (fDict [labelNum] == "label");
    	QC_ASSERT (fDict [distNum] == "dist");
    	
    	QC_ASSERT (fDict. size () == features. size ());
    	
    	
    	Steiner* parentNode = nullptr;
    	if (parent != no_id)
    	{
    	  parentNode = id2steiner [parent];
    	  if (! parentNode)
    	    throw runtime_error ("Parent id " + toString (parent) + " is not found");
    	}
    	
    	Real len = NaN;
    	if (parentNode && ! features [distNum]. empty ())
    	{
    	  len = max (0.0, str2real (features [distNum]));
    	  ASSERT (len >= 0.0);
    	}    	
    	
   	  const string name (features [labelNum]);
    	if (name. empty ())
    	  id2steiner [id] = new Steiner (tree, parentNode, len);
    	else
    	{
    	  tree. name2leaf [name] = new Leaf (tree, parentNode, len, name);
    	  if (featureF)
    	  {
          if (! featureHeader)
          {
      	    *featureF << '#';
      	    save (*featureF, fDict, '\t');
      	    *featureF << endl;
      	    featureHeader = true;
          }
          save (*featureF, features, '\t');
    	    *featureF << endl;
    	  }
    	}
    }
};




struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Convert a tree in textual ASN.1 format to internal format")
	{
    version = VERSION;
	  // Input
	  addPositional ("input_tree", "Tree in textual ASN.1 format (BioTreeContainer)");
	  addKey ("feature", "Output file with features");
	}



	void body () const final
  {
		const string input_tree   = getArg ("input_tree");    
		const string featureFName = getArg ("feature");


    DistTree tree;
    if (! featureFName. empty ())
      featureF = new OFStream (featureFName);
    {
      TreeAsn asn (input_tree, tree);
      asn. asnRead ();
      asn. qc ();
      const string title ("BioTreeContainer");
      if (asn. title != title)
        throw runtime_error ("Title " + strQuote (title) + " is expected");
    }
    tree. qc ();         
    tree. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


