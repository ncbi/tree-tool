// asn.cpp

#undef NDEBUG
#include "common.inc"

#include "asn.hpp"
using namespace Common_sp;

 

namespace Asn_sp
{


struct Brace : Nocopy
{
private:
	Asn &asn;
public:
	explicit Brace (Asn &asn_arg)
	  : asn (asn_arg)
	  { asn. expect ("{"); }
 ~Brace ()
	  { if (! uncaught_exception ())
	  	  asn. expect ("}"); 
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
	  	const Token id (asn. in, Token::eInteger);
	  	if (id. n < 0)
	  	  asn. in. error ("Non-negative integer");
	  	const size_t n = (size_t) id. n;
	    expect (",");
	    expect ("name");
	    const Token text (asn. in, Token::eText);
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
	  	const Token id (asn. in, Token::eInteger);
	  	if (id. n < 0)
	  	  asn. in. error ("Non-negative integer");
	  	const size_t n = (size_t) id. n;
	    expect (",");
	    expect ("value");
	  	const Token value (asn. in, Token::eText);
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
	  { prog ();
	  //if (prog. n >= 10)
	  	//exit (1);  
	  	expect ("id");
	  	{	const Token id (asn. in, Token::eInteger);
  	  	if (id. n < 0)
  	  	  asn. in. error ("Non-negative integer");
	  	  asn. id = (uint) id. n;
	  	}
      asn. parent = Asn::no_id;
      asn. features. clear ();
      asn. features. resize (asn. fDict. size ());
	    if (expectTry (","))
	    { if (expectTry ("parent"))
		    { const Token parent (asn. in, Token::eInteger);
    	  	if (parent. n < 0)
    	  	  asn. in. error ("Non-negative integer");
	  	  	asn. parent = (uint) parent. n;
	  	  }
		    if (expectTry (","))
		    { expect ("features");
			    Features features (asn);
			    features. asnRead ();
			  }
		  }
	    asn. processNode ();
	  }
};




// Asn

void Asn::asnRead ()
{
	const Token titleToken (in, Token::eName);
	title = titleToken. name;
	expect (":");
	expect (":");
	expect ("=");
	
	Brace br1 (*this);
	
	expect ("fdict");
  FeatureDictionary fDictList (*this);
  fDictList. asnRead ();

	expect (",");    	
	expect ("nodes");
	Nodes nodes (*this);
	nodes. asnRead ();
}



bool Asn::expectTry (const string& text)
{ 
	ASSERT (! text. empty ());
	
	if (last. empty ()) 
		last = Token (in);
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



void Asn::printFeatures (ostream &os) const
{
	ASSERT (fDict. size () == features. size ());
	FOR (size_t, i, fDict. size ())
	  if (! features [i]. empty ())
	  	os << fDict [i] << ": " << features [i] << endl;
}



}
