// pairs2tsv.cpp

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
*   Convert a directory with attribute-value files to a tsv-table
*
*/

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"
#include "tsv.hpp"



namespace
{
  
  
struct Attr
{
  string name;
  size_t num;
  
  void saveText (ostream &os) const
    { os << name << '\t' << num << endl; }
  bool operator< (const Attr& other) const
    { LESS_PART (other, *this, num);
      LESS_PART (*this, other, name);
      return false;
    }
};

  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert a directory with attribute-value files to a tsv-table")
  	{
      version = VERSION;
  	  addPositional ("dir", "directory with attribute-value files; separators are ':', '=' or tab");
  	  addKey ("file_col", "Column name for the file name");
  	}
  	
  	
 
	void body () const final
	{
		const string dirName = getArg ("dir");
		const string fileCol = getArg ("file_col");
		
		
    constexpr const char* sep (":=\t");

    map<string,size_t> attr2num;
    {
      DirItemGenerator dir (100, dirName, false);  // PAR
      string fName;
      while (dir. next (fName))
      {
        LineInput f (dirName + "/" + fName);
        while (f. nextLine ())
        {
          const string err ("Error at " + f. lineStr ());
          trim (f. line);
          const size_t pos = f. line. find_first_of (sep);
          if (pos == string::npos)
            throw runtime_error (err + ": no attribute-value separator " + strQuote (sep));
          string attr (f. line. substr (0, pos));
          trim (attr);
          if (attr. empty ())
            throw runtime_error (err + ": empty attribute");
          attr2num [attr] ++;
        }
      }
    }
    if (attr2num. empty ())
      throw runtime_error ("No attributes");
    
    TextTable tt;
    tt. pound = true;
    map<string,TextTable::ColNum> attr2col;
    if (! fileCol. empty ())
    {
      attr2col [fileCol] = 0;
      tt. header << TextTable::Header (fileCol);
    }
    {
      Vector<Attr> attrs;  attrs. reserve (attr2num. size ());
      for (const auto& it : attr2num)
        attrs << Attr {it. first, it. second};
      attrs. sort ();
      for (const Attr& attr: attrs)
      {
        if (contains (attr2col, attr.name))
          throw runtime_error ("Duplicate column name " + attr. name);
        attr2col [attr. name] = tt. header. size ();
        tt. header << TextTable::Header (attr. name);
      }
    }
    
    {
      DirItemGenerator dir (100, dirName, false);  // PAR
      string fName;
      while (dir. next (fName))
      {
        tt. rows << StringVector (tt. header. size ());
        StringVector& row = tt. rows. back ();
        LineInput f (dirName + "/" + fName);
        while (f. nextLine ())
        {
          if (! fileCol. empty ())
            row [0] = fName;
          trim (f. line);
          const size_t pos = f. line. find_first_of (sep);
          ASSERT (pos != string::npos);
          string attr (f. line. substr (0, pos));
          trim (attr);
          ASSERT (! attr. empty ());
          ASSERT (contains (attr2col, attr));
          const TextTable::ColNum col = attr2col [attr];          
        //if (! row [col]. empty ()) 
          //throw runtime_error ("Error at line " + to_string (f. lineNum + 1) + ": multiple attribute " + strQuote (attr));
          string val (f. line. substr (pos + 1));
          trim (val);
          if (val. size () >= 2)
            if (   val. front () == '"' 
                && val. back  () == '"' 
               )
            {
              val. erase (val. size () - 1);
              val. erase (0, 1);
            }              
          if (! row [col]. empty ()) 
            row [col] += "; ";
          row [col] = val;
        }
      }
    }
    tt. qc ();
    
    tt. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



