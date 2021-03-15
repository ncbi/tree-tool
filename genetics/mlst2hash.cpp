// mlst2hash.cpp

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
*   Print 64-bit MLST hashes
*
*/


#undef NDEBUG 
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



namespace 
{


struct Blastn
{
  string seq;
  size_t len {0};
  
  Blastn () = default;
  explicit Blastn (string &&seq_arg)
    : seq (move (seq_arg))
    , len (seq. size ())
    {}
};




// ThisApplication

constexpr const char* noType = "NONE";



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print 64-bit MLST hashes")
    {
      version = VERSION;
      addPositional ("loci", "Number of MLST loci");
      addPositional ("blastn", "BLASTN output in the format: qseqid sseq");  
      addPositional ("out", "Output file with 64-bit MLST hash, or " + strQuote (noType) + " if type cannot be computed");
    }



  void body () const final
  {
    const size_t loci        = str2<size_t> (getArg ("loci"));
    const string blastnFName = getArg ("blastn");
    const string outFName    = getArg ("out");
    

    static const string dnaAlphabet ("ACGT");
    map <string/*locus*/,Blastn> locus2blastn;
    {
      LineInput li (blastnFName);
      Istringstream iss;
      string locus;
      string seq;
      while (li. nextLine ())
      {
        iss. reset (li. line);
        locus. clear ();
        seq. clear ();
        iss >> locus >> seq;
        if (seq. empty ())
          throw runtime_error ("End of file");
        replaceStr (seq, "-", ""); 
        if (stringInSet (seq, dnaAlphabet) != seq. end ())
          continue;
        if (locus2blastn [locus]. len < seq. size ())
          locus2blastn [locus] = move (Blastn (move (seq)));
      }
    }
    if (locus2blastn. size () > loci)
      throw runtime_error ("Too many loci");

      
    OFStream out (outFName);
    if (locus2blastn. size () < loci)
      out << noType << endl;
    else
    {
      string s;
      for (const auto& it : locus2blastn)
        s += it. second. seq + "#";
      hash<string> hashFunc;
      out << hashFunc (s) << endl;
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



