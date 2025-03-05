// objHash_find.cpp

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
*   Print closest objects
*
*/


#undef NDEBUG 

#include "common.hpp"
using namespace Common_sp;
#include "version.inc"

#include "common.inc"



namespace 
{


struct Object : Named
{
  size_t freq;
  
  Object (const string &name_arg,
          size_t freq_arg)
    : Named (name_arg)
    , freq (freq_arg)
    { 
      ASSERT (freq); 
    }
};



bool operator< (const Object& g1,
                const Object& g2)
{
  LESS_PART (g2, g1, freq);
  LESS_PART (g1, g2, name);
  return false;
}



// ThisApplication

struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Print closest objects")
    {
      version = VERSION;
      addPositional ("obj_hash", "File with lines: <object> <hash>");
      addPositional ("target", "File with unique hashes");
      addPositional ("common_min", "Min. number of common hashes");
      addPositional ("out_size", "Number of output objects");
      addKey ("superset", "Superset of <object>'s");
    }



  void body () const final
  {
    const string objectHashFName = getArg ("obj_hash");
    const string targetFName     = getArg ("target");
    const size_t common_min      = (size_t) arg2uint ("common_min");
    const size_t outSize         = (size_t) arg2uint ("out_size");
    const string supersetFName   = getArg ("superset");
    
    
    if (! outSize)
      return;


    StringVector superset;
    if (! supersetFName. empty ())
    {
      LineInput f (supersetFName);  
      while (f. nextLine ())
      {
        trim (f. line);
        if (f. line. empty ())
          continue;
        superset << f. line;
      }
    }    
    superset. sort ();
    QC_ASSERT (superset. isUniq ());        

    Vector<size_t> hashes;
    {
      LineInput f (targetFName);  
      while (f. nextLine ())
      {
        const size_t hash = str2<size_t> (f. line);
        QC_ASSERT (hash);
        hashes << hash;
      }
    }    
    hashes. sort ();
    ASSERT (hashes. isUniq ());

    map<string,size_t> object2freq;   
    {
      LineInput f (objectHashFName, 100000);  
      Istringstream iss;
      string object;
      size_t hash;
      while (f. nextLine ())
      {
        iss. reset (f. line);
        hash = 0;
        iss >> object >> hash;
        QC_ASSERT (hash);
        if (   hashes. containsFast (hash)
            && (supersetFName. empty () || superset. containsFast (object))
           )
          object2freq [object] ++;
      }
    }    
    
    Vector<Object> objects;  objects. reserve (object2freq. size ());
    for (const auto& it : object2freq)
      if (it. second >= common_min)
        objects << std::move (Object (it. first, it. second));
        
    objects. sort ();  
    FFOR (size_t, i, min (objects. size (), outSize))
    {
      const Object& g = objects [i];
      cout << g. name;
      if (verbose ())
        cout << '\t' << g. freq;
      cout << '\n';
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



