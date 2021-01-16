// mutation2feature.cpp

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
*   Convert DNA mutations to 3-values Boolean attributes for makeFeatureTree
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



namespace 
{
  
  
 
struct Mutation : Root
{
  size_t pos {no_index};
    // In reference
    // < no_index
  string ref;
  string allele;
  
  bool ambig {false};
    // Function of allele
  

  explicit Mutation (const string& mut)
    { size_t posStart = no_index;
      size_t alleleStart = no_index;
      FFOR (size_t, i, mut. size ())
        if (isDigit (mut [i]))
        {
          if (posStart == no_index)
            posStart = i;
        }
        else
        {
          if (posStart != no_index)
            if (alleleStart == no_index)
              alleleStart = i;
        }
      QC_ASSERT (posStart < alleleStart);
      QC_ASSERT (alleleStart < no_index);
      ref    =                mut. substr (0, posStart);
      pos    = (size_t) stoi (mut. substr (posStart, alleleStart - posStart));
      allele =                mut. substr (alleleStart);
      ASSERT (pos < no_index);
      QC_ASSERT (pos);
      pos--;
      if (ref == "INS")
        ref. clear ();
      if (allele == "DEL")
        allele. clear ();
      QC_ASSERT (! (ref. empty () && allele. empty ()));    
      static const string acgt ("acgt");
      for (const char c : ref)
        QC_ASSERT (charInSet (c, acgt));

      // ambig        
      for (const char c : allele)
        if (! charInSet (c, acgt))
          ambig = true;
    }
  Mutation () = default;
  Mutation (const Mutation& ) = default;
  void saveText (ostream &os) const override
    { os << nvl (ref, "INS") << pos + 1 << nvl (allele, "DEL"); }
    

  bool operator== (const Mutation& other) const
    { return    pos    == other. pos
             && ref    == other. ref
             && allele == other. allele;
    }
  bool operator< (const Mutation& other) const
    { LESS_PART (*this, other, pos);
      LESS_PART (*this, other, ref);
      LESS_PART (*this, other, allele);
      return false;
    }
  struct Hash
  {
    size_t operator() (const Mutation &mut) const
      { static hash<string> strHash;
        return mut. pos ^ strHash (mut. ref) ^ strHash (mut. allele);
      }
  };  
  size_t stop () const
    { return pos + ref. size (); }    
};

  
  
void save (const Vector<Mutation> &muts, 
           const Vector<Vector<Mutation>> &pos2muts,
           ostream &os) 
{ 
  for (const Mutation& mut : muts)
    if (mut. ambig)
    {
      FFOR_START (size_t, i, mut. pos, min (mut. stop (), pos2muts. size ()))
        for (const Mutation& other : pos2muts [i])
        {
          ASSERT (! other. ambig);
          ASSERT (other. pos == i);
          ASSERT (! muts. containsFast (other));
          if (other. stop () <= mut. stop ())  // other.ref is inside mut.ref
          {
            other. saveText (os);
            os << " 1" << endl;
          }
        }
    }
    else
    {
      mut. saveText (os);
      os << " 0" << endl;
    }
}

  


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert DNA mutations to 3-values Boolean attributes for makeFeatureTree")
	  {
	  	addPositional ("in", "Input directory with DNA mutations in the format <ref seq><ref pos><target seq>");
	  	addPositional ("feature_dir", "Output directory for 3-valued Boolean atrtributes for makeFeatrueTree");
	  	addFlag ("ambig5end", "Replace 5' end deletions by N-substitutions");
	  }



  void body () const final
  {
    const string inDirName  = getArg ("in");
    const string outDirName = getArg ("feature_dir");
    const bool   ambig5end  = getFlag ("ambig5end");


    unordered_map<string,Vector<Mutation>> obj2muts;  obj2muts. rehash (100000);  // PAR
    Vector<Vector<Mutation>> pos2muts;  // !ambig
    {
      size_t pos_max = 0;
      unordered_set<Mutation,Mutation::Hash> allMuts;  allMuts. rehash (10000);  // PAR
      {
        FileItemGenerator fig (1000, true, inDirName);  // PAR
        string fName;
    	  while (fig. next (fName))
        {
          Vector<Mutation> muts;  muts. reserve (16);  // PAR
          {
            LineInput li (inDirName + '/' + fName);
            while (li. nextLine ())
            {
              Mutation mut (li. line);
              if (ambig5end && ! mut. pos)
              {
                if (mut. ref. empty ())
                  continue;
                mut. allele = mut. ref;
                for (char& c : mut. allele)
                  c = 'n';
                mut. ambig = true;
              }
              if (! mut. ambig)
              {
                allMuts. insert (mut);
                maximize (pos_max, mut. pos);
              }
              muts << move (mut);
            }
          }
          muts. sort ();
          muts. uniq ();
          obj2muts [fName] = move (muts);        
        }
      }
      cout << "Objects: " << obj2muts. size () << endl;
      cout << "Mutations: " << allMuts. size () << endl;    
      cout << "Max. position: " << pos_max + 1 << endl;
      
      pos2muts. resize (pos_max + 1);
      for (const Mutation& mut : allMuts) 
        pos2muts [mut. pos] << mut;
    }
    
    size_t pos_used = 0;
    for (const auto& muts : pos2muts)
      if (! muts. empty ())
        pos_used++;
    cout << "Positions: " << pos_used << endl;
    
    {
      Progress prog (obj2muts. size (), 100);  // PAR
      for (const auto& it : obj2muts)
      {
        prog (it. first);
        OFStream of (outDirName + "/" + it. first);
        save (it. second, pos2muts, of);
      }
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



