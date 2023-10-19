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
#include "../genetics/seq.hpp"
using namespace Seq_sp;
#include "../version.inc"



namespace 
{
  
  
 
map<string/*Mutation::geneName*/,Vector<Vector<Mutation>>> gene2pos2muts;  
  // !Mutation::ambig



void save (const Vector<Mutation> &muts, 
           ostream &os) 
{ 
  Set<Mutation> realMuts;
  Set<Mutation> addedMuts;
  for (const Mutation& mut : muts)
    if (mut. ambig)
    {
      const Vector<Vector<Mutation>>& pos2muts = gene2pos2muts [mut. geneName];      
      FFOR_START (size_t, i, mut. pos, min (mut. stop (), pos2muts. size ()))
        for (const Mutation& refMut : pos2muts [i])
        {
          ASSERT (! refMut. ambig);
          ASSERT (refMut. pos == i);
          if (muts. containsFast (refMut))
          {
            ASSERT (mut. prot);
            continue;
          }
          if (refMut. stop () <= mut. stop ())  // refMut.ref is inside mut.ref
            addedMuts << refMut;
        }
    }
    else
      realMuts << mut;
    
  // !Mutation::ambig
  for (const Mutation& mut : realMuts)
    os << mut << " 0" << endl;    
  for (const Mutation& mut : addedMuts)
    if (! realMuts. contains (mut))
      os << mut << " 1" << endl;    
}

  


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert mutations to 3-valued Boolean attributes of non-ambiguous mutations for makeFeatureTree")
	  {
      version = VERSION;
	  	addPositional ("in", "Input directory with DNA mutations in the format <ref seq><ref pos><target seq>");
	  	addPositional ("feature_dir", "Output directory for 3-valued Boolean atrtributes for makeFeatrueTree");
	  	addFlag ("aa", "Mutations of protein sequences, otherwise of DNA sequences");
	  	addFlag ("append", "Append existing files in <feature_dir>");
  	  addFlag ("large", "Create files in subdirectories \"0\" .. \"" + to_string (hash_class_max - 1) + "\" which are the hashes of file names");
	  }



  void body () const final
  {
    const string inDirName  = getArg ("in");
    const string outDirName = getArg ("feature_dir");
    const bool   aa         = getFlag ("aa");
    const bool   append     = getFlag ("append");
		const bool   large      = getFlag ("large");


    unordered_map<string,Vector<Mutation>> obj2muts;  obj2muts. rehash (100000);  // PAR
    {
      map<string/*Mutation::geneName*/,unordered_set<Mutation,Mutation::Hash>> gene2muts;  
        // !Mutation::ambig
      {
        DirItemGenerator dig (1000, inDirName, large);  // PAR
        string fName;
    	  while (dig. next (fName))
        {
          Vector<Mutation> muts;  muts. reserve (16);  // PAR
          {
            LineInput li (inDirName + (large ? "/" + to_string (str2hash_class (fName)) : "") + "/" + fName);
            while (li. nextLine ())
            {
              Mutation mut (aa, li. line);
              try { mut. qc (); }
                catch (const exception &e)
                  { throw runtime_error (fName + "\n" + li. line + "\n" + e. what ()); }
              if (! mut. ambig)
              {
                unordered_set<Mutation,Mutation::Hash>& gene_muts = gene2muts [mut. geneName];
                if (gene_muts. empty ())
                  gene_muts. rehash (10000);  // PAR
                gene_muts. insert (mut);
              }
              muts << std::move (mut);
            }
          }
          muts. sort ();
          muts. uniq ();
          obj2muts [fName] = std::move (muts);        
        }
      }
      cout << "Objects: " << obj2muts. size () << endl;
      
      // gene2pos2muts
      for (const auto& it : gene2muts)
      {
        const unordered_set<Mutation,Mutation::Hash>& gene_muts = it. second;
        ASSERT (! gene_muts. empty ());
        size_t pos_max = 0;
        const bool prot = (* gene_muts. begin ()). prot;
        for (const Mutation& mut : gene_muts)
        {
          ASSERT (mut. prot == prot);
          maximize (pos_max, mut. pos);
        }
        Vector<Vector<Mutation>>& pos2muts = gene2pos2muts [it. first];
        pos2muts. resize (pos_max + 1);
        for (const Mutation& mut : gene_muts) 
          pos2muts [mut. pos] << mut;
      }
    }
    
    
    {
      Progress prog (obj2muts. size (), 100);  // PAR
      for (const auto& it : obj2muts)
      {
        prog (it. first);
        string dir (outDirName);
        if (large)
        {
          dir += "/" + to_string (str2hash_class (it. first));
          Dir (dir). create ();
        }
        const ios_base::openmode mode = append ? ios_base::app : ios_base::out;
        ofstream of (dir + "/" + it. first, mode);
        try { save (it. second, of); }
          catch (const exception &e)
            { throw runtime_error (it. first + "\n" + e. what ()); }
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



