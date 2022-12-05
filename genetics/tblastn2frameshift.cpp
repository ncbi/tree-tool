// tblastn2frameshift.cpp

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
*   Find frame shuifts
*
*/
   
   
#undef NDEBUG 
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



namespace
{
  
  
struct Hsp : Root
{
  string sseqid;
  string qseqid;
  // bp
  size_t sstart {0};
  size_t send {0};
  bool sstrand {true};
  // aa
  size_t qstart {0};
  size_t qend {0};

  
  Hsp (Istringstream &iss)
    { iss >> sseqid >> qseqid >> sstart >> send >> qstart >> qend;
      QC_ASSERT (qend);
      //
      QC_ASSERT (sstart != send);
      sstrand = (sstart < send);
      if (! sstrand)
        swap (sstart, send);
      ASSERT (sstart < send);
      QC_ASSERT (sstart);
      sstart--;
      //
      QC_ASSERT (qstart < qend);
      QC_ASSERT (qstart);
      qstart--;
      //
      ASSERT (! empty ());
    }
  Hsp () = default;
  bool empty () const override
    { return sseqid. empty (); }
  void saveText (ostream &os) const override
    { os         << sseqid 
         << '\t' << qseqid
         << '\t' << sstart
         << '\t' << send
         << '\t' << sstrand
         << '\t' << qstart
         << '\t' << qend
         << endl;
    }
    
    
private:
  size_t sstart_orig () const
    { return sstrand ? sstart + 1 : send; }
public:
  bool operator<= (const Hsp &other) const
    { LESS_PART (*this, other, sseqid);
      LESS_PART (*this, other, qseqid);
      LESS_PART (*this, other, sstart_orig ());
      return true;
    }
  long global_start () const
    { return sstrand
               ? (long) sstart - 3 * (long) qstart
               : (long) send   + 3 * (long) qstart; 
    }
};



} // namespace
 



// ThisApplication

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find frame shifts using tblastn output.\nOutput: sseqid qseqid <aa pos>{ins|del}<bp length>bp")
    {
      version = VERSION;
      // Input
      addPositional ("tblastn", "tblastn output in the format: sseqid qseqid sstart send qstart qend, ordered by: sseqid, qseqid, sstart"); 
    }



  void body () const final
  {
    const string blastFName = getArg ("tblastn");


    constexpr size_t diff_max_aa = 30;  // PAR

    LineInput f (blastFName);
    Istringstream iss;
    Hsp hsp_prev;
	  while (f. nextLine ())
	  {
	    iss. reset (f. line);
	    Hsp hsp (iss);
	    QC_ASSERT (hsp_prev <= hsp);
	    if (   ! hsp_prev. empty ()
	        && hsp_prev. sseqid  == hsp. sseqid
	        && hsp_prev. qseqid  == hsp. qseqid
	        && hsp_prev. sstrand == hsp. sstrand
	        && difference (hsp_prev. send, hsp. sstart) <= 3 * diff_max_aa
	        && (   (  hsp. sstrand && difference (hsp_prev. qend,   hsp. qstart) <= diff_max_aa)
  	          || (! hsp. sstrand && difference (hsp_prev. qstart, hsp. qend)   <= diff_max_aa)
  	         )
  	     )
      {
        const long diff = hsp. global_start () - hsp_prev. global_start ();
        QC_ASSERT (diff % 3);  // Otherwise no frame shift
        cout         << hsp. sseqid 
             << '\t' << hsp. qseqid
             << '\t' << (hsp. sstrand ? hsp_prev. qend : hsp. qend) + 1 << (diff > 0 ? "ins" : "del") << abs (diff) << "bp" 
              << endl;
      }
	    hsp_prev = move (hsp);
	  }
  }
};




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



