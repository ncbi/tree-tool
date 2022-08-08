// fsrteam_test.cpp

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
*   fstream test
*
*/

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "version.inc"



namespace
{


constexpr size_t BUFFER_SIZE = 1024 * 1024;  // PAR
constexpr size_t ITER_COUNT = 100;



template <typename TType>
void test_write (const TType data [], 
                 char const* const file, 
                 char const* const message)
{
	const Chronometer_OnePass cop (message, cerr, false);  
  {
    std::basic_ofstream<TType, std::char_traits<TType>> ofile (file, std::ios::binary);
    FOR (size_t, i, ITER_COUNT)
      ofile. write (data, BUFFER_SIZE);
  }
}



#if 0
template <typename TType>
void test_read (TType data [], 
                char const* const file,
                char const* const message)
{
	const Chronometer_OnePass cop (message, cerr, false);  
	{
    std::basic_ifstream<TType, std::char_traits<TType>> ifile (file /*, std::ios::binary*/);
    FOR (size_t, i, ITER_COUNT)
    {
      ifile. read (data, BUFFER_SIZE);
      ifile. seekg (BUFFER_SIZE, ios_base::cur);
    }
  }
}
#endif



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("fstream test")
  	{
      version = VERSION;
      addPositional ("go", "Go");
  	}
  	
  	
 
	void body () const final
	{
	  // From: https://github.com/microsoft/STL/issues/817
	  
	  {
      auto data_char = new char [BUFFER_SIZE];      
      FOR (size_t, i, BUFFER_SIZE)
        data_char [i] = '\n';
      test_write<char>    (data_char,    "io_char.dat",    "Write to ofstream<char>");
      cout << endl;
    //test_read<char>    (data_char,    "io_char.dat",    "Read from ifstream<char>");
    }

    constexpr size_t buf_size = BUFFER_SIZE * 50;  // PAR
    auto buf = new char [buf_size];   

    {
      const Chronometer_OnePass cop ("basic_ifstream", cerr, false);
      std::basic_ifstream<char, std::char_traits<char>> f ("io_char.dat");
      while (f. get () != EOF);
      ASSERT (f. eof ());
    }
    
    {
      const Chronometer_OnePass cop ("FILE", cerr, false);
      FILE* f = fopen ("io_char.dat", "r");
      while (fgetc (f) != EOF);
      ASSERT (feof (f));
    }    

    {
      const Chronometer_OnePass cop ("FILE: buf", cerr, false);
      FILE* f = fopen ("io_char.dat", "r");
      setbuffer (f, buf, buf_size);
      while (fgetc (f) != EOF);
      ASSERT (feof (f));
    }    

    {
      const Chronometer_OnePass cop ("ifstream", cerr, false);
      ifstream f ("io_char.dat");
      while (f. get () != EOF);
      ASSERT (f. eof ());
    }    

    {
      const Chronometer_OnePass cop ("ifstream: buf", cerr, false);
      ifstream f ("io_char.dat");
      f. rdbuf () -> pubsetbuf (buf, buf_size);
      while (f. get () != EOF);
      ASSERT (f. eof ());
    }    
    
    {
      const Chronometer_OnePass cop ("ifstream: getline()", cerr, false);
      ifstream f ("io_char.dat");
      string line;
      do
        getline (f, line);
      while (! f. eof ());
    }    

    {
      const Chronometer_OnePass cop ("istream", cerr, false);
      ifstream f ("io_char.dat");
      istream* is = & f;
      string line;
      do
        getline (*is, line);
      while (! is->eof ());
    }    

    cout << endl;
    {
      const Chronometer_OnePass cop ("LineInput", cerr, false);
      LineInput f ("io_char.dat", buf_size);
      while (f. nextLine ())
        ;
      PRINT (f. lineNum);
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



