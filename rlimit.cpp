// rlimit.cpp

#undef NDEBUG

#include "common.hpp"
using namespace Common_sp;

#include <sys/resource.h>

#include "common.inc"



namespace
{
  
  
string lim2str (rlim_t rlim)
{
  if (rlim == RLIM_INFINITY)
    return "INF";
  return to_string (rlim);
}



void getRLimit (int resource,
                const string &name)
{
  rlimit rlim;
  EXEC_ASSERT (getrlimit (resource, & rlim) == 0);
  cout << name << '\t' << lim2str (rlim. rlim_cur) << "\t" << lim2str (rlim. rlim_max) << '\n';
}



struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Print tsv-table with resource limits")
  	{
  	  addPositional ("go", "Argument");
  	}
  


 	void body () const final
	{
  //const string arg = getArg ("arg");
	  
	  cout << "#Resource\tCurLimit\tMaxLimit\n";
    #define GET_RLIMIT(S)   getRLimit (RLIMIT_ ## S, #S); 
    GET_RLIMIT (CPU);
    GET_RLIMIT (DATA);
    GET_RLIMIT (FSIZE);
    GET_RLIMIT (LOCKS);
    GET_RLIMIT (MEMLOCK);
    GET_RLIMIT (NOFILE);
    GET_RLIMIT (NPROC);
    GET_RLIMIT (SIGPENDING);
    GET_RLIMIT (STACK);
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


