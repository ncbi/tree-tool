// curl-easy_test.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;

#include "curl_easy.hpp"



namespace
{
  
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Test curl_easy")
  	{
  	  addPositional ("url", "URL to read");
  	}
  


 	void body () const final
	{
    const string url = getArg ("url");


    const SoftwareVersion ver (CURL_sp::getLibVersion ());
    cout << "libcurl version: " << ver << endl << endl;
    
    CURL_sp::Curl curl;
    cout << curl. read (url) << endl;
  }
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


