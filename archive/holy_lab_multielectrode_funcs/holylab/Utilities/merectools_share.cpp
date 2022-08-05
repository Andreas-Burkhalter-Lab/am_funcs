#include <iostream>
#include <string>

using namespace std;

//-------------------------------------------------------------------

bool getValue(string& lines, string key, string& outValue)
{
   string tKey=key+'=';//doesn't handle spaces between key and =
   int pos=lines.find(tKey);
   if(pos==string::npos)return false;
   pos+=tKey.size();
   outValue="";
   while(lines[pos]!='\n')outValue+=lines[pos++];
   return true;

}
//-------------------------------------------------------------------
int getHeadSize(char buf[1024])
{
   buf[1023]='\0';
   string ts=buf;
   string result;
   if(!getValue(ts, "header size", result))return -1;

   return atoi(result.c_str());
}

//-------------------------------------------------------------------


