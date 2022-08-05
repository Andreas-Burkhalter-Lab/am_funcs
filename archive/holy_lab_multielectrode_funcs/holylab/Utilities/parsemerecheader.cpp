
#define MEREC_HEADER_AS_USRHDR //if defined, use full merec header as usrhdr field in old header format
                               //else use comment
#include "../commong/lfs_common.h"

#include "parsemerecheader.h"

#include <iostream>
#include <string>

using namespace std;


//#include "miscpp.h" //this include directive causes problem, so I just copy splitStringSkipEmptyItem() from miscpp.cpp to here. The problem is matlab complains that "undefined symbol: __dynamic_cast_2" in the mex file. Without this include directive, you can use " mex ReadMerecHeader.cpp parsemerecheader.cpp " to create ReadMerecHeader.mexglx
//---------------------------------------------------------------------------
//another implementation of splitString() that uses strtok().
//but this one has problem that it skips empty items:
//      e.g. splitString2("a,,c",",") output "a" "c" instead of "a" "" "c"
vector<string> splitStringSkipEmptyItem(string const & rStr, string const delimiter)
{
   char* p=new char[rStr.size()+1];
   strcpy(p, rStr.c_str());
   vector<string> result;
   p=strtok(p, delimiter.c_str());
   while(p){
      result.push_back(p);
      p = strtok(NULL, delimiter.c_str());
   }

   return result;
}
//--------------------------------------------------

#include "merectools_share.cpp" //for such small shared func sets, don't bother use .h file

//--------------------------------------------------
bool parseMerecHeader(int fd, TMerecHeader& rHeader)
{
   char bufMinHeader[1024];

   off_t posNow=lseek(fd, 0, SEEK_SET);
   int nRead=read(fd, bufMinHeader, 1024);
   if(nRead<1024){perror("read error"); return false;}
   int headsize=getHeadSize(bufMinHeader);
   char *buf=new char[headsize];

   //reread the whole header:
   posNow=lseek(fd, 0, SEEK_SET);
   nRead=read(fd, buf, headsize);
   if(nRead<headsize){perror("read error"); return false;}
   string header=buf;

   //close(fd);

   string ts;

   rHeader.headersize=headsize;

   getValue(header, "nscans", ts);
   rHeader.nscans=atol(ts.c_str());

   getValue(header, "channel list",ts);
   vector<string> channels=splitStringSkipEmptyItem(ts, " ");
   rHeader.numch=channels.size();

   rHeader.channels.clear();
   for(int i=0; i<rHeader.numch; ++i)rHeader.channels.push_back(atoi(channels[i].c_str()));

   getValue(header, "scan rate",ts);
   rHeader.scanrate=atol(ts.c_str());

   int minSample, maxSample;
   getValue(header, "min sample",ts);
   minSample=atol(ts.c_str());
   getValue(header, "max sample",ts);
   maxSample=atol(ts.c_str());
   getValue(header, "min input",ts);
   rHeader.voltageMin=atof(ts.c_str());
   getValue(header, "max input",ts);
   rHeader.voltageMax=atof(ts.c_str());
   rHeader.scalemult=(rHeader.voltageMax-rHeader.voltageMin)/(maxSample-minSample);
   rHeader.scaleoff=rHeader.voltageMin;

   getValue(header, "datetime",ts);
   int pos=ts.find(',');
   if(pos==string::npos){
      pos=ts.find(' '); //assume date and time are separated by space
   }//if, date time is not separated by comma
   rHeader.date=ts.substr(0, pos);
   rHeader.time=ts.substr(pos+1);
   
#if defined(MEREC_HEADER_AS_USRHDR)
   rHeader.usrhdr=header;
#else
   getValue(header, "comment",ts);
   rHeader.usrhdr=ts;
#endif

   return true; //success
}
