#ifndef parse_merec_header_h
#define parse_merec_header_h

#include <string>
#include <vector>

using namespace std;

struct TMerecHeader {
   int headersize;      //0
   long nscans;         //1
   int numch;           //2
   vector<int> channels;//3
   int scanrate;        //4
   double scalemult;    //5
   double scaleoff;     //6
   double voltageMin;   //7
   double voltageMax;   //8
   string date;         //9
   string time;         //10
   string usrhdr;       //11

};

static const char * fields[]={
   "headersize",
   "nscans",
   "numch",
   "channels",
   "scanrate",
   "scalemult",
   "scaleoff",
   "voltageMin",
   "voltageMax",
   "date",
   "time",
   "usrhdr"

};


bool parseMerecHeader(int fd, TMerecHeader& rHeader);


#endif //parse_merec_header_h

/*
     headersize: 263 
         nscans: 20000
          numch: 64
       channels: [1x64 double]
       scanrate: 10000
      scalemult: 0.0024
       scaleoff: -5
     voltageMin: -5
     voltageMax: 5
           date: '25-Feb-2003'
           time: ' 3:35 PM'
         usrhdr: [1x66 char]
*/

