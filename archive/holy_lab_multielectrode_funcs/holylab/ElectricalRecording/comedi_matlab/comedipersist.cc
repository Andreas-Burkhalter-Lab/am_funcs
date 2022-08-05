/* 
    comedipersist: a layer on top of comedilib that facilitates
    binding for languages like MATLAB

    Copyright (c) 2002 Tim Holy <holy@pcg.wustl.edu>

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation, version 2.1
    of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
    USA.

*/
/*
   This uses a little bit of C++, so I can use the map class, but
   the interface to the rest of the world looks like C.
*/
#include <map>
#include <string>
#include <unistd.h>
#include <comedilib.h>
#include "comedipersist.h"

using namespace std;

typedef string      keytype;
typedef comedi_t*   valtype;

map<keytype,valtype> comedi_dev_map;

/* 
   These functions implement a shared library
   allowing access to the devices via their name rather than
   a comedi_t pointer.
   This allows one to use the comedilib functions even if you don't
   or can't keep track of comedi_t pointers.
   dev = open_comedi_device("/dev/comedi0")
       Opens the device and returns the comedi_t pointer
   comedi_device_isopen("/dev/comedi0")
       Tests to see if this device is already open
   dev = get_comedi_device("/dev/comedi0")
       Returns the pointer for an already-open device
   dev = close_comedi_device("/dev/comedi0")
       Closes the comedi device
*/
   

comedi_t* open_comedi_device(const char *cstrfilename)
{
  map<keytype,valtype>::iterator comedi_devp;
  valtype dev;
  const string filename(cstrfilename);

  // Is this device supposedly open?
  comedi_devp = comedi_dev_map.find(filename);
  if (comedi_devp != comedi_dev_map.end()) {
    fprintf(stderr,"device %s is already registered as being open.\n",
	    cstrfilename);
    exit(1);
  }
  dev = comedi_open(cstrfilename);
  if (dev == NULL) {
    comedi_perror(cstrfilename);
    fprintf(stderr,"Error opening device %s\n",cstrfilename);
    exit(1);
  }
  comedi_dev_map[filename] = dev;
  return comedi_dev_map[filename];
}

bool comedi_device_isopen(const char *cstrfilename)
{
  map<keytype,valtype>::iterator comedi_devp;
  const string filename(cstrfilename);

  // Is this device already open?
  comedi_devp = comedi_dev_map.find(filename);
  return (comedi_devp != comedi_dev_map.end());
}

comedi_t* get_comedi_device(const char *cstrfilename)
{
  map<keytype,valtype>::iterator comedi_devp;
  const string filename(cstrfilename);

  // Is this device already open?
  comedi_devp = comedi_dev_map.find(filename);
  if (comedi_devp == comedi_dev_map.end()) {
    fprintf(stderr,"device %s is not registered as being open.\n",
	    cstrfilename);
    exit(1);
  }
  return comedi_devp->second;
}

void close_comedi_device(const char *cstrfilename)
{
  map<keytype,valtype>::iterator comedi_devp;
  const string filename(cstrfilename);

  // Is this device really open?
  comedi_devp = comedi_dev_map.find(filename);
  if (comedi_devp == comedi_dev_map.end()) {
    fprintf(stderr,"device %s is not registered as being open.\n",
	    cstrfilename);
    exit(1);
  }
  comedi_close(comedi_devp->second);
  comedi_dev_map.erase(comedi_devp);
}

