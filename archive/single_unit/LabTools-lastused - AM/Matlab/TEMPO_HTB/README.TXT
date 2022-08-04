                         TEMPO Technical Support
                                1 Dec 99

This directory (c:\tempo\matlab) contains Matlab (v5.0) functions for
manipulating TEMPO HTB database files.  They provide the basic routines
necessary for using TEMPO's HTB databases in Matlab.

The following functions are available. Each .M file contains
documentation on how to use the function. Alternatively, you can use
Matlab's HELP facility to display information about each function.

There is some error checking in these routines but it is not exhaustive.
For example, htbWrite() does not check the validity of the HTB header
or that the epochData that you provide to it is of the correct size.
htbGetHd() does not check the validity of the database header fields.

htbtest  m      Example script that uses other htbXXX() functions
htbshow  m      Example function displays an HTB database epoch graphically
htbinit  m      Returns an empty HTB header

htbopen  m      Opens an HTB file for read-only access
htbopenw m      Opens an HTB file for read-append access

htbcount m      Returns number of databases in HTB file

htbgethd m      Returns an HTB database header as a structure

htbgetda m      Returns all epochs of a database as a matrix
htbgetep m      Returns a particular epoch of a database as a matrix

htbwrite m      Append an HTB header and data to an HTB file

htbclose m      Close HTB file

readme   txt    This file
test     htb    A test HTB file with 6 databases

While we have made every attempt to test these routines, there may be
errors in them. It is your responsibility to insure that they do what
you require.

These are provided as a courtesy by Reflective Computing to TEMPO users.
We do not provide ongoing support for Matlab or these routines.  You may
use or modify these for your use provided you have purchased at least
one TEMPO license from Reflective Computing. These routines are provided
"as is" and without warrantee of any kind.

Note: The standard C source code to the HTB.EXE utility is available.
Please contact Reflective Computing for a copy of it.

Copyright (c) 2000 Reflective Computing.  All Rights Reserved.