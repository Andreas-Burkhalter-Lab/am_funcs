//header shared by all file dealing with LFS
//@note: 1. this file MUST be the first file to be included.
//       2. if LFS_MACROS_ONLY is defined, include files for file access 
//          are not included in this file.
//       3. In kylix, you should avoid to use lseek/lseek64, use lseek_int64_g() 
//          or lseek64_g() instead. lseek_int64_g() is prefered in Kylix.
//       4. If you use lseek64_g() and/or lseek_int64_g() in kylix, you need add  
//          lseek64_g.o (compiled by gcc from lseek64_g.c) into project manager.
//       5. when call lseek...() functions and argment offset > file's size,
//          return value may ==offset and !=-1: for file opened w/ readonly mode,
//          after such call, the actual file position != offset; for file opened
//          w/ writeonly mode, after such call, file is enlarged. @see: man lseek.
//       6. you can use lseek...(fd, 0, SEEK_CUR) to get file position (may beyond EOF).
//       

#ifndef LFS_COMMON_H
#define LFS_COMMON_H

//1. for using ftello() and fseeko() function:
#define _LARGEFILE_SOURCE
//2. for using ...64() functions and off64_t:
#define _LARGEFILE64_SOURCE
//3. for mapping old 32-bit interface (open(), read(), fopen(), ... and off_t)
//   to ...64() functions and off64_t. I.e. old interface is just aliases of
//   new interface:
#define _FILE_OFFSET_BITS 64

#if !defined(LFS_MACROS_ONLY)

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#if defined(__BORLANDC__) || defined(__BCPLUSPLUS__)
#include "lseek64_g.h"
#endif

#endif // #if not defined LFS_MACROS_ONLY, include necessary files for file accesses.


#endif //LFS_COMMON_H
