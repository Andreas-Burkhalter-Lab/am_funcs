function headerfile=imdatafilename2headerfilename(datafile)
   headerfile=[datafile '.txt'];
   if(fileexist(headerfile))
      return
   end
   headerfile=replace_extension(datafile, '.txt');
   
