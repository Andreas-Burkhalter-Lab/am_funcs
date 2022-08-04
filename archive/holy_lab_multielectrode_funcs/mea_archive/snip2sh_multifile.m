function sh = snip2sh_multifile
% perform snipfile2sortheader on all .ssnp files in the directory
% by AM

snip_files = dir('*.ssnp');
thisdir = pwd;

sh = snipfile2sortheader({snip_files.name})
