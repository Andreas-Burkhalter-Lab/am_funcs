%make sound clip
%clear all



cd  'd:/data/baby_song/pups_isol_tab_ko/day8/ko';



options.todisk=true;
options.window=[1 3];
whisshowplay('a147.sng',options);
%now put in break at twhis line. substitute the times you want (in a coumn
%vector)
%%now print out the matching snglot
sngparms.plot = false;sngparms.nfreq = 256; sngparms.freqrange = [25000 110000];sngparms.threshold = .23;
filenames2= dir ('*.bin');
if 1%~isequal(size(filenames),size(filenames2))
    filenames = dir( '*.bin');
    for f = 2:size(filenames,1)
        a=filenames(f).name;
        outfile=sprintf('%s.sng',a(1:end-4));
        sound2sng(filenames(f).name,sngparms,outfile);
    end
end

