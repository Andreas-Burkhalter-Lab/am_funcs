merecfiles=UIGetFiles('*.merec', ...
                    'please select merec data files' ...
                    );

for idx=1:length(merecfiles)
   suc=merec2berrydata(merecfiles{idx}, replace_extension(merecfiles{idx}, '_berry.merec'));
   if(~suc)
      error(['cannot process file ' merecfiles{idx}]);
   end
end
