function test_keypress(sender, event)
% USAGE:
%  figure; set(gcf, 'keyPressFcn', @test_keypress)
   if(event.Key=='a') keyboard; end
   
   if(isempty(event.Character)) 
      disp('no ch');
   else
      disp(event.Character);
   end
   
   if(isempty(event.Modifier)) 
      disp('no mod');
   else
      disp(event.Modifier);
   end
   
   disp(event.Key);
   
