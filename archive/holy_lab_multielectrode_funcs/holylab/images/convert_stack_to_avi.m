function movie = convert_stack_to_avi(filename, instack, map, fps, quality)
% convert_stack_to_movie : converts a basic stack to an avi movie
%    filename is a string.  If .avi is not added, it will be added for you
%    instack is a 3-d image stack (x, y, z)
%    map must be a matlab colormap (predefined or self defined)
%    fps is the frames per second
%    quality is a 1-100 integer, 1 is worst, 100 is best, 100 = no compression

scaled_min = min(min(min(instack,[],3),[],2),[],1);
scaled_max = max(max(max(instack,[],3),[],2),[],1);

for index = 1: size(instack,3)
   mov(index) = im2frame(uint8(instack(:,:,index)-scaled_min)/(scaled_max-scaled_min),map);
end

if isempty(regexp(filename, '.avi'))
    filename = [filename '.avi'];  % add '.avi' if not already put in
end

if quality == 100
    movie2avi(mov, filename,...
                          'colormap', map,...
                          'fps', fps);
else
    movie2avi(mov, filename,...
                          'colormap', map,...
                          'fps', fps,...
                          'compression', 'indeo5',...
                          'quality', quality);
end

movie = mov;

end
