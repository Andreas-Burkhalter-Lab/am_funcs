function owner=id2owner(id)

   pos=find(id=='_');
   if(isempty(pos))
      error('invalid id');
   else
      owner=id(1:pos-1);
   end
