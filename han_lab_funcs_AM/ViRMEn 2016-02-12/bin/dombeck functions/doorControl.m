function vr = doorControl(vr,oldWorld)

vr.doorSpeed = 1; %in game units per second
%{
  vr.doorVertice -- List of all vertices for the door
  vr.doorigin -- origin of all the vertices the door
    
%}
if(vr.hasDoor)            %%MASTER LOOP FOR DOOR BAR JPB FOR MARK 
        if vr.position(2) >50 %%BAR and JPB added
            vr.worlds{oldWorld}.surface.vertices(1,vr.doorVertices(1):vr.doorVertices(2)) = vr.worlds{oldWorld}.surface.vertices(1,vr.doorVertices(1):vr.doorVertices(2))+vr.doorSpeed*vr.dt;
            vr.worlds{oldWorld}.edges.endpoints(vr.doorEdgeIndexes,[1 3]) =  vr.worlds{oldWorld}.edges.endpoints(vr.doorEdgeIndexes,[1 3]) +  vr.doorSpeed*vr.dt;

        else
            vr.worlds{oldWorld}.surface.vertices(:,vr.doorVertices(1):vr.doorVertices(2))  = vr.doorigin;
            vr.worlds{oldWorld}.edges.endpoints(vr.doorEdgeIndexes,:) = vr.doorEdgeOrigin;
        end
end
    