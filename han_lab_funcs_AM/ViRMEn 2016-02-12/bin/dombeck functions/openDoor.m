function vr = doorControl(vr)

if(vr.hasDoor)            %%MASTER LOOP FOR DOOR BAR JPB FOR MARK 
        if vr.position(2) >50 %%BAR and JPB added
            vr.worlds{oldWorld}.surface.vertices(1,vr.doorVertices(1):vr.doorVertices(2)) = vr.worlds{oldWorld}.surface.vertices(1,vr.doorVertices(1):vr.doorVertices(2))+vr.doorSpeed;
        else
              vr.worlds{oldWorld}.surface.vertices(1,vr.doorVertices(1):vr.doorVertices(2))  = vr.doorigin;
        end
end
    