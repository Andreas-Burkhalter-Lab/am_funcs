function [xtemp,ytemp]=find_wall(angle,loc,t)

if sin(angle)<0
    %Hits either A or B
    if cos(angle)>0
        rb=(180-loc(t,2))/cos(angle);
        ra=-loc(t,1)/sin(angle);
        if ra<rb %hits A
            xtemp=0;
            ytemp=ra*cos(angle)+loc(t,2);
        else %hits B
            ytemp=180;
            xtemp=rb*sin(angle)+loc(t,1);
        end
        
    else %Hits either A or D
        ra=-loc(t,1)/sin(angle);
        rd=-loc(t,2)/cos(angle);
        if ra<rd %hits A
            xtemp=0;
            ytemp=ra*cos(angle)+loc(t,2);
        else %hits D
            ytemp=0;
            xtemp=rd*sin(angle)+loc(t,1);
        end
    end
    
else %looking right
    if cos(angle)>0
        rb=(180-loc(t,2))/cos(angle);
        rc=(14-loc(t,1))/sin(angle);
        if rb<rc %hits B
            ytemp=180;
            xtemp=rb*sin(angle)+loc(t,1);            
        else %hits C
            xtemp=14;
            ytemp=rc*cos(angle)+loc(t,2);
        end
    else
        rc=(14-loc(t,1))/sin(angle);
        rd=-loc(t,2)/cos(angle);
        if rc<rd %hits C
            xtemp=14;
            ytemp=rc*cos(angle)+loc(t,2);
        else %hits D
            ytemp=0;
            xtemp=rd*sin(angle)+loc(t,1);
        end
    end
end