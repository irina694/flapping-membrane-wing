%function [word,dir,thick,markposition,markletter,markthick] = preparesRules(rule,dir0,thick0)
function Prules = preparesRules(Prules)

    Nrules = length(Prules.letters);
    
    for iR=1:Nrules
        letters = Prules.letters{iR};
        direction = Prules.direction{iR};
        thickRule = Prules.thickRule{iR};
        elastRule = Prules.elastRule{iR};
        pressRule = Prules.pressRule{iR};
        
        [directRule.letters,directRule.thick,directRule.elast,directRule.press,directRule.dir,...
         directRule.markindex,directRule.markposition,directRule.markletter,...
         directRule.markthick] = removesBrackets(letters,thickRule,elastRule,pressRule,direction);
        
        [letters,direction,thickRule]=invertrule(letters,direction,thickRule,elastRule,pressRule);
        [inverseRule.letters,inverseRule.thick,inverseRule.elast,inverseRule.press,inverseRule.dir,...
         inverseRule.markindex,inverseRule.markposition,inverseRule.markletter,...
         inverseRule.markthick] = removesBrackets(letters,thickRule,elastRule,pressRule,direction);

        Prules.directRule{iR} = directRule;
        Prules.inverseRule{iR} = inverseRule;
    end

end



function [letters,thick,elast,press,dir,markindex,markposition,markletter,markthick] = removesBrackets(letters0,thick0,elast0,press0,dir0)

        letters0 = letters0';
        thick0 = thick0';
        elast0 = elast0';
        press0 = press0';
        dir0 = dir0';
        
        N = length(letters0);
        inbracket = false;        
        letters = [];
        dir = [];
        thick = [];
        elast = [];
        press = [];
        position = 0;
        markposition = [];
        markletter = [];
        markthick = [];
        markindex = [];
        mkindex = 0;
        for i=1:N
            if(letters0(i)=='[')
                mkindex = mkindex +1;
                markindex(mkindex) = mkindex;
                markposition = [markposition,position];
                markletter = [markletter,100*double(letters0(i+1))+double(letters0(i+2))];
                markthick = [markthick,thick0(i+2)];
                inbracket = true;
            elseif(letters0(i)==']')
                inbracket = false;
            end

            if(~inbracket && letters0(i)~=']')
                position = position + 1;
                letters = [letters,letters0(i)];
                thick = [thick,thick0(i)]; 
                elast = [elast,elast0(i)]; 
                press = [press,press0(i)]; 
                dir = [dir,dir0(i)];
            end
        end
                            
end





function [invrule,invdir,invthick,invelast,invpress]=invertrule(rule,dir,thick,elast,press)

    nn=length(rule);
    
    invrule=[];
    invdir=[];
    invthick=[];
    invelast=[];
    invpress=[];
    i=nn;
    while (i>0)
        if(rule(i)==']')
            if(rule(i-2)=='+')    
                invrule=[invrule,'[-',rule(i-1),']'];
            else
                invrule=[invrule,'[+',rule(i-1),']'];
            end
            
            if(dir(i-1)=='+')    
                invdir=[invdir,'00-0'];
            else
                invdir=[invdir,'00+0'];
            end

            invthick=[invthick,[0 0 thick(i-1) 0]];
            invelast=[invelast,[0 0 elast(i-1) 0]];
            invpress=[invpress,[0 0 press(i-1) 0]];
            i=i-4;
        else
            invrule=[invrule,rule(i)];
            if(dir(i)=='+')    
                invdir=[invdir,'-'];
            else
                invdir=[invdir,'+'];
            end
            invthick=[invthick,thick(i)];
            invelast=[invelast,elast(i)];
            invpress=[invpress,press(i)];
            i=i-1;
        end
    end

end