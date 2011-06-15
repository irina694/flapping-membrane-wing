function ruleStr = chromo2rule(Chromossome)

limBracket = 0.7;
limAngle = 0.5;
limLength = 0.5;
limThick = 0.5;

ruleStr = '';

chromoGene = Chromossome(1);
if(chromoGene>limBracket)
        ruleStr = [ruleStr,'['];
end;

chromoGene = Chromossome(2);
if(chromoGene > limAngle)
        ruleStr = [ruleStr,'+'];
elseif(chromoGene < -limAngle)
        ruleStr = [ruleStr,'-'];
end;

chromoGene = Chromossome(3);
if(chromoGene > limLength)
        ruleStr = [ruleStr,'*'];
elseif(chromoGene < -limLength)
        ruleStr = [ruleStr,'/'];
end;

chromoGene = Chromossome(4);
if(chromoGene > limThick)
        ruleStr = [ruleStr,'#'];
elseif(chromoGene < -limThick)
        ruleStr = [ruleStr,'|'];
end;

chromoGene = Chromossome(5);
if(chromoGene>=1)
ruleStr = [ruleStr,char(chromoGene+64)];
end;

chromoGene = Chromossome(6);
if(chromoGene>=1)
ruleStr = [ruleStr,char(chromoGene+64)];
end;

chromoGene = Chromossome(1);
if(chromoGene>limBracket)
        ruleStr = [ruleStr,']'];
end;

chromoGene = Chromossome(7);
if(chromoGene>0)
        ruleStr = [ruleStr,'['];
end;

chromoGene = Chromossome(8);
if(chromoGene > limAngle)
        ruleStr = [ruleStr,'+'];
elseif(chromoGene < -limAngle)
        ruleStr = [ruleStr,'-'];
end;

chromoGene = Chromossome(9);
if(chromoGene > limLength)
        ruleStr = [ruleStr,'*'];
elseif(chromoGene < -limLength)
        ruleStr = [ruleStr,'/'];
end;

chromoGene = Chromossome(10);
if(chromoGene > limThick)
        ruleStr = [ruleStr,'#'];
elseif(chromoGene < -limThick)
        ruleStr = [ruleStr,'|'];
end;

chromoGene = Chromossome(11);
if(chromoGene>=1)
ruleStr = [ruleStr,char(chromoGene+64)];
end;

chromoGene = Chromossome(12);
if(chromoGene>=1)
ruleStr = [ruleStr,char(chromoGene+64)];
end;

chromoGene = Chromossome(7);
if(chromoGene>0)
        ruleStr = [ruleStr,']'];
end;

