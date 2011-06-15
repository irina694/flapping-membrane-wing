function [rules,direct,thickness] = setRules

    alphabet = 'ABx';
    
    %double spiral
    %rules{1} = 'xB[+A]x';
    %rules{2} = 'A';

    
    rules{1} = 'B[-A]x[+A]B';
    direct{1} = '+00+0+00+0+';
    thickness{1} = [2 0 0 1 0 1 0 0 1 0 0.5];

    rules{2} = 'A';
    direct{2} = '+';
    thickness{2} = '1'; 
    
%     rules{3} = 'x';
%     direct{3}= '+';
    
%     rules{4} = 'x[+A]BA';
%     direct{4}= '+00+0++';
%     
%     rules{5} = 'D';
%     direct{5}= '-';
%     
%     rules{6} = 'B[+x]B';
%     direct{6}= '+00+0+';
    
    rules{end+1} = 'x';
    Nrules = length(rules);

end