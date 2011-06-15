function strOut = createStr(axiom,productions,nDev);

global nChromo nAlpha nSym;

strIn = axiom;
for i=1:nDev
    strOut = '';
    strLength = length(strIn);
    for j=1:strLength
        letter = strIn(j);
        if(j==1)
            letterLeft = strIn(end);
        else
            letterLeft = strIn(j-1);
        end;
        
        if(j==strLength)
            letterRight = strIn(1);
        else
            letterRight = strIn(j+1);
        end;
        
        numLetter = abs(double(letter)-64);
        if(numLetter < nAlpha+1)
            switch (letterLeft)
                case '['
                    numLetterLeft = nAlpha+1;
                case ']'
                    numLetterLeft = nAlpha+2;
                case '+'
                    numLetterLeft = nAlpha+3;
                case '-'
                    numLetterLeft = nAlpha+4;
                case '*'
                    numLetterLeft = nAlpha+5;
                case '/'
                    numLetterLeft = nAlpha+6;
                case '#'
                    numLetterLeft = nAlpha+7;
                case '|'
                    numLetterLeft = nAlpha+8;
                otherwise
                    numLetterLeft = abs(double(letterLeft)-64);
            end;
        
            switch (letterRight)
                case '['
                    numLetterRight = nAlpha+1;
                case ']'
                    numLetterRight = nAlpha+2;
                case '+'
                    numLetterRight = nAlpha+3;
                case '-'
                    numLetterRight = nAlpha+4;
                case '*'
                    numLetterRight = nAlpha+5;
                case '/'
                    numLetterRight = nAlpha+6;
                case '#'
                    numLetterRight = nAlpha+7;
                case '|'
                    numLetterRight = nAlpha+8;
                otherwise
                    numLetterRight = abs(double(letterRight)-64);
            end;
            strOut=[strOut,char(productions(numLetterLeft,numLetter,numLetterRight))];
        else
            strOut=[strOut,strIn(j)];
        end;
    end;
    strIn = strOut;
end;
