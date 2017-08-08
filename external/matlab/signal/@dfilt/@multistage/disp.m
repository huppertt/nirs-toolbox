function disp(this)
%DISP Display method of discrete-time filter.
%
%   See also DFILT.   
  
%   Author: V. Pellissier
%   Copyright 1988-2005 The MathWorks, Inc.

if length(this) > 1
    vectordisp(this);
    return;
end

disp(['     FilterStructure: ' this.FilterStructure]);

lvlTwooffset = ['        '];
lvlThreeoffset = [lvlTwooffset lvlTwooffset];

for lvlOne=1:length(this.Stage)
    lvlOnesect = this.Stage(lvlOne);
    filtstruct = lvlOnesect.FilterStructure;
    disp(['            Stage(',num2str(lvlOne),'): ' filtstruct]);
    if ~isempty(strmatch(filtstruct,{'Cascade', 'Parallel'})),
        for lvlTwo=1:length(lvlOnesect.Stage)
            lvlTwosect = lvlOnesect.Stage(lvlTwo);
            filtstruct = lvlTwosect.FilterStructure;
            disp([lvlTwooffset, '              Stage(',num2str(lvlTwo),'): ' filtstruct]);
            if ~isempty(strmatch(filtstruct,{'Cascade', 'Parallel'})),
                for lvlThree=1:length(lvlTwosect.Stage),
                    filtstruct = lvlTwosect.Stage(lvlThree).FilterStructure;
                    disp([lvlThreeoffset, '            Stage(',num2str(lvlThree),'): ' filtstruct]);
                end
            end
        end
    end
end
siguddutils('dispstr', this, {'PersistentMemory'}, 20);

% [EOF]
