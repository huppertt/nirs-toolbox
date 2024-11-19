function tbl=add_conditional_table_items(tbl,test)
% This function is used in the nirs.modules.MixedEffects code to all you to
% add additional demographics to a table based on a test condition
%
%The test should be a string (or cell) with the format:
% Example 1:  Add a new column named "like" with the value 'like' if cond ismember('A') and 'dislike' otherwise
%   test=['cond ismember(''A'') -> like = ''like'' else ''dislike'''];
% Example 2 using "and" in the condition test  (also allows "or")
%  test=['cond ismember(''A'') and subject ismember({''S1'',''S2'',''S3''}) -> like = true else false'];


if(~iscell(test))
    test={test};
end
for ii=1:length(test)
    lhs = test{ii}(1:strfind(test{ii},'->')-1);
    rhs = test{ii}(strfind(test{ii},'->')+2:end);
    splt=[2; strfind(lhs,'and'); strfind(lhs,'or'); length(lhs)];
    teststatement='';

    for idx=2:length(splt)
        lhss{idx-1}=lhs(splt(idx-1)-1:splt(idx)-1);
        if(contains(lhss{idx-1},'and'));
            opp=' & ';
        elseif(contains(lhss{idx-1},'or'));
            opp=' | ';
        else;
            opp='';
        end
        lhss{idx-1}=strrep(lhss{idx-1},'or','');
        lhss{idx-1}=strrep(lhss{idx-1},'and','');
        lhs_var{idx-1} = strtrim(lhss{idx-1}(1:strfind(lhss{idx-1},'ismember')-1));
        lhs_ismember{idx-1}=strtrim(lhss{idx-1}(strfind(lhss{idx-1},'ismember'):end));
        teststatement = [ teststatement opp strrep(lhs_ismember{idx-1},'ismember(',['ismember(tbl.' lhs_var{idx-1} ','] )];
    end
    setvar = strtrim(rhs(1:strfind(rhs,'=')-1));
    truevar = strtrim(rhs(strfind(rhs,'=')+1:strfind(rhs,'else')-1));
    falsevar = strtrim(rhs(strfind(rhs,'else')+4:end));
    bool = eval(teststatement);
    for idx=1:length(bool)
        if(bool(idx))
            var{idx,1}=eval(truevar);
        else
            var{idx,1}=eval(falsevar);
        end
    end
    if(any(contains(tbl.Properties.VariableNames,setvar)))
        tbl.(setvar)=[];
    end
    tbl=[tbl table(var,'variableNames',{setvar})];
    disp(['Variable ' setvar ' added to demographics table']);
end