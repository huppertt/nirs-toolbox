function converted = convertScoreTransform(st,to,K)
% Convert score transform ST to TO. ST is either a function handle to be
% converted to a string or a string to be converted to a function handle. K
% is the number of classes in the data.

switch to
    case 'handle' % convert to handle
        if     isa(st,'function_handle')
            converted = st;
        elseif ischar(st)
            if strcmpi(st,'none')
                converted = @classreg.learning.transform.identity;
            else
                converted = str2func(['classreg.learning.transform.' st(:)']);
            end
        else
            error(message('stats:classreg:learning:internal:convertScoreTransform:BadType'));
        end
        try
            converted(zeros(1,K));
        catch me
            error(message('stats:classreg:learning:internal:convertScoreTransform:BadTransform', me.message));
        end
    case 'string' % convert to string
        converted = func2str(st);
        if strcmp(converted,'classreg.learning.transform.identity')
            converted = 'none';
        end
        idx = strfind(converted,'classreg.learning.transform.');
        if ~isempty(idx)
            converted = converted(1+length('classreg.learning.transform.'):end);
        end
end
end
