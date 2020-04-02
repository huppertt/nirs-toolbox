function ROItbl = threshold2ROI( GroupStats , threshold , type )
% Create an ROI entry based on what survives at the group level.
% Can be used with nirs.modules.ApplyROI() or advanced.nirs.modules.ApplyROIMFX()
%
% Usage: ROItbl = advanced.nirs.util.threshold2ROI( GroupStats , threshold , type );
%
% Example: ROItbl = advanced.nirs.util.threshold2ROI( GroupStats.ttest('2back-1back') , 'q<.05' , 'hbo' );
if ~exist('GroupStats','var') || isempty(GroupStats)
    error('Please specify GroupStats variable');
end
if ~exist('threshold','var') || isempty(threshold)
    error('Please specify significance threshold');
end
if ~exist('type','var') || isempty(type)
    error('Please specify type (''hbo'' or ''hbr'')');
end
if length(GroupStats.conditions)>1
    error('Stats must have only 1 condition. Please use the ttest() method to submit the contrast of interest.');
end

threshold = strrep(threshold,' ','');
tparts = strsplit(threshold,'<');
thresh_field = tparts{1};
thresh_value = str2double(tparts{2});
is_good = (GroupStats.(thresh_field) < thresh_value) & strcmpi(GroupStats.variables.type,type);

if ~any(is_good)
    warning('Nothing survived. Exiting.');
    ROItbl=table();
    return;
end

vars = GroupStats.variables(is_good,:);
source = {vars.source'};
detector = {vars.detector'};
name = {sprintf('%s (%s)',GroupStats.conditions{1},threshold)};

ROItbl = table(source,detector,name);

end