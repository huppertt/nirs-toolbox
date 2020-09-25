function [J,meas] = jacobian( obj ,type)
%JACOBIAN Summary of this function goes here
%   Detailed explanation goes here



if nargin < 2
    isSpectral = false;
elseif strcmpi( type,'standard' )
    isSpectral = false;
elseif strcmpi( type,'spectral' );
    isSpectral = true;
else
    error('Jacobian can either be ''standard'' or ''spectral''.')
end


% TODO change fft to properly due a dtft
% for now, choose proper timeStep
fbin = obj.Fm*1e6 * obj.timeStep * obj.nTimeGates;
assert( fbin - fix(fbin) < 1e-9 )
fbin = fix(fbin);

obj.saveFluence();

meas = load([obj.directory filesep 'measurement.mat']);
meas = meas.meas;

types = unique( obj.probe.link.type );
assert( isnumeric( types ) );

for iLink = 1:size(obj.probe.link,1)
    iSrc = obj.probe.link.source(iLink);
    iDet = obj.probe.link.detector(iLink);
    type=obj.probe.link.type(iLink);
    src = load([obj.directory filesep 'src' num2str(iSrc) '_' num2str(type) 'nm.mat'],'fluence');
    det = load([obj.directory filesep 'det' num2str(iDet) '_' num2str(type) 'nm.mat'],'fluence');
    
    Jmua(iLink,:) = src.fluence .* det.fluence;
    
end




[~,~,iType] = unique(obj.probe.link.type );
if ~isSpectral
    J.mua = Jmua;
else
    % convert jacobian to conc
    ext = nirs.media.getspectra( types );
    
    ehbo = ext(iType,1);
    ehbr = ext(iType,2);
    
    J.hbo = bsxfun(@times,ehbo,Jmua);
    J.hbr = bsxfun(@times,ehbr,Jmua);
    
end

if obj.cleanup == 1
    rmdir(obj.directory,'s');
end
end


