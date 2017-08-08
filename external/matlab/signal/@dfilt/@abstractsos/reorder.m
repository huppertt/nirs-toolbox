function varargout = reorder(this, numorder, denorder, svorder)
%REORDER   Reorder the sections.
% Help in dfilt/reorder.

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.

error(nargchk(2,4,nargin,'struct'));

if nargout>0,
    Hd = copy(this);
else
    Hd = this;
end

if ischar(numorder),
    switch numorder,
        case {'lowpass','highpass','bandpass','bandstop'},
            fd = feval(['fdesign.',numorder]);
            sosreorder(fd,Hd);
        case 'auto',
            fd = getfdesign(Hd);
            % Reorder if fdesign found, otherwise no-op
            if ~isempty(fd),
                sosreorder(fd,Hd);
            end
        case 'none',
            % No-op
        otherwise,
            % Make sure we use the reference coefficients.
            newsv  = get(Hd, 'refScaleValues');
            oldsos = get(Hd, 'refsosMatrix');
            
            % Make sure leading coeffs are equal to one
            newsv(1:size(oldsos,1)) = newsv(1:size(oldsos,1)).*oldsos(:,1);
            oldsos(:,1:3) = oldsos(:,1:3)./repmat(oldsos(:,1),1,3);

            % Find all the poles and zeros.
            z = []; p = [];
            for indx = 1:size(oldsos, 1)
                [z1, p1, k] = sos2zp(oldsos(indx,:));
                z = [z; z1];
                p = [p; p1];
                newsv(indx) = newsv(indx)*k;
            end

            newsos = zp2sos(z, p, 1, numorder,'none',true);

            if nargin < 3,

                % The default scale values movement is dependent on whether the
                % structure has secondary scaling points that shift the numerator.
                if shiftsecondary(Hd)
                    denorder = 'poles';
                else
                    denorder = 'zeros';
                end
            end

            rindx = [];
            if strcmpi(denorder, 'poles'),
                sindx = [4 5 6];
            else,
                sindx = [1 2 3];
            end

            % Find the reorder index used.
            for indx = 1:size(newsos, 1)
                newindx = find(all([abs(newsos(indx,sindx(1))-oldsos(:,sindx(1))) < sqrt(eps), ...
                    abs(newsos(indx,sindx(2))-oldsos(:,sindx(2))) < sqrt(eps), ...
                    abs(newsos(indx,sindx(3))-oldsos(:,sindx(3))) < sqrt(eps)], 2));
                if length(newindx) == 1
                    rindx = [rindx; newindx];
                else
                    for jndx = 1:length(newindx)
                        if isempty(find(rindx == newindx(jndx))),
                            rindx = [rindx; newindx(jndx)];
                            break;
                        end
                    end
                end
            end
            newsv = [newsv(rindx);newsv(end)];
            
            Hd.ScaleValues = newsv;
            Hd.sosMatrix   = newsos;
    end
else
    
    if islogical(numorder), numorder = find(numorder); end
    
    if nargin < 3, denorder = numorder; end
    if nargin < 4, svorder  = [numorder max(numorder)+1]; end
    
    % Change the logicals into indices to make the error checking simpler.
    if islogical(denorder), denorder = find(denorder); end
    if islogical(svorder),  svorder  = find(svorder);  end
    
    if max(numorder) > nsections(Hd),
        error(message('signal:dfilt:abstractsos:reorder:exceedsdims'));
    end
    
    if any(numorder < 1),
        error(message('signal:dfilt:abstractsos:reorder:badsubscript'));
    end
    
    if length(numorder) ~= length(denorder) || length(svorder)-length(numorder) > 1,
        error(message('signal:dfilt:abstractsos:reorder:arrayMismatch'));
    end
    
    w    = warning('off', 'signal:dfilt:basefilter:warnifreset:PropWillBeReset');
    
    % Get the scale values that the user added.
    privsv = Hd.refScaleValues(1:end-Hd.NumAddedSV);
    
    % Reorder all the scale values, including our added ones.
    newsv  = Hd.refScaleValues(svorder);
    
    % Remove any scale values that occur after those that were set by the user.
    if length(privsv) <= max(svorder),
        newsv = newsv(1:max(find(svorder <= length(privsv))));
    end
    
    newsos = [Hd.refsosMatrix(numorder, 1:3) Hd.refsosMatrix(denorder, 4:6)];
    
    warning(w);
    
    Hd.ScaleValues = newsv;
    Hd.sosMatrix   = newsos;
end

if nargout>0,
    varargout{1} = Hd;
end


% [EOF]
