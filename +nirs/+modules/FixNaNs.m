classdef FixNaNs < nirs.modules.AbstractModule
%% FixNaNs - Attempts to fix NaN values by interpolation.
%
% Options:
%    ifFailReplaceWith - value to replace NaNs with if interpolation fails
%
% Notes:
%     1 is a good value for raw data; 0 otherwise

    properties
        ifFailReplaceWith = 1; % value to replace NaNs with if interpolation fails
<<<<<<< HEAD
=======
        type='linear';
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f
    end
    
    methods

        function obj = FixNaNs( prevJob )
           obj.name = 'Fix NaNs';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)

                d = data(i).data;
                t = data(i).time;
                
                lst = isnan(d);
                
                if any(lst(:))
                    try
                        for j = 1:size(d,2)
                            if(all(lst(:,j)))
                                 d(lst(:,j),j) = obj.ifFailReplaceWith;
                            elseif any(lst(:,j))
                                % interpolation
                                l = lst(:,j);
<<<<<<< HEAD
                                d(l,j) = interp1(t(~l), d(~l,j), t(l),'linear','extrap');
=======
                                d(l,j) = interp1(t(~l), d(~l,j), t(l),obj.type);
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f
                                data(i).data = d;
                            end
                        end
                        
                    catch
                        % just replace with white noise
                        d(lst) = obj.ifFailReplaceWith;
                    end

                end
<<<<<<< HEAD
=======
                
                % repeat to get the edges using nearest
                lst = isnan(d);
                
                if any(lst(:))
                    try
                        for j = 1:size(d,2)
                            if(all(lst(:,j)))
                                 d(lst(:,j),j) = obj.ifFailReplaceWith;
                            elseif any(lst(:,j))
                                % interpolation
                                l = lst(:,j);
                                d(l,j) = interp1(t(~l), d(~l,j), t(l),'nearest','extrap');
                                data(i).data = d;
                            end
                        end
                        
                    catch
                        % just replace with white noise
                        d(lst) = obj.ifFailReplaceWith;
                    end

                end
                
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f
                data(i).data=d;
            end
        end
    end
    
end

