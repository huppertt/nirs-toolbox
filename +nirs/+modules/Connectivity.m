classdef Connectivity < nirs.modules.AbstractModule
    %% CONNECTIVITY - Computes all-to-all connectivity model.
    % Outputs nirs.core.ConnStats object

    properties
        corrfcn;  % function to use to compute correlation (see +nirs/+sFC for options)
        divide_events;  % if true will parse into multiple conditions
        min_event_duration;  % minimum duration of events
        ignore;  % time at transitions (on/off) to ignore (only valid if dividing events)
        AddShortSepRegressors = false;
        verbose=true;
    end
    methods
        function obj = Connectivity( prevJob )
            obj.name = 'Connectivity';
            obj.corrfcn = @(data)nirs.sFC.ar_corr(data,'4xFs',true);  %default to use AR-whitened robust correlation
            obj.divide_events=false;
            obj.min_event_duration=30;
            obj.ignore=10;
            if nargin > 0
                obj.prevJob = prevJob;
            end
            obj.citation{1}='Santosa, H., Aarabi, A., Perlman, S. B., & Huppert, T. J. (2017). Characterization and correction of the false-discovery rates in resting state connectivity using functional near-infrared spectroscopy. Journal of Biomedical Optics, 22(5), 055002-055002.';
            obj.citation{2}='Lanka, Pradyumna, Heather Bortfeld, and Theodore J. Huppert. "Correction of global physiology in resting-state functional near-infrared spectroscopy." Neurophotonics 9.3 (2022): 035003-035003.';
        end

        function connStats = runThis( obj, data )
            if(obj.AddShortSepRegressors)
                if(obj.verbose)
                    disp('removing short-seperation noise');
                end
                job=advanced.nirs.modules.ShortDistanceFilter;
                data=job.run(data);
            end


            nFiles = numel(data);
            connStats = repmat(nirs.core.sFCStats(), 1, nFiles);

            % Cache loop-invariant values for parfor
            corrfcn_  = obj.corrfcn;
            divide_   = obj.divide_events;
            ignore_   = obj.ignore;
            minDur_   = obj.min_event_duration;
            verbose_  = obj.verbose;

            parfor i = 1:nFiles

                cs = nirs.core.sFCStats();
                cs.type = corrfcn_;
                cs.description = ['Connectivity model of ' data(i).description];
                cs.probe = data(i).probe;
                cs.demographics = data(i).demographics;

                if divide_
                    stim = data(i).stimulus;
                    cnt = 1;
                    for idx = 1:length(stim.keys)
                        s = stim(stim.keys{idx});
                        lst = find(s.dur - 2*ignore_ > minDur_);
                        if ~isempty(lst)
                            s.onset = s.onset(lst);
                            s.dur = s.dur(lst);

                            nCh = size(data(i).data, 2);
                            r = zeros(nCh, nCh, length(s.onset));
                            dfe = zeros(length(s.onset), 1);

                            itime = data(i).time;
                            idata = data(i).data;
                            for j = 1:length(s.onset)
                                tmp = nirs.core.Data;
                                lstpts = find(itime > s.onset(j)+ignore_ & ...
                                    itime < s.onset(j)+s.dur(j)-ignore_);
                                tmp.data = idata(lstpts,:);
                                tmp.time = itime(lstpts);
                                [r(:,:,j), ~, dfe(j)] = corrfcn_(tmp);
                            end

                            cs.dfe(cnt) = sum(dfe);
                            cs.R(:,:,cnt) = tanh(mean(atanh(r), 3));
                            cs.conditions{cnt} = stim.keys{idx};
                            cnt = cnt + 1;
                        end
                    end

                else
                    tmp = data(i);
                    lst = find(tmp.time < ignore_ | tmp.time > tmp.time(end) - ignore_);
                    tmp.data(lst,:) = [];
                    tmp.time(lst) = [];
                    [r, ~, dfe] = corrfcn_(tmp);

                    cs.dfe = dfe;
                    cs.R = r;
                    cs.conditions = cellstr('rest');
                end

                connStats(i) = cs;

            end
        end

    end

end
