classdef TestPreprocessing < matlab.unittest.TestCase
    
    properties
        od
        od_resample
        hb
    end
    
    methods
        function obj = TestPreprocessing()
            
            %% Create probe
            srcPos = [0 0 0;
                30 30 0];
            
            detPos = [ 0 30 0;
                30 0 0];
            
            link = [1 1 690;
                1 1 830;
                1 2 690;
                1 2 830;
                2 1 690;
                2 1 830;
                2 2 690;
                2 2 830 ];
            
            link = table(link(:,1), link(:,2), link(:,3), ...
                'VariableNames', {'source', 'detector', 'type'});
            
            probe = nirs.core.Probe(srcPos, detPos, link);
            
            %% Create task design
            time = (0:2000) * (1/10);
            
            stim = Dictionary();
            
            stim('Slow') = nirs.design.StimulusEvents('Low',[25 97 149],[17 6 12],[1 1 1]);
            stim('Fast') = nirs.design.StimulusEvents('Medium',[5 80 112],[12 15 16],[1 1 1]);
            stim('Extreme') = nirs.design.StimulusEvents('High',[50 60 134],[6 10 7],[1 1 1]);
            
            beta = [20 60 100]';
            
            %% Generate data
            rng('default');
            noise = nirs.testing.simARNoise(probe,time);
            raw = nirs.testing.simData(noise,stim,beta);
            
            %% Randomize probe link sequence
            rand_inds = randperm(height(raw.probe.link));
            raw.probe.link = raw.probe.link(rand_inds,:);
            raw.data = raw.data(:,rand_inds);
            link = raw.probe.link;
            
            %% Preprocess w/ toolbox
            job = nirs.modules.OpticalDensity();
            obj.od = job.run( raw );
            
            job = nirs.modules.Resample();
            job.Fs = 4;
            obj.od_resample = job.run( obj.od );
            
            job = nirs.modules.BeerLambertLaw();
            obj.hb = job.run( obj.od_resample );
            
            %% Preprocess manually
            d = raw.data;
            
            % OD
            d = -log( bsxfun( @rdivide , abs(d) , mean(abs(d)) ) );
            obj.od(2) = raw;
            obj.od(2).data = d;
            
            % Resample
           if(isempty(ver('signal')) || verLessThan('signal','7.0'));
               t=obj.od(2).time;
               ord = floor( length(t) / 10 );
               Fc = 4/obj.od(2).Fs;
               
               b = fir1(ord,Fc);
               
               % zero phase filtering
               d = filtfilt(b,1,d);
               
               % interpolation
               N = floor((t(end)-t(1)) *4)+1;
               new_t = t(1) + (0:N-1)' /4;
               d = interp1(t,d,new_t,'linear','extrap');
               time=new_t;
           else
                [d,time]=resample(d,obj.od(2).time,4);
           end
            obj.od_resample(2) = raw;
            obj.od_resample(2).time = time;
            obj.od_resample(2).data = d;
            
            % Beer-Lambert law
            ext = nirs.media.getspectra( [690 830] );
            E = ext(:,[1 2]);
            EL = E * 30 * .1;
            iEL = pinv(EL);
            obj.hb(2) = raw;
            obj.hb(2).probe.link.type = repmat({'hbo';'hbr'},4,1);
            obj.hb(2).time = time;
            obj.hb(2).data = zeros(size(d,1),8);
            for src = 1:2
                for det = 1:2
                    ind1 = find( link.source==src & link.detector==det & link.type==690 );
                    ind2 = find( link.source==src & link.detector==det & link.type==830 );
                    tmp = d(:,[ind1 ind2]);
                    outind = ((src-1)*4+(det-1)*2+1);
                    obj.hb(2).data(:,[outind outind+1]) = (tmp * iEL') * 1e6;
                    obj.hb(2).probe.link.source([outind outind+1]) = src;
                    obj.hb(2).probe.link.detector([outind outind+1]) = det;
                end
            end
            
        end
    end
    methods (Test)
        
        function testOpticalDensity( obj )
            MAD = median( abs(obj.od(1).data(:) - obj.od(2).data(:)) );
            obj.verifyLessThan( MAD , 1e-12 );
        end
        
        function testResample( obj )
            obj.verifyEqual( obj.od_resample(1).time , obj.od_resample(2).time );
            MAD = median( abs(obj.od_resample(1).data(:) - obj.od_resample(2).data(:)) );
            obj.verifyLessThan( MAD , 1e-3);
        end
                
        function testBeerLambertLaw( obj )
            obj.verifyEqual( obj.hb(1).probe , obj.hb(2).probe );
            MAD = median( abs(obj.hb(1).data(:) - obj.hb(2).data(:)) );
            obj.verifyLessThan( MAD , 1 );
        end
    end
end