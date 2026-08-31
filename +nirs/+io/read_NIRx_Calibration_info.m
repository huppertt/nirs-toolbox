function [info,DetectorDarkNoise]=read_NIRx_Calibration_info(filename)

calibinfo = nirs.util.json2table(filename);

info.amplitudes=calibinfo.amplitudes{1};
info.coefficient_of_variation = calibinfo.coefficient_of_variation{1};
info.signal_quality_level = vertcat(calibinfo.signal_quality{1}.level);
info.signal_quality_issaturated=vertcat(calibinfo.signal_quality{1}.saturated);
info.signal_quality_xtalk = vertcat(calibinfo.signal_quality{1}.xtalk);

chan=fields(calibinfo.dark_noise{1});
DetectorDarkNoise = zeros(length(chan),4);
for i=0:length(chan)-1
    x=['x' num2str(i)];
    dn = calibinfo.dark_noise{1}.(x);
    DetectorDarkNoise(i+1,:)=dn;

end

info = struct2table(info)
