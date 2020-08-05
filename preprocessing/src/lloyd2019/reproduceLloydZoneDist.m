format compact; format long;
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);
addpath(genpath('../../../../../libmatlab/')) % lib codes

Radio = readLloydRadioData(kfacType);
outPath = "../../out/lloyd2019/";
fileType = ".pdf";
hold on; box on;
for irtype = 1:Radio.Type.count
    rtype=Radio.Type.Name{irtype};
    cdfplot( Radio.(rtype).Zone.Val )
end
hold off;
legend({"Dark","Bright"}...
        ,'Location','southwest');
xlabel("Z + 1");
ylabel("CDF");
export_fig(outPath + "ZoneCDF" + fileType,'-m4 -transparent')
