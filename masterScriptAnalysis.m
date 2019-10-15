
%%
addpath(genpath('~/associativeLearningSingleUnit/'));

%%
[ pId ] = {'P02','P04','P05','P07','P08','P09','P03ERL','P22AMS','P23AMS'};%
[ expMode ] = {'fVSpEM','cnEM'};
[ spkMode ] = {'Sorting','noSorting'};% can be either noSorting or Sorting
[ spk2LFPmode ] = {'plv','ppc'};% can be either noSorting or Sorting
[ stratMode ] = {'On'};%{'On','Off'}

%%
[ savePath ] = '/media/rouxf/rds-share/resultsAUG2019/';
[ rdsPath ] = '/media/rouxf/rds-share/iEEG_DATA/MICRO/';

%% compute spk parameters
%computeSPKparamsEMtask( pId, expMode, spkMode, savePath, rdsPath  );

%% compute spectral parameters
%timeFreqAnalysisScriptEMtask( pId, expMode, savePath, rdsPath );

%% compute spk 2 lfp coupling
nRand = 400;
alpha = 0.05;
%spk2LFPcouplingScriptEMtask( pId, expMode, spkMode, spk2LFPmode, nRand, alpha, savePath, rdsPath );

%% compute spk 2 lfp coupling for hits vs misses
spk2lfCoupling4HitsANDmisses( pId, expMode, spkMode, spk2LFPmode, rdsPath, savePath, alpha, nRand , stratMode);
