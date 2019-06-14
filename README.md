# HGCElectronicsValidation
A set of analyzers/analysis scripts to help with the validation of the electronics simulation in CMSSW

```
cmsrel CMSSW_11_0_X_2019-06-13-1100
cd CMSSW_11_0_X_2019-06-13-1100/src
cmsenv
git cms-init
git cms-checkout-topic PFCal-dev:hgc_eolt_11_0_X
git clone https://github.com/PFCal-dev/HGCElectronicsValidation.git UserCode/HGCElectronicsValidation
scram b -j 8
```

## Radiation map analysis

A set of control plots for the noise, charge collection efficiency, fluence is produced based on a radiation map file
and on a geometry.

## Occupancy analysis

A set of control plots is produced globally per layer and per wafer-equivalent position in each layer.
The number of plots is very large and an additional script is used to loop over them
and produce a summary in terms of quantiles of the distributions.
A third histogram is used to build final profiles (function of radius, layer, etc.)
or to dump an ascii file with the summarized information per wafer equivalent and layer.

```
cmsRun test/hgcoccupancyanalysis_cfg.py input=input_dir
```

## Generator-level jet profiles

May be useful to profile the generator level jets which are being thrown on the HGCal surface area
to get an idea of its characteristics (pT, eta, em fraction spectra etc.).
The first script fills the histograms using the genJets collection, the second one dumps the plots

```
python test/scripts/genJetAnalyzer.py input_dir output_name
python test/scripts/compareJetProfiles.py name1=file1.root ...
```
