# GenOnTheFly
CMS EDAnalyzers, python configurations, and submission tools for performing gen-level studies "on the fly" by analyzing events as they are generated and then discarding them

The `QCDAnalyzer` class looks at collections of `GenParticles` and `GenJets` in 
each event to produce spectra of particles, charged particles, postitively 
charged particles, and jets. It also matches particles with the  
most energetic jet within a given radius to analyze the fragmentation of the 
`GenJets`. Particles may be selected based on their species, and jets on their flavor. 
Note the jet flavor identification is currently not working in PYTHIA8.

To run a small 100 event test:

```
scram p CMSSW CMSSW_7_3_1
cd CMSSW_7_3_1/src
cmsenv

git clone https://github.com/rkunnawa/GenLevelStudy.git
cd GenOnTheFly
scram b
cd test 
cmsRun anaQCD_PYTHIA6_cfg.py
```

The analyzer output from different pt-hat ranges can be combined using the 
ROOT macro `makeCombinedPtHatSample.C` which takes the cross sections and 
input files as arguments, and weights all of the histograms appropriately.
An example usage of this macro is `test/examplePtHatCombination.C`.
