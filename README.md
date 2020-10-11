This analyzer was written to read the miniAOD tree for checking the GEN level information for ISR jet removing.


# ISR_miniAOD

Set up analyzer::
```
cmsrel CMSSW_10_2_9
cd CMSSW_10_2_9/src
cmsenv
git clone git@github.com:StealthStop/ISR_miniAOD.git
scram b -j 4
```

Run the analyzer:
```
cd ISR_miniAOD/miniAOD/python/
cmsRun ConfFile_cfg.py 
```
