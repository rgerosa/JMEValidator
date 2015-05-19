JMEValidator
======

Common framework for CMS JEC / JER analyzes.

### Dependencies

This package depends on [TreeWrapper](https://github.com/blinkseb/TreeWrapper). Inside the `src` folder of you CMSSW release, do:

```sh
git clone git@github.com:blinkseb/TreeWrapper.git JMEAnalysis/TreeWrapper

# E/Gamma ID
git cms-merge-topic ikrav:egm_id_74X_v0

scram b -j8
```
