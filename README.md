# dti_stats
DTI stats pipeline used for quality assurance.

# Installation instructions:
1) Install [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
2) Install [camino](http://camino.cs.ucl.ac.uk/)
3) Download repos and (optional) example data:
```
git clone https://github.com/justinblaber/system_utils.git
git clone https://github.com/justinblaber/nifti_utils.git
git clone https://github.com/justinblaber/dwmri_visualizer.git
git clone https://github.com/justinblaber/dti_stats.git

# Optionally download example data
wget https://justinblaber.org/downloads/github/dti_stats/PREPROCESSED.zip
unzip PREPROCESSED.zip
```
4) In MATLAB:
```
>> addpath('system_utils');
>> addpath(genpath('nifti_utils'));
>> addpath(genpath('dwmri_visualizer'));
>> addpath('dti_stats');
```
If you've downloaded the example data, then edit the test script (only the `fsl_path` and `camino_path` variables should have to be changed) and run it:

```
>> edit test_dti_stats
```
The output PDF should look like:

<a href="https://justinblaber.org/downloads/github/dti_stats/dti_stats.pdf">
<p align="center">
  <img width="769" height="995" src="https://i.imgur.com/P6vLxCP.png">
</p>
<p align="center">
  <img width="768" height="994" src="https://i.imgur.com/2iuPADP.png">
</p>
</a>
