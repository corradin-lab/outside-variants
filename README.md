
# Outside Variants Repository
In this repository is the code for the analysis pipeline accompanying the paper:
**Cell type specificity of intralocus interactions reveals oligodendrocyte intrinsic mechanisms for multiple sclerosis (Corradin et al, 2019)**

The main component is a script to predict loci that significantly alters genetic risk when acting in conjunction with the putative GWAS loci for a trait.

In order to download the scripts, you should clone this repository via the commands

```
git clone https://github.com/corradin-lab/outside-variants.git
cd outside-variants
```
In order to install the Python dependencies, you will need the [Anaconda](https://store.continuum.io/cshop/anaconda/) Python distribution and package manager. After installing Anaconda, run the following commands to create an environment with dependencies:
```
conda env create —file environment.yml
source activate OVP_env
```

Once the above has completed, you can run:
```
python_scripts/OVP_script.py -h
```
to print a list of all command-line options. If these commands fail with an error, then something as gone wrong during the installation process.

To do a walk-through tutorial on how to use the scripts and the format of inputs/outputs, click [here](https://mybinder.org/v2/gh/corradin-lab/outside-variants/master?urlpath=lab/tree/OVP_tutorial.ipynb)

To learn about the theoretical framework of our approach, click [here].


