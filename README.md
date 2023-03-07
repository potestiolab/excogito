# README #

EXCOGITO is the program to investigate the mapping problem in coarse-grained modelling of biomolecules.

If you use EXCOGITO in your research please cite:

**EXCOGITO, an EXtensible COarse-GraIning TOol**, M Giulini, R Fiorentini, L Tubiana, R Potestio, *in preparation*

**An Information-Theory-Based Approach for Optimal Model Reduction of Biomolecules**, M Giulini, R Menichetti, MS Shell, R Potestio, *Journal of chemical theory and computation 16 (11), 6795-6813*

**A journey through mapping space: characterising the statistical and metric properties of reduced representations of macromolecules**, R Menichetti, M Giulini, R Potestio, *The European Physical Journal B 94 (10), 1-26*

# 1. Installation #

## 1.1 General requirements on Linux systems ##

The only requirement is to have [Openmp](https://www.openmp.org/) installed on your machine.

## 1.2 Additional requirements on MAC OS ##

* Install [argp](https://formulae.brew.sh/formula/argp-standalone) by using homebrew. At the terminal, run this command: 
```
brew install argp-standalone 
```

* Install [xcode](https://apps.apple.com/it/app/xcode/id497799835?mt=12) if your version is higher than MacOs 10.7. You are not required to install the Xcode App from AppStore. 
At the terminal, just run this command (about 15 Gb are required free on your disk even though, at the end of installation, only 2 Gb will be consumed)
```
xcode-select --install 
```

In order to have access to OpenMP libraries you can install [libomp](https://formulae.brew.sh/formula/libomp) by using homebrew. At the terminal, run this command: 
```
brew install libomp
```

## 1.3 Compiling ## 

The code can be compiled using [CMake](https://cmake.org/). A minimal installation is obtained following these steps:

1.  create a directory in *excogito*, such as *build*
```bash
mkdir build
cd build
```
2. run cmake from *build*, calling the outer directory
```bash
cmake ..
```
3. run make
```bash
make
```

### 1.3.2 Compilation options ###
Cmake allows to specify several options, such as the C compiler, compilation links and compilation flags. For instance, if the optimized Intel C compiler (icc) is available, step 2 may be substitued by:
```bash
cmake .. -DCMAKE_C_COMPILER=icc -DCMAKE_C_FLAGS="-Ofast -fopenmp -I./include -mkl -xSSE4.2 -parallel -ipo -mcpu=native"
```

On MacOs, the C compiler identification should be AppleClang (check the first line printed on terminal after launching the command ```cmake ..```).

# 2. Running #

The typical usage of the program consists in a call to *excogito* with one of the following options:
- **optimize**: to optimize the coarse-grained mapping by minimising its mapping entropy;
- **random**: to randomly generate coarse-grained representations and measure the associated mapping entropies;
- **measure**: to measure the mapping entropy of a mapping provided by the user (in the form of a .txt file);
- **norm**: to calculate the norm of a mapping (provided by the user) throughout a trajectory;
- **cosine**: to calculate pairwise distance and cosine between a pair of mappings (provided by the user) throughout a trajectory;
- **distance**: to calculate the distance matrix between a data set of mappings (provided by the user) over a single conformation;
- **optimize_kl**: to optimize the coarse-grained mapping by minimising its mapping entropy, calculated using the original Kullback-Leibler divergence;
- **random_kl**: to randomly generate coarse-grained representations and measure the associated mapping entropies, calculated using the original Kullback-Leibler divergence;
- **measure_kl**: to measure the mapping entropy of a mapping provided by the user (in the form of a .txt file), calculated using the original Kullback-Leibler divergence.

Each task can require different input files, which are provided to the program in the form of command-line options. 

For further information, please type on terminal ```./excogito --help``` or ```./excogito -h```

Alternatively, for printing a short usage message, please type: ```./excogito --usage``` or ```./excogito -u```

After selecting which task is suitable for your purposes, read carefully the documentation below according to your choice.

## 2.1. Optimize Task ##

The **optimize** task requires the _protein code_ string and three input files: _parameter_, _trajectory_, and _energy_. 

In order to launch the **optimize** task follow this syntax: 

```bash 
./excogito optimize -p $parameter_file.ini -t $trajectory_file.xyz -e $energy_file.txt -c $prot_code

or

./excogito optimize --p $parameter_file.ini --t $trajectory_file.xyz --e $energy_file.txt --code $prot_code
```

For further information, please type on terminal ```./excogito optimize```

## 2.2. Random Task ##

The **random** task requires the _protein code_ string and three input files: _parameter_, _trajectory_, and _energy_. 

In order to launch the **random** task follow this syntax: 

```bash   
./excogito random -p $parameter_file.ini -t $trajectory_file.xyz -e $energy_file.txt -c $prot_code

or

./excogito random --p $parameter_file.ini --t $trajectory_file.xyz --e $energy_file.txt --code $prot_code
```

For further information, please type on terminal ```./excogito random```

## 2.3. Measure Task ##

The **measure** task requires the _protein code_ string and four input files: _parameter_, _trajectory_, _energy_, and _mapping_. 

In order to launch the **measure** task follow this syntax: 


```bash
./excogito measure -p $parameter_file.ini -t $trajectory_file.xyz -e $energy_file.txt -c $prot_code -m $mapping_file.txt

or 

./excogito measure --p $parameter_file.ini --t $trajectory_file.xyz --e $energy_file.txt --prot_code $prot_code --m1 $mapping_file.txt
```

For further information, please type on terminal ```./excogito measure```

## 2.4. Norm Task ##

The **norm** task requires the _protein code_ string and three input files: _parameter_, _trajectory_, and _mapping_. 

In order to launch the **norm** task follow this syntax: 

```bash 
./excogito norm -p $parameter_file.ini -t $trajectory_file.xyz -c $prot_code -m $mapping_file.txt

or 

./excogito norm --p $parameter_file.ini --t $trajectory_file.xyz --prot_code $prot_code --m1 $mapping_file.txt
```

For further information, please type on terminal ```./excogito norm```

## 2.5. Cosine Task ##

The **cosine** task requires the _protein code_ string and four input files: _parameter_, _trajectory_, _1st mapping_, and _2nd mapping_.

In order to launch the **cosine** task follow this syntax: 


```bash 
./excogito cosine -p $parameter_file.ini -t $trajectory_file.xyz -c $prot_code -m $mapping_file.txt -n $mapping_file2.txt

or 

./excogito cosine --p $parameter_file.ini --t $trajectory_file.xyz --prot_code $prot_code --m1 $mapping_file.txt --m2 $mapping_file2.txt
```

For further information, please type on terminal ```./excogito cosine```

## 2.6. Distance Task ##

The **distance** task requires the _protein code_ string and thre input files: _parameter_, _trajectory_, _mapping matrix_.

In order to launch the **distance** task follow this syntax:

```bash 
./excogito distance -p $parameter_file.ini -t $trajectory_file.xyz -c $prot_code -x $mapping_matrix_file.txt

or

./excogito distance --p $parameter_file.ini --t $trajectory_file.xyz --prot_code $prot_code --matrix $mapping_matrix_file.txt

```

For further information, please type on terminal ```./excogito distance```

## 2.7. Optimize_kl Task ##

The **optimize_kl** task requires the _protein code_ string and three input files: _parameter_, _trajectory_, and _probability_. 

In order to launch the **optimize_kl** task follow this syntax: 

```bash 
./excogito optimize -p $parameter_file.ini -t $trajectory_file.xyz -r $probability_file.txt -c $prot_code

or

./excogito optimize --p $parameter_file.ini --t $trajectory_file.xyz --probs $probability_file.txt --code $prot_code
```

For further information, please type on terminal ```./excogito optimize_kl```

## 2.8. Random_kl Task ##

The **random_kl** task requires the _protein code_ string and three input files: _parameter_, _trajectory_, and _probability_. 

In order to launch the **random_kl** task follow this syntax: 

```bash   
./excogito random_kl -p $parameter_file.ini -t $trajectory_file.xyz -r $probability_file.txt -c $prot_code

or

./excogito random_kl --p $parameter_file.ini --t $trajectory_file.xyz --probs $probability_file.txt --code $prot_code
```

For further information, please type on terminal ```./excogito random_kl```

## 2.9. Measure_kl Task ##

The **measure_kl** task requires the _protein code_ string and four input files: _parameter_, _trajectory_, _probability_, and _mapping_. 

In order to launch the **measure_kl** task follow this syntax: 


```bash
./excogito measure_kl -p $parameter_file.ini -t $trajectory_file.xyz -r $probability_file.txt -c $prot_code -m $mapping_file.txt

or 

./excogito measure_kl --p $parameter_file.ini --t $trajectory_file.xyz --probs $probability_file.txt --prot_code $prot_code --m1 $mapping_file.txt
```

For further information, please type on terminal ```./excogito measure_kl```


# 3. Which arguments are mandatory? A short explanation #

As shown in **Section 2.x**, the _protein code_ string  and two files are always mandatory, namely the _parameter file_ and the _xyz trajectory file_. The other files can be mandatory, depending on the chosen task.

What are these files? 

- **$parameter_file.ini**       → Set of parameters in _ini_ format for the algorithm (see 3.1). Examples are present in */examples/parameters*;
- **$trajectory_file.xyz**      → Trajectory in _xyz_ format (see the Section 3.2). An example is present in */examples/trajectories*;
- **$energy_file.txt**          → File with the energies corresponding to each configuration in the trajectory (see the Section 3.3). An example is present in */examples/energies*;
- **$prot_code**                → Unique string that identifies the structure (see 3.4). It will be used to generate the output files;
- **$mapping_file.txt**         → Mapping file, containing the indices of the retained atoms (see 3.5). An example is present in */examples/mappings*;
- **$mapping_file2.txt**        → 2<sup>nd</sup> Mapping file, containing the indices of the retained atoms (see 3.5). An example is present in */examples/mappings*;
- **$mapping_matrix_file.txt**  → Matrix with _n_mappings_ CG mappings (see 3.6).
- **$probability_file.txt**     → File with the probabilities corresponding to each configuration in the trajectory (see 3.7). They must sum to 1.0. An example is present in */examples/probabilities*;

## 3.1. Parameter FILE ##

The core element of EXCOGITO is the parameter file, which is employed to define the constants used in the different tasks. 

A sample parameter file for each task can be found in /examples/parameters.

There exist 16 parameters, but only few of them are mandatory for the selected task. They are illustrated in the following table:

| Parameter | Description | Type | Mandatory | Suggested value |
| ----------- | ----------- | ---- | --------- | ---------------|
| atomnum     | number of atoms in the system | int | all | |
| frames   | number of frames in the trajectory | int | all | < 5000 on laptops, < 15000 if criterion != 1 |
| cgnum | number of CG sites | int | all | between atomnum/20 and atomnum/2|
| criterion | criterion for clustering | int | O-R-M |  0, 1, 2, 3|
| nclust | number of CG macrostates | int | C0 - C3 | between frames/500, and frames/100|
| n_mappings | number of mappings in tasks **random** and **distance** | int | R-D | | 
| MC_steps | number of MC step in task **optimize** | int | O |  > 5000|
| rotmats_period | MC steps between two full alignments in task **optimize** | int | O | |
| t_zero | starting temperature in task **optimize** | double | O | |
| distance |  cophenetic distance threshold | double | C3 | |
| max_nclust | upper number of clusters | int | C2 | between frames/100 and frames/50|
| min_nclust | lower number of clusters | int | C2 | between frames/1000 and frames/500 (must be < max\_nclust)|
| Ncores | number of cores | int | no | |
| decay_time | governs temperature decay in task **optimize** | double | O | |
| rsd | use rsd (if 1) instead of rmsd (if 0) | int | no | |
| stride | number of structures between two pivot configurations | int | C1 | ~ 10 if frames between 1000 and 10000|

O-R-M-D refer to the tasks (optimize/optimize_kl, random/random_kl, measure/measure_kl, distance) in which the parameter is mandatory. C0 .. C3 indicates that the parameter is mandatory if the clustering criterion is equal to 0 .. 3, respectively.

## Clustering

Four criteria for hierarchical clustering:

- **0** *Maxclust* clustering: configurations are lumped into *Nclust* macrostates;
- **1** *Fast clustering*: as in criterion **0**, but applied to a set of pivot configurations. Labels of intermediate structures are assigned to the closer pivot;
- **2** *Multiple maxclust*: as described in *Giulini et al.* (JCTC, 2020);
- **3** *Maxdist* clustering: clustering with the cophenetic distance;

## 3.2. Trajectory FILE ##

The trajectory should be provided in the xyz format. The first line of each frame indicates the number of atoms, while the second can contain an arbitrary string. As an example, a trajectory with 2 frames and 3 atoms should resemble the following string:

```
3

X   2.53  2.09   3.55
X   2.57  1.95   3.51
X   2.45  1.87   3.46
3

X   2.69  1.96   3.40
X   2.80  1.91   3.43
X   2.67  2.03   3.28
```

In the *python* subdirectory there is a script that helps with the conversion from GROMACS XTC to the XYZ format.

## 3.3. Energy FILE ##

Energy files, mandatory for tasks **optimize**, **random**, and **measure**, should contain one value for each frame in the trajectory.

## 3.4 Protein Code ##

The protein code is a string that is used to create output files. Don't insert spaces or special characters in this string

## 3.5 Mapping FILES ##

A mapping file, mandatory for tasks **measure**, **norm**, and **cosine** is a file with an integer per line. The value correspond to the index of the atom in the xyz trajectory. As an example, a mapping with 8 sites on a peptide of 50 sites should respect the following format:

```
3
7
19
21
26
34
40
47
```

## 3.6. Mapping Matrix FILES ##

A mapping matrix is mandatory for task **distance**. It is simply a series of transposed mappings. If we aim at computing the distance matrix between three mappings with 8 sites on a peptide of 50 sites, we must respect the following syntax:

```
3 7 19 21 26 34 40 47
2 8 19 24 25 38 41 44
0 10 12 20 29 31 35 49
```

## 3.7. Probability FILE ##

Probability files, mandatory for tasks **optimize_kl**, **random_kl**, and **measure_kl**, must contain one value for each frame in the trajectory and should be properly normalized to 1. For a trajectory of 5 frames, the following file is acceptable:

```
0.1
0.15
0.6
0.05
0.1
```

# 4. Examples #

Inside the directory examples there are example files for the *6d93* protein, allowing the user to try all the different tasks:

- **optimize**: ``` ./build/excogito optimize -p examples/parameters/parameters_optimize_6d93_N31_small.ini -t examples/trajectories/6d93_100frames.xyz -e examples/energies/6d93_energies_100frames.txt -c 6d93 ```

- **random**: ```./build/excogito random -p examples/parameters/parameters_random_6d93_N31_small.ini -t examples/trajectories/6d93_100frames.xyz -e examples/energies/6d93_energies_100frames.txt -c 6d93 ```

- **measure**: ```./build/excogito measure -p examples/parameters/parameters_loadca_6d93_N31.ini -t examples/trajectories/6d93_1000frames.xyz -e examples/energies/6d93_energies_1000frames.txt -c 6d93 -m examples/mappings/tamapin_ca_mapping.txt```

- **norm**: ```./build/excogito norm -p examples/parameters/parameters_norm_6d93_N31.ini -t examples/trajectories/6d93_1000frames.xyz -e examples/energies/6d93_energies_1000frames.txt -c 6d93 -m examples/mappings/tamapin_ca_mapping.txt```

- **cosine**: ```./build/excogito cosine -p ./examples/parameters/parameters_cosine_6d93_N31.ini -t ./examples/trajectories/6d93_1000frames.xyz -e ./examples/energies/6d93_energies_1000frames.txt -c 6d93 -m ./examples/mappings/tamapin_ca_mapping.txt --m2 ./examples/mappings/tamapin_nextca_mapping.txt```

- **distance**: ```./build/excogito distance -p examples/parameters/parameters_distance_6d93_N31.ini -t ./examples/trajectories/6d93_1frame.xyz -x examples/mappings/6d93_mapping_matrix.txt -c 6d93```

- **optimize**: ``` ./build/excogito optimize_kl -p examples/parameters/parameters_optimizekl_6d93_N31_notemp.ini -t examples/trajectories/6d93_100frames.xyz -r examples/probabilities/6d93_probs_100frames.txt -c 6d93 ```

- **random_kl**: ```./build/excogito random_kl -p examples/parameters/parameters_randomkl_6d93_N31.ini -t examples/trajectories/6d93_100frames.xyz -r examples/probabilities/6d93_probs_100frames.txt -c 6d93```

- **measure_kl**: ```./build/excogito measure_kl -p examples/parameters/parameters_measurekl_6d93_N31.ini -t examples/trajectories/6d93_100frames.xyz -r examples/probabilities/6d93_probs_100frames.txt -c 6d93 -m examples/mappings/tamapin_ca_mapping.txt```

# 5. Scaling values #

The approximated mapping entropy is calculated (tasks **optimize**, **random** and **measure**) without the scaling factor $` \frac{k_B \beta^2}{2}  `$ (see. [Giulini et al.](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.0c00676)). This factor should be computed by the user according to the temperature employed to simulate the system.

# 6. Documentation #

File ```refman.pdf``` in the ```docs``` directory contains detailed documentation authomatically generated with doxygen version 1.8.5.

A custom documentation can be generated in ```html``` and ```tex``` format by running ```doxygen excogito_doxygen.conf```.

# 7. Contacts #

Marco Giulini (mrcgiulini@gmail.com)
Raffaello Potestio (raffaello.potestio@unitn.it)
