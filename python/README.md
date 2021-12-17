# python scripts #

This folder contains a minimal conda environment that allows the user to:
- convert a GROMACS xtc file to an XYZ file (script ```sample_convert_xtc_to_xyz.py```)
- generate a custom .ini parameter file (script ```setup_parfile.py```)
- test the software (more info inside the *tests* folder)

### Installation ###

In order to install the software you must have [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) installed.

Then, it is sufficient to run the following command:

```
conda env create --file conda_env_excogito.yml
```

to create the **excogito** environment. Once the packages are downloaded, **excogito** can be activated via:

```
conda activate excogito
```

### XTC to XYZ conversion

In order to convert a GROMACS XTC to XYZ you just need to run:

```
python3 sample_convert_xtc_to_xyz.py
```

Once you provide GROMACS XTC and GRO files and a reasonable name for your output, the script will perform the conversion making use of the [MDTraj 1.9.5](https://www.mdtraj.org/1.9.5/index.html) software 

### Parameter file setup

Running

```
python3 setup_parfile.py
```

will help you with the setup of the *ini* parameter file needed by EXCOGITO.

### Contacts ###

Marco Giulini (mrcgiulini@gmail.com)
