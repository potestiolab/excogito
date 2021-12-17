# script to generate a valid parameter file
from datetime import datetime

print("guided generation of the parameter file for EXCOGITO")

tasks = ["optimize", "random", "measure", "norm", "cosine", "distance", "optimize_kl", "measure_kl"]

mandatory_pars = {
    "optimize" : [ "atomnum", "frames", "cgnum", "MC_steps", "rotmats_period", "Ncores" ],
    "random" : [ "atomnum", "frames", "cgnum", "n_mappings" ],
    "measure" : [ "atomnum", "frames", "cgnum" ],
    "norm": [ "atomnum", "frames", "cgnum" ],
    "cosine" : [ "atomnum", "frames", "cgnum" ],
    "distance" : [ "atomnum", "frames", "cgnum", "n_mappings" ],
    "optimize_kl" : ["atomnum", "frames", "cgnum", "MC_steps", "Ncores"],
    "measure_kl" : [ "atomnum", "frames", "cgnum" ]
}

optional_pars = {
    "optimize" : [ "criterion", "t_zero", "decay_time"],
    "random" : [ "criterion"],
    "measure" : [ "criterion"],
    "norm": [],
    "cosine" : [],
    "distance" : [],
    "optimize_kl" : ["criterion", "t_zero", "decay_time"],
    "measure_kl" : [ "criterion"]
}

pars_description = {
    "atomnum" : "number of atoms in your structure",
    "frames" : "number of frames in your trajectory",
    "cgnum" : "number of coarse-grained sites",
    "n_mappings": "number of coarse-grained mappings",
    "MC_steps" : "number of Monte Carlo steps",
    "rotmats_period" : "steps between two calculations of rotation matrices ",
    "criterion" : "criterion for hierarchical clustering:\n0 : fixed number of clusters\n1 : cophenetic distance\n2 : multiple numbers of clusters\n3 : fast clustering on a subset of configurations (only for continuous trajectories)",
    "Ncores" : "number of cores to use",
    "t_zero" : "starting temperature for Simulated Annealing",
    "nclust" : "number of clusters",
    "distance" : "cophenetic distance",
    "max_nclust" : "upper limit to the number of clusters",
    "min_nclust" : "lower limit to the number of clusters",
    "decay_time" : "temperature decay for Simulated Annealing",
    "rsd" : "use RSD instead of RMSD",
    "stride" : "number of structures between two pivot configurations"    
}

pars_type = {
    "atomnum" : int,
    "frames" : int,
    "cgnum" : int,
    "n_mappings": int,
    "MC_steps" : int,
    "rotmats_period" : int,
    "criterion" : int,
    "Ncores" : int,
    "t_zero" : float,
    "criterion" : int, 
    "nclust" : int,
    "distance" : float,
    "max_nclust" : int,
    "min_nclust" : int,
    "Ncores" : int,
    "decay_time" : float,
    "rsd" : int,
    "stride" : int
}

clustering_pars = {
    0: [ "nclust" ],
    1: [ "distance" ],
    2: [ "min_nclust", "max_nclust" ],
    3: [ "stride", "nclust"]
}

# mandatory parameters for every task

task=input("Insert the task you would like to perform among the following: " + str(tasks) + "\n")
if task not in tasks:
    raise ValueError("The inserted value is not in the list of available tasks, try again\n")

my_pars = {}

def retrieve_parameter(par_name, par_type, par_desc):
    if par_type == int:
        while True:
            par = input("Insert the " + par_desc + "(type " + str(par_type) + "):")
            try:
                par=int(par)
                break
            except ValueError:
                print("The inserted value for ", par_name, " is not an integer number, try again\n")
    elif par_type == float:
        while True:
            par = input("Insert the " + par_desc + "(type " + str(par_type) + "):")
            try:
                par=float(par)
                break
            except ValueError:
                print("The inserted value for ", par_name, " is not a float number, try again\n")
    else:
        raise Exception("Unrecognized data type ", par_type)
    return par

def get_mandatory_parameters(task):
    for par_name in mandatory_pars[task]:
        print(par_name)
        param = retrieve_parameter(par_name,pars_type[par_name],pars_description[par_name])
        print(par_name, " = ", param)
        my_pars[par_name] = param

def get_optional_parameters(task):
    for par_name in optional_pars[task]:
        opt = input("Define optional parameter " + par_name + "? (y/n)")
        if opt == "y" or opt == "yes":
            param = retrieve_parameter(par_name,pars_type[par_name],pars_description[par_name])
            my_pars[par_name] = param
            if par_name == "criterion":
                if param not in clustering_pars.keys():
                    raise ValueError("invalid criterion for clustering")
                else:
                    for cl_par in clustering_pars[param]:
                        param = retrieve_parameter(cl_par,pars_type[cl_par],pars_description[cl_par])
                        my_pars[cl_par] = param

def write_parameters(task, pars_dict):
    filename = "parameters_sample_" + task + "_N" + str(pars_dict["cgnum"]) + ".ini"
    parfile = open(filename,"w")
    parfile.write("# Sample parameter file for Excogito \n")
    parfile.write("#\n# Generated on " + str(datetime.now())+ "\n")
    for key in pars_dict.keys():
        key_spaces = (18 - len(key))*" "
        val_spaces = (12 - len(str(pars_dict[key])))*" "
        line = key + key_spaces + "= "+ str(pars_dict[key]) + val_spaces + "; " + pars_description[key].split("\n")[0] + "\n"
        parfile.write(line)
    parfile.close()

get_mandatory_parameters(task)
print(my_pars)
opt = input("Insert optional parameters? (y/n)")
if opt == "y" or opt == "yes":
    get_optional_parameters(task)
print(my_pars)
write_parameters(task, my_pars)

