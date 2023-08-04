""" \class test_suite

The file contains several python tests to check the correct installation of METool package. 
"""

import unittest
import pathlib as pl
import os
import subprocess
import datetime as dt
import numpy as np

t_start = dt.datetime.now()
bash_script = subprocess.Popen("./test_suite.sh", shell=True, stdout=subprocess.PIPE)
bash_script.wait()
print(bash_script.returncode)
print("test scripts run in ", dt.datetime.now() - t_start)

class test0(unittest.TestCase):
    """
    class that checks a three-cores Simulated Annealing run
    """
    def test0_exist(self):
        """check existence of .dat files"""
        count = len([name for name in os.listdir('./test0/') if name.endswith(".dat")])
        self.assertEqual(count , 1)
    def test0_SA(self):
        """open files and check that the optimizations correctly finished"""
        count = 0
        list_files = [name for name in os.listdir('./test0/') if name.endswith(".dat")]
        for el in list_files:
            with open(os.path.join("test0",el),"r") as sa_run:
                for line in sa_run:
                    if (line == "END OF SIMULATED ANNEALING OPTIMISATION\n"):
                        count += 1
        self.assertEqual(count , 1)
    def test0_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test0/log_test0.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test0_head_tail(self):
        """check log's head, tail and size"""
        rfile = open("test0/log_test0.log", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line = lines[0]
        tail_line = lines[-2]
        self.assertEqual(head_line, "parameter file is ../files/parameters/parameters_optimize_6d93_N31_small.ini") 
        self.assertEqual(tail_line, "finished mapping optimisation") # finished SA
        self.assertTrue(len(lines) < 10000)

class test1(unittest.TestCase):
    """
    class that checks the correct generation of random mappings and the measurement of their mapping entropy
    """
    def test1_exist(self):
        """check existence of log file"""
        path = pl.Path("test1/log_test1.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test1_exist_dat(self):
        """check existence of output file"""
        path = pl.Path("test1/6d93_random_N31.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test1_count_mapping(self):
        """check for consistency between the declared and effective number of random mappings (and the number of calculations of smap)"""
        nrd_declared = -1
        nrd_effective = 0
        nrd_smaps = 0
        with open("test1/log_test1.log","r") as rd_log:
            for line in rd_log:
                if line.startswith("Parameter: n_mappings"):
                    nrd_declared = int(line.split()[-1])
        with open("test1/6d93_random_N31.dat","r") as rd_run:
            for line in rd_run:
                if line.startswith("generating random mapping number"):
                    nrd_effective += 1
                elif line.startswith("random_smap = "):
                    nrd_smaps += 1
        self.assertEqual(nrd_declared,nrd_effective)
        self.assertEqual(nrd_declared,nrd_smaps)
    def test1_head_tail(self):
        """check log's head, tail and size"""
        rfile = open("test1/log_test1.log", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line = lines[0]
        tail_line = lines[-2]
        self.assertEqual(head_line, "parameter file is ../files/parameters/parameters_random_6d93_N31_small.ini")  # ini invece di Dat 
        self.assertEqual(tail_line, "finished random mappings generation")
        self.assertTrue(len(lines) < 1000)

class test2(unittest.TestCase):
    """
    class that checks loading mapping task
    """
    # check existence of file
    def test2_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test2/log_test2.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test2_output_exist(self):
        """check existence of output file"""
        path = pl.Path("test2/6d93_loaded_N31.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test2_head_tail(self):
        """check parameter file is correct and that file length is not excessive"""
        rfile = open("test2/log_test2.log", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line = lines[0]
        self.assertEqual(head_line, "parameter file is ../files/parameters/parameters_measure_6d93_N31.ini")
        self.assertTrue(len(lines) < 1000)
    def test2_correct_smap(self):
        """check that the mapping entropy is correct """
        rfile = open("test2/6d93_loaded_N31.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        tail_line = lines[-1]
        self.assertEqual(tail_line,"smap for mapping ../files/mappings/tamapin_ca_mapping.txt : 7387.674178")

class test3(unittest.TestCase):
    """
    class that checks the estimation of T_start inside optimisation
    """
    def test3_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test3/log_test3.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test3_head_tail(self):
        """check log's head and tail"""
        rfile = open("test3/log_test3.log", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line = lines[0]
        tail_line = lines[-2]
        self.assertEqual(head_line, "parameter file is ../files/parameters/parameters_optimize_6d93_N31_notemp.ini")
        self.assertEqual(tail_line[:29], "finished mapping optimisation")
    def test3_deltas(self):
        """check delta dat files"""
        ncores = None
        rfile = open("test3/log_test3.log", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        for line in lines:
            if line.startswith("Parameter: Ncores ="):
                ncores = int(line.split()[-1])
        ls = [el for el in os.listdir("./test3/") if el.startswith("6d93_fast_delta")]
        self.assertEqual(len(ls), ncores)
        for el in ls:
            rfile = open("test3/" + el, "r")
            lines = rfile.read().split("\n")
            rfile.close()
            avg_delta = float(lines[-2].split()[-1])
            self.assertTrue(50.0 <= avg_delta <= 1000.0)
            temp = float(lines[-1].split()[-1])
            self.assertTrue(-avg_delta/np.log(0.75) - 0.1 <= temp <= -avg_delta/np.log(0.75) +0.1)
            
class test4(unittest.TestCase):
    """
    class that checks an invalid task ID (ex. optimizer)
    """
    def test4_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test4/log_test4.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test4_err_exist(self):
        """check existence of error.dat"""
        path = pl.Path("test4/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test4_err_correct(self):
        """check error file"""
        rfile = open("test4/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line = lines[0]
        self.assertEqual(head_line,"Error. optimizer not in the list of accepted tasks")

class test5(unittest.TestCase):
    """
    class that checks the expected number of alignments
    """
    def test5_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test5/log_test5.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test5_count_alignments(self):
        """count the effective (from dat file) and expected (from log/parameter file) number of alignments. See if they match"""
        n_alignment = 0
        ncores = None
        rotmats_period = None
        MC_steps = None
        list_files = [name for name in os.listdir('./test5/') if name.endswith(".dat")]
        for el in list_files:
            with open(os.path.join("test5",el),"r") as sa_run:
                for line in sa_run:
                    if line.startswith("Updating rotation matrices"):
                        n_alignment += 1
        with open("test5/log_test5.log","r") as sa_log:
            for line in sa_log:
                if line.startswith("Parameter: Ncores"):
                    ncores = int(line.split()[-1])
                elif line.startswith("Parameter: rotmats_period"):
                    rotmats_period = int(line.split()[-1])
                elif line.startswith("Parameter: MC_steps"):
                    MC_steps = int(line.split()[-1])
        prev_alignments = (ncores*MC_steps/rotmats_period)-1
        self.assertEqual(n_alignment , prev_alignments)

class test6(unittest.TestCase):
    """
    class that checks the error output if the trajectory does not respect the declared number of frames
    """
    def test6_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test6/log_test6.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test6_err_exist(self):
        """check existence of error.dat"""
        path = pl.Path("test6/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test6_err_correct(self):
        """check error file"""
        rfile = open("test6/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line_start = lines[0][:23]
        self.assertEqual(head_line_start,"Error. Frames completed") 

class test7(unittest.TestCase):
    """
    class that checks the error output if the trajectory is cut at the last frame
    """
    # check existence of log
    def test7_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test7/log_test7.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test7_err_exist(self):
        """check existence of error.dat"""
        path = pl.Path("test7/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test7_err_correct(self):
        """check error file"""
        rfile = open("test7/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line_start = lines[0][:23]
        self.assertEqual(head_line_start,"Error. Frames completed")

# 
class test8(unittest.TestCase):
    """
    class that checks the error output if the energy file is not complete
    """
    def test8_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test8/log_test8.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test8_err_exist(self):
        """check existence of error.dat"""
        path = pl.Path("test8/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test8_err_correct(self):
        """check error file"""
        rfile = open("test8/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line_start = lines[0][:18]
        self.assertEqual(head_line_start,"Error. Energy file")


class test9(unittest.TestCase):
    """
    class that checks the error output if the mapping file is not complete (shorter than n_cg beads) 
    """
    def test9_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test9/log_test9.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test9_err_exist(self):
        """check existence of error file"""
        path = pl.Path("test9/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test9_err_correct(self):
        """check error file"""
        rfile = open("test9/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line_start = lines[0][:27]
        self.assertEqual(head_line_start,"Error. Mapping file shorter")

class test10(unittest.TestCase):
    """
    class that checks the error output if the mapping file is not complete (longer than n_cg beads)
    """
    def test10_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test10/log_test10.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test10_err_exist(self):
        """check existence of error file"""
        path = pl.Path("test10/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test10_err_correct(self):
        """check error file"""
        rfile = open("test10/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line_start = lines[0][:26]
        self.assertEqual(head_line_start,"Error. Mapping file longer")

class test11(unittest.TestCase):
    """
    class that checks the error output if the mapping file is not complete (value not between [0 ;n_at) ) 
    """
    def test11_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test11/log_test11.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test11_err_exist(self):
        """check existence of error file"""
        path = pl.Path("test11/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test11_err_correct(self):
        """check error file"""
        rfile = open("test11/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line_start = lines[0][:33]
        self.assertEqual(head_line_start,"Error. Mapping file contains atom")

class test12(unittest.TestCase):
    """
    class that checks the error output if the mapping file is not complete (each value must be INT)
    """
    def test12_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test12/log_test12.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test12_err_exist(self):
        """check existence of error file"""
        path = pl.Path("test12/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test12_err_correct(self):
        """check error file"""
        rfile = open("test12/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line_start = lines[0][:31]
        self.assertEqual(head_line_start,"Error. The line must be integer") 

class test13(unittest.TestCase):
    """
    class that checks the error output if the mapping file is not complete (it contains duplicates)
    """
    def test13_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test13/log_test13.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test13_err_exist(self):
        """check existence of error file"""
        path = pl.Path("test13/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test13_err_correct(self):
        """check error file"""
        rfile = open("test13/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line_start = lines[0][:38]
        self.assertEqual(head_line_start,"Error. Mapping file contains duplicate") 

### START ELIO CLASSES 

class test14(unittest.TestCase):
    """
    class that checks the error output if the parameter file contains, at least, a string VALUE instead of integer number
    """
    def test14_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test14/log_test14.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test14_err_exist(self):
        """check existence of error file"""
        path = pl.Path("test14/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test14_err_correct(self):
        """check error file"""
        rfile = open("test14/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line_start = lines[0][:32]
        self.assertEqual(head_line_start,"Error. The line must be integer.")

class test15(unittest.TestCase):
    """
    class that checks the error output if the energy file contains, at least, an integer value instead of float
    """
    def test15_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test15/log_test15.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test15_err_exist(self):
        """check existence of error file"""
        path = pl.Path("test15/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test15_err_correct(self):
        """check error file"""
        rfile = open("test15/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line_start = lines[0][:31]
        self.assertEqual(head_line_start,"Error. Each line must be FLOAT.")

class test16(unittest.TestCase):
    """
    class that checks the error output if the energy file contains, at least, one row containing more than one column
    """
    def test16_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test16/log_test16.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test16_err_exist(self):
        """check existence of error file"""
        path = pl.Path("test16/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test16_err_correct(self):
        """check error file"""
        rfile = open("test16/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line_start = lines[0][:45]
        self.assertEqual(head_line_start,"Error. Each row must contain ONLY ONE column.")

class test17(unittest.TestCase):
    """
    class that checks the error output if the trajectory file contains an integer number != n_atoms when n_column == 1
    """
    def test17_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test17/log_test17.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test17_err_exist(self):
        """check existence of error file"""
        path = pl.Path("test17/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test17_err_correct(self):
        """check error file"""
        rfile = open("test17/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line_start = lines[0][:26]
        self.assertEqual(head_line_start,"Error. The number of atoms")

class test18(unittest.TestCase):
    """
    class that checks the error output if the trajectory file contains, at least, one letter, instead of float at 2nd column when n_columns == 3
    """
    def test18_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test18/log_test18.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test18_err_exist(self):
        """check existence of error file"""
        path = pl.Path("test18/error.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test18_err_correct(self):
        """check error file"""
        rfile = open("test18/error.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line_start = lines[0][:31]
        self.assertEqual(head_line_start,"Error. Each line must be FLOAT.")

class test19(unittest.TestCase):
    """
    class that checks the correct calculation of the norm of the mapping for 4AKE
    """
    def test19_output_exist(self):
        """check existence of output file"""
        path = pl.Path("test19/4ake_norms_N214.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test19_correct_coord_number(self):
        """check correct atomistic coordination number"""
        exp_coord_num = 13.716035
        rfile = open("test19/4ake_norms_N214.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        for ln in lines:
            if ln.startswith("normalized n_coord_at"):
                true_coord_num = float(ln.split()[-1])
        self.assertEqual(true_coord_num,exp_coord_num)

class test20(unittest.TestCase):
    """
    class that checks calculation of cosines and distances between two identical mappings
    """
    def test20_output_exist(self):
        """check existence of output file"""
        path = pl.Path("test20/6d93_cosine_N31.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test20_cosines(self):
        """check that all cosines are calculated and equal to one"""
        with open("test20/log_test20.log","r") as dist_log:
            for line in dist_log:
                if line.startswith("Parameter: frames"):
                    frames = int(line.split()[-1])
        expected_cosines = np.ones(frames)
        true_cosines = []
        rfile = open("test20/6d93_cosine_N31.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        for ln in lines:
            if ln.startswith("cosine"):
                true_cosines.append(float(ln.split()[1]))
        self.assertTrue((true_cosines == expected_cosines).all())
    def test20_distances(self):
        """check that all distances are calculated and equal to zero"""
        with open("test20/log_test20.log","r") as dist_log:
            for line in dist_log:
                if line.startswith("Parameter: frames"):
                    frames = int(line.split()[-1])
        expected_distances = np.zeros(frames)
        true_distances = []
        rfile = open("test20/6d93_cosine_N31.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        for ln in lines:
            if ln.startswith("distance"):
                true_distances.append(float(ln.split()[1]))
        self.assertTrue((true_distances == expected_distances).all())

class test21(unittest.TestCase):
    """
    class that checks the correct calculation of the norm of the Calpha mapping for 6D93
    """
    def test21_output_exist(self):
        """check existence of output file"""
        path = pl.Path("test21/6d93_norms_N31.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test21_consistent_norms(self):
        """check consistency with calculated norms"""
        prev_norms = np.loadtxt("./files/6d93_norms_test21.csv")
        rfile = open("test21/6d93_norms_N31.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        curr_norms = []
        for ln in lines:
            if ln.startswith("normalized_norm"):
                curr_norms.append(float(ln.split()[-1]))
        curr_norms = np.array(curr_norms)
        self.assertTrue((curr_norms == prev_norms).all())

class test22(unittest.TestCase):
    """
    class that checks the correct calculation of distance matrix between mappings
    """
    def test22_log_exist(self):
        """check existence of log file"""
        path = pl.Path("test22/log_test22.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test22_output_exist(self):
        """check existence of output file"""
        path = pl.Path("test22/4ake_distmat_N214.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test22_distmat_exist(self):
        """check existence of distance matrix file"""
        path = pl.Path("test22/4ake_distmat_N214.csv")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test22_distmat_shape(self):
        """check shape of distance matrix"""
        dmat = np.loadtxt("test22/4ake_distmat_N214.csv")
        prev_shape = (10,10)
        self.assertTrue(prev_shape == dmat.shape)

class test23(unittest.TestCase):
    """
    class that checks the correct functioning of criterion 3 (fast clustering)
    """
    def test23_output_exist(self):
        """check existence of output file"""
        path = pl.Path("test23/6d93_loaded_N31.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test23_correct_smap(self):
        """check that the mapping entropy is correct """
        rfile = open("test23/6d93_loaded_N31.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        tail_line = lines[-1]
        self.assertEqual(tail_line,"smap for mapping ../files/mappings/tamapin_ca_mapping.txt : 7535.859930")
    def test23_check_pairs(self):
        """check existence of log file and that the number of pairs matches its expected value"""
        path = pl.Path("test23/log_test23.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
        rlogfile = open(path, "r")
        log_lines = rlogfile.read().split("\n")
        rlogfile.close()
        for ln in log_lines:
            if ln.startswith("Parameter: frames"):
                frames = int(ln.split()[-1])
            if ln.startswith("Parameter: stride"):
                stride = int(ln.split()[-1])
            if ln.startswith("effective frames"):
                pivots = int(ln.split()[-1])
        exp_pivots = frames/stride + 1
        self.assertEqual(exp_pivots , pivots)

class test24(unittest.TestCase):
    """
    class that checks the correct functioning of task optimize_kl
    """
    def test24_output_exist(self):
        """check existence of output file"""
        path = pl.Path("test24/6d93kl_SA_N31_0.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test24_check_probabilities(self):
        """check that the program is checking probabilities"""
        path = pl.Path("test24/log_test24.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
        rfile = open(path, "r")
        lines = rfile.read().split("\n")
        rfile.close()
        checks = [0, 0, 0] # checking stuff
        for ln in lines:
            if ln == "reading probabilities":
                checks[0] = 1
            if ln.startswith("sum of probs 1.00000"):
                checks[1] = 1
            if ln == "properly normalized array of probabilities!":
                checks[2] = 1
        self.assertTrue(checks == [1,1,1])
    def test24_use_probabilities(self):
        """check that the program is using probabilities"""
        rfile = open("test24/6d93kl_SA_N31_0.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        head_line = lines[1]
        self.assertEqual(head_line,"kl_flag is 1: explicitly calculating Kullback-Leibler divergences")

class test25(unittest.TestCase):
    """
    class that checks the correct functioning of task measure_kl
    """
    def test25_output_exist(self):
        """check existence of output file"""
        path = pl.Path("test25/6d93_measurekl_N31.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test25_check_probabilities(self):
        """check that the program is checking probabilities"""
        path = pl.Path("test25/log_test25.log")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
        rfile = open(path, "r")
        lines = rfile.read().split("\n")
        rfile.close()
        checks = [0, 0, 0] # checking stuff
        for ln in lines:
            if ln == "reading probabilities":
                checks[0] = 1
            if ln.startswith("sum of probs 1.00000"):
                checks[1] = 1
            if ln == "properly normalized array of probabilities!":
                checks[2] = 1
        self.assertTrue(checks == [1,1,1])
    def test25_use_probabilities_correct_smap(self):
        """check that the program is using probabilities and that the mapping entropy is correct"""
        rfile = open("test25/6d93_measurekl_N31.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        kl_line = lines[-2]
        self.assertEqual(kl_line,"explicitly calculating Kullback-Leibler divergences")
        tail_line = lines[-1]
        self.assertEqual(tail_line,"smap for mapping ../files/mappings/tamapin_ca_mapping.txt : 0.158291")

class test26(unittest.TestCase):
    """
    class that checks the correct functioning of task random_kl
    """
    def test26_output_exist(self):
        """check existence of output file"""
        path = pl.Path("test26/6d93_randomkl_N31.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test26_use_probabilities(self):
        """check that the program is using probabilities"""
        rfile = open("test26/6d93_randomkl_N31.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        kl_line = lines[1]
        self.assertEqual(kl_line,"explicitly calculating Kullback-Leibler divergences")

class test27(unittest.TestCase):
    """
    class that checks the correct functioning of task random_kl
    """
    def test27_output_exist(self):
        """check existence of output files"""
        path = pl.Path("test27/spinkls_fast_delta_N2_0.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
        path = pl.Path("test27/spinkls_SA_N2_0.dat")
        self.assertEqual((str(path), path.is_file()), (str(path), True))
    def test27_smaps(self):
        """check that all 3 possible values of smap are calculated"""
        rfile = open("test27/spinkls_SA_N2_0.dat", "r")
        lines = rfile.read().split("\n")
        rfile.close()
        smaps = [float(ln.split()[-1]) for ln in lines if ln.startswith("new_smap")]
        assert 0.007167 in smaps
        assert 0.0 in smaps



if __name__ == '__main__':
    unittest.main()
