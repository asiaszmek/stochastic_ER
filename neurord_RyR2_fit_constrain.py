#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''trying to fit RyR2 dynamics to Xu and Meissner 2004 with  CaM'''

import sys
import os
homepath = os.path.expanduser("~")
sys.path.insert(0, os.path.join(homepath, "ajustador"))

import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,nrd_fitness
from ajustador.helpers import converge,save_params
import os
from pathlib import Path

dirname = 'RyR2_CaM'
#name of model xml file for optimization
model_set = 'Model_KS_fit'
#name of experimental data, a simulation file in this case
exp_name = 'ryr_cam_oconc_no_cam'
#molecule to compare between 'experiments' and simulations
mol = {'RyRO':['RyRO1', "RyRO2", "RyRCaMO1", "RyRCaMO2"]}
#directory to store output during optimization

start_stim = 20 # time stim start sec 
#norm_method  = 'percent'


# number of iterations, use 1 for testing
iterations = 250
# default popsize = 8, use 3 for testing
popsize = 8
test_size = 20#
rootdir = os.getcwd()
if not dirname in os.listdir(rootdir):
    os.mkdir(os.path.join(rootdir,dirname))
tmpdir = os.path.join(rootdir, dirname, 'tmp')
#use #2 exp since loading exp 
#exp  =  aju.xml.NeurordResult(exp_set)
#this command indicates that experimental data are concentration in csv formatted files
os.chdir(dirname)
exp = loadconc.CSV_conc_set(exp_name,
                            stim_time=start_stim)
#print('***********************************',exp.data[0].waves['Cof'].wave)
#specify parameters to vary, either from ReactionScheme or InitialConditions
P  =  aju.xml.XMLParam
#Double check. #2
params  =  aju.optimize.ParamSet(
    P('FirstCa_fwd_rate',5e-11, min = 5e-13, max = 5e-10, xpath = '//Reaction[@id="RyRe"]/forwardRate'),
    P('FirstCa_bckd_rate',0.0001658, fixed='FirstCa_fwd_rate', constant=1000000000.0, xpath='//Reaction[@id="RyRe"]/reverseRate'),
    P('O1_fwd_rate',5e-10, min=5.91e-12, max=5.91e-9, xpath='//Reaction[@id="RyRa"]/forwardRate'),
    P('O1_bcw_rate',9.6, min=9.6e-3, max=9.6e2, xpath='//Reaction[@id="RyRa"]/reverseRate'),
    P('O2_fwd_rate', 5e-11, min=5e-15, max=5e-10, xpath='//Reaction[@id="RyRb"]/forwardRate'),
    P('O2_bcw_rate',13, min=13e-2, max=130, xpath='//Reaction[@id="RyRb"]/reverseRate'),
    P('C1C2_fwd_rate', 5e-12, fixed='O1_fwd_rate',constant=1e-3, xpath='//Reaction[@id="RyRd"]/forwardRate'),
    P('C1C2_bcw_rate', 1.235e-3, fixed="O1_bcw_rate", constant=9.5e-5,  xpath='//Reaction[@id="RyRd"]/reverseRate'),
    P('C2O2_fwd_rate',3.3333e-15, fixed="O2_fwd_rate",constant=6.6666e-05, xpath='//Reaction[@id="RyRc"]/forwardRate'),
    P('C2O2_bcw_rate',66.667e-3,fixed='O2_bcw_rate',constant= 0.005128, xpath='//Reaction[@id="RyRc"]/reverseRate'),
    P('C2C3_fwd_rate',1e-14, min=1e-17, max=1e-10, xpath='//Reaction[@id="RyRf"]/forwardRate'),
    P('C2C3_bcw_rate', 3, min=3e-5, max=30, xpath='//Reaction[@id="RyRf"]/reverseRate'))

#this command indicates that experiments are from a previous simulation
###################### END CUSTOMIZATION #######################################

fitness = nrd_fitness.specie_concentration_fitness(species_list=mol)

fit = aju.optimize.Fit(tmpdir, exp, model_set, None, fitness, params,
                       _make_simulation=aju.xml.NeurordSimulation.make,
                       _result_constructor=aju.xml.NeurordResult)
fit.load()
fit.do_fit(iterations, sigma=0.3)
mean_dict,std_dict,CV=converge.iterate_fit(fit,test_size,popsize)
########################################### Done with fitting

#to look at centroid [0] or stdev [6] of cloud of good results:
if callable(fit.optimizer.result):
    result = fit.optimizer.result()
else:
    result = fit.optimizer.result
for i,p in enumerate(fit.params.unscale(result[0])):
        print(fit.param_names()[i],'=',p, '+/-', fit.params.unscale(result[6])[i])

#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)

#to save
save_params.save_params(fit,0,1)
