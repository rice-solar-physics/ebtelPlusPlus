#ebtel2fl_dem_analysis_main.py

#Will Barnes
#15 May 2015

#Import needed modules
import sys
import os
import argparse
import numpy as np
sys.path.append('../../../bin/')
import ebtel2fl_dem as ebd
import ebtel2fl_plot_em as ebpe
import ebtel2fl_plot as ebp

#set root directory
root_dir = '/data/datadrive2/EBTEL-2fluid_runs/'
#set figure file format
figdir = '%s_heating_runs/alpha%.1f/'
figname = 'ebtel2fl_L%.1f_tpulse%.1f_alpha%.1f_%s_heating'

#make parameter vectors
loop_length = [20.0,40.0,60.0,80.0,100.0,120.0]
alpha = [1.5,2.0,2.5]
#make waiting time vector
Tn = np.arange(250,5250,250)

#set static parameters
tpulse = 100.0
solver = 'rka4'
mc = 100

#parse species argument
parser = argparse.ArgumentParser(description='Script that performs DEM analysis for EBTEL-2fluid runs.')
#Add arguments to parser
parser.add_argument("-s","--species",help="Species to which the heating was applied for particular run.")
#Declare the parser dictionary
args = parser.parse_args()

#declare instance of Plotter class
surf_plot = ebp.Plotter()

#iterate over variable parameters
for i in range(len(alpha)):
    temp_max_save = []
    em_max_save = []
    for j in  range(len(loop_length)):
        #print status
        print "Processing L = %.1f, alpha = %.1f"%(loop_length[j],alpha[i])
        #get data
        dema = ebd.DEMAnalyzer(root_dir,args.species,alpha[i],loop_length[j],tpulse,solver,mc=mc,Tn=Tn)
        dema.process_raw()
        dema.many_slopes()
        dema.em_max()
        temp_max_save.append(np.mean(dema.temp_max,axis=1))
        em_max_save.append(np.mean(dema.em_max,axis=1))
        #plot data
        figname_temp = figdir%(args.species,alpha[i])+figname%(loop_length[j],tpulse,alpha[i],args.species)
        demp = ebpe.DEMPlotter(dema.temp_em,dema.em,alpha[i],Tn=Tn)
        demp.plot_em_max(dema.temp_max,dema.em_max,print_fig_filename=root_dir+figname_temp+'_TmaxVTn')
        demp.plot_em_slopes(dema.a_cool,dema.a_hot,print_fig_filename=root_dir+figname_temp+'_hs_compare')
        demp.plot_em_curves(print_fig_filename=root_dir+figname_temp+'_dem')
        #plot all em curves for given tn
        if not os.path.exists(root_dir+figname_temp+'_dem_mc/'):
            os.makedirs(root_dir+figname_temp+'_dem_mc/')
        for k in range(len(Tn)):
            demp.plot_em_curve(k,print_fig_filename=root_dir+figname_temp+'_dem_mc/'+figname%(loop_length[j],tpulse,alpha[i],args.species)+'_'+str(k)+'_dem')
    #build surface plot
    surf_plot.plot_surface(Tn,loop_length,temp_max_save,xlab=r'$L$ (Mm)',ylab=r'$T_n$ (s)',plot_title=r'$T_{max}$ Surface, $\alpha=$'+str(alpha[i]),print_fig_filename=root_dir+figdir+'t_max_surface_'+args.species+'_alpha'+str(alpha[i])+'_tpulse'+str(tpulse)+solver)
    surf_plot.plot_surface(Tn,loop_length,em_max_save,xlab=r'$L$ (Mm)',ylab=r'$T_n$ (s)',plot_title=r'EM$_{max}$ Surface, $\alpha=$'+str(alpha[i]),print_fig_filename=root_dir+figdir+'em_max_surface_'+args.species+'_alpha'+str(alpha[i])+'_tpulse'+str(tpulse)+solver)