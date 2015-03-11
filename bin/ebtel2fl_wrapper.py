#Name: ebtel_wrapper.py
#Author: Will Barnes
#Date: 1 February 2015

#Description: Wrapper for the EBTEL-C program

#Import modules
import os,os.path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import astropy.constants as aconst

def plot_collisional_exchange(data_directory,data_file,**kwargs):
    """Plot the collisional exchange term and collisional frequency as a function of time.
    
    Arguments:
    data_directory -- directory that contains EBTEL-C output files
    data_file -- file that contains EBTEL-C output to be plotted
    
    Optional keyword arguments
    print_fig_filename -- specify filename to print figure to (default: show on screen)
    
    """
    #Load the data
    data = np.loadtxt(data_directory+data_file)
    #Slice the array to get the vectors we want
    time = data[:,0]
    temp_e = data[:,1]
    temp_i = data[:,2]
    dens = data[:,3]

    #Define adiabatic index
    gamma = 5.0/3.0

    #Define coloumb logarithm
    lnlambda = 23.0
    #Define electric charge
    qe = aconst.e.gauss.value
    #Define electron and proton masses
    me = aconst.m_e.cgs.value
    mi = aconst.m_p.cgs.value
    #Define Boltzmann constant
    kb = aconst.k_B.cgs.value
    
    #Calculate collisional frequency
    def nu_ei(te,n):
        return 16.0*(np.pi)**(1/2)/3.0*qe**4/(me*mi)*(2*kb*te/me)**(-3/2)*n*lnlambda
        
    #Make the collisional exchange vector
    col = np.empty(len(temp_e))
    nu = np.empty(len(temp_e))
    for i in range(len(temp_e)):
        nu[i] = nu_ei(temp_e[i],dens[i])
        col[i] = kb*dens[i]/(gamma-1)*nu[i]*(temp_e[i]-temp_i[i])

    #Set up the figure
    fig = plt.figure(figsize=(10,10))
    ax = fig.gca()
    fs = 18
    line_col = ax.plot(time,col,label=r'col')
    ax.set_ylabel(r'collisional exchange (erg cm$^{-3}$ s$^{-1}$)',fontsize=fs)
    #Set up second axis for collisional frequency
    ax_nu = ax.twinx()
    line_nu = ax_nu.plot(time,nu,'k',label=r'$\nu$')
    ax_nu.set_ylabel(r'$\nu$ (s$^{-1}$)',fontsize=fs)
    ax_nu.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax_nu.locator_params(nbins=5)
    ax_nu.ticklabel_format(axis='y', style='sci', scilimits=(-2,2) )
    #configure legend
    lines = line_col + line_nu 
    labels = [l.get_label() for l in lines]
    ax.legend(lines,labels,loc=4)
    #set x-axis and title parameters
    ax.set_xlabel(r'$t$ (s)',fontsize=fs)
    ax.set_title(r'EBTEL Two-fluid Collisions',fontsize=fs)
    ax.set_xlim([time[0],time[-1]])
    
    #Check if output filename is specified
    if 'print_fig_filename' in kwargs:
        plt.savefig(kwargs['print_fig_filename'],format='eps',dpi=1000)
    else:
        plt.show()
    

def plot_ebtel_dem(data_directory,data_file,**kwargs):
    """Plot differential emission measure parameters generated by the EBTEL-C model.
    
    Arguments:
    data_directory -- directory that contains EBTEL-C output files
    data_file -- file that contains EBTEL-C output to be plotted
    
    Optional keyword arguments
    print_fig_filename -- specify filename to print figure to (default: show on screen)
    
    """
    
    #Load DEM data
    data = np.loadtxt(data_directory+data_file)
    #Slice array to get necessary vectors
    temp = data[:,0]
    dem_tr = data[:,1]
    dem_cor = data[:,2]
    dem_tot = data[:,3]
    
    #Set up the figure
    fig = plt.figure(figsize=(10,10))
    ax = fig.gca()
    fs = 18
    ax.plot(temp,dem_tr,label=r'TR')
    ax.plot(temp,dem_cor,label=r'corona')
    ax.plot(temp,dem_tot,label=r'total')
    ax.legend()
    ax.set_title(r'EBTEL Two-fluid DEM',fontsize=fs)
    ax.set_xlabel(r'$\log(T_{DEM})$ (K)',fontsize=fs)
    ax.set_ylabel(r'$\log($DEM$)$ (cm$^{-5}$ K$^{-1}$)',fontsize=fs)
    ax.set_xlim([5.5,7.5])
    
    #Check if output filename is specified
    if 'print_fig_filename' in kwargs:
        plt.savefig(kwargs['print_fig_filename'],format='eps',dpi=1000)
    else:
        plt.show()
    

def plot_ebtel(data_directory,data_file,**kwargs):
    """Plot plasma parameters generated by the EBTEL-C model.
    
    Arguments:
    data_directory -- directory that contains EBTEL-C output files
    data_file -- file that contains EBTEL-C output to be plotted
    
    Optional keyword arguments
    print_fig_filename -- specify filename to print figure to (default: show on screen)
    
    """
    
    #Load the data
    data = np.loadtxt(data_directory+data_file)
    #Slice the array to get the vectors we want
    time = data[:,0]
    temp_e = data[:,1]
    temp_i = data[:,2]
    dens = data[:,3]
    temp_apex_e = data[:,7]
    temp_apex_i = data[:,8]
    dens_apex = data[:,9]
    heat = data[:,15]
    
    #Set up the figure
    fig,axes = plt.subplots(3,1,figsize=(15,10))
    #Set the fontsize
    fs = 16
    #Plot the heating
    axes[0].plot(time,heat)
    axes[0].set_ylabel(r'$h$ (erg cm$^{-3}$ s$^{-1}$)',fontsize=fs)
    axes[0].set_title(r'EBTEL Two-fluid Plasma Parameters',fontsize=fs)
    axes[0].set_xlim([time[0],time[-1]])
    axes[0].locator_params(nbins=5)
    axes[0].ticklabel_format(axis='y', style='sci', scilimits=(-2,2) )
    
    #Plot the temperatures
    line_te = axes[1].plot(time,temp_e/10**6,label=r'$T_e$')
    line_ti = axes[1].plot(time,temp_i/10**6,'r--',label=r'$T_i$')
    axes[1].set_ylabel(r'$T$ (MK)',fontsize=fs)
    axes[1].yaxis.set_major_locator(MaxNLocator(prune='lower'))
    axes[1].locator_params(nbins=5)
    axes[1].ticklabel_format(axis='y', style='sci', scilimits=(-2,2) )
    #Set up second axis for density
    ax_n = axes[1].twinx()
    line_n = ax_n.plot(time,dens/10**8,'k',label=r'$n$')
    ax_n.set_ylabel(r'$n$ (10$^8$ cm$^{-3}$)',fontsize=fs)
    ax_n.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax_n.locator_params(nbins=5)
    ax_n.ticklabel_format(axis='y', style='sci', scilimits=(-2,2) )
    #configure legend
    lines = line_te + line_ti + line_n
    labels = [l.get_label() for l in lines]
    axes[1].legend(lines,labels,loc=1)
    #set limits on x axes
    axes[1].set_xlim([time[0],time[-1]])
    
    #Plot the apex temperature
    axes[2].plot(time,temp_apex_e/10**6)
    axes[2].plot(time,temp_apex_i/10**6,'r--')
    axes[2].set_ylabel(r'$T_a$ (MK)',fontsize=fs)
    axes[2].yaxis.set_major_locator(MaxNLocator(prune='lower'))
    axes[2].locator_params(nbins=5)
    axes[2].ticklabel_format(axis='y', style='sci', scilimits=(-2,2) )
    #Set up the second axis for the apex density
    ax_na = axes[2].twinx()
    ax_na.plot(time,dens_apex/10**8,'k')
    ax_na.set_ylabel(r'$n_a$ (10$^8$ cm$^{-3}$)',fontsize=fs)
    ax_na.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax_na.locator_params(nbins=5)
    ax_na.ticklabel_format(axis='y', style='sci', scilimits=(-2,2) )
    #Set the x-axis label and limit
    axes[2].set_xlim([time[0],time[-1]])
    axes[2].set_xlabel(r'$t$ (s)',fontsize=fs)
    
    #Check if output filename is specified
    if 'print_fig_filename' in kwargs:
        plt.savefig(kwargs['print_fig_filename'],format='eps',dpi=1000)
    else:
        plt.show()
         

def run_ebtel(exec_directory,config_directory,**kwargs):
    """Run ebtel executable for a single configuration or a whole directory of configurations
    
    Arguments:
    exec_directory -- path to directory that contains executable
    config_directory -- path to config file directory
    
    Optional keyword arguments:
    config_file -- specific config file (run only one instance of EBTEL)
    
    """
    #Change to executable directory
    subprocess.call(['cd',exec_directory],shell=True)
    
    #Check if we want to run a single file or a whole directory
    if 'config_file' in kwargs:
        #Get full output when running a single config file
        output = subprocess.check_output([exec_directory+'ebtel-2fl',config_directory+kwargs['config_file']])
    else:
        for name in os.listdir(config_directory):
            if os.path.isfile(config_directory+name):
                #Only get exit code when running many configurations
                output = subprocess.call([exec_directory+'ebtel-2fl',config_directory+name,'quiet'])
    
    #Print the output of the subprocess call
    print output

def print_xml_config(config_dictionary,**kwargs):
    """Print XML configuration file for the EBTEL-C model
    
    Arguments:
    config_dictionary -- dictionary that holds config file inputs
    
    Optional keyword arguments:
    config_file -- specify configuration filename (default: 'ebtel_config.xml')
    
    """
    
    #Check if we have passed a filename
    #If not, pass a default filename
    if 'config_file' in kwargs:
        config_file = kwargs['config_file']
    else:
        config_file = 'ebtel_config.xml'
        
    #Open the file
    f = open(config_file,'w')
    
    #Print necessary header info
    f.write('<?xml version="1.0" ?>\n')
    f.write('<input>\n')

    #Loop through dictionary and print to xml file
    for key in config_dictionary:
        #Print tab delimiter, brackets and keyword
        f.write('\t<')
        f.write(key)
        f.write('>')
        #Check if entry is a list
        #If so print it as a list
        if isinstance(config_dictionary[key],list) or type(config_dictionary[key]).__name__ == 'ndarray':
            #Make temporary list
            temp = config_dictionary[key]
            #Skip to new line
            f.write('\n')
            #Begin loop over list
            for i in range(len(config_dictionary[key])):
                f.write('\t\t<')
                f.write(key+str(i))
                f.write('>')
                f.write(str(temp[i]))
                f.write('</')
                f.write(key+str(i))
                f.write('>\n')
            #Print additional tab to preserve alignment
            f.write('\t')
        else:
            #Print value
            f.write(str(config_dictionary[key]))
        #Close the brackets and print newline
        f.write('</')
        f.write(key)
        f.write('>\n')
    
    #Close the main node of the file
    f.write('</input>')
    
    #Close the file
    f.close()