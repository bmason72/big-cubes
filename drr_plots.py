# import libraries

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import re
import math
from astropy.table import Table, QTable, join
import ipdb

array_color = {'12m': 'darkslateblue',
               '7m': 'darkorange',
               'TP': 'seagreen'}

def create_per_mous_db(mydb):
    '''
    Purpose: create a per mous version of Andres' table.

    Input: data base with eb per line

    Output: data base with mous per line

    NOTE: code not longer needed (8/13/2025)

    Date        Programmer      Description of Changes
    ---------------------------------------------------
    8/7/2025    A.A. Kepley     Original Code
    '''

    # get groups
    mydb_by_mous = mydb.group_by('mous_status_uid')

    # create output dictionary
    newdb_dict = {}
    for mykey in mydb_by_mous.keys():
        newdb_dict[mykey] = []

    # add variable to turn off messages after first round
    keymsg = True

    # iterate over groups and calculate aggregate values
    for mygroup in mydb_by_mous.groups:
        for mykey in mygroup.keys():

            # take sum
            if (mykey in ['dq_if_non-ima_m5_per_eb_mb','m5_wsu_asdm_size_gb']):
                newval = np.sum(mygroup[mykey])
                newdb_dict[mykey].append(newval)

            else:
                newdb_dict[mykey].append(mygroup[mykey][0])


    mous_db = Table(newdb_dict)
    
    return mous_db

    

def plot_asdm_vs_dq(mous_db, dq_val = 'dq_if_non_ima_gb'):
    '''
    Purpose: Plot ASDM size vs. non-imaging DQ

    Input: mous data base

    Ouput: plot

    Date        Programmer      Description of Changes
    --------------------------------------------------
    8/7/2025    A.a. Kepley     Original Code
    '''

  
    xaxis_vals = mous_db['m5_size_gb']
    yaxis_vals = mous_db[dq_val]
    
    if dq_val == 'm5_dq_if_non_ima_gb':
        ylabel = 'Log WSU MS5 Non-Imaging DQ (GB)'
    elif dq_val == 'm5_dq_if_ima_gb':
        ylabel = 'Log WSU MS5 Imaging DQ (GB)'
    elif dq_val == 'm5_dq_gb':
        ylabel = 'Log WSU M5 DQ (GB)'
    else:
        ylabel = dq_val

    xaxis_bins = np.arange(np.nanmin(np.log10(xaxis_vals)), np.nanmax(np.log10(xaxis_vals)), 0.2)
    yaxis_bins = np.arange(np.nanmin(np.log10(yaxis_vals)), np.nanmax(np.log10(yaxis_vals)), 0.2)
    
    # set up plot region
    fig,ax = plt.subplots(1)
    
    ## maybe make this a 2d colorized plot
    ax.hist2d(np.log10(xaxis_vals),np.log10(yaxis_vals),
              bins=[xaxis_bins, yaxis_bins],
              cmap='magma_r')

    ax.set_aspect('equal')
    ax.grid()
    ax.plot(xaxis_bins,xaxis_bins+0,color='black',linewidth=2)
    ax.text(2.7,2.2+0.35,'1x',color='black',size=12,weight='bold',horizontalalignment='left')

    ax.plot(xaxis_bins,xaxis_bins-1,color='black',linewidth=2)
    ax.text(3.7,2.2+0.35,'0.1x',color='black',size=12,weight='bold',horizontalalignment='left')

    ax.plot(xaxis_bins,xaxis_bins-2,color='black',linewidth=2)
    ax.text(4.7,2.2+0.35,'0.01x',color='black',size=12,weight='bold',horizontalalignment='left')

    
    ax.set_xlabel('log WSU MS5 ASDM Size (GB)')
    ax.set_ylabel(ylabel)

    ax.text(0.05,0.9,'per MOUS',transform=ax.transAxes,weight='bold',color='black')

def plot_asdm_vs_dq_points(mous_db, dq_val = 'm5_dq_gb', filename=None,
                           xlim = None, ylim = None):
    '''
    Purpose: Plot ASDM size vs. non-imaging DQ

    Input: mous data base

    Ouput: plot

    TODO:
    * Adapt to new variable names -- DONE
    * distinguish points by array -- DONE

    Date        Programmer      Description of Changes
    --------------------------------------------------
    8/7/2025    A.a. Kepley     Original Code
    '''

      
    if dq_val == 'm5_dq_if_non_ima_gb':
        ylabel = 'Log WSU MS5 Non-Imaging DQ (GB)'
        array_list = ['12m','7m']
        factor_list = [10, 1,0.1,0.01]
    elif dq_val == 'm5_dq_if_ima_gb':
        ylabel = 'Log WSU MS5 Imaging DQ (GB)'
        array_list = ['12m','7m']
        factor_list = [10, 1,0.1,0.01]
    elif dq_val == 'm5_dq_gb':
        ylabel = 'Log WSU M5 DQ (GB)'
        array_list = ['12m','7m','TP']
        factor_list = [10, 1,0.1,0.01]
    else:
        ylabel = dq_val
        array_list = ['12m','7m','TP']

    # set up plot region
    fig,ax = plt.subplots(1, figsize=(8.5,8))
    fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)
    
    ## maybe make this a 2d colorized plot


    for myarray in array_list:
        idx = mous_db['array'] == myarray
        ax.loglog(mous_db['m5_size_gb'][idx],mous_db[dq_val][idx],
                  marker='.',alpha=0.2,linestyle='none',markeredgecolor='none',
                  label=myarray, markerfacecolor=array_color[myarray])

        


    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    ax.set_aspect('equal')
    ax.grid()
        
    (xmin, xmax) = ax.get_xlim()

    for myfactor in factor_list:

        
        if myfactor >= 1:
            myline = '--'
            mylinecolor='black'
        else:
            myline = ':'
            mylinecolor = 'black'
            
        model_xvals = np.logspace(np.log10(xmin),np.log10(xmax),1000)
        model_yvals = model_xvals * myfactor
        
        ax.loglog(model_xvals, model_yvals, color=mylinecolor,linestyle=myline)
        ax.text(0.09,0.06*myfactor,'{:4s}'.format(str(myfactor)),
                color='black',size=10,weight='bold',horizontalalignment='right')

        
    ax.set_xlabel('WSU MS5 Total ASDM Size (GB)')
    ax.set_ylabel(ylabel)


    ax.text(0.05,0.9,'per MOUS',transform=ax.transAxes,weight='bold',color='black')


    mylegend = ax.legend(loc='lower right')

    for handle in mylegend.legendHandles:
        handle.set_alpha(1.0)

    
    if filename:
        fig.savefig(filename,facecolor='white',edgecolor='white',dpi=600)


def plot_db1_vs_db2_sample_datarate(db1, db2,
                                    quant1='wsu_datarate_goal_stepped2_typical',
                                    quant2='m5_data_rate_gb_per_s',
                                    if_only=True,
                                    filename=None):
    '''
    
    Purpose: Plot cumulative distribution of sample 1 vs. sample 2. If they are the
    same the curve should be the same.

    Input:

    db1: first database

    db2: second database

    quant1: quantity to plot for database 1

    quant1: quantity to plot for database 2

    Date        Programmer      Description of Changes
    ----------------------------------------------------
    8/11/2025   A.A. Kepley     Original Code       
    '''

    mybins = 500
    
    (vals1, bins1, patches1) = plt.hist(np.log10(db1[quant1].value),
                                        cumulative=-1,histtype='step',
                                        bins=mybins,
                                        log=True,
                                        density=True,
                                        linewidth=2,
                                        label='Cycle 7 and 8 sample')


    (vals2, bins2, patches2) = plt.hist(np.log10(db2[quant2]),
                                        cumulative=-1,histtype='step',
                                        bins=bins1,
                                        log=True,
                                        density=True,
                                        linewidth=2,
                                        label='Cycle 9 and 10 sample')


    plt.axhline(0.1,color='gray',linestyle=':')
    plt.text(-3,0.1,'10% larger')

    plt.axhline(0.05,color='gray',linestyle=':')
    plt.text(-3,0.05,'5% larger')

    plt.axhline(0.01,color='gray',linestyle=':')
    plt.text(-3,0.01,'1% larger')

    
    plt.grid(which='both',axis='both',linestyle=':')


    #plt.axvline(np.log10(1.98),color='gray', linestyle='-.', linewidth=2)
    #plt.text(np.log10(1.98)-0.05, 0.6, '1.98\nGB/s',color='gray',horizontalalignment='right')

    plt.axvline(np.log10(3.95),color='gray', linestyle='-', linewidth=2)
    plt.text(np.log10(3.95)+0.03, 0.15, '3.95\nGB/s',color='gray')

    plt.xlabel('Data Rate (GB/s)')
    plt.ylabel('Fraction of Larger Data')

    if if_only:
        plt.title('IF Only')
    
    locs, labels = plt.xticks()

    newlabels = ['10$^{{ {:2.0f} }}$'.format(val) for val in locs]
    plt.xticks(locs[1:],newlabels[1:])

    
    plt.legend(loc='lower left')
    
    if filename:
        plt.savefig(filename,facecolor='white',edgecolor='white',dpi=600)

def plot_db1_vs_db2_sample_datavol(db1, db2,
                                   quant1='wsu_datavol_goal_stepped2_typical_total',
                                   quant2='m5_size_gb',
                                   if_only = True,
                                   filename=None):
    '''
    
    Purpose: Plot cumulative distribution of sample 1 vs. sample 2. If they are the
    same the curve should be the same.

    Input:

    db1: first database

    db2: second database

    quant1: quantity to plot for database 1

    quant1: quantity to plot for database 2

    Date        Programmer      Description of Changes
    ----------------------------------------------------
    8/11/2025   A.A. Kepley     Original Code       
    '''

    mybins = 500
    
    (vals1, bins1, patches1) = plt.hist(np.log10(db1[quant1].to('TB').value),
                                        cumulative=-1,histtype='step',
                                        bins=mybins,
                                        log=True,
                                        density=True,
                                        linewidth=2,
                                        label='Cycle 7 and 8')


    (vals2, bins2, patches2) = plt.hist(np.log10(db2[quant2]/1000.0),
                                        cumulative=-1,histtype='step',
                                        bins=bins1,
                                        log=True,
                                        density=True,
                                        linewidth=2,
                                        label='Cycle 9 and 10')


    plt.axhline(0.1,color='gray',linestyle=':')
    plt.text(-3,0.1,'10% larger')

    plt.axhline(0.05,color='gray',linestyle=':')
    plt.text(-3,0.05,'5% larger')

    plt.axhline(0.01,color='gray',linestyle=':')
    plt.text(-3,0.01,'1% larger')

    
    plt.grid(which='both',axis='both',linestyle=':')


    plt.xlabel('Total ASDM Size per MOUS (TB)')
    plt.ylabel('Fraction of Larger Data')

    if if_only:
        plt.title('IF Only')
    
    locs, labels = plt.xticks()

    newlabels = ['10$^{{ {:2.0f} }}$'.format(val) for val in locs]
    plt.xticks(locs[1:],newlabels[1:])

    
    plt.legend()

    #plt.clf()
    
    #plt.plot(vals1,vals2)

    #xvals = np.linspace(np.min(vals1),np.max(vals1),10)
    #yvals = xvals
    
    #plt.plot(xvals, yvals, color='black',linestyle=':')

    if filename:
        plt.savefig(filename,facecolor='white',edgecolor='white',dpi=600)



def plot_l80_vs_nsrc(mous_db, dq_val='m5_dq_ratio'):
    '''
    Purpose: Plot L80 vs. Nsource to show where larger points lie

    Input: mous database

    Output: plot

    Date        Programmer      Description of Changes
    ----------------------------------------------------
    8/13/2025   A.A. Kepley     Original Code
    '''


    array_list = ['12m','7m']
    color_list =['blue','orange']
    
    # set up plot region
    fig,ax = plt.subplots(1, figsize=(8,6))
    
    ## maybe make this a 2d colorized plot
    
    for (myarray ,mycolor) in zip(array_list,color_list):
        idx_lt1 = (mous_db['array'] == myarray) & (mous_db[dq_val] < 1.0)
        idx_gt1 = (mous_db['array'] == myarray) & (mous_db[dq_val] >= 1.0)
        
        ax.loglog(mous_db['l80_avg_m'][idx_lt1],mous_db['n_sources'][idx_lt1],
                  marker='.',alpha=0.15,markeredgecolor='none', linestyle='none',
                  label=myarray)


        ax.loglog(mous_db['l80_avg_m'][idx_gt1],mous_db['n_sources'][idx_gt1],
                  marker='o',alpha=0.3,markeredgecolor='none',linestyle='none',
                  label=myarray,color=mycolor)

    ax.legend()
    ax.set_xlabel('L80 (m)')
    ax.set_ylabel('Number of sources')

def plot_l80_vs_nchan(mous_db, dq_val='m5_dq_ratio',
                      array_list=['12m','7m','TP'],
                      filename=None):
    '''
    Purpose: Plot L80 vs. Number of channels to show where larger points lie

    Input: mous database

    Output: plot

    Date        Programmer      Description of Changes
    ----------------------------------------------------
    8/13/2025   A.A. Kepley     Original Code
    '''
    
    # set up plot region
    fig,ax = plt.subplots(1, figsize=(8,6))
    
    ## maybe make this a 2d colorized plot
    
    for myarray in array_list:
        idx_lt = (mous_db['array'] == myarray) & (mous_db[dq_val] < 1.0)
        idx_gt = (mous_db['array'] == myarray) & (mous_db[dq_val] >= 1.0)

    
        ax.scatter(mous_db['l80_avg_m'][idx_lt],mous_db['m5_channel_number'][idx_lt],
                   marker='o',alpha=0.5,edgecolor='none',
                   label=myarray + ', DQ/ASDM < 1', s=mous_db['n_sources'][idx_lt]+5.0,
                   color=array_color[myarray])

        
        ax.scatter(mous_db['l80_avg_m'][idx_gt],mous_db['m5_channel_number'][idx_gt],
                   marker='o',alpha=0.5,edgecolor='none',
                   label=myarray+ ', DQ/ASDM >= 1', s=mous_db['n_sources'][idx_gt]+5.0,
                   color='firebrick')


    ax.set_xscale('log')
    ax.set_yscale('log')
        
    mylegend = ax.legend()

    for handle in mylegend.legendHandles:
        handle.set_alpha(1.0)
    
    ax.set_xlabel('L80 (m)')
    ax.set_ylabel('$N_{chan}$ (WSU)')
    

    ax.text(0.1,0.05,'Larger points = \nmore sources per MOUS',transform=ax.transAxes,weight='bold',color='black')

    if filename:
        plt.savefig(filename,facecolor='white',edgecolor='white',dpi=600)


    


def plot_dq_ratio_vs_l80(mous_db, dq_val='m5_dq_ratio',
                         array_list=['12m','7m','TP'],                        
                         filename=None):
    '''
    Purpose: Plot DQ ratio vs. l80


    Input: mous database


    Output: plot

    Date        Programmer      Description of Changes
    ---------------------------------------------------
    8/13/2025   A.A. Kepley     Original Code
    '''

    if dq_val == 'm5_dq_ratio':
        ylabel = 'WSU DQ / WSU ASDM'
    elif dq_val == 'm5_dq_if_ima_ratio':
        ylabel = 'WSU Imaging DQ / WSU ASDM'
    elif dq_val == 'm5_dq_gb':
        ylabel = 'WSU DQ Size (GB)'
    else:
        ylabel = dq_val
    
    
    # set up plot region
    fig,ax = plt.subplots(1, figsize=(8,6))

    
    for myarray in array_list:
        idx = (mous_db['array'] == myarray) 

        #if myarray == 'TP':
        #    parts = ax.violinplot([mous_db[dq_val][idx]], positions=[12],
        #                  widths=[10])
        #    for pc in parts['bodies']:
        #        pc.set_facecolor(array_color[myarray])

        #else:
        ax.scatter(mous_db['l80_avg_m'][idx],mous_db[dq_val][idx],
                   marker='o',alpha=0.5,edgecolor='none', 
                   label=myarray, color=array_color[myarray],
                   s=mous_db['n_sources'][idx]+2.0)
        


    ax.set_xscale('log')
    ax.set_yscale('log')

    mylegend = ax.legend()

    for handle in mylegend.legendHandles:
        handle.set_alpha(1.0)

    ax.set_xlabel('L80 (m)')
    ax.set_ylabel(ylabel)

    if dq_val != 'm5_dq_gb':        
        ax.axhline(1.0, color='black',linestyle='--')

    ax.text(0.67,0.05,'Larger points = \nmore sources per MOUS',transform=ax.transAxes,weight='bold',color='black')

    if filename:
        fig.savefig(filename,facecolor='white',edgecolor='white',dpi=600)

    

def plot_dq_ratio_vs_nsrc(mous_db, dq_val='m5_dq_ratio',
                          filename=None):
    '''
    Purpose: Plot DQ ratio vs. l80


    Input: mous database


    Output: plot

    Date        Programmer      Description of Changes
    ---------------------------------------------------
    8/13/2025   A.A. Kepley     Original Code
    '''


    
    array_list = ['7m']
    color_list =['orange']
    
    # set up plot region
    fig,ax = plt.subplots(1, figsize=(8,6))
    
    ## maybe make this a 2d colorized plot
    
    for (myarray,mycolor) in zip(array_list,color_list):
        idx_lt = (mous_db['array'] == myarray) & (mous_db['m5_channel_width_kms'] < 10)
        idx_gt = (mous_db['array'] == myarray) & (mous_db['m5_channel_width_kms'] >= 10)

        
        
        ax.loglog(mous_db['n_sources'][idx_lt],mous_db[dq_val][idx_lt],
                  marker='.',alpha=0.5,markeredgecolor='none', linestyle='none',
                  label=myarray+ ', < 10km/s',color=mycolor)


        ax.loglog(mous_db['n_sources'][idx_gt],mous_db[dq_val][idx_gt],
                  marker='o',alpha=0.5,markeredgecolor='none', linestyle='none',
                  label=myarray+', >= 10km/s',color=mycolor)
                


    ax.legend()

    ax.set_xlabel('Number of Sources')
    ax.set_ylabel(dq_val)

    ax.axhline(1.0, color='red', linestyle='--')


def plot_dq_ratio_vs_nchan(mous_db, dq_val='m5_dq_ratio', array_list=['7m'],
                           filename=None):
    '''
    Purpose: Plot DQ ratio vs. l80


    Input: mous database


    Output: plot

    Date        Programmer      Description of Changes
    ---------------------------------------------------
    8/13/2025   A.A. Kepley     Original Code
    '''

    if dq_val == 'm5_dq_ratio':
        ylabel = 'WSU DQ / WSU ASDM'
    elif dq_val == 'm5_dq_if_ima_ratio':
        ylabel = 'WSU Imaging DQ / WSU ASDM'
    else:
        ylabel = dq_val
    
    # set up plot region
    fig,ax = plt.subplots(1, figsize=(8,6))
    
    ## maybe make this a 2d colorized plot
    
    for myarray in array_list:
        idx = mous_db['array'] == myarray
        
        ax.scatter(mous_db['m5_channel_number'][idx],mous_db[dq_val][idx],
                  marker='o',alpha=0.5,edgecolor='none', 
                   label=myarray,color=array_color[myarray],
                   s=mous_db['n_sources'][idx]+2.0)



    ax.set_xscale('log')
    ax.set_yscale('log')

    if array_list == ['TP']:
        mylegend = ax.legend(loc='upper left')
    else:
        mylegend = ax.legend(loc='upper right')
        
    for handle in mylegend.legendHandles:
        handle.set_alpha(1.0)
 
    ax.set_xlabel('Number of channels')
    ax.set_ylabel(ylabel)

    ax.axhline(1.0, color='black', linestyle='--')

    ax.text(0.1,0.05,'Larger points = \nmore sources per MOUS',transform=ax.transAxes,weight='bold',color='black')

    if filename:
        fig.savefig(filename,facecolor='white',edgecolor='white',dpi=600)
