import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mlines
import matplotlib.font_manager as font_manager
import pandas as pd
import argparse
import sys

#global parameters
font_prop_title = font_manager.FontProperties(size=40)
font_prop_axis = font_manager.FontProperties(size=30)
font_prop_legend = font_manager.FontProperties(size=26)

def produce_plot():
    
    #argument parser
    parser = argparse.ArgumentParser(description='plot single node thread scaling')
    parser.add_argument('--inputfile', type=str, nargs=1, help='specify the input file in csv format')
    parser.add_argument('--outputfile', type=str, nargs=1, default='single_node_thread_scaling.eps', help='specify the output file in eps format')
    args = parser.parse_args()
    
    if not args.inputfile:
        parser.print_help()
        sys.exit()

    #get the filename
    filename=args.inputfile[0]

    #load file into dataframe
    singlenodedf=pd.DataFrame.from_csv(filename,index_col=False)
    singlenodedf.sort_values(by="ncores",inplace=True)
    
    #color list
    colors=['r','b','g','y']
    
    #plot
    fig=plt.figure(num=None, figsize=(20, 34), facecolor='w', edgecolor='k')
    ax=plt.subplot(2,1,1)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontproperties(font_prop_axis)
        label.set_fontsize(26)
    
    #set padding
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(15)
    for tick in ax.yaxis.get_major_ticks():
        tick.set_pad(10)
    
    #loop over hyperthreads
    for nhtid,nht in enumerate(singlenodedf["nht"].unique()):
        #ranges:
        Xhsw=list(singlenodedf["ncores"].ix[ (singlenodedf.arch=="hsw") & (singlenodedf.nht==nht) ].unique())
        Yhsw=list(singlenodedf["time"].ix[ (singlenodedf.arch=="hsw") & (singlenodedf.nht==nht) ])
        Xknl=list(singlenodedf["ncores"].ix[ (singlenodedf.arch=="knl") & (singlenodedf.nht==nht) ].unique())
        Yknl=list(singlenodedf["time"].ix[ (singlenodedf.arch=="knl") & (singlenodedf.nht==nht) ])
        #plt.xlim(0,np.max(Xhsw+Xknl)+1)
        #plt.ylim(np.min(Yknl+Yhsw)*0.9,np.max(Yknl+Yhsw)*1.1)
        
        #connecting line
        lineknl,=plt.plot(Xknl,Yknl, '-', linewidth=6,color=colors[nhtid])
        linehsw,=plt.plot(Xhsw,Yhsw, '-', linewidth=6,color=colors[nhtid+2])
        
        #datapoints
        pointsknl=plt.scatter(Xknl,Yknl,marker='o',s=400,color=colors[nhtid],label="KNL, "+str(nht)+"ht")
        pointshsw=plt.scatter(Xhsw,Yhsw,marker='s',s=400,color=colors[nhtid+2],label="HSW, "+str(nht)+"ht")
    
    #logscale
    plt.xscale('log',basex=2)
    plt.yscale('log',basey=2)
    
    #label the plot
    plt.legend(loc='upper right', prop=font_prop_legend, scatterpoints=1)
    plt.xlabel('#cores',fontproperties=font_prop_axis)
    plt.ylabel('time to solution [s]',fontproperties=font_prop_axis)
    fig.savefig(args.outputfile[0],bbox_inches='tight',format='eps')

def main():
    produce_plot()

if __name__ == "__main__":
    main()
