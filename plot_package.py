"""
Package for organizing and maintaining plotting funtions using some customized format
"""
import numpy as np

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt


def return_max(cmax, nmax):
    if cmax < nmax:
        return nmax
    else:
        return cmax

def plot_simple(x, y, label, title, xaxis_label, yaxis_label, axis=None, save_file=None, show=True, reference=False):
    """plot_simple is for a simple x-y plot with the dots connected.  """
    #set the save_file name to the same as the title name
    if save_file == None:
        save_file = title
    #specify generic color order for use
    if reference:
        colors = ["k","b","g","r","c","m","y","b","g","r","c","m","y","b","g","r","c","m","y","b","g","r","c","m","y"]
    else:
        colors = ["b","g","r","c","m","y","b","g","r","c","m","y","b","g","r","c","m","y","b","g","r","c","m","y"]
    linetype = ["-", "--",":"]

    plt.figure()
    maxvalue = 0.0
    maxcenter = 0.0

    #plot it!
    for i in range(np.shape(label)[0]):
        if reference and i==0:
            alpha_value = 1.0
        else:
            alpha_value = 0.75
        plt.plot(x[i], y[i], alpha=alpha_value, linewidth=2, linestyle=linetype[i/6], color=colors[i], label="%s"%label[i], marker="o")
        maxvalue = return_max(maxvalue, np.max(y[i]))
        maxcenter = return_max(maxcenter, np.max(x[i]))

    #specify the axis size, labels, etc.
    if axis == None:
        plt.axis([0,int(maxcenter*1.5)+1, 0, maxvalue*1.2],fontsize=20)
    else:
        plt.axis(axis)
    plt.legend()
    plt.xlabel(xaxis_label, fontsize=25)
    plt.ylabel(yaxis_label,fontsize=25)
    plt.title(title, fontsize=25)
    plt.tight_layout()

    plt.savefig("%s.png"%save_file)
    if show:
        plt.show()



def plot_error_bar(x, y, label, title, xaxis_label, yaxis_label, xerr=None, yerr=None, axis=None):
    """plot_simple is for a simple x-y plot with the dots connected.  """

    #specify generic color order for use
    colors = ["b","g","r","c","m","y","b","g","r","c","m","y","b","g","r","c","m","y","b","g","r","c","m","y"]
    linetype = ["-", "--",":"]

    #initialize stuff, specify starting maxvalue and maxcenter for determining axis size.
    plt.figure()
    maxvalue = 0.0
    maxcenter = 0.0

    #plot it!
    for i in range(np.shape(label)[0]):
        #plt.errorbar(x[i], y[i], xerr=xerr, yerr=yerr, linewidth=2, alpha=0.75, linestyle=linetype[i/6], color=colors[i], label="%s"%label[i], marker="o")
        plt.errorbar(x[i], y[i], xerr=xerr, yerr=yerr[i], linewidth=2, alpha=0.75, linestyle=linetype[i/6], color=colors[i], label="%s"%label[i])
        maxvalue = return_max(maxvalue, np.max(y[i]))
        maxcenter = return_max(maxcenter, np.max(x[i]))

    #specify the axis size, labels, etc.
    if axis == None:
        plt.axis([0,int(maxcenter*1.5)+1, 0, maxvalue*1.2],fontsize=20)
    else:
        plt.axis(axis)
    plt.legend()
    plt.xlabel(xaxis_label, fontsize=25)
    plt.ylabel(yaxis_label,fontsize=25)
    plt.title(title, fontsize=25)

    plt.show()

    plt.savefig("%s.png"%title)


def multiplot_stacked(x, y, labels, xaxis_label, yaxis_label, axis=None, save_file=None, fig_size=(4,8)):
    nplots = len(x)

    fig, axarr = plt.subplots(nplots, sharex=True)
    fig.set_size_inches(fig_size[0],fig_size[1])
    xaxis_max = 0
    xaxis_min = 0
    for idx in range(nplots):
        if np.max(x[idx]) > xaxis_max:
            xaxis_max = np.max(x[idx])
        if np.min(x[idx]) < xaxis_min:
            xaxis_min = np.min(x[idx])

    for idx in range(nplots):
        this_y = y[idx]
        this_x = x[idx]
        maxval = np.max(this_y)
        minval = 0
        if np.min(this_y) < minval:
            minval = np.min(this_y)
        intervals = (maxval - minval) / 4.
        yticks = np.arange(minval, maxval+intervals, intervals)
        axarr[idx].axis([xaxis_min, xaxis_max, minval, maxval*1.1])
        axarr[idx].yaxis.set_ticks(yticks)
        axarr[idx].plot(this_x, this_y)
        axarr[idx].set_ylabel(yaxis_label)
        axarr[idx].set_title(labels[idx])

    axarr[-1].set_xlabel(xaxis_label)

    fig.tight_layout()

    if not save_file is None:
        plt.savefig(save_file)
