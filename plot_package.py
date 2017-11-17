"""
Package for organizing and maintaining plotting funtions using some customized format
"""
import numpy as np

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

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
        try:
            plt.plot(x[i], y[i], alpha=alpha_value, linewidth=2, linestyle=linetype[i/6], color=colors[i], label="%s"%label[i], marker="o")
        except:
            print "FAILURE:"
            print i
            print np.shape(x[i])
            print np.shape(y[i])
            print x[i]
            print y[i]
            raise
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

def determine_number_states(check_array):
    try:
        retain = np.shape(check_array)[1]
    except:
        retain = 1 # thsi is terrible code, but it works

    return retain

def plot_timescales(lags_list, timescales_list, labels, save_name, time_step=1., mark_lines=None, units=None, manual_x_max=None, confidence_lags=None, confidence_interval=None):
    """ Plot timescales versus lag time with the nice gray area at the bottom

    """

    colors = ["blue", "red", "orange", "green", "magenta", "yellow", "cyan"]

    fig = plt.figure()
    assert len(lags_list) == len(timescales_list)
    assert len(lags_list) == len(labels)

    all_linestyles = ["-", "--", "."]

    if confidence_interval is not None:
        this_line = all_linestyles[0]
        lower = confidence_interval[0]
        upper = confidence_interval[1]
        lags = lags_list[count]
        retain = determine_number_states(lower)
        if retain == 1:
            plt.fill_between(confidence_lags, lower, upper, where=yhigh>=ylow, interpolate=True, color=colors[0],alpha=0.5)
        else:
            for i in range(retain):
                plt.fill_between(confidence_lags, lower[:,i], upper[:,i], where=yhigh>=ylow, interpolate=True, color=colors[i], alpha=0.5)

    label_lines = []
    max_time = np.max(timescales_list[0])
    min_time = np.min(timescales_list[0])
    for count in range(len(lags_list)):
        all_timescales = timescales_list[count]
        lags = lags_list[count]
        this_line = all_linestyles[count%len(all_linestyles)]
        retain = determine_number_states(all_timescales)
        if retain == 1:
            plt.semilogy(lags*time_step, all_timescales*time_step, marker="o", linestyle=this_line, linewidth=2.0, color=colors[0])
        else:
            for i in range(retain):
                plt.semilogy(lags*time_step, all_timescales[:,i]*time_step, marker="o", linestyle=this_line, linewidth=2.0, color=colors[i])
        max_timescales = np.max(all_timescales)
        min_timescales = np.min(all_timescales)
        if max_time < max_timescales:
            max_time = max_timescales
        if min_time > min_timescales:
            min_time = min_timescales

        label_lines.append(mlines.Line2D([0],[0], linewidth=2, color="black", linestyle=this_line))

    miny = min_time * time_step * 0.1
    maxy = max_time * time_step * 3.

    if manual_x_max is not None:
        maxx = manual_x_max
    else:
        maxx = np.max(np.array(lags)) * time_step

    if mark_lines is not None:
        yrange = [miny, maxy]
        for step in mark_lines:
            plt.plot([step*time_step, step*time_step], yrange, label="%d Time Steps" % step, linestyle="--")

    # make the legend
    fig.axes[0].legend(label_lines, labels, loc="lower right")

    xspacing = np.arange(0, maxx*1.1, maxx / 1000. )
    num_fill = np.shape(xspacing)[0]
    yhigh = xspacing
    ylow = np.ones(num_fill) * miny
    plt.fill_between(xspacing, ylow, yhigh, where=yhigh>=ylow, interpolate=True, color="grey")

    plt.axis([0, maxx, miny, maxy])

    if units is None:
        units = "steps"
    plt.xlabel("Lag Time (%s)" % units)
    plt.ylabel("Timescales (%s)" % units)
    plt.savefig(save_name)
