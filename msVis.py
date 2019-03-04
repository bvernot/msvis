from __future__ import division
import sys, os, random, itertools
import fileinput
from operator import itemgetter
import argparse
#from numpy.random import binomial
#import tables
import re
import math, matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy

###############################################################################################################################
##  QUESTIONS FOR HUDSON             ##########################################################################################
###############################################################################################################################
#
# 1) Does growth start from a size of N or N0? In the manual, it looks like it's N0
# 2) Is the growth equation exp(-alpha*time) or exp(-alpha*time_from_start_of_growth)?  In the manual, looks like former, but 
#    that doesn't make a lot of sense to me.
# 3) Is there a difference between -en .1 1 10 -eg .1 1 5 and -eg .1 1 5 -en .1 1 10 ?  i.e., does the order matter?
# 4) Do you have advice for simulating a simultaneous split of one population into many populations (pastward), such that
#    the present day population is made up of equal parts of each of the many past populations?

## this is the EA AA ms command
# ms 5 5 -I 2 1 1 -t 5 -r 5 10000 -n 1 58.002735978 -n 2 70.041039672 -eg 0 1 482.46 -eg 0 2 570.18 -em 0 1 2 0.7310 -em 0 2 1 0.7310 -eg 0.006997264 1 0 -eg 0.006997264 2 89.7668 -en 0.006997264 1 1.98002736 -en 0.031463748 2 0.141176471 -en 0.03146375 2 0.254582763 -em 0.03146375 1 2 4.386 -em 0.03146375 2 1 4.386 -em 0.069767442 1 2 0 -em 0.069767442 2 1 0 -ej 0.069767442 2 1 -en 0.069767442 1 1.98002736 -en 0.20246238 1 1

###############################################################################################################################
##  SET UP COMMAND LINE ARGUMENTS  ############################################################################################
###############################################################################################################################

class ValidateTimedParam(argparse.Action):
    def __call__(self, parser, args, values, option_string = None):
        if 'timed_count' not in args: setattr(args, 'timed_count', 0)
        ## get previous timed command of this type
        if getattr(args, self.dest) == None: setattr(args, self.dest, [])
        previous = getattr(args, self.dest)
        
        if self.dest == 'em':
            vals = (int(values[1]), int(values[2]), float(values[3]))
        elif self.dest == 'es':
            vals = (int(values[1]), float(values[2]))
        elif self.dest == 'ej':
            vals = (int(values[1]), int(values[2]))
        elif self.dest == 'en':
            vals = (int(values[1]), float(values[2]))
        elif self.dest == 'eN':
            vals = float(values[1])
        elif self.dest == 'eG':
            vals = float(values[1])
        elif self.dest == 'eg':
            vals = (int(values[1]), float(values[2]))
        elif self.dest == 'ema':
            vals = values[1:]
        else:
            print 'unsupported timed command', self.dest
            sys.exit(-1)
            pass
            
        args.timed_count += 1
        new_cmd = {'card' : args.timed_count,
                   'cmd' : self.dest,
                   'time' : float(values[0]),
                   'vals' : vals}
        
        previous.append(new_cmd)
        setattr(args, self.dest, previous)

        ## get all previous timed commands
        if 'ecmd' not in args: setattr(args, 'ecmd', [])
        args.ecmd.append(new_cmd)

        return
    pass


parser = argparse.ArgumentParser(description='Visualize ms command line arguments.')

ms_group = parser.add_argument_group('ms arguments')

ms_group.add_argument('nsam', type=int)
ms_group.add_argument('nreps', type=int)
theta_group = ms_group.add_mutually_exclusive_group(required=True)
theta_group.add_argument('-t', metavar='theta', dest='theta', type=float)
theta_group.add_argument('-s', metavar='segsites', dest='segsites', type=int)

ms_group.add_argument('-I', type = float, nargs='+', default = None)

ms_group.add_argument('-ej', nargs = 3, action = ValidateTimedParam, default = [])
ms_group.add_argument('-es', nargs = 3, action = ValidateTimedParam, default = [])
ms_group.add_argument('-em', nargs = 4, action = ValidateTimedParam, default = [])
ms_group.add_argument('-ema', nargs = '+', action = ValidateTimedParam, default = [])
ms_group.add_argument('-en', nargs = 3, action = ValidateTimedParam, default = [])
ms_group.add_argument('-eN', nargs = 2, action = ValidateTimedParam, default = [])
ms_group.add_argument('-eg', nargs = 3, action = ValidateTimedParam, default = [])
ms_group.add_argument('-eG', nargs = 2, action = ValidateTimedParam, default = [])

ms_group.add_argument('-m', nargs = 3, type = float, action = 'append', default = [])
ms_group.add_argument('-ma', nargs = '+', default = None)
ms_group.add_argument('-n', nargs = 2, type = float, action = 'append', default = [])
ms_group.add_argument('-r', nargs = 2, type = float)
ms_group.add_argument('-c', nargs = 2, type = float)
ms_group.add_argument('-G', type = float, default = None)
ms_group.add_argument('-g', nargs = 2, type = float, action = 'append', default = [])

ms_group.add_argument('-T', action = 'store_true')
ms_group.add_argument('-L', action = 'store_true')

ms_group.add_argument('-seeds', nargs = '+', type = int)

display_group = parser.add_argument_group('display arguments')

display_group.add_argument('-N0', type = int, default = 10000, help = 'N0 (used to calculate population sizes and time in years).')
display_group.add_argument('-gen', type = int, default = 25, help = 'Generation time (used to calculate time in years).')
display_group.add_argument('-curve', '--curve-points', type = int, default = 100, help = 'The number of points plotted for each curve (a higher number makes smoother curves, but takes longer to plot).')
display_group.add_argument('-range', type = float, default = None, nargs = 2, help = 'Range of time to plot (in scaled time; real time if -years flag is used).')
display_group.add_argument('-random', action = 'store_const', metavar = ['start_pops', 'pop_add', 'pop_subtract', 'pop_change_prob', 'n_add', 'n_freq'], const = [5, .1, .05, .1, -.1, .1], help = 'Make a random ms figure (doesn\'t work yet).')
display_group.add_argument('-indices', type = int, default = None, nargs = '+', help = 'Plot order of populations.')
display_group.add_argument('-pops', '--population-names', type = str, default = None, nargs = '+', help = 'Population names (in plot order).')
display_group.add_argument('-bw', action = 'store_true', help = 'Plot in black and white.')
display_group.add_argument('-logn', action = 'store_true', help = 'Plot population sizes as log of real size (uses N0, which is set by default to 10k).')
display_group.add_argument('-logt', action = 'store_true', help = 'Plot y axis as log of years (not yet working).')
display_group.add_argument('-nomig', action = 'store_true', help = 'Do not plot migration.')
display_group.add_argument('-years', action = 'store_true', help = 'Plot y axis as years instead of scaled time.')
display_group.add_argument('-colors', default = [ '#bd2309', '#bbb12d', '#1480fa', '#14fa2f',
                                                  '#faf214', '#2edfea', '#ea2ec4', '#ea2e40', '#cdcdcd',
                                                  '#577a4d', '#2e46c0', '#f59422', '#219774', '#8086d9' ], nargs = '+', help = 'A list of colors to use when plotting populations.  If the number of populations is larger than the number of colors, colors will be repeated.')
display_group.add_argument('-v', '--version', action = 'store_true', help = 'Print the version and exit.')
display_group.add_argument('-o', '--output-file', help = 'Output file name.', default = 'ms_plot.pdf')
display_group.add_argument('-width', type = float, default = 15, help = 'Width of plot (in inches).')
display_group.add_argument('-height', type = float, default = 15, help = 'Height of plot (in inches).')

args = parser.parse_args()
print args

if args.version:
    print
    print 'msVis version 0.1.1'
    print 'Benjamin Vernot'
    print
    sys.exit(-1)
    pass

## make sure that if no timed commands are given, we at least have an empty list of timed commands (e commands)
if 'ecmd' not in args: setattr(args, 'ecmd', [])

if args.years and args.range != None:
    args.range = [t / 4 / args.N0 / args.gen for t in args.range]
    pass

###############################################################################################################################
##  RANDOMLY GENERATE AN MS COMMAND  ##########################################################################################
###############################################################################################################################

## if we're generating a random ms command, do it now, then re-call parse_args()
if args.random != None:
    print 'random', args.random
    cmd = '1 1 -t 5 -I %d ' % args.random[0]
    cmd += ' '.join('1' * args.random[0])
    max_pops = len(args.random[0])
    for t in [i/100 for i in range(100)]:
        pass
    print cmd
    args = parser.parse_args(cmd)
    sys.exit()
    pass


###############################################################################################################################
##  SET UP INITIAL POPULATION STRUCTURE AND GLOBAL MIGRATION RATE (INTERPRET -I ARGUMENT)  ####################################
###############################################################################################################################

## set up population structure
initial_pops = [1]
indices = [1]
if args.I != None:
    if args.I[0] + 1 != len(args.I) and args.I[0] + 2 != len(args.I):
        print 'incorrect number of args for -I'
        print args.I
        print args.I[0], len(args.I)
        sys.exit(-1)
        pass
    
    if int(args.I[0]) != args.I[0]:
        print 'First argument to -I must be an integer'
        print args.I
        sys.exit(-1)
        pass

    if args.I[0] + 1 == len(args.I):
        args.I = [int(i) for i in args.I]
    elif args.I[0] + 2 == len(args.I):
        args.I = [int(i) for i in args.I[:-1]] + [args.I[-1]]
        pass

    print 'I', args.I

    initial_pops = [i+1 for i in range(args.I[0])]
    pop_chrs = args.I[1:len(initial_pops)+1] # (not needed yet)
    if sum(pop_chrs) != args.nsam:
        print "Sample counts for -I argument must equal the number of samples generated!"
        print sum(pop_chrs), args.nsam
        sys.exit(-1)
        pass

    if args.I[0] + 2 == len(args.I):
        global_mig = args.I[-1]
    else:
        global_mig = 0
        pass
    ## allow the user to set a population order (or just plot them in numeric order)
    if args.indices == None:
        indices = range(1, args.I[0]+1)
    else:
        indices = args.indices
        pass

    pass


seq_pops = {0:set(initial_pops)}
print seq_pops


###############################################################################################################################
##  SET UP MIGRATION DATA STRUCTURE  ##########################################################################################
###############################################################################################################################

# times are in units of 4N0
mig_matrices = {0:[['x' if i == j else global_mig / 4 / args.N0 / (len(initial_pops)-1) for i in initial_pops] for j in initial_pops]}
if args.ma != None and len(args.ma) != len(initial_pops) ** 2:
    print 'Migration matrix given by -ma is incorrect shape:'
    print 'Contains %d elements, should contain %d' % (len(args.ma), len(initial_pops) ** 2)
    print args.ma
    sys.exit(-1)
    pass
if args.ma != None:
    dim = len(initial_pops)
    print [[i*dim+j for j in range(dim)] for i in range(dim)]
    print args.ma
    mig_matrices = {0: [[float(args.ma[i*dim+j]) / 4 / args.N0 if i != j else args.ma[i*dim+j] for j in range(dim)] for i in range(dim)]}
    pass

print mig_matrices
for i,j,m in args.m:
    mig_matrices[0][int(i)-1][int(j)-1] = m / 4 / args.N0
    pass
print mig_matrices

###############################################################################################################################
##  SET UP POPULATION SIZE / GROWTH DATASTRUCTURES  ###########################################################################
###############################################################################################################################

# pop sizes are in units of N0
# times are in units of 4N0
## CHECK: that no -g precedes a -G
pop_sizes = {0:[1.0 for i in initial_pops]}
pop_growths = {0:[(args.G,0) if args.G != None else None for i in initial_pops]}
for i,n in args.n:
    pop_sizes[0][int(i)-1] = n
    pass
for i,g in args.g:
    pop_growths[0][int(i)-1] = (g,0)
    pass

###############################################################################################################################
##  FUNCTIONS FOR ACCESSING THOSE DATASTRUCTURES  #############################################################################
###############################################################################################################################

def get_pop_list_at_time(t, add_pop = False, modify = True):
    if t not in seq_pops:
        ts  = [k for k in sorted(seq_pops) if k <= t]
        pops = set(seq_pops[ts[-1]])
        if modify: seq_pops[t] = pops
        pass
    else:
        pops = seq_pops[t]
        pass
    if add_pop:
        pops.add(1.0)
        pass
    return pops

def get_pop_sizes_at_time(t, add_pop = False):
    if t not in pop_sizes:
        ts  = [k for k in sorted(pop_sizes) if k <= t]
        pop_sizes[t] = [None for i in pop_sizes[ts[-1]]]
        pop_growths[t] = [None for i in pop_sizes[t]]
        pass
    if add_pop:
        pop_sizes[t].append(1.0)
        pop_growths[t].append(None)
        pass
    return pop_sizes[t]

def get_pop_growths_at_time(t):
    if t not in pop_growths:
        get_pop_sizes_at_time(t)
        pass
    return pop_growths[t]

def get_mig_matrix_at_time(t, add_pop = False):
    if t not in mig_matrices:
        ts  = [k for k in sorted(mig_matrices) if k < t]
        mig_matrices[t] = [list(r) for r in mig_matrices[ts[-1]]]
        pass
    if add_pop:
        j = len(mig_matrices[t]) + 1
        mig2 = [[mig_matrices[t][jj][ii] if jj != j-1 and ii != j-1 else 0 for ii in range(j)] for jj in range(j)]
        mig2[j-1][j-1] = 'x'
        mig_matrices[t] = mig2
        pass
    return mig_matrices[t]

def get_next_size_change(t, i):
    ts  = [k for k in sorted(pop_sizes) if k > t and pop_sizes[k][i-1] != None]
    #print 'getting next size change', t, i, ts, pop_sizes
    if len(ts) == 0:
        return max(pop_sizes)
    return ts[0]

def get_const_size_at_time(t, i, allow_t = True):
    print 'getting previous size change', t, i, pop_sizes
    ts  = [k for k in sorted(pop_sizes) if ((allow_t and k <= t) or k < t) \
               and i <= len(pop_sizes[k]) \
               and pop_sizes[k][i-1] != None]
    pt = ts[-1]
    print '        previous size change', t, i, pop_growths[pt][i-1], pt, ts, pop_sizes
    if pt != t and pop_growths[pt][i-1] != None:
        alpha, start_time = pop_growths[pt][i-1]
        new_size = pop_sizes[pt][i-1] * math.exp(-alpha * (t - start_time))
        print 'calculating const start time from previous growth curve:', alpha, start_time, t, pop_sizes[pt][i-1], new_size
        return new_size
    return pop_sizes[pt][i-1]

def get_next_mig_event(t):
    ts  = [k for k in sorted(mig_matrices) if k > t]
    if len(ts) > 0:
        return ts[0]
    return max(seq_pops)

def get_start_of_mig_event(t, i, j, m):
    #print 'getting start of mig event:', t, i, j, m
    # we're looking for the lowest time k such that the migrations between pops i and j (to and from)
    #  are the same as at time t.  then, time k is when this "migration event" began.
    for k in reversed(sorted(mig_matrices)):
        if k > t: continue
        # don't bother looking if either pop didn't exist at time  k
        if len(mig_matrices[k]) < i or len(mig_matrices[k]) < j: continue
        #print 'k:', k, mig_matrices[k][i-1][j-1], mig_matrices[t][i-1][j-1]
        #print 'k:', k, mig_matrices[k][j-1][i-1], mig_matrices[t][j-1][i-1]
        if mig_matrices[k][i-1][j-1] == mig_matrices[t][i-1][j-1] \
                and mig_matrices[k][j-1][i-1] == mig_matrices[t][j-1][i-1]:
            prior_k = k
            continue
        break
    return prior_k



###############################################################################################################################
##  BUILD UP DATA STRUCTURES TO REPRESENT TIMED EVENTS (MIGRATION, SPLITTING, JOINING, POPSIZE CHANGES, GROWTH, ETC)  #########
###############################################################################################################################

## keep track of the largest population number seen (for adding pops via splits)
max_pop = max(initial_pops)

## progress through timed events in order of 1) time, and 2) cardinality
for event in sorted(args.ecmd, key = itemgetter('time', 'card')):
    print 'starting timed event'
    
    ## every timed event has a time (in  units of 4N0) and a cardinality (the order it was given on the command line)
    ## events are processed starting at time 0, with events of equal time ranked by cardinality

    ## joins
    if event['cmd'] == 'ej':
        print '-ej', event
        t = event['time']
        i,j = event['vals']
        pops = get_pop_list_at_time(t)

        print 'pre join', pops, i
        if i not in pops or j not in pops:
            print 'cannot join populations that do not exist at time %f: %d, %d' % (t, i, j)
            print 'current pops: %s' % pops
            sys.exit(-1)
            pass
        pops.remove(i)
        print 'post join', pops

        # also modify migration matrix to stop migration between i and j
        mig = get_mig_matrix_at_time(t)
        print 'previous mig matrix', mig, mig[i-1][j-1], mig[j-1][i-1]
        mig[i-1][j-1] = 0
        print 'MIGRATION MATRICES MIGHT BE A LITTLE MESSED UP - THE ITH ROW IS MIGRATION *TO* POPULATION I, BUT I\'M SETTING THAT ROW TO ZERO WHEN TWO POPS JOIN (SHOULD BE SETTING THE ITH ITEM IN EVERY ROW?)'
        mig[i-1][:] = [0 if x != i-1 else 'x' for x in range(len(mig[i-1]))]
        mig[j-1][i-1] = 0
        print 'new mig matrix', mig

        # also set the size of pop i to zero
        sizes = get_pop_sizes_at_time(t)
        sizes[i-1] = 0

    elif event['cmd'] == 'es':
        print '-es', event
        t = event['time']
        i,p = event['vals']
        pops = get_pop_list_at_time(t)

        max_pop += 1
        j = max_pop
        event['vals'] = (i,j,p) # save the new pop (j)

        print 'pre split', pops, i,  j
        if i not in pops:
            print 'cannot split a population that does not exist at time %f: %d' % (t, i)
            print 'current pops: %s' % pops
            sys.exit(-1)
            pass

        pops.add(j)
        # keep a list of all pops we've ever seen, so we know their indices
        if args.indices == None:
            indices.insert(indices.index(i)+1, j)
            pass
        print 'post split', pops

        # also modify migration matrix to include this new population
        get_mig_matrix_at_time(t, add_pop = True)
        
        # also add a population to pop_sizes
        get_pop_sizes_at_time(t, add_pop = True)

    elif event['cmd'] == 'en':
        print '-en', event
        t = event['time']
        i,n = event['vals']
        pops = get_pop_sizes_at_time(t)
        pops[i-1] = n
        # also set a seq_pops change point here, just so it looks up the size value
        get_pop_list_at_time(t)

    elif event['cmd'] == 'eg':
        # growth
        print '-eg', event
        t = event['time']
        i,alpha = event['vals']
        growths = get_pop_growths_at_time(t)
        growths[i-1] = (alpha, t)
        pops = get_pop_sizes_at_time(t)
        pops[i-1] = get_const_size_at_time(t,i)
        # also set a seq_pops change point here, just so it looks up the size value
        get_pop_list_at_time(t)

    elif event['cmd'] == 'em':
        print '-em', event
        t = event['time']
        i,j,m = event['vals']
        mig = get_mig_matrix_at_time(t)
        print 'previous mig matrix', mig
        mig[i-1][j-1] =  m / 4 / args.N0

    elif event['cmd'] == 'ema':
        print '-ema', event
        t = event['time']
        mat = event['vals']
        mig = get_mig_matrix_at_time(t)
        print 'previous mig matrix', mig
        mig[i-1][j-1] =  m / 4 / args.N0

    elif event['cmd'] == 'eN':
        # size change for all
        print '-en', event
        t = event['time']
        n = event['vals']
        pops = get_pop_sizes_at_time(t)
        for i in range(len(pops)):
            pops[i] = n
            pass
        # also set a seq_pops change point here, just so it looks up the size value
        get_pop_list_at_time(t)

    elif event['cmd'] == 'eG':
        # growth for all
        print '-eg', event
        t = event['time']
        alpha = event['vals']
        growths = get_pop_growths_at_time(t)
        pops = get_pop_sizes_at_time(t)
        for i in range(len(growths)):
            growths[i] = (alpha, t)
            pops[i] = get_const_size_at_time(t,i)
            pass
        # also set a seq_pops change point here, just so it looks up the size value
        get_pop_list_at_time(t)

    else:
        print 'unsupported timed command', self.dest
        sys.exit(-1)
        pass
    
    print
    pass

# make a last set of pops and migrations, for the top of the graph
the_end_time = max([t for t in seq_pops] + [t for t in mig_matrices]) * 1.1 if args.range == None else args.range[1]
if the_end_time == 0: the_end_time = 1
get_pop_list_at_time(the_end_time)
get_mig_matrix_at_time(the_end_time)
get_pop_sizes_at_time(the_end_time)


###############################################################################################################################
##  SOME BASIC FIGURE DRAWING CODE (COLORS, POPULATION POSITIONS, ETC)  #######################################################
###############################################################################################################################

if args.indices != None and len(args.indices) != max_pop:
    print 'Incorrect number of indices given!  Must match the final number of populations:', max_pop
    print 'indices:', args.indices
    sys.exit(-1)
    pass

if args.population_names != None and len(args.population_names) != max_pop:
    print 'Incorrect number of population names given!  Must match the final number of populations:', max_pop
    print 'Names:', args.population_names
    sys.exit(-1)
    pass

## allow the user to set population names
if args.population_names == None:
    args.population_names = list(indices)
    pass

## set population colors
if args.bw:
    args.colors = ['#000000']
args.colors = args.colors * int(math.ceil(max_pop / len(args.colors)))
print 'cols', args.colors

tot_time = max([t for t in seq_pops])
if args.range == None:
    tot_plot_time = tot_time
else:
    tot_plot_time = args.range[1] - args.range[0]
    pass

def get_x(*args):
    r = [indices.index(x)+1 if type(x) is int else indices.index(int(x))+1 for x in args]
    if len(r) > 1: return r
    return r[0]

###############################################################################################################################
##  WE ARE ABOUT TO START PLOTTING - FIRST PRINT THE FINAL DATA STATES (MOSTLY FOR DEBUGING)  #################################
###############################################################################################################################

print
print "final pop sequence, max time:", tot_time
for t in sorted(seq_pops):
    print t, seq_pops[t]
    pass
print

print "final mig sequence"
for t in sorted(mig_matrices):
    print t, mig_matrices[t]
    pass
print

print "final pop size sequence"
max_pop_size = 0
max_pop_size2 = 0
min_pop_size = sys.maxint
for t in sorted(pop_sizes):
    print t, pop_sizes[t], pop_growths[t]
    max_pop_size = max(max_pop_size, max([k for k in pop_sizes[t]]))
    tmp_sizes = [k for k in pop_sizes[t] if k != None]
    if len(tmp_sizes) > 0:
        min_pop_size = min(min_pop_size, min(tmp_sizes))
        pass
    max_pop_size2 = max([get_const_size_at_time(t, i) for i in range(1, len(pop_sizes[t])+1)])
    pass
print 'max pop size:', max_pop_size, max_pop_size2
print

#print 'debug size from growth at time 1:', get_const_size_at_time(1, 1)
print

print "plotting indices"
print indices
print

###############################################################################################################################
##  MAKE THE FIGURE OBJECT, SET RANGE (YLIM, XLIM), ETC  ######################################################################
###############################################################################################################################

#fig = plt.figure(figsize=(max(initial_pops),max(seq_pops.keys())))
fig = plt.figure(figsize=(args.width, args.height))
#ax = plt.axes([0,0,max(initial_pops),max(seq_pops.keys())])
#ax = plt.axes([0,0,1,1])
#ax = plt.axes([.10,.10,.8,.8])
ax = plt.axes()
ax.set_xlim(0, max_pop + 1)
if args.range != None:
    ax.set_ylim(args.range[0] - tot_plot_time * .1, args.range[1] + tot_plot_time * .1)
else:
    print 'setting ylim:', -max([t for t in seq_pops]) * .1, max([t for t in seq_pops]) * 1.1
    ax.set_ylim(-max([t for t in seq_pops]) * .1, max([t for t in seq_pops]) * 1.1)
    pass
ax.set_xticks(range(1,len(indices)+1))
ax.set_xticklabels(args.population_names, size = 'xx-large')
ax.set_xlabel('Population IDs', size = 'x-large')

if args.years:
    yticks = ax.get_yticks()
    yticks = yticks * args.gen * 4 * args.N0
    print 'scaling yticks:', yticks
    ax.set_yticklabels(['%d' % int(y) for y in yticks])
    ax.set_ylabel('Years (N0 = %.1f, gen = %.1f)' % (args.N0, args.gen), size = 'x-large')
else:
    print 'not scaling yticks'
    ax.set_ylabel('Scaled Time', size = 'x-large')
    pass

print [0,0,max(initial_pops),max(seq_pops.keys())]


###############################################################################################################################
##  PLOT MIGRATION BOXES FIRST (SO THEY AREN'T PRINTED OVER THE POPULATIONS)  #################################################
###############################################################################################################################

## plot migration boxes
if max_pop > 1:
    max_mig = max([max([max([m[i][j] for j in xrange(len(m[i])) if i != j]) for i in xrange(len(m))]) for m in mig_matrices.values()])
else:
    max_mig = 0
    pass
arrow_plot_fns = []
def arrow_fn(j, i, t, nmt, m):
    i = get_x(i)
    j = get_x(j)
    if i < j:
        dir = 1
        va = 'bottom'
    else:
        dir = -1
        va = 'top'
        pass
    def fn():
        ax.arrow(i + .2 * dir, t + (nmt - t) / 2, j - i - .4 * dir, 0, head_width = tot_plot_time / 40, head_length = .1, color = 'k', alpha = .3, length_includes_head=True, lw=0, shape = 'right', width = tot_plot_time / 100)
        ax.text(i + (j-i)/2, t + (nmt - t) / 2 + tot_plot_time / 100 * dir, '%g' % m, horizontalalignment='center', verticalalignment=va)
        pass
    arrow_plot_fns.append(fn)
    return
print
print 'plotting migrations, max mig amount:', max_mig
for t in sorted(mig_matrices):
    if args.nomig: continue
    nmt = get_next_mig_event(t)
    print 'considering migration from times %s to %s' % (t, nmt)
    #print mig_matrices[t]
    pops = get_pop_list_at_time(t, modify = False)
    for i in pops:
        for j in pops:
            if i == j: continue
            m = mig_matrices[t][i-1][j-1]
            #print 'comparing migs for pops %d and %d between times %s and %s (mig rate is %f):' % (i, j, t, nmt, m)
            ## plot migration if:
            ## - it's between two different pops (i != j)
            ## - it exists (m>0)
            ## - we're at the end of time, so you have to plot whatever migrations you haven't plotted yet
            ## - migration between i and j is going to change at the next time point
            if m > 0 and (t == the_end_time or mig_matrices[t][i-1][j-1] != mig_matrices[nmt][i-1][j-1] or mig_matrices[t][j-1][i-1] != mig_matrices[nmt][j-1][i-1]):
                pt = get_start_of_mig_event(t, i, j, m)
                print 'plotting mig for pops %d to %d at rate %f from time %s to (%s, %s)' % (i, j, m, pt, nmt, t)
                # ax.bar(i, j-i, pt, nmt-pt)
                ax.fill(get_x(i, i, j, j), (pt, nmt, nmt, pt), color = 'black', alpha = .1, zorder = -1)
                arrow_fn(i, j, pt, nmt, m)
                #ax.plot(get_x(i, i, j, j, i), (pt, nmt, nmt, pt, pt), color = 'black')
                pass
            pass
        pass
    pass

###############################################################################################################################
##  PLOT THE BASIC POPULATIONS, WITH GROWTH / SIZE CHANGES  ###################################################################
###############################################################################################################################

## plot population stages
left_x = {p:[] for p in range(1,max_pop+1)}
right_x = {p:[] for p in range(1,max_pop+1)}
left_y = {p:[] for p in range(1,max_pop+1)}
right_y = {p:[] for p in range(1,max_pop+1)}

for t in sorted(seq_pops):
    print 'ploting pops at time', t
    for p in seq_pops[t]:
        size = pop_sizes[t][p-1]
        growth = pop_growths[t][p-1]

        if size != None and growth == None:
            # plot constant size from time t forward (pastward) (until the next size change)
            nt = get_next_size_change(t, p)
            if args.logn:
                w = math.log10(size * args.N0) / math.log10(max_pop_size * args.N0) * .4
            else:
                #print size, max_pop_size
                w = size / max_pop_size * .4
                pass
            #print 'width to adjust size:', p, t, w
            x = get_x(p)
            #ax.fill((x - w, x - w, x + w, x + w), (t, nt, nt, t), color = args.colors[x-1], lw = 1, ec = args.colors[x-1])
            ax.arrow(x+w+.11, t, -.1, 0, width = tot_plot_time / 300, lw = 0, head_width = tot_plot_time / 100, head_length = .02, length_includes_head = True)
            left_x[p] += [x - w, x - w]
            right_x[p] += [x + w, x + w]
            left_y[p] += [t, nt]
            right_y[p] += [t, nt]
            
        elif growth != None:
            # plot growth from time t forward (pastward) (until the next size change)
            alpha, start_time = growth
            nt = get_next_size_change(t, p)
            #print 'plotting growth for pop %d from time %s to %s, growth param %s from start time %f and start size %f' % (p, t, nt, alpha, start_time, size)
            t_points = [t + i/args.curve_points * (nt-t) for i in range(args.curve_points + 1)]
            if args.logn:
                w_points = [math.log10(max(size * math.exp(-alpha * (i-start_time)) * args.N0, 1)) / math.log10(max_pop_size * args.N0) * .4 for i in t_points]
                w_s = math.log10(size * args.N0) / math.log10(max_pop_size * args.N0) * .4
            else:
                w_points = [size * math.exp(-alpha * (i-start_time)) / max_pop_size * .4 for i in t_points]
                w_s = size / max_pop_size * .4
                pass
            #print w_points[:10]
            #print t_points[:10]
            if w_points[0] < w_points[1]: y2 = nt
            else: y2 = t
            x = get_x(p)
            w = min(w_points)
            ## replace these with patches, so that you can draw the entire population in one image (and avoid both lw = 1 and slivers of gaps)
            #ax.fill((x - w, x - w, x + w, x + w), (t, nt, nt, t), color = args.colors[x-1], lw = 1, ec = args.colors[x-1])
            #ax.fill_between([x + w for w in w_points], t_points, y2 = y2, color = args.colors[x-1], lw = 1, edgecolor = args.colors[x-1])
            #ax.fill_between([x - w for w in w_points], t_points, y2 = y2, color = args.colors[x-1], lw = 1, edgecolor = args.colors[x-1])

            ax.arrow(x+w_s+.11, t, -.1, 0, width = tot_plot_time / 300, lw = 0, head_width = tot_plot_time / 100, head_length = .02, length_includes_head = True)

            left_x[p] += [x - w for w in w_points]
            right_x[p] += [x + w for w in w_points]
            left_y[p] += t_points
            right_y[p] += t_points

            pass
        pass
    pass

for p in left_x:
    right_x[p].reverse()
    right_y[p].reverse()
    # print 'left coords for pop', p
    # for i,x in enumerate(left_x[p]):
    #     print x, left_y[p][i]
    #     pass
    # print 'right coords for pop', p
    # for i,x in enumerate(right_x[p]):
    #     print x, right_y[p][i]
    #     pass
    ax.fill(left_x[p] + right_x[p], left_y[p] + right_y[p], lw = 0, ec = 'k', color = args.colors[get_x(p)-1])
    pass

###############################################################################################################################
##  PLOT JOINS, SPLITS  #######################################################################################################
###############################################################################################################################

## plot population joins
for event in sorted(args.ej, key = itemgetter('time', 'card')):
    t = event['time']
    i,j = event['vals']
    #ax.plot(get_x(i,j), (t, t), '--', c = args.colors[get_x(i)-1], linewidth = 3)
    dt = tot_plot_time / 1000
    t_points = numpy.array([x * dt for x in range(20)])
    print "plotting population joins for event", event
    print t
    print t_points
    print t-t_points
    print t+t_points
    wi = numpy.array([get_const_size_at_time(t,i, allow_t = False) / max_pop_size * .4])
    wj = numpy.array([get_const_size_at_time(t,j) / max_pop_size * .4])
    #wi = numpy.array([get_const_size_at_time(k,i, allow_t = False) / max_pop_size * .4 for k in t - t_points])
    #wj = numpy.array([get_const_size_at_time(k,j) / max_pop_size * .4 for k in t + t_points])
    print 'plotting join for pops %d to %d, at time %f' % (i,j,t)
    xi = get_x(i)
    xj = get_x(j)
    print 'i', i, xi, wi, xi-wi[0]
    print 'j', j, xj, wj, xj-wj[0]
    if xi > xj:
        wi = -wi
        wj = -wj
        pass
    #join_dist = ((xj - wj[0]) - (xi + wi[0]))
    #x_points = list(xi + wi) + [join_dist / 3 + xi + wi[0]] + list(xj - wj) + [join_dist / 3 * 2 + xi + wi[0]]
    #y_points = list(t - t_points) + [t] + list(t + t_points) + [t]
    #print x_points
    #print y_points
    # ax.fill(x_points, y_points, color = args.colors[xi-1], lw = 1, ec = args.colors[xi-1])
    ax.plot((xi - wi[0], xj - wj[0]), (t, t), lw = 2, solid_capstyle = 'butt', color = args.colors[xi-1])
    #ax.fill((xi + wi[0], xi + wi[0], xj - wj[0], xj) - wj[0]), (t, t-dt*5, t-dt*5, t), lw = 0, color = args.colors[xi-1])
    #c = 2
    #ax.plot(xi + max(wi) + numpy.sqrt(c - t_points**2), t + t_points)
    pass

## plot population splits
for event in sorted(args.es, key = itemgetter('time', 'card')):
    t = event['time']
    i,j,p = event['vals']
    wi = get_const_size_at_time(t,i) / max_pop_size * .4
    wj = get_const_size_at_time(t,j) / max_pop_size * .4
    #print 'plotting join for pops %d to %d, at time %f' % (i,j,t)
    #print wi,wj
    xi = get_x(i)
    xj = get_x(j)
    if xi < xj:
        wi = -wi
        wj = -wj
        pass
    # ax.fill(x_points, y_points, color = args.colors[xi-1], lw = 1, ec = args.colors[xi-1])
    ax.plot((xi - wi, xj - wj), (t, t), lw = 2, solid_capstyle = 'butt', color = args.colors[xj-1])
    ax.text(xj - (xj-xi)/2, t, '%.2f' % (1-p), horizontalalignment='center', verticalalignment='top')
    #ax.plot(get_x(i,j), (t, t), '--', c = args.colors[get_x(i)-1], linewidth = 3)
    #ax.arrow(get_x(i), t, get_x(j) - get_x(i), 0, ls = 'dashed', color = args.colors[get_x(i)-1], shape='full', head_width = .01, head_length = .1, length_includes_head=True, lw=.5)
    pass

###############################################################################################################################
##  PLOT MIGRATION ARROWS AND RATES (SO THEY GO OVER THE POPULATIONS)  ########################################################
###############################################################################################################################

## plot migration arrows
for fn in arrow_plot_fns:
    fn()
    pass

###############################################################################################################################
##  SAVE THE DAMN FIGURE  #####################################################################################################
###############################################################################################################################

fig.savefig(args.output_file)

