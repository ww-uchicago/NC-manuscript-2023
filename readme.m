%% analysis summary of Vm fast events in SAC in K internal, Cs internal and 1 mM TEA
first run Vm event detection for every abf file, double check event detection, 
then run summary Vm events to align start time
then run get_all_fast_event_all_cells to get prominence during baseline (bin -3 and before, and bin 20 and after) for each cell
then run cdf_plots to plot cummulative distribution of event amplitude across drug conditions.

then do the above for cs + LY and LY + TEA

to get average bin events during baseline and bar, run get all fast event all cells poolCounts
