import os, copy, scipy, argparse
from scipy import signal
import numpy as np
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt
from matplotlib import interactive

def get_parser():
    parser = argparse.ArgumentParser(description='Process breathing data (from Smartlab software).')
    parser.add_argument("-smartlab",
                        help="Filename of the Smartlab breathing data (with partial CO2 measurements)",
                        type=str,
                        dest="fname_smartlab",
                        required=True)
    parser.add_argument("-path",
                        help="Path to the data.",
                        type=str,
                        dest="path_data",
                        required=False,
                        default='./')
    parser.add_argument("-v",
                        help="Verbose: 2 = plot data, 1 = no plots",
                        type=int,
                        dest="verbose",
                        choices=[1, 2],
                        default=2)

    return parser 


def get_sec(time_str):
    h, m, s = time_str.split(':')
    return int(h) * 3600 + int(m) * 60 + round(float(s), 4)

def get_data(fname_smartlab):
	nlines_header = 7
	# df = pd.read_csv(os.path.join(path_data, fname_smartlab), header=nlines_header, sep=',')
	#
	f = open(fname_smartlab)
	lines = f.readlines()
	f.close()
	#
	dic_data = {'time': [], 'flow': [], 'concentration': [], 'partial_pressure': []}
	dic_data_small = copy.deepcopy(dic_data)
	#
	for i, l in enumerate(lines[nlines_header:]):
		dat = l.split(',')
		if len(dat) > 4:
			t, flow, concentration, partial_pressure, alert = dat
		else:
			t, flow, concentration, partial_pressure = dat
		#
		t = get_sec(t)
		flow = float(flow)
		concentration = float(concentration)
		partial_pressure = float(partial_pressure)
		#
		dic_data['time'].append(t)
		dic_data['flow'].append(flow)
		dic_data['concentration'].append(concentration)
		dic_data['partial_pressure'].append(partial_pressure)
		# 
		if i%10 == 0:
			dic_data_small['time'].append(t)
			dic_data_small['flow'].append(flow)
			dic_data_small['concentration'].append(concentration)
			dic_data_small['partial_pressure'].append(partial_pressure)
	#
	df = pd.DataFrame(dic_data)
	df_small = pd.DataFrame(dic_data_small)
	#
	return df, df_small


def get_subset(df, start=0, end='max'):
	end = int(max(df.time)+1) if end == 'max' else end
	df_subset = df[df.time>start]
	df_subset = df_subset[df_subset.time<end]
	return df_subset

def plot_data(df, fname='plot_all', todo='save'):
	# interactive(True)
	if todo == 'save':
		figsize = (100, 15)
	else:
		figsize = (30, 5)
	#
	fig, ax = plt.subplots(3, 1, figsize = figsize)
	#
	start = int(min(df.time))
	end= int(max(df.time))
	#
	ax[0].plot(df.time, df.flow, '--', color='b')
	ax[0].set_ylabel('Flow [LPM]')
	ax[1].plot(df.time, df.concentration, '--', color='r')
	ax[1].set_ylabel('CO2 Concentration [%]')
	ax[2].plot(df.time, df.partial_pressure, '--', color='y')
	ax[2].set_ylabel('CO2 Partial Pressure [mmHg]')
	#
	ax[0].set_xlabel('time [sec]')
	#
	for a in ax:
		a.set_xlim(start, end)
	#
	plt.setp(ax, xticks=range(start, end, 10))
	# plt.show()
	#
	if todo == 'save':
		plt.savefig(fname+'.png', bbox_inches='tight',pad_inches=1)
		plt.close()
	else:
		plt.show()
	

def smooth(a, window=10):
	b = np.ones(window)
	conv = np.convolve(b/b.sum(), a)
	remove_end = int(window/2)
	remove_init = remove_end-1 if window%2==0 else remove_end
	return conv[remove_init:-remove_end]

# def fit_square_wave(df, val='partial_pressure'):
# 	scipy.optimize.curve_fit(signal.square, df.time, df[val])
# 	# scipy.optimize.curve_fit(signal.square, df_pre.time, df_pre[val])

def get_end_tidal_pp(df, max_range_peak=100, smooth_window=5, lag=0.6, v=2, prefix=''):
	#
	smoothed_pp = smooth(df.partial_pressure, smooth_window)
	grad = np.gradient(np.asarray(df.partial_pressure))
	grad_from_smoothed = np.gradient(np.asarray(smoothed_pp))
	# grad_smoothed = smooth(grad)
	#
	peaks_ind = signal.find_peaks_cwt(-grad_from_smoothed, np.arange(1,max_range_peak))
	#
	# peaks_ind_pp = signal.find_peaks_cwt(df.partial_pressure, np.arange(1,max_range_peak))
	#
	time_corr = [round(df.time.iloc[ind]-lag, 2) for ind in peaks_ind]
	end_tidal = [round(list(df.partial_pressure[df.time == t])[0] if len(list(df.partial_pressure[df.time == t])) == 1 else np.nan, 3) for t in time_corr]
	# plot
	if v==2:
		fig, ax = plt.subplots(2,1, figsize = (80, 10))
		ax[0].plot(df.time, df.partial_pressure, '--', color='r')
		ax[0].plot(df.time, smoothed_pp, '--', color='k')
		ax[0].set_ylabel('CO2 Partial Pressure [mmHg]')
		ax[0].set_xlabel('Time [sec]')
		#
		ax[1].plot(df.time, grad, '--', color='b')
		ax[1].plot(df.time, grad_from_smoothed, '--', color='k')
		ax[1].set_ylabel('Gradient of the partial pressure')
		ax[1].set_xlabel('Time [sec]')
		# ax[1].plot(df_post.time, grad_smoothed, '--', color='b')
		#
		for ind in peaks_ind:
			ax[0].axvline(x=df.time.iloc[ind], color='r', linewidth=0.5)
			ax[1].axvline(x=df.time.iloc[ind], color='r', linewidth=0.5)
		#
		for t in time_corr:
			ax[0].axvline(x=t, linestyle='--', color='r', linewidth=0.5)
			ax[1].axvline(x=t, linestyle='--', color='r', linewidth=0.5)			
		#
		# for ind in peaks_ind_pp:
		# 	ax[0].axvline(x=df.time.iloc[ind], color='k', linewidth=0.2)
		# 	ax[1].axvline(x=df.time.iloc[ind], color='k', linewidth=0.2)
		#
		plt.setp(ax, xticks=range(int(min(df.time)), int(max(df.time)), 5))
		fname = prefix+'end_tidal_values_smooth'+str(smooth_window)+'_range'+str(max_range_peak)+'.png'
		plt.savefig(fname, bbox_inches='tight',pad_inches=1)
		plt.close()
		print 'Breathing data and found end-tidal times saved in: '+fname
		print '(End tidal times are displayed with the red dotted lines)'
	#
	return end_tidal

def get_end_tidal_pre_post(fname_smartlab, ofolder, start_pre=None, end_pre=None, start_post=None, end_post=None, v=2):
	# get breathing data
	# use df_small if data were registered every 10ms
	# use df if data were registered every 100ms
	df, df_small = get_data(fname_smartlab)
	df_to_use = df #_small
	#
	if None in [start_pre, end_pre, start_post, end_post]:
		# show data
		print 'Note the start and end times (in sec) on the scans BEFORE and DURING the CO2 challenge: '
		plot_data(df_to_use, todo='show')
		#
		# get start and end times from user input
		start_pre = int(raw_input('Start time of scan PRE-challenge ? (in sec)\n'))
		end_pre = int(raw_input('End time of scan PRE-challenge ? (in sec)\n'))
		start_post = int(raw_input('\nStart time of scan  DURING challenge ? (in sec)\n'))
		end_post = int(raw_input('End time of scan DURING challenge ? (in sec)\n'))
	# JARED:
	# start_pre = 370
	# end_pre = 660
	# start_post = 790
	# end_post = 900
	# BEN:
	# start_pre = 1450
	# end_pre = 1750
	# start_post = 1880
	# end_post = 2180
	#
	# select pre-challenge data and during challenge (post) data
	df_pre = get_subset(df_to_use, start=start_pre, end=end_pre)
	df_post = get_subset(df_to_use, start=start_post, end=end_post)
	#
	if v == 2:
		plot_data(df_to_use, fname=os.path.join(ofolder,'plot_al'))
		plot_data(df_pre, fname=os.path.join(ofolder,'plot_pre'))
		plot_data(df_post, fname=os.path.join(ofolder,'plot_post'))

	end_tidal_pre = get_end_tidal_pp(df_pre, prefix=os.path.join(ofolder, 'pre_'), v=v)
	end_tidal_post = get_end_tidal_pp(df_post, prefix=os.path.join(ofolder, 'post_'), v=v)
	#
	if v != 0:
		print 'mean ETCO2 - basal', np.mean(end_tidal_pre)
		print 'mean ETCO2 - hyper', np.mean(end_tidal_post) 
	#
	return end_tidal_pre, end_tidal_post, start_pre, end_pre, start_post, end_post

if __name__ == "__main__":
	parser = get_parser()
	param = parser.parse_args()

	# path_data = '/Users/sdupont/Documents/HIV/vasoreactivity/2018-02-12-Ben_CO2_challenge'
	# #
	# # fname_smartlab = 'SmartLab_1-22-2018_test_ASL_jared_hdr.txt'
	# # fname_smartlab = 'SmartLab_1-24-2018_test_save_smartlab_after_replay.txt'
	# fname_smartlab = 'SmartLab_2-9-2018_testBen.txt'
	# #
	# param = parser.parse_args(['-smartlab', fname_smartlab, '-path', path_data])

	os.chdir(param.path_data)

	get_end_tidal_pre_post(param.fname_smartlab, './', v=param.verbose)

