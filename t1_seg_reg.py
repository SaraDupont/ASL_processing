import os, commands
import pandas as pd
import numpy as np
import xml.etree.ElementTree
from utils import *


def segment_brain(fname_t1, centered_fov=False, ofolder=None):
	# extract brain with bet
	fname_t1_brain = add_suffix(fname_t1, '_brain')
	#
	if ofolder is not None and os.path.isdir(ofolder):
		fname_t1_brain = os.path.join(ofolder, fname_t1_brain.split('/')[-1])
	#
	cmd_bet = 'bet '+fname_t1+' '+fname_t1_brain+' -m '
	if not centered_fov : 
		im_t1 = Image(fname_t1)
		i_dim_z = im_t1.orientation.index('I') if 'I' in im_t1.orientation else im_t1.orientation.index('S')
		factor_get_center = [0.5, 0.5, 0.5]
		factor_get_center[i_dim_z] = 2.0/3
		center = ' '.join(list((im_t1.dim * np.asarray(factor_get_center)).astype(int).astype(str)))
		cmd_bet += ' -c '+center+' '
	#
	print 'Segmenting brain in T1: '
	print cmd_bet
	s, o = commands.getstatusoutput(cmd_bet)
	#
	fname_t1_brain_mask = add_suffix(fname_t1_brain, '_mask')
	if s == 0:
		return fname_t1_brain, fname_t1_brain_mask
	else:
		print 'ERROR while extracting brain with command: '+cmd_bet
		return None, None


def segment_gm_wm(fname_t1_brain, ofolder=None):
	'''
	:param fname_t1_brain: T1 image * brain mask = T1 inside the brain, zeros outside
	'''
	fname_gm = add_suffix(fname_t1_brain, '_pve_1')
	fname_wm = add_suffix(fname_t1_brain, '_pve_2')
	cmd_fast = 'fast '+fname_t1_brain
	if not os.path.isfile(fname_gm) or not os.path.isfile(fname_wm):
		print 'Segmenting GM/WM in T1: '
		print cmd_fast
		s, o = commands.getstatusoutput(cmd_fast)
	else:
		s=0
	#
	if s == 0: 
		if ofolder is not None and ofolder != '/'.join(fname_t1_brain.split('/')[:-1]): 
			new_fname_gm = os.path.join(ofolder, fname_gm.split('/')[-1])
			new_fname_wm = os.path.join(ofolder, fname_wm.split('/')[-1])
			#
			shutil.move(fname_gm, new_fname_gm)
			shutil.move(fname_wm, new_fname_wm)
			# 
			fname_gm = new_fname_gm
			fname_wm = new_fname_wm
		#
		return fname_gm, fname_wm
	else:
		print "ERROR while segmenting brain with command: "+cmd_fast
		return None, None


def register_image(fname_in, fname_dest, init_mat=None, init_nl_transfo=None, fname_mask_in=None, fname_mask_dest=None, linear_reg=True, non_linear_reg=True):
	#
	name_in = fname_in.split('/')[-1].split('.')[0]
	# fname_transfo = name_in+'_to_'+fname_dest_reorient.split('.')[0]+'_affine'
	fname_im_out = name_in+'_to_'+fname_dest.split('/')[-1].split('.')[0]+'_affine'
	fname_transfo_lin = fname_im_out+'_transform.mat'
	if not os.path.isfile(fname_transfo_lin) and linear_reg:
		# MNI152_T1_0.5mm.nii.gz is also available
		print '--- registering '+name_in+' to image '+fname_dest
		# list_param_cmd = ['flirt', '-in', fname_in, '-ref', fname_dest_reorient, '-omat', fname_transfo_lin, '-out', fname_transfo]
		list_param_cmd = ['flirt', '-in', fname_in, '-ref', fname_dest, '-omat', fname_transfo_lin, '-out', fname_im_out]
		if init_mat is not None:
			list_param_cmd.append('-init')
			list_param_cmd.append(init_mat)
		# 
		cmd = ' '.join(list_param_cmd)
		print cmd
		status_lin, output_lin = commands.getstatusoutput(cmd) 
		#
		# fname_return = fname_transfo_lin if status == 0 else None
		# print '-> done !'
		#
	else: 
		status_lin = 0
		if not linear_reg:
			print 'Not doing linear registration'
			fname_transfo_lin = None
	#
	if status_lin == 0:
		if non_linear_reg:
			fname_out_nl = name_in+'_to_'+fname_dest.split('.')[0]+'_non_linear'
			fname_transfo_nl = fname_out_nl+'_transform'
			if not os.path.isfile(fname_transfo_nl+'.nii.gz'):
				cmd = 'fnirt --in='+fname_in+' --ref='+fname_dest+' --iout='+fname_out_nl+' --cout='+fname_transfo_nl
				cmd += ' --aff='+fname_transfo_lin if linear_reg else ''
				cmd += ' --aff='+init_mat if not linear_reg and init_mat is not None else ''
				cmd += ' --inwarp='+init_nl_transfo if not linear_reg and init_nl_transfo is not None else ''
				cmd += ' --inmask='+fname_mask_in if fname_mask_in is not None else ''
				cmd += ' --refmask='+fname_mask_dest if fname_mask_dest is not None else ''
				cmd += ' --infwhm=0,0,0,0 --reffwhm=0,0,0,0'
				# cmd += ' --warpres=1.5,1.5,3.3'
				print cmd
				status_nl, output_nl = commands.getstatusoutput(cmd)
				fname_out_nl += '.nii.gz'
				#
				# fname_return = fname_transfo_nl+'.nii.gz' if status_nl==0 else fname_return
			else: 
				status_nl = 0
			fname_transfo_nl += '.nii.gz'
		else: 
			print 'Not doing non linear registration'
			status_nl = -1
		#
		if status_nl != 0:
			print '\n\n Status non linear reg IS NOT 0'
			message = 'Error while doing non linear spline registration: \n'+output_nl if status_nl != -1 else 'Not doing non linear registration'
			print message
			fname_transfo_nl = None
	else: 
		print '\n\n Status linear reg IS NOT 0'
		print 'Error while doing linear affine registration: \n'+output_lin
		fname_transfo_lin, fname_transfo_nl = None, None
	#
	fname_im_out += '.nii.gz'
	# fname_out_nl += '.nii.gz'
	#
	s_fsl = 'fslview '+fname_dest+' '
	s_fsl += fname_im_out+' ' if os.path.isfile(fname_im_out) else ''
	if non_linear_reg:
		s_fsl += fname_out_nl+' ' if os.path.isfile(fname_out_nl) else ''
	s_fsl += '&\n'
	print 'To checkout results, type: \n', s_fsl
	#
	return fname_transfo_lin, fname_transfo_nl	
	#		


def apply_transfo(fname_in, fname_dest, fname_transfo, interpolation='trilinear',  name_transformed='', ofolder='./'): #'nearestneighbour'
	# flirt -in "$path_atlases""$fname_atlas""$ext" -ref "$fname_metric_img""$suffix_reorient""$ext" -out "$fname_atlas_reg"  -init "$fname_mat_metric_img" -applyxfm -interp "$interp"
	name_transformed = os.path.basename(fname_in).split('.')[0] if name_transformed == '' else name_transformed
	name_dest = fname_dest.split('/')[-1]
	fname_out = os.path.join(ofolder, name_transformed+'_reg_to_'+name_dest)
	# list_param_cmd = ['flirt', '-in', fname_in, '-ref', fname_dest, '-init', fname_transfo, '-out', fname_out, '-applyxfm', '-interp', interpolation]
	#
	arg_transfo = '--premat' if '.mat' in fname_transfo else '--warp'
	cmd='applywarp --in='+fname_in+' --ref='+fname_dest+' --out='+fname_out+' '+arg_transfo+'='+fname_transfo+' --interp='+interpolation
	#
	status, output = commands.getstatusoutput(cmd) 
	#
	s_fsl = 'fslview '+fname_dest+' '
	s_fsl += fname_out+' ' if os.path.isfile(fname_out) else ''
	s_fsl += '&\n'
	print 'To checkout results, type: \n', s_fsl
	#
	return fname_out


def get_val_from_ROI(im_metric, im_roi, val_roi=1, list_measures=["mean", "median", "std", "min", "max"], perc_keep=95):
	# get metric data inside ROI
	if val_roi == 'pve':
		# taking partial volume effect into account (single ROI per file i.e GM or WM)
		data_roi = im_metric.data*im_roi.data
		# data_roi = data_roi[np.abs(data_roi)>0.001]
	else:
		# taking data equal to the specified label for that roi (=data from atlas)
		data_roi = im_metric.data[im_roi.data==val_roi]
	data_roi = data_roi[np.abs(data_roi)>0.02]
	# 
	if perc_keep<100:
		n_count = len(data_roi) 
		data_roi = np.sort(data_roi)
		#
		perc_end_tail = (100.0-perc_keep)/2.0
		n_end_tail = int((n_count * perc_end_tail)/100.0)
		#
		data_roi = data_roi[n_end_tail:-n_end_tail]
	## get relevant values from that data
	# med = np.median(data_roi)
	# mean = np.mean(data_roi)
	# std = np.std(data_roi)
	# mi = np.min(data_roi)
	# ma = np.max(data_roi)
	# #
	# dict_val = {'median': med, 'mean': mean, 'std': std, 'min': mi, 'max': ma}
	dict_val = {}
	#
	for measure in list_measures:
		try: 
			func = getattr(np, measure)
			dict_val[measure] = func(data_roi)
		except: 
			'ERROR: function %s is not defined in numpy' %measure
		#
	return dict_val


def get_atlas_info(fname_xml):
	f = xml.etree.ElementTree.parse(fname_xml)
	root = f.getroot()
	data = root.findall('data')[0]
	list_labels = data.findall('label')
	dic_atlas_labels = {}
	for label in list_labels:
		region_name = label.text
		region_index = int(label.attrib['index'])
		## add 1 because the label indices start at 0 and the label values in the nifti file start at 1
		dic_atlas_labels[region_index+1]=region_name
	#
	return dic_atlas_labels


class Registration:
	def __init__(self, fname_t1, fname_t1_brain_mask, fname_cvr, fname_metric_mask, fname_cbf_pre=None, fname_cbf_post=None, dict_param=None, fname_gm=None, fname_wm=None, ofolder=None):
		#
		self.fname_t1 = fname_t1
		self.fname_t1_brain_mask = fname_t1_brain_mask
		self.fname_cvr = fname_cvr
		self.fname_cbf_pre = fname_cbf_pre
		self.fname_cbf_post = fname_cbf_post
		self.fname_metric_mask = fname_metric_mask
		#
		s, fsl_dir = commands.getstatusoutput("echo $FSL_DIR")
		assert fsl_dir != '', 'ERROR $FSL_DIR not defined...' 
		self.fsl_dir = fsl_dir
		self.path_atlas_im = os.path.join(fsl_dir, 'data/standard/MNI152_T1_1mm.nii.gz')
		self.path_atlas_im_brain = os.path.join(fsl_dir, 'data/standard/MNI152_T1_1mm_brain.nii.gz')
		self.path_atlas_im_brain_mask = os.path.join(fsl_dir, 'data/standard/MNI152_T1_1mm_brain_mask.nii.gz')
		#
		self.path_atlas_cort = os.path.join(fsl_dir, 'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz')
		self.path_atlas_subcort = os.path.join(fsl_dir, 'data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr0-1mm.nii.gz')
		#
		self.path_atlas_info_cort = os.path.join(fsl_dir, 'data/atlases/HarvardOxford-Cortical.xml')
		self.path_atlas_info_subcort = os.path.join(fsl_dir, 'data/atlases/HarvardOxford-Subcortical.xml')
		#
		self.fname_gm = fname_gm
		self.fname_wm = fname_wm
		#
		self.dict_param = {"atlas": ["cortical", "subcortical", "WMGM"],
		"metric": ["CBF_pre", "CBF_post", "CVR"],
		"measures": ["mean", "median", "std", "min", "max"]} if dict_param is None else dict_param
		#
		self.ofolder = ofolder if ofolder is not None else os.path.abspath('.')
	#
	def pipeline(self):
		# go to output directory
		os.chdir(self.ofolder)
		# register MNI to t1 image
		fname_mat_reg_t1_lin, fname_mat_reg_t1_nl = register_image(fname_in=self.path_atlas_im, fname_dest=self.fname_t1, fname_mask_in=self.path_atlas_im_brain_mask, fname_mask_dest=self.fname_t1_brain_mask, linear_reg=self.dict_param['registration']['linear'], non_linear_reg=self.dict_param['registration']['non-linear'])
		#
		#
		# register template to 1 metric
		fname_metric = self.fname_cbf_pre if self.fname_cbf_pre is not None else self.fname_cvr
		## segment brain on metric if mask is None: 
		if self.fname_metric_mask is None: 
			fname_metric_brain, fname_metric_brain_mask = segment_brain(fname_metric, centered_fov=True)
			self.fname_metric_mask = fname_metric_brain_mask
			fname_metric = fname_metric_brain
		# do registration
		fname_mat_reg_metric_lin, fname_mat_reg_metric_nl = register_image(fname_in=self.path_atlas_im_brain, fname_dest=fname_metric, init_mat=fname_mat_reg_t1_lin, fname_mask_in=self.path_atlas_im_brain_mask, fname_mask_dest=self.fname_metric_mask, linear_reg=self.dict_param['registration']['linear'], non_linear_reg=self.dict_param['registration']['non-linear'])
		fname_mat_reg_metric = fname_mat_reg_metric_nl if fname_mat_reg_metric_nl is not None else fname_mat_reg_metric_lin
		#
		# register t1 to metric (for GM WM extraction)
		if 'WMGM' in self.dict_param['atlas']:
			fname_t1_brain = mask_image(self.fname_t1, self.fname_t1_brain_mask, suffix='_brain', ofolder=self.ofolder)
			#
			fname_mat_t1_to_metric_lin, fname_mat_t1_to_metric_nl = register_image(fname_t1_brain, fname_metric, fname_mask_in=self.fname_t1_brain_mask, fname_mask_dest=self.fname_metric_mask, linear_reg=self.dict_param['registration']['linear'], non_linear_reg=self.dict_param['registration']['non-linear'])
			fname_mat_t1_to_metric = fname_mat_t1_to_metric_nl if fname_mat_t1_to_metric_nl is not None else fname_mat_t1_to_metric_lin
		#
		for metric, fname_metric in zip(["CBF_pre", "CBF_post", "CVR"], [self.fname_cbf_pre, self.fname_cbf_post, self.fname_cvr]):
			if metric in self.dict_param['metric'] and fname_metric is not None:
				# extract values from atlas ROIs: 
				fname_excel_out = 'atlas_values_'+metric+'.xls'
				# excel writer to output the pandas
				writer = pd.ExcelWriter(fname_excel_out) 
				#
				if 'cortical' in self.dict_param['atlas']:
					name_atlas = 'cortical'
					data_atlas = self.extract_values_atlas(fname_metric, fname_mat_reg_metric, name_atlas)
					data_atlas.to_excel(writer, name_atlas)
				#
				if 'subcortical' in self.dict_param['atlas']:
					name_atlas = 'subcortical'
					data_atlas = self.extract_values_atlas(fname_metric, fname_mat_reg_metric, name_atlas)
					data_atlas.to_excel(writer, name_atlas)
				#
				if 'WMGM' in self.dict_param['atlas']:
					# register T1 to metric and use transfo to bring GM/WM in metric space + extarct metric values
					data_res_gm_wm = self.extract_values_gm_wm(fname_metric, fname_mat_t1_to_metric)
					data_res_gm_wm.to_excel(writer, 'GM_WM')
				#
				writer.save()
		#
	#
	#
	def extract_values_atlas(self, fname_metric, fname_mat_reg, name_atlas):
		# load metric and ROI images
		im_metric = Image(fname_metric)
		#
		dict_atlases = {'cortical': (self.path_atlas_cort, self.path_atlas_info_cort),
						'subcortical': (self.path_atlas_subcort, self.path_atlas_info_subcort)
						}
		#
		print 'Extracting Atlas values from '+fname_metric+' using registration '+fname_mat_reg
		#
		path_atlas, path_atlas_info = dict_atlases[name_atlas]
		#
		fname_atlas_reg = apply_transfo(path_atlas, fname_metric, fname_mat_reg, interpolation='nn', name_transformed='MNI'+name_atlas, ofolder=self.ofolder)
		# if fname_atlas_reg.split('.')[-1] == 'nii':
		#     fname_atlas_reg += '.gz'
		im_atlas = Image(fname_atlas_reg)
		#
		# mask atlas with brain mask to keep only values inside the brain
		im_mask = Image(self.fname_metric_mask)
		im_atlas.data = im_atlas.data * im_mask.data
		#
		dict_val_by_roi = {}
		#
		list_labels = list(np.unique(im_atlas.data))
		list_labels.pop(list_labels.index(0))
		for l in list_labels:
			dict_val = get_val_from_ROI(im_metric, im_atlas, val_roi=l, list_measures=self.dict_param['measures'])
			dict_val['label'] = l
			dict_val_by_roi[l] = dict_val
		#
		data_atlas = pd.DataFrame(dict_val_by_roi).transpose()
		#
		dic_atlas_info = get_atlas_info(path_atlas_info)
		col_names = [dic_atlas_info[l] for l in data_atlas.label]
		data_atlas['region_name'] = col_names
		# 
		#
		return data_atlas
	#
	#
	def extract_values_gm_wm(self, fname_metric, fname_mat_t1_to_metric):
		dict_res = {}
		#
		print 'Extracting GM and WM values from '+fname_metric+' using registration '+fname_mat_t1_to_metric
		#
		fname_gm_reg = apply_transfo(self.fname_gm, fname_metric, fname_mat_t1_to_metric,  name_transformed='gm', ofolder=self.ofolder)
		fname_wm_reg = apply_transfo(self.fname_wm, fname_metric, fname_mat_t1_to_metric,  name_transformed='wm', ofolder=self.ofolder)
		# extract values from ROI on registered GM and WM
		im_metric = Image(fname_metric)
		im_gm = Image(fname_gm_reg)
		im_wm = Image(fname_wm_reg)
		#
		dict_res['GM'] = get_val_from_ROI(im_metric, im_gm, val_roi='pve', list_measures=self.dict_param['measures'])
		dict_res['WM'] = get_val_from_ROI(im_metric, im_wm, val_roi='pve', list_measures=self.dict_param['measures'])
		#
		return pd.DataFrame(dict_res).transpose()
	#

