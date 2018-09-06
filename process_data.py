import argparse, json, os, commands
from json_minify import json_minify

from utils import * 
from t1_seg_reg import *
from assess_labeling_efficiency import get_corrected_alpha
from process_asl import *
from process_breathing import * 
from get_cvr import *


def get_parser():
    parser = argparse.ArgumentParser(description='Process ASL data to (1) get CBF, (2) get CVR, (3) extract values within antomical ROI. This script relies on FSL tools and home made methods.')
    parser.add_argument("-c", "--config", 
                        help="Filename of the JSON config file.",
                        type=input_file,
                        dest="fname_json",
                        required=True)
    parser.add_argument("-v",
                        help="verbose",
                        type=int,
                        dest="verbose", 
                        default=1)
    return parser


class Config:
	def __init__(self, command_line_args):
		self.command_line_args = command_line_args
		#
		self.config_params = self.read_json_config()
		#
		self.processing = self.config_params['process']
		self.roi = self.config_params['roi']
		self.path_data = self.config_params['data']
		# 
		self.processing = self.format_bool_dict(self.processing)
		self.format_path_dict()
		self.check_input_match_processing()
		#
		self.config_params['t1'] = self.format_bool_dict(self.config_params['t1'])
		self.roi['registration'] = self.format_bool_dict(self.roi['registration'])
	#
	def read_json_config(self):
		# read JSON
		with open(self.command_line_args.fname_json) as config_file:
			json_str = config_file.read()
			config_params = json.loads(json_minify(json_str))
		#
		return config_params
	#
	def format_bool_dict(self, dictionary):
		for k, v in dictionary.items():
			dictionary[k] = bool(v)
		return dictionary
	#
	def format_path_dict(self):
		path_wd = self.path_data['path']
		path_ofolder = self.path_data['output_folder']
		#
		if not os.path.isdir(path_wd) and path_wd != '':
			raise IOError('input folder does not exist: %s' %path_wd)
		#
		path_ofolder = create_folder(path_ofolder)
		self.path_data['output_folder'] = path_ofolder
		#
		for k, fname in self.path_data.items():
			if k != 'path' and k != 'output_folder' and fname != '':
				if not os.path.isabs(fname):
					fname = os.path.join(path_wd, fname)
				if not os.path.isfile(fname):
					if os.path.isdir(fname) and '.dcm' in os.listdir(fname)[1]:
						fname = self.convert_dicom_to_nifti(fname, k)
					else:
						raise IOError('input file does not exist: %s' %fname)
				#
				self.path_data[k] = fname
			if fname == '':
				self.path_data[k] = None
	#
	def convert_dicom_to_nifti(self, path_dcm, file_out):
		default_ext = '.nii.gz'
		fname_out = os.path.join(self.path_data['output_folder'], file_out+default_ext)
		if not os.path.isfile(fname_out):
			#
			path_dcm = replace_char(path_dcm)
			#
			cmd = "dcm2niix -o "+self.path_data['output_folder']+" -z y -b n -f "+file_out+" "+path_dcm+"/ "
			#
			print 'Converting dicom to nifti: ', cmd
			status, output = commands.getstatusoutput(cmd) 
			#
			if status != 0:
				print ' --> OUTPUT DCM2NIIX: ', output
				raise IOError('Error during dicom to nifti conversion of folder %s ...' %path_dcm)
			else: 
				fname_other_ext = os.path.join(self.path_data['output_folder'], file_out+'.nii')
				if not os.path.isfile(fname_out) and os.path.isfile(fname_other_ext):
					default_ext = '.nii'
					s, o = commands.getstatusoutput('gzip -6 '+fname_other_ext)
					# shutil.move(os.path.join(self.path_data['output_folder'], file_out+'.nii'), fname_out)
				#
				if os.path.isfile(os.path.join(self.path_data['output_folder'], file_out+'a'+default_ext)):
					im = Image(fname_out)
					for c in 'abcdef':
						fname_sup = os.path.join(self.path_data['output_folder'], file_out+c+default_ext)
						if os.path.isfile(fname_sup):
							im_sup = Image(fname_sup)
							new_shape = list(im_sup.data.shape)
							new_shape.append(1)
							dat = im.data.reshape(new_shape) if len(im.data.shape)==3 else im.data
							dat_sup = im_sup.data.reshape(new_shape) if len(im_sup.data.shape)==3 else im_sup.data
							#
							im.data = np.concatenate((dat, dat_sup), axis=3)
							im.save()
							#
							os.remove(im_sup.fname)
		#
		return fname_out

	#
	def check_input_match_processing(self):
		self.requirement('Brain_seg_t1', 't1')
		self.requirement('WMGM_seg_t1', 't1')
		self.requirement('CBF_pre', 'asl_pre')
		self.requirement('CBF_post', 'asl_post')
		self.requirement('Breathing', 'breathing')
		self.requirement('CVR', 'asl_pre')
		self.requirement('CVR', 'asl_post')
		self.requirement('Get_values_roi', 't1')
		#
	#	
	def requirement(self, key_processing, key_file):
		if self.processing[key_processing]:
			assert self.path_data[key_file] is not None, "ERROR: Processing %s requires file %s. \nPlease fill in config file accordingly: \n--> %s " %(key_processing, key_file, self.command_line_args.fname_json)
	#
	def save(self, suffix='_processed'):
		file_in = self.command_line_args.fname_json.split('/')[-1]
		fname_out = os.path.join(self.path_data['output_folder'], add_suffix(file_in, suffix))
		#
		self.config_params['process'] = self.processing
		self.config_params['data'] = self.path_data
		self.config_params['roi'] = self.roi
		#
		j = json.dumps(self.config_params, indent=4)
		f = open(fname_out, 'w')
		print >> f, j
		f.close()
		# with open(fname_out, 'w') as outfile:
		# 	json.dump(self.config_params, outfile)


class Pipeline:
	def __init__(self, command_line_args):
		self.config = Config(command_line_args)
	#
	def process(self):
		# 
		## PROCESS ASL PRE TO GET CBF
		if self.config.processing['CBF_pre'] and self.config.path_data['cbf_pre'] is None:
			#
			if int(self.config.config_params['labeling_efficiency']['tCBF_pre']) != 0:
				## correct labeling efficiency 
				alpha_pre = get_corrected_alpha(self.config.config_params['labeling_efficiency']['tCBF_pre'], self.config.config_params['asl'], self.config.path_data['asl_pre'], fname_asl_mask=self.config.path_data['asl_pre_brain_mask'], fname_t1_mask=self.config.path_data['t1_brain_mask'], fname_t1=self.config.path_data['t1'], ofolder=self.config.path_data['output_folder'], fname_alpha='alpha_pre.txt')
				#
				self.config.config_params['asl']['alpha'] = alpha_pre
				#
			#
			asl_pre = AslData(self.config.config_params['asl'], self.config.path_data['asl_pre'], fname_mask=self.config.path_data['asl_pre_brain_mask'], ofolder=self.config.path_data['output_folder'])
			asl_pre.process()
			self.config.path_data['cbf_pre'] = asl_pre.im_cbf.fname
			self.config.path_data['asl_pre_brain_mask'] = asl_pre.im_mask.fname
		#
		## PROCESS ASL PRE TO GET CBF
		if self.config.processing['CBF_post'] and self.config.path_data['cbf_post'] is None:
			# 
			if int(self.config.config_params['labeling_efficiency']['tCBF_post']) != 0:
				## correct labeling efficiency 
				alpha_post = get_corrected_alpha(self.config.config_params['labeling_efficiency']['tCBF_post'], self.config.config_params['asl'], self.config.path_data['asl_post'], fname_asl_mask=self.config.path_data['asl_post_brain_mask'], fname_t1_mask=self.config.path_data['t1_brain_mask'], fname_t1=self.config.path_data['t1'], ofolder=self.config.path_data['output_folder'], fname_alpha='alpha_post.txt')
				self.config.config_params['asl']['alpha'] = alpha_post
				#
			#
			asl_post = AslData(self.config.config_params['asl'], self.config.path_data['asl_post'], fname_mask=self.config.path_data['asl_post_brain_mask'], ofolder=self.config.path_data['output_folder'])
			asl_post.process()				
			#
			self.config.path_data['cbf_post'] = asl_post.im_cbf.fname
			self.config.path_data['asl_post_brain_mask'] = asl_post.im_mask.fname
			#
			if self.config.path_data['cbf_pre'] is not None and self.config.path_data['asl_pre_brain_mask'] is not None:
				# register CBF post to CBF pre 
				fname_cbf_post_out = add_suffix( asl_post.im_cbf.fname, '_flirt_reg')
				cmd_moco = 'flirt -in '+asl_post.im_cbf.fname+' -ref '+self.config.path_data['cbf_pre']+' -out '+fname_cbf_post_out
				s, o = commands.getstatusoutput(cmd_moco)
				#
				self.config.path_data['cbf_post'] = fname_cbf_post_out
				self.config.path_data['asl_post_brain_mask'] = self.config.path_data['asl_pre_brain_mask'] 
				#
		#
		## PROCESS SMARTLAB BREATHING DATA
		if self.config.processing['Breathing']:
			# process breathing
			for k, v in self.config.config_params['challenge'].items():
				if v == 0:
					self.config.config_params['challenge'][k] = None
			#
			# add verbose with v=x
			end_tidal_pre, end_tidal_post = get_end_tidal_pre_post(self.config.path_data['breathing'], self.config.path_data['output_folder'], start_pre=self.config.config_params['challenge']['start_pre'], end_pre=self.config.config_params['challenge']['end_pre'], start_post=self.config.config_params['challenge']['start_post'], end_post=self.config.config_params['challenge']['end_post'])
			# save figures only if verbose 2 ?
			self.config.config_params['challenge']['end_tidal_pre'] = np.nanmean(end_tidal_pre)
			self.config.config_params['challenge']['end_tidal_post'] = np.nanmean(end_tidal_post)
		#
		## GET CVR
		if self.config.processing['CVR'] and self.config.path_data['cvr'] is None:
			# process CVR using ETCO2 values, CBF pre and CBF post
			if self.config.config_params['challenge']['type'] == 'CO2':
				assert self.config.config_params['challenge']['end_tidal_pre'] is not None, "ERROR: End Tidal CO2 value for the pre-challenge scan is empty"
				assert self.config.config_params['challenge']['end_tidal_post'] is not None, "ERROR: End Tidal CO2 value for the post-challenge scan is empty"
			elif self.config.config_params['challenge']['type'] == 'Diamox':
				self.config.config_params['challenge']['end_tidal_pre'] = None 
				self.config.config_params['challenge']['end_tidal_post'] = None
			else: 
				raise IOError('Type of CVR challenge not recognized: %s. Accepted types are: ["CO2", "Diamox"]' %self.config.config_params['challenge']['type'])
			#
			assert self.config.path_data['cbf_pre'] is not None, "ERROR: No CBF map pre challenge"
			assert self.config.path_data['cbf_post'] is not None, "ERROR: No CBF map post challenge"
			#
			# get CVR
			im_cbf_pre = Image(self.config.path_data['cbf_pre'])
			# if self.config.path_data['cbf_post'].split('.')[-1] == 'nii':
	  #   		    self.config.path_data['cbf_post'] += '.gz'
			im_cbf_post = Image(self.config.path_data['cbf_post'])
			#
			im_mask = Image(self.config.path_data['asl_pre_brain_mask']) if self.config.path_data['asl_pre_brain_mask'] is not None else None
			#
			im_cvr = get_cvr(im_cbf_pre, im_cbf_post, self.config.config_params['challenge']['end_tidal_pre'], self.config.config_params['challenge']['end_tidal_post'], im_mask=im_mask)
			self.config.path_data['cvr'] = im_cvr.fname
		#
		#
		## SEGMENT BRAIN ON T1
		if self.config.processing['Brain_seg_t1']:
			fname_t1_brain, fname_t1_brain_mask = segment_brain(self.config.path_data['t1'], bool(self.config.config_params['t1']['centered_fov']), ofolder=self.config.path_data['output_folder'])
			#
			self.config.path_data['t1_brain_mask'] = fname_t1_brain_mask
		else:
			fname_t1_brain, fname_t1_brain_mask = None, None
		##
		## SEGMENT PARTIAL GM AND WN VOLUMES
		if self.config.processing['WMGM_seg_t1'] or (self.config.processing['Get_values_roi'] and 'WMGM' in self.config.roi['atlas'] and (self.config.path_data['t1_wm_mask'] is None or self.config.path_data['t1_gm_mask'] is None)):
			if fname_t1_brain is None: 
				if self.config.path_data['t1_brain_mask'] is not None:
					fname_t1_brain = mask_image(self.config.path_data['t1'], self.config.path_data['t1_brain_mask'], suffix='_brain', ofolder=self.config.path_data['output_folder'])
				else:
					fname_t1_brain, fname_t1_brain_mask = segment_brain(self.config.path_data['t1'], bool(self.config.config_params['t1']['centered_fov']), ofolder=self.config.path_data['output_folder'])
			##
			fname_gm, fname_wm = segment_gm_wm(fname_t1_brain, ofolder=self.config.path_data['output_folder'])
			#
			self.config.path_data['t1_gm_mask'] = fname_gm
			self.config.path_data['t1_wm_mask'] = fname_wm
			#
		## 
		# GET VALUES IN ROI
		if self.config.processing['Get_values_roi']:
			# Register template to T1 and metrics
			reg = Registration(self.config.path_data['t1'], self.config.path_data['t1_brain_mask'], self.config.path_data['cvr'], self.config.path_data['asl_pre_brain_mask'], fname_cbf_pre=self.config.path_data['cbf_pre'], fname_cbf_post=self.config.path_data['cbf_post'], dict_param=self.config.roi, fname_gm=self.config.path_data['t1_gm_mask'], fname_wm=self.config.path_data['t1_wm_mask'], ofolder=self.config.path_data['output_folder'])
			reg.pipeline()
			self.config.path_data['asl_pre_brain_mask'] = reg.fname_metric_mask
		#
		##
		# SAVE CONFIG FILE TO STORE OUT FILES INFO 
		self.config.save()

if __name__ == "__main__":
	parser = get_parser()
	command_line_args = parser.parse_args()
	#
	# command_line_args = parser.parse_args(['--config', '/Users/sdupont/Documents/HIV/vasoreactivity/2018-07-24-test_processing_andy/config_test.json'])	#_no_comment
	#
	pipeline = Pipeline(command_line_args)
	pipeline.process()


