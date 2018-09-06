import os, argparse, commands
from process_breathing import *
import process_asl #import *
from process_asl import Image 

## TODO: ADD ASL ACQUISITION PARAMS IN THE PARSER (AND FEED IT TO THE CBF FUNCTIONS)
def get_parser():
    parser = argparse.ArgumentParser(description='Process ASL data to get cerebral blood flow (CBF) data (or)')
    parser.add_argument("-smartlab",
                        help="Filename of the Smartlab breathing data (with partial CO2 measurements). If no Smartlab file is provided, type of challenge will be assumed to be diamox (i.e. CVR will not be normalized by the difference in ETCO2).",
                        type=str,
                        dest="fname_smartlab",
                        required=False)
    parser.add_argument("-asl-pre",
                        help="Filename of the ASL raw data BEFORE CO2 challenge (nifti file)",
                        type=str,
                        dest="fname_asl_pre",
                        required=True)
    parser.add_argument("-asl-post",
                        help="Filename of the ASL raw data DURING CO2 challenge (nifti file)",
                        type=str,
                        dest="fname_asl_post",
                        required=True)
    parser.add_argument("-type-asl",
                        help="type of ASL data: pulsed ASL (pasl), pseudo-continous ASL (pcasl) or pre-computed CBF maps (cbf)",
                        type=str,
                        dest="type_asl",
                        choices=['pasl', 'pcasl', 'cbf'],
                        default='pasl')
    parser.add_argument("-ofolder",
                        help="Folder to output results",
                        type=create_folder,
                        dest="ofolder",
                        default="./")
    #
    return parser



def get_mask(fname_im, radius=60):
	#bet raw_asl.nii.gz raw_asl_brain -m -n -r 65
	file_name = fname_im.split('.')[0]
	ext = '.'.join(fname_im.split('.')[1:])
	file_name_brain = file_name+'_brain'
	file_name_brain_mask = file_name+'_brain_mask'
	#
	if not os.path.isfile(file_name_brain_mask+'.'+ext):
		print 'Getting brain mask with BET'
		cmd = 'bet '+fname_im+' '+file_name_brain+' -m -n -r '+str(radius)
		s, o = commands.getstatusoutput(cmd)
	else:
		print 'Loading existing brain mask'
		s=0
	#
	if s == 0:
		im_mask = Image(file_name_brain_mask+'.'+ext)
	else:
		print 'ERROR during brain extraction'
		im_mask = None
	return im_mask

def set_path(im, path):
	im.set_fname(os.path.join(path, im.file+im.ext))
	return im


def get_cbf_pre_post(param, moco=True, order='ct'):
	im_asl_pre = Image(param.fname_asl_pre)
	im_asl_post = Image(param.fname_asl_post)
	#
	im_mask_pre = get_mask(param.fname_asl_pre)
	im_mask_post = get_mask(param.fname_asl_post)
	#
	im_asl_pre = set_path(im_asl_pre, param.ofolder)
	im_asl_post = set_path(im_asl_post, param.ofolder)
	im_mask_pre = set_path(im_mask_pre, param.ofolder)
	im_mask_post = set_path(im_mask_post, param.ofolder)
	#
	ti = 2.1
	#
	parser_asl = process_asl.get_parser()
	param_asl = parser_asl.parse_args(['-asl', im_asl_pre.fname, '-mask', im_mask_pre.fname, '-ti', str(ti), '-type-asl', param.type_asl, '-order', order_asl])
	#
	asl_pre = process_asl.AslData(param_asl)
	asl_pre.process()
	im_cbf_pre = asl_pre.im_cbf 
	#
	param_asl.fname_asl = im_asl_post.fname
	param_asl.fname_mask = im_mask_post.fname
	asl_post = process_asl.AslData(param_asl)
	asl_post.process()
	im_cbf_post = asl_post.im_cbf
	# if param.type_asl == 'pasl':
	# 	im_cbf_pre = compute_cbf_pasl(im_asl_pre, order=order, moco=True, im_mask=im_mask_pre, lambd=0.9, ti=ti, t1_blood=1.65, alpha=0.98, bolus_dur=0.8, tr=4.6)
	# 	im_cbf_post = compute_cbf_pasl(im_asl_post, order=order, moco=True, im_mask=im_mask_post, lambd=0.9, ti=ti, t1_blood=1.65, alpha=0.98, bolus_dur=0.8, tr=4.6)
	# elif param.type_asl == 'pcasl':
	# 	im_cbf_pre = compute_cbf_pcasl(im_asl_pre, order=order, moco=True, im_mask=im_mask_pre, lambd=0.9, t1_blood=1.650, alpha=0.85, tau=1.800, pld=2.000)
	# 	im_cbf_post = compute_cbf_pcasl(im_asl_post, order=order, moco=True, im_mask=im_mask_post, lambd=0.9, t1_blood=1.650, alpha=0.85, tau=1.800, pld=2.000)
	# elif param.type_asl == 'cbf':
	# 	im_cbf_pre = im_asl_pre
	# 	im_cbf_post = im_asl_post
	#
	if moco: 
		fname_cbf_post_out = process_asl.add_suffix(im_cbf_post.fname, '_flirt_reg')
		cmd_moco = 'flirt -in '+im_cbf_post.fname+' -ref '+im_cbf_pre.fname+' -out '+fname_cbf_post_out
		s, o = commands.getstatusoutput(cmd_moco)
		if s==0:
			print 'Using motion corrected CBF post map: ', fname_cbf_post_out
			im_cbf_post = Image(fname_cbf_post_out)
		else: 
			print 'ERROR during motion corection --> keeping original CBF maps.'
	return im_cbf_pre, im_cbf_post, im_mask_pre


def get_cvr(im_cbf_pre, im_cbf_post, ave_et_pre, ave_et_post, im_mask=None, v=1):
	im_cvr = im_cbf_pre.copy()
	fname_cvr = im_cbf_pre.path + 'CVR' + im_cbf_pre.ext
	im_cvr.set_fname(fname_cvr)
	#
	if v != 0:
		print 'computing CVR ...'
	#
	dat_diff = im_cbf_post.data - im_cbf_pre.data
	#
	# pre0 = (im_cbf_pre.data == 0.0)
	pre0 = (abs(im_cbf_pre.data) < 0.1)
	#
	dat_pre_denom = copy.deepcopy(im_cbf_pre.data)
	dat_pre_denom[pre0] = 1.0
	dat_diff[pre0] = 0.0
	#
	# im_diff = im_cbf_pre.copy()
	# im_diff.set_fname(im_cbf_pre.path + 'DIFF_CBF' + im_cbf_pre.ext)
	# im_diff.set_data(dat_diff)
	# im_diff.save() #im_diff.fname)
	# #
	# im_asl_pre_denom =  im_cbf_pre.copy()
	# im_asl_pre_denom.set_fname(im_cbf_pre.path + 'CBF_PRE_DENOM' + im_cbf_pre.ext) 
	# im_asl_pre_denom.data = dat_pre_denom
	# im_asl_pre_denom.save()
	# formula from (Leoni 2012)
	if ave_et_pre is not None and ave_et_post is not None:
		diff_et_co2 = float(ave_et_post - ave_et_pre)
	else:
		diff_et_co2 = 1.0
	#
	dat_cvr = (100.0/diff_et_co2)*np.divide(dat_diff, dat_pre_denom)
	#
	if v!= 0 and diff_et_co2 != 1.0:
		print 'diff etco2: ', diff_et_co2
	#
	# if im_mask is  None: ## DO BET ON IM_CBF_PRE TO GET MASK AND MASK DATA ??
	if im_mask is not None:
		dat_cvr = dat_cvr * im_mask.data
	#
	#
	im_cvr.set_data(dat_cvr)
	im_cvr.data = dat_cvr
	im_cvr.save()
	return im_cvr

if __name__ == "__main__":

	parser = get_parser()
	param = parser.parse_args()

	# fname_smartlab = 'SmartLab_2-9-2018_testBen.txt'
	# fname_asl_pre = 'nifti_files/ASL_3D_ax_pre.nii.gz'
	# fname_asl_post = 'nifti_files/ASL_3D_ax_challenge.nii.gz'

	# param = parser.parse_args(['-smartlab', fname_smartlab, '-asl-pre', fname_asl_pre, '-asl-post', fname_asl_post])
	if param.fname_smartlab is not None: 
		end_tidal_pre, end_tidal_post = get_end_tidal_pre_post(param.fname_smartlab, param.ofolder)
		mean_et_pre = np.nanmean(end_tidal_pre)
		mean_et_post = np.nanmean(end_tidal_post)
	else:
		mean_et_pre = None
		mean_et_post = None		
	# #
	order_asl = 'tc'
	# order_asl = 'ct'
	im_cbf_pre, im_cbf_post, im_mask = get_cbf_pre_post(param, moco=True, order=order_asl)
	#
	im_cvr = get_cvr(im_cbf_pre, im_cbf_post, mean_et_pre, mean_et_post) #, im_mask = im_mask)




