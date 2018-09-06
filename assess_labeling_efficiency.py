import os, argparse, commands
import process_asl
from utils import *

def get_parser():
    parser_asl_data = process_asl.get_parser(add_help=False)
    #
    parser = argparse.ArgumentParser(description="Assess labeling efficiency of ASL from tCBF computed on Phase Contrast data and ASL data", parents=[parser_asl_data])
    parser.add_argument("-tot-cbf",
                    help="Total CBF value computed from the PC data (in mL/min).",
                    type=float,
                    dest="tCBF_PC",
                    required=True)
    parser.add_argument("-brain-density",
                    help="Brain density (in g/mL).",
                    type=float,
                    dest="rho",
                    required=False,
                    default = 1.06)
    parser.add_argument("-t1",
                    help="File name of the T1 image for brain volume estimation. (Provide a nifti file)",
                    type=str,
                    dest="fname_t1",
                    required=False)
    parser.add_argument("-t1-mask",
                    help="File name of the brain mask for the T1. (Provide a nifti file)",
                    type=str,
                    dest="fname_t1_mask",
                    required=False)
    #
    parser.add_argument("-ofile",
                help="Save the labeling efficiency to a text file ",
                type=str,
                dest="fname_alpha",
                required=False,
                default="alpha.txt")
    return parser


def get_volume_image(fname_im, unit='cm3'):
	im = Image(fname_im)
	# vol in mm3
	vol=np.sum(im.data)*np.prod(im.pixdim)
	# convert vol in cm3
	vol = vol/1000.0 if unit in ['cm3', 'cc', 'ml', 'mL'] else vol
	return round(vol, 2)


def get_aveCBF(tCBF_PC, fname_t1_mask=None, fname_t1=None, rho=1.06, ofolder='./'):
	#
	# get brain mask from T1 file if not provided
	if fname_t1_mask is None:
		assert os.path.isfile(str(fname_t1)), 'ERROR: cannot find T1 file: %s' %fname_t1
		# segment brain on t1 image 
		fname_t1_mask = os.path.join(ofolder, 't1_brain.nii.gz')
		cmd_bet = 'bet '+fname_t1+' '+fname_t1_mask+' -m -n '
		fname_t1_mask = add_suffix(fname_t1_mask, '_mask')
		#
		s, o = commands.getstatusoutput(cmd_bet)
		#
		if not os.path.isfile(fname_t1_mask):
			raise IOError('ERROR while segmenting brain on T1 file (%s):\n %s' %(fname_t1, o))
	else:
		assert os.path.isfile(fname_t1_mask), 'ERROR: cannot find T1 mask file: %s' %fname_t1_mask
	#
	# get brain volume from T1 scan (intracranial volume)
	brain_vol = get_volume_image(fname_t1_mask, unit='mL') # in mL
	#
	# gte brain mass
	brain_mass = brain_vol * rho # in g (rho = 1.06g/mL)
	#
	# get average whole brain blood flow
	tCBF_ave_PC = 100.0*tCBF_PC/brain_mass # in mL/min/100g
	#
	return tCBF_ave_PC


def get_corrected_alpha(tCBF_PC, dict_param_asl, fname_asl, fname_asl_mask=None, fname_t1_mask=None, fname_t1=None, rho=1.06, ofolder='./', fname_alpha='alpha.txt', v=1):
	# 
	tCBF_ave_PC = get_aveCBF(tCBF_PC, fname_t1_mask=fname_t1_mask, fname_t1=fname_t1, rho=rho, ofolder=ofolder)
	if v != 0:
		print "Total averaged CBF from PC is %s mL/min/100g" %tCBF_ave_PC
	#
	# process asl
	dict_param_asl['alpha'] = 1.0 ## assume labeling efficiency of 1.0 to get the ratio from PC data
	#
	asl_data = process_asl.AslData(dict_param_asl, fname_asl, fname_mask=fname_asl_mask, ofolder=ofolder)
	asl_data.process()
	#
	tCBF_ave_ASL = np.mean(asl_data.im_cbf.data[np.nonzero(asl_data.im_cbf.data)])
	#
	# get labeling efficiency
	alpha = tCBF_ave_ASL/tCBF_ave_PC
	#
	if v != 0:
		print 'Labeling efficiency: ', round(alpha, 3)
	#
	f = open(os.path.join(ofolder, fname_alpha), 'w')
	f.write(str(alpha))
	f.close()
	#
	return alpha

	

if __name__ == "__main__":
	parser = get_parser()
	param=parser.parse_args()
	#
	if param.fname_t1 is None and param.fname_t1_mask is None: 
		raise IOError('Please provide at least one of the two parameters -t1 or -t1-mask for brain volume estiamtion.')
	#
	# convert param into a dictionary
	dict_param_asl = process_asl.get_dict_param(param)
	# 
	alpha = get_corrected_alpha(param.tCBF_PC, dict_param_asl, param.fname_asl, fname_asl_mask=param.fname_mask, fname_t1_mask=param.fname_t1_mask, fname_t1=param.fname_t1, rho=param.rho, ofolder=param.ofolder, fname_alpha=param.fname_alpha)
	#



