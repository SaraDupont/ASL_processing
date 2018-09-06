#!/usr/bin/env python
# Script for ASL MRI data processing
# Author: Sara Dupont <sara.m.dupont@gmail.com>
import os, shutil, copy, sys, argparse
import nibabel as nib
import numpy as np 
from utils import * 


##############################################################################################
# parser & parser functions
def get_parser(add_help=True):
    parser = argparse.ArgumentParser(description='Process ASL data to get cerebral blood flow (CBF)', add_help=add_help)
    parser.add_argument("-asl",
                        help="Filename of the ASL raw data (nifti file)",
                        type=str,
                        dest="fname_asl",
                        required=True)
    parser.add_argument("-mask",
                        help="Filename of the brain mask for the ASL data (nifti file)",
                        type=str,
                        dest="fname_mask",
                        required=False)
    parser.add_argument("-pd",
                        help="Filename of the proton density image acquired with the same param as the ASL (M0 image) (nifti file). If not provided, the first control frame is used instead.",
                        type=str,
                        dest="fname_pd",
                        required=False)
    parser.add_argument("-type-asl",
                        help="Type of ASL scan: Pulsed ASL (pasl), pseudo-continuous ASL (pcasl)",
                        type=str,
                        dest="type_asl",
                        choices=['pasl', 'pcasl', 'ge-pcasl'],
                        default='pasl')
    parser.add_argument("-order",
                        help="Order of the tag and control frames in the ASL data, first control then tag (ct) or first tag then control (tc)",
                        type=str,
                        dest="order",
                        choices=['ct', 'tc'],
                        default='ct')
    # processing param
    parser.add_argument("-moco",
                        help="Do Motion correction for CBF computation. Input 0 or 1. ",
                        type=bool_param,
                        dest="moco",
                        required=False,
                        default=True)
    parser.add_argument("-smooth",
                        help="Do gaussian smoothing on each frame of the ASL data. Specify here the value of the 3D kernel used for smoothing. Input 0 to skip this step.",
                        type=float,
                        dest="smooth",
                        required=False,
                        default=1.5)
    parser.add_argument("-int-norm",
                        help="Do intensity normalization on each frame of the ASL data based on intensity of the first control frame. Specify here the goal mean intensity value (intensity averaged within brain). Input 0 to skip this step.",
                        type=float,
                        dest="int_norm",
                        required=False,
                        default=100.0)
    parser.add_argument("-rad",
                        help="Estimated radius of the brain in mm for BET initialization (see param -r in bet documentation)",
                        type=int,
                        dest="radius",
                        required=False,
                        default=55)
    #
    # acquisition param of ASL
    # pasl: ti=1.8, t1_blood=1.65, alpha=0.98, bolus_dur=0.8, tr=4.6
    # pcasl: lambd=0.9, t1_blood=1.650, alpha=0.85, tau=1.800, pld=2.000
    parser.add_argument("-lambda",
                        help="Blood brain partition coefficient (lambda value).",
                        type=float,
                        dest="lambd",
                        required=False,
                        default=0.9)
    parser.add_argument("-t1-blood",
                        help="T1 relaxation time of blood.",
                        type=float,
                        dest="t1_blood",
                        required=False,
                        default=1.65)
    parser.add_argument("-alpha",
                        help="Labelling efficiency.",
                        type=float,
                        dest="alpha",
                        required=False,
                        default=0.98)
    parser.add_argument("-ti",
                        help="Inversion time (for PASL only).",
                        type=float,
                        dest="ti",
                        required=False,
                        default=1.8)
    parser.add_argument("-bol-dur",
                        help="Bolus duration (for PASL only).",
                        type=float,
                        dest="bolus_dur",
                        required=False,
                        default=0.80)
    parser.add_argument("-tr",
                        help="Repetition time (for PASL only).",
                        type=float,
                        dest="tr",
                        required=False,
                        default=4.6)
    parser.add_argument("-pld",
                        help="Post labelling delay (for pcASL only).",
                        type=float,
                        dest="pld",
                        required=False,
                        default=2.0)
    parser.add_argument("-tau",
                        help="Label duration (for pcASL only).",
                        type=float,
                        dest="tau",
                        required=False,
                        default=1.8)
    # 
    # other arguments
    parser.add_argument("-ofolder",
                        help="Folder to output the created images.",
                        type=create_folder,
                        dest="ofolder",
                        default="./")
    parser.add_argument("-v",
                        help="verbose",
                        type=bool_param,
                        dest="verbose",
                        default=True)
    #
    return parser


def get_dict_param(param):
    dict_param = {
    "type_asl": param.type_asl, 
    "order": param.order, 
    "moco": param.moco, 
    "smooth": param.smooth, 
    "int-norm": param.int_norm, 
    "lambda": param.lambd,
    "t1-blood": param.t1_blood,
    "alpha": param.alpha, 
    "ti": param.ti, 
    "bolus-duration": param.bolus_dur, 
    "tr": param.tr, 
    "pld": param.pld, 
    "tau": param.tau,
    "bet_radius": param.radius 
    }
    #
    return dict_param

##############################################################################################\
# class ASL dataset
class AslData():
    def __init__(self, dict_param, fname_asl, fname_mask=None, ofolder=None, verbose=1):
        self.param = dict_param
        self.verbose = verbose
        #
        self.im_asl = Image(fname_asl)
        self.im_mask = Image(fname_mask) if fname_mask is not None else None
        self.im_pd = None
        #
        self.ofolder = self.im_asl.path if (ofolder is None or not os.path.isdir(ofolder)) else ofolder
        #
        print '\n\nOFOLDER: ', self.ofolder
        #
        self.im_asl.set_path(self.ofolder)
        if self.im_mask is not None:
            self.im_mask.set_path(self.ofolder)
        #
        self.path_tmp_files = make_tmp_dir(self.im_asl.path, 'tmp_dir_'+self.im_asl.file)

        self.list_indiv_im = []
        # get correct indices for control and tag frames
        self.i_control, self.i_tag = self.get_indices() if self.param['type_asl'] != 'ge-pcasl' else (None, None)
        #
        self.im_diff = None
        self.im_cbf = self.im_asl.copy()
        self.im_cbf.set_fname(add_suffix(self.im_asl.fname, '_cbf'))
        # 
    #
    # ##### main processing pipeline
    def process(self):
        # load list of tag and control images
        self.get_indiv_images()
        # get proton density image from 1st control if not in the inputs 
        self.im_pd = self.list_indiv_im[1] if self.im_pd is None else self.im_pd
        #
        # preprocess images
        self.preprocess()
        # 
        # compute cbf
        self.compute_cbf()
        #
    #
    #
    # ##### sub functions
    def get_indices(self):
        n_frames = self.im_asl.dim[-1]
        if self.param['order'] == 'ct':
            i_control = range(0, n_frames, 2)
            i_tag = range(1, n_frames, 2)
        elif self.param['order'] == "tc":
            i_control = range(1, n_frames, 2)
            i_tag = range(0, n_frames, 2)
        #
        return i_control, i_tag
    #
    #
    def get_indiv_images(self):
        if self.param['type_asl'] != 'ge-pcasl':
            if self.verbose:
                print '\n Saving individual frames'
            #
            for i in range(self.im_asl.data.shape[-1]):
                im_i = self.im_asl.copy()
                im_i.set_data(self.im_asl.data[:,:,:,i])
                #
                fname_i = add_suffix(self.im_asl.fname, '_frame'+str(i))
                im_i.set_fname(fname_i)
                #
                # set path to tmp output path:
                im_i.set_path(self.path_tmp_files) 
                # print im_i.fname
                # save indiv image
                im_i.save()
                #
                self.list_indiv_im.append(im_i)
        else:
            #
            im_diff = self.im_asl.copy()
	    print self.im_asl.data.shape
            im_diff.set_data(self.im_asl.data[:,:,:,0])
            im_diff.set_path(self.path_tmp_files) 
            im_diff.set_fname('im_diff.nii.gz')
            im_diff.save()
            self.im_diff = im_diff
            #
            im_pd = self.im_asl.copy()
            im_pd.set_data(self.im_asl.data[:,:,:,1])
            im_pd.set_path(self.path_tmp_files) 
            im_pd.set_fname('im_proton_density_m0.nii.gz')
            im_pd.save()
            self.im_pd = im_pd
        #
    #
    # 
    def preprocess(self):
        if self.verbose:
            print '\n Starting preprocessing'
        # Do motion correction
        if self.param['moco'] and self.param['type_asl'] != 'ge-pcasl':
            self.motion_correction()
        #
        # create mask if not given as an input
        if self.im_mask is None:
            self.get_mask()
        #
        # mask images to keep only signal inside the brain
        self.mask_images()
        #
        # perform gaussian smoothing on each frame 
        if self.param['smooth'] != 0 and self.param['type_asl'] != 'ge-pcasl':
            self.smooth_images()
        #
        # Intensity normalization (based on the intesnity of the first control frame)
        if self.param['int-norm'] != 0:
            self.intensity_normalization()
        #
    #
    #
    def compute_cbf(self):
        if self.verbose:
            print '\n Getting CBF'        
        # get the mean control-tag difference
        if self.im_diff is None:
            self.get_ave_diff()
        #
        # correct intensity of proton density image if the repetition time is too short
        if self.param['tr'] < 5:
            if self.verbose:
                print '... Correct proton density intensity for short TR' 
            # ref for t1_gm http://onlinelibrary.wiley.com/doi/10.1002/mrm.20605/epdf
            t1_gm = 1.82
            self.im_pd.data = self.im_pd.data /(1-np.exp(-self.param['tr']/t1_gm))
        #
        if self.param['type_asl'] == 'pasl':
            self.get_cbf_pasl()
        elif 'pcasl' in self.param['type_asl']:
            self.get_cbf_pcasl()
        #
        if self.verbose:
            print '... Saving CBF image' 
        self.im_cbf.save()
    #
    #
    # ##### pre processing functions
    def motion_correction(self):
        if self.verbose:
            print '... Motion correction' 
        # 
        list_im_moco = []
        i_ref = 1
        #
        for i, im in enumerate(self.list_indiv_im):
            if i != i_ref:
                fname_ref = self.list_indiv_im[i_ref].fname 
                fname_in = self.list_indiv_im[i].fname 
                fname_out = add_suffix(fname_in, '_reg_to_frame'+str(i_ref))
                # 
                if not os.path.isfile(fname_out):
					print 'flirt -in '+fname_in+' -ref '+fname_ref+' -out '+fname_out
					os.system('flirt -in '+fname_in+' -ref '+fname_ref+' -out '+fname_out)
                # 
                im_out = Image(fname_out)
                list_im_moco.append(im_out)
            else: 
                list_im_moco.append(self.list_indiv_im[i_ref])
        #
        self.list_indiv_im = list_im_moco
    #
    #
    def get_mask(self):
        if self.verbose:
            print '... Getting brain mask with BET on 1st control frame' 
        #
        fname_in = self.im_pd.fname # 
        fname_mask = add_suffix(fname_in, '_brain')
		# if fname_mask.split('.')[-1] == 'nii':
		#     fname_mask += '.gz'
        os.system('bet '+fname_in+' '+fname_mask+' -m -n -r '+str(self.param['bet_radius']))
        fname_mask = add_suffix(fname_mask, '_mask')
        self.im_mask = Image(fname_mask)
    #
    #
    def mask_images(self, update_name=True):
        if self.verbose:
            print '... Masking all frames'
        list_im = self.list_indiv_im if self.list_indiv_im != [] else [self.im_pd, self.im_diff]
        for im in list_im:
            im.data = im.data * self.im_mask.data
            if update_name:
                im.set_fname(add_suffix(im.fname, '_m'))
            im.save()
    #
    #
    def smooth_images(self):
        if self.verbose:
            print '... Gaussian smoothing with kernel of '+str(self.param['smooth'])+' mm' 
        list_im_smooth = []
        for im in self.list_indiv_im: 
            fname_in = im.fname
            fname_out = add_suffix(fname_in, '_smooth')
		    # if fname_out.split('.')[-1] == 'nii':
		    #     fname_out += '.gz'
            #
            os.system('fslmaths '+fname_in+' -s '+str(self.param['smooth'])+' '+fname_out)
            #
            im_out = Image(fname_out)
            list_im_smooth.append(im_out)
        self.list_indiv_im = list_im_smooth
        self.mask_images(update_name=False)
    #
    #
    def intensity_normalization(self):
        if self.verbose:
            print '... Intensity normalization based on 1st control frame intensity. Goal mean intensity within brain: '+str(self.param['int-norm']) 
        #
        data_ref = self.im_pd.data * self.im_mask.data
        mean_ref = np.mean(data_ref[np.nonzero(data_ref)])
        #
        list_im_norm = []
        #
        list_im = self.list_indiv_im if self.list_indiv_im != [] else [self.im_diff]
        #
        for im in list_im:
            fname_norm = add_suffix(im.fname, '_norm')
            im_norm = im.copy()
            im_norm.data = im.data * self.param['int-norm'] / mean_ref
            im_norm.set_fname(fname_norm)
            im_norm.save()
            #
            list_im_norm.append(im_norm)
        #
        self.list_indiv_im = list_im_norm if self.list_indiv_im != [] else self.list_indiv_im
        self.mask_images(update_name=False)
        self.im_diff = list_im_norm[0] if self.list_indiv_im == [] else self.im_diff 
    #
    #
    # ##### get CBF functions
    def get_ave_diff(self):
        if self.verbose:
            print '... Get mean control-tag difference' 
        self.im_diff = self.list_indiv_im[1].copy()
        self.im_diff.set_fname(add_suffix(self.im_asl.fname, '_control_tag_diff'))
        self.im_diff.set_path(self.path_tmp_files)
        #
        list_data_control = np.asarray([im.data for im in np.asarray(self.list_indiv_im)[self.i_control]])
        list_data_tag = np.asarray([im.data for im in np.asarray(self.list_indiv_im)[self.i_tag]])
        # average controls and tags
        data_control_mean = np.mean(list_data_control, axis=0)
        data_tag_mean = np.mean(list_data_tag, axis=0)
        # compute diff
        data_diff_mean = data_control_mean - data_tag_mean
        #
        self.im_diff.data = data_diff_mean
        self.im_diff.save()
    #
    #
    def get_cbf_pasl(self):
        if self.verbose:
            print '... Get CBF from PASL data' 
        #
        data_cbf = np.divide((6000*self.param['lambda']*self.im_diff.data*np.exp(self.param['ti']/self.param['t1-blood'])), (2*self.param['alpha']*self.param['bolus-duration']*self.im_pd.data))
        data_cbf = np.nan_to_num(data_cbf)
        #
        self.im_cbf.data = data_cbf
    #
    #
    def get_cbf_pcasl(self):
        if self.verbose:
            print '... Get CBF from PC-ASL data' 
        #
        data_cbf = np.divide(6000*self.param['lambda']*self.im_diff.data*np.exp(self.param['pld']/self.param['t1-blood']), 2*self.param['alpha']*self.param['t1-blood']*self.im_pd.data*(1-np.exp(-self.param['tau']/self.param['t1-blood'])))
        data_cbf = np.nan_to_num(data_cbf)
        data_cbf[data_cbf<-2000] = -2000
        data_cbf[data_cbf>2000] = 2000
        #
        self.im_cbf.data = data_cbf




##############################################################################################
# MAIN
def main():
    parser = get_parser()
    param=parser.parse_args()
    # param=parser.parse_args('-asl asl_pre.nii.gz'.split(' '))
    # param = parser.parse_args('-asl /Users/sdupont/Documents/HIV/vasoreactivity/2018-05-08-test_retest_jianxun/scan3_10%_new_cbf/asl_post.nii.gz -order tc'.split(' '))
    #
    dict_param = get_dict_param(param)
    #
    asl_data = AslData(dict_param, param.fname_asl, fname_mask=param.fname_mask, ofolder=param.ofolder)
    asl_data.process()



if __name__ == "__main__":
    main()




