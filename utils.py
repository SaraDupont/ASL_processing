#!/usr/bin/env python
# Script for MRI data processing
# Author: Sara Dupont <sara.m.dupont@gmail.com>
import os, shutil, copy, sys, argparse
import nibabel as nib
import numpy as np 

def create_folder(path_folder):
    if not os.path.isdir(path_folder):
        try:
            os.mkdir(path_folder)
        except:
            pass
    path_folder = os.path.abspath(path_folder)
    #
    return path_folder

def input_file(path_file):
    '''Checks if input file exists. Raises an error if it doesn't. To be used as a parser type.
    :param path_file: a str corresponding to the path to the input file
    :return: the path to the existing input file
    '''
    # convert path to ouput file to an absolute path
    path_file = os.path.abspath(path_file)
    # check if file exists
    if not os.path.isfile(path_file):
        raise IOError('input file does not exist: %s' %path_file)
    else:
        return path_file


def list_param(s):
    return s.split(',')

def bool_param(b):
    return bool(int(b))


##############################################################################################
# class image & image static functions
class Image():
    def __init__(self, fname):
        im = nib.load(fname)
        self.fname = fname
        self.path = '/'.join(fname.split('/')[:-1]) + '/' if fname.split('/')[:-1] != [] else './'
        self.file = fname.split('/')[-1].split('.')[0]
        self.ext = '.'+'.'.join(fname.split('/')[-1].split('.')[1:])
        self.data = im.get_data()
        self.hdr = im.get_header()
        self.orientation = get_orientation(self.hdr)
        self.dim = self.data.shape
        self.pixdim = self.hdr.get_zooms()
    #
    def set_data(self, data):
        self.data = data
        self.dim = self.data.shape
        self.hdr.set_data_shape(self.dim)
    #
    def set_fname(self, fname):
        self.fname = fname
        self.path = '/'.join(fname.split('/')[:-1]) + '/'
        self.file = fname.split('/')[-1].split('.')[0]
        self.ext = '.'+'.'.join(fname.split('/')[-1].split('.')[1:])
    #
    def set_path(self, path):
        self.path = path
        self.fname = os.path.join(self.path, self.file+self.ext)
    #
    def copy(self):
        return copy.deepcopy(self)
    #
    def save(self):
        im_to_save = nib.Nifti1Image(self.data, None, self.hdr)
        print 'Saving ', self.fname
        nib.save(im_to_save, self.fname)

def get_orientation(hdr):
    from nibabel import orientations
    orientation_dic = {
        (0, 1): 'L',
        (0, -1): 'R',
        (1, 1): 'P',
        (1, -1): 'A',
        (2, 1): 'I',
        (2, -1): 'S',
    }
    orientation_matrix = orientations.io_orientation(hdr.get_best_affine())
    ori = "".join([orientation_dic[tuple(i)] for i in orientation_matrix])
    return ori

def add_suffix(path, suffix):
	list_path = path.split('/')
	fname = list_path[-1]
	list_fname = fname.split('.')
	list_fname[0]+= suffix
	fname_suffix = '.'.join(list_fname)
	list_path[-1] = fname_suffix
	path_suffix = '/'.join(list_path)
	return path_suffix

def mask_image(fname_im, fname_mask, suffix='_mask', ofolder=None):
    im = Image(fname_im)
    mask = Image(fname_mask)
    #
    im.data = im.data*mask.data
    #
    fname_im_masked = add_suffix(im.fname, suffix)
    im.set_fname(fname_im_masked)
    if ofolder is not None and os.path.isdir(ofolder):
        im.set_path(ofolder)
        fname_im_masked = im.fname
    im.save()
    #
    return fname_im_masked

def make_tmp_dir(path, tmp_dir):
    path_tmp_dir = os.path.join(path, tmp_dir)+'/'
    print 'Creating folder: ', path_tmp_dir
    if tmp_dir not in os.listdir(path):
        os.mkdir(path_tmp_dir)
    else:
        print "Folder "+tmp_dir+" already exists in "+path
    # 
    return path_tmp_dir

def replace_char(path, char_list="+'( )%"):
    for c in char_list+'"':
        path = path.replace(c, "\\"+c)
    #
    return path


if __name__ == "__main__":
    print "This is a utils script, it is not supposed to be called in the terminal"




