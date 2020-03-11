#!/usr/bin/env python3
# 
# -*- coding: utf-8 -*-
# title           : kd.py
# description     : [description]
# author          : Adebayo B. Braimah
# e-mail          : adebayo.braimah@cchmc.org
# date            : 2020 03 11 19:07:33
# version         : 0.0.1
# usage           : kd.py [-h,--help]
# notes           : [notes]
# python_version  : 3.7.3
#==============================================================================

# Define usage
"""

Applies pre-computed (affine) linear and non-linear transforms to 4D EP (bold) image data in parallel.
Additional options include the ability to resample the image back to its native voxel size and/or native
space dimensions (in the case a native space image is provided).

This script at minimum requires that FSL be installed and added to system path.
Image resampling options required that MIRTK be installed and added to the system path.

Usage:
  reg4D.py [options | expert options] [required arguments]

Required arguments
    -i,--in INPUT       Input 4D file
    -o,--out PREFIX     Output prefix for transformed 4D file
    -r,--ref IMAGE      Reference image that corresponds to the transform(s) or
                        target space to be resampled/interpolated to

Options:
    --ref-tar IMAGE     Reference (target) image to resample to post-transform
                        (usually the native space image) - behaves best with
                        linear transfroms.
    -TR,--TR FLOAT      Repetition time (TR) of the input 4D EPI. If this value is not
                        specified, it will then be read from the nifti header [default: infer]
    -n,--num-jobs INT   Number of jobs to run in parallel. Default behavior is to use the
                        max number of cores available [default: infer]

Expert Options:
    -d,--dim CMD        (Command) Dimension to split and concatenate along the 4D image
                        (valid options: "-x", "-y", "-z", "t", "-tr") [default: -tr]
    -w,--warp IMAGE     FSL-style non-linear warp field file to reference image
    --warp-app CMD      (Command) Warp field treatment and application
                        (valid options: "relative", "absolute") [default: "relative"]
    --premat MAT        FSL-style pre-transform linear transformation matrix
    --postmat MAT       FSL-style post-transform linear transformation matrix
    --resamp-vox        Resample voxel-size to reference target image (if '--ref-tar' option is not
                        used then image voxel-size is resampled back to native voxel-size)
    --resamp-dim        Resample image dimension to reference target image (if '--ref-tar' option is not
                        used then image dimensions are resampled back to native dimensions)
    -m,--mask IMAGE     Binary mask image file in reference space to use for applying transform(s)
    --interp CMD        Interpolation method, options include: "nn","trilinear","sinc","spline"
    --padding-size INT  Extrapolates outside original volume by n voxels
    --use-qform         Use s/qforms of ref_vol and nii_file images 
                        - NOTE: no other transorms can be applied with this option
    --data-type CMD     Force output data type (valid options: "char" "short" "int" "float" "double")
    --super-sampling    Intermediary supersampling of output. [default: False]
    --super-level CMD   Level of intermediary supersampling, a for 'automatic' or integer level.
                        Only used when '--super-sampling' option is enabled. [default: 2]
    --no-parallel       Do not apply linear/non-linear transforms in parallel.
                        NOT RECOMMENDED for 4D neuroimage EPIs greater than 300 frames (TRs) or
                        for use with non-linear transforms as FSL's applywarp has to read the entire
                        timeseries into memory.

    -v,--verbose        Enable verbose output [default: False]
    -h,--help           Prints help message, then exits.
    --version           Prints version, then exits.

NOTE: '--resamp-dim' and '--resamp-vox' options require MIRTK to be installed if these options are specified.
"""

# Import modules and packages
import nibabel as nib
import os
import sys
import subprocess
import glob
import random
import shutil
import multiprocessing

# Import packages and modules for argument parsing
from docopt import docopt

# Define classes

class File():
    '''
    File class doc-string

    Attributes:
        file_name (string): Filename with path to file. File does not need to exist.
    '''

    def __init__(self, file_name):
        '''
        Return file_name string
        '''
        self.file_name = file_name

    def touch(file_name):
        '''
        Analagous to UNIX's touch command. Creates an empty file.

        Arguments:
            file_name (string): Filename (with path)

        Returns:
            file_name (string): Path to file. This now file exists and can now be used as inputs to other functions.
        '''

        if not os.path.exists(file_name):
            with open(file_name, "w"):
                pass
            file_name = os.path.abspath(file_name)
        else:
            print(f"file: {file_name} already exists. File will not be overwritten.")

        return file_name

    def make_filename_path(file_name):
        '''
        Analagous to UNIX's touch command. Creates an empty file, and then removes it - preserving the
        absolute path to the file.

        Arguments:
            file_name (string): Filename (with path)

        Returns:
            file_name (string): Path to file. This now file exists and can now be used as inputs to other functions.
        '''

        if not os.path.exists(file_name):
            with open(file_name, "w"):
                pass
            file_name = os.path.abspath(file_name)
            os.remove(file_name)
        else:
            print(f"file: {file_name} already exists. File will not be overwritten.")

        return file_name


class Command():
    '''
    Creates a command and an empty command list for UNIX command line programs/applications. Primary use and
    use-cases are intended for the subprocess module and its associated classes (i.e. run).

    Attributes:
        command: Command to be performed on the command line
    '''

    def __init__(self):
        '''
        Init doc-string for Command class.
        '''
        pass

    def init_cmd(self, command):
        '''
        Init command function for initializing commands to be used on UNIX command line.

        Arguments:
            command (string): Command to be used. Note: command used must be in system path

        Returns:
            cmd_list (list): Mutable list that can be appended to.
        '''
        self.command = command
        self.cmd_list = [f"{self.command}"]
        return self.cmd_list

# Define functions

def split_4D(nii_4d, out_prefix, dim="-t", verbose=False):
    '''
    Wrapper function for FSL's fslsplit. Splits 4D nifti timeseries data and returns a list of
    the split timeseries files.

    Arugments:
        nii_4d (nifti file): Input 4D file to be separated
        out_prefix (string): Output file prefix (no '.nii.gz', preferably to separate directory)
        dim (string): Dimension to separate along [default: '-t']
        verbose (bool, optional): Enable verbose output [default: False]

    Returns:
        nii_list (list): Sorted list of the zeropadded filenames of the split 4D timeseries.
    '''

    # Check dimension argument
    if dim != "-t" and dim != "-x" and dim != "-y" and dim != "-z":
        print("Unrecognized option for dimension. Using default of time.")
        dim = "-t"

    # Initialize command
    split = Command().init_cmd("fslsplit")

    # Add input/output arguments
    split.append(nii_4d)
    split.append(out_prefix)
    split.append(dim)

    if verbose:
        print(' '.join(split))

    # Perform/Execute command
    subprocess.call(split)

    # Create list of files
    nii_list = sorted(glob.glob(out_prefix + "*.nii*"))

    return nii_list


def merge_4D(nii_list, out_prefix, dim="-tr", tr=2.000, verbose=False):
    '''
    Wrapper function for FSL's fslmerge. Merges a list of 3D nifti files and returns a merged file in some
    arbitrary dimension (t (time, with TR), x, y, or, z).

    Arugments:
        nii_list (list): List of the zeropadded filenames of the split 4D timeseries.
        out_prefix (string): Output file prefix (no '.nii.gz')
        dim (string): Dimension to separate along [default: '-tr']
        tr (float): Repetition Time (TR, in sec.) [default: 2.000]
        verbose (bool, optional): Enable verbose output [default: False]

    Returns:
        out_file (string): Merged output nifti file.
    '''

    # Check dimension argument
    if dim != "-tr" and dim != "-t" and dim != "-x" and dim != "-y" and dim != "-z":
        print("Unrecognized option for dimension. Using default of time with defined TR.")
        dim = "-tr"

    # Create empty file and path to retrieve
    out_prefix = File.make_filename_path(out_prefix)

    # Initialize command
    merge = Command().init_cmd("fslmerge")

    # Add input/output arguments
    merge.append(dim)
    merge.append(out_prefix)
    merge.extend(nii_list)

    if dim == "-tr":
        merge.append(str(tr))

    if verbose:
        print(' '.join(merge))

    # Perform/Execute command
    subprocess.call(merge)

    out_file = out_prefix + ".nii.gz"

    return out_file


def apply_xfm_3D(nii_file, out_prefix, ref_vol, warp="", warp_app="relative", premat="", postmat="",
                 mask="", interp="", padding_size="", use_qform=False, data_type="", super_sampling=False,
                 super_level=2, verbose=False):
    '''
    Wrapper function for FSL's applywarp. Applies linear and non-linear transforms (xfms) in a one-step manner on
    3D images. Generally, this reduces image blurring as a result of reduced image resampling accross the application of
    multiple transforms.

    Note:
    - If wanting to apply more than two linear transforms, then the first two transformation matrices should be concatenated (via FSL's convert_xfm)
    - If wanting to apply more than one non-linear warp field, then the two warps should be combined (via FSL's convertwarp)

    Arugments:
        nii_file (nifti file): Input nifti file
        out_prefix (string): Output file (no '.nii.gz')
        ref_vol (nifti file): Reference volume
        warp (nifti file, optional): Warp file to reference
        warp_app (nifti file, optional): Warp field treatment. Valid options are: "relative" (default), and "absolute"
        premat (FSL mat file, optional): Pre-transform linear transformation matrix
        postmat (FSL mat file, optional): Post-transform linear transformation matrix
        mask (nifti file, optional): Mask, in reference space
        interp (string, optional): Interpolation method, options include: "nn","trilinear","sinc","spline"
        padding_size (int, optional): Extrapolates outside original volume by n voxels
        use_qform (bool, optional): Use s/qforms of ref_vol and nii_file images
        data_type (string,optional): Force output data type ["char" "short" "int" "float" "double"].
        super_sampling (bool, optional): Intermediary supersampling of output [default: False]
        super_level (int, optional): Level of intermediary supersampling, a for 'automatic' or integer level. [default: 2]
        verbose (bool, optional): Enable verbose output [default: False]

    Returns:
        out_file (string): Output file with the applied transform(s).
    '''

    # Create empty file and path to retrieve
    out_prefix = File.make_filename_path(out_prefix)

    # Initialize command
    xfm = Command().init_cmd("applywarp")

    # Add required input/output arguments
    xfm.append(f"--in={nii_file}")
    xfm.append(f"--out={out_prefix}")
    xfm.append(f"--ref={ref_vol}")

    # Add optional input arguments
    if warp:
        xfm.append(f"--warp={warp}")

    if warp_app == "relative":
        xfm.append("--rel")
    elif warp_app == "absolute":
        xfm.append("--abs")

    if premat:
        xfm.append(f"--premat={premat}")

    if postmat:
        xfm.append(f"--postmat={postmat}")

    if data_type:
        if data_type == "char" or data_type == "short" or data_type == "int" or data_type == "float" or data_type == "double":
            xfm.append(f"--datatype={data_type}")
        else:
            print(f"Unrecognized option: {data_type}. Using input data type.")

    if mask:
        xfm.append(f"--mask={mask}")

    if super_sampling:
        xfm.append("--super")
        if super_level:
            xfm.append(f"--superlevel={super_level}")

    if interp:
        if interp == "nn" or interp == "trilinear" or interp == "sinc" or interp == "spline":
            xfm.append(f"--interp={interp}")

    if padding_size:
        xfm.append(f"--paddingsize={padding_size}")

    if use_qform:
        xfm.append("--usesqform")

    if verbose:
        xfm.append("--verbose")
        print(' '.join(xfm))

    # Perform/Execute command
    subprocess.call(xfm)

    out_file = out_prefix + ".nii.gz"

    return out_file


def get_img_dims(nii_file):
    '''
    Reads and stores nifti image dimensions from input nifti file.

    Arugments:
        nii_file (nifti file): Input nifti file

    Returns:
        n_dim (int): N-dimensions of image (i.e. 2, for 2D, 3 for 3D, and 4 for 4D)
        x_dim (int): X-dimension length
        y_dim (int): Y-dimension length
        z_dim (int): Z-dimension length
        n_vol (int): Number of volumes in image (1 for 2D and 3D images, N>1 for 4D images)
        x_vox (float): X-dimension voxel size
        y_vox (float): Y-dimension voxel size
        z_vox (float): Z-dimension voxel size
        tr (float): Repetition time from nifti header
    '''

    # Load image
    img = nib.load(nii_file)

    # Image dimensions
    n_dim = int(img.header['dim'][0])
    x_dim = int(img.header['dim'][1])
    y_dim = int(img.header['dim'][2])
    z_dim = int(img.header['dim'][3])
    n_vol = int(img.header['dim'][4])

    # Image voxel size
    x_vox = float(img.header['pixdim'][1])
    y_vox = float(img.header['pixdim'][2])
    z_vox = float(img.header['pixdim'][3])
    tr = float(img.header['pixdim'][4])

    return n_dim, x_dim, y_dim, z_dim, n_vol, x_vox, y_vox, z_vox, tr


def img_res(nii_file, out_prefix, ref, resamp_vox=True, resamp_dim=True, verbose=False):
    '''
    Wrapper function for MIRTK's image resample-image binary. Performs image resampling of an image in referance
    to another image. This function is primarily intended for resampling a transformed image back to its
    original native space - implying that the reference (ref) image should be the original native space image.

    Note:
    - Not all of MIRTK's resample-image options are present. This is to ensure consistent behavior.
    - Only intended for 2D or 3D images. If the image is 4D, this function will exit as this will likely consume too much memory.

    Arugments:
        nii_file (nifti file): List of the zeropadded filenames of the split 4D timeseries.
        out_prefix (string): Output file prefix (no '.nii.gz')
        ref (nifti file): (Native space) Reference volume
        resamp_vox (bool, optional): Resample voxel size to reference image [default: True]
        resamp_dim (bool, optional): Resample image dimension to reference image [default: True]
        verbose (bool, optional): Enable verbose output [default: False]

    Returns:
        out_file (string): Merged output nifti file.
    '''

    # Create empty file and path to retrieve
    out_prefix = File.make_filename_path(out_prefix)

    # Initialize command
    rsm = Command().init_cmd("mirtk")
    rsm.append("resample-image")  # append mirtk sub-command

    # Get image dimensions (input)
    [n_dim, x_dim, y_dim, z_dim, n_vol, x_vox, y_vox, z_vox, tr] = get_img_dims(nii_file=nii_file)

    if n_dim > 3:
        print(f"Image {nii_file} is greater than 3 dimensions. Exiting")
        sys.exit(1)

    # Get image dimensions (reference)
    [n_dim, x_dim, y_dim, z_dim, n_vol, x_vox, y_vox, z_vox, tr] = get_img_dims(nii_file=ref)

    if n_dim > 3:
        print(f"Image {ref} is greater than 3 dimensions. Exiting")
        sys.exit(1)

    # Add required input/output arguments
    if '.nii.gz' not in nii_file:
        # print(True)
        nii_file = nii_file + '.nii.gz'
    rsm.append(nii_file)

    out_prefix = out_prefix + ".nii.gz"
    rsm.append(out_prefix)

    if resamp_vox:
        rsm.append("-size")
        rsm.append(str(x_vox))
        rsm.append(str(y_vox))
        rsm.append(str(z_vox))

    if resamp_dim:
        rsm.append("-imsize")
        rsm.append(str(x_dim))
        rsm.append(str(y_dim))
        rsm.append(str(z_dim))

    if verbose:
        rsm.append("-verbose")
        print(' '.join(rsm))

    # Perform/Execute command
    subprocess.call(rsm)

    out_file = out_prefix

    return out_file


def multiproc_cmd_list_tuple(in_list, out_list, **kwargs):
    '''
    Creates a tuple of lists that contains the input arguments in the format that multiprocessing's Pool uses.
    This function is intended for use with multiprocessing's Pool class and reg4D's apply_xfm_4D function.

    Arguments:
        in_list (list): Input list of filenames
        out_list (list): Output list of filenames
        **kwargs (string,dict, key,value pairs): key=value pairs - Most useful when function arguments are passed as a dictionary

    Returns:
        cmd_tup (tuple): Tuple of (command) lists that can be used as input to  multiprocessing's Pool class.
    '''

    # Empty list
    cmd_list = list()

    # for idx,file in enumerate(in_list):
    for idx in range(0, len(in_list), 1):
        tmp_list = [in_list[idx], out_list[idx]]
        for key, item in kwargs.items():
            tmp_list.append(f"{item}")
        cmd_list.append(tmp_list)

    cmd_tup = tuple(cmd_list)
    return cmd_tup


def make_out_list(out_name, max_val):
    '''
    Creates an output name list provided an output name and a maximum value (not inclusive) to iterate to.

    Arguments:
        out_name (string): Output name prefix
        max_val (int): Maximum value (not inclusive) to iterate to

    Returns:
        name_list (list): Output names list with zeropadded numbers added to the end.
    '''

    name_list = list()

    for idx in range(0, max_val, 1):
        tmp_num = '{:04}'.format(int(idx))
        tmp_name = out_name + f"_{tmp_num}"
        name_list.append(tmp_name)

    return name_list


def cp_img(file, dest_name):
    '''
    Copies a gzipped nifti image file. Primarily intended for copying single file image data.

    Arguments:
        file (nifit file): File path to source (image) file
        dest_name (string): Absolute path output file name.

    Returns:
        out_file (string): Absolute path to output file.
    '''

    if '.nii.gz' not in file:
        file = file + '.nii.gz'

    if '.nii.gz' not in dest_name:
        out_file = dest_name + '.nii.gz'

    shutil.copy(file, out_file)

    return out_file


def cp_dir_imgs(in_list, out_list, max_val):
    '''
    Copies a directory of image data provided an input, and output list of filenames in addition to a maximum
    value (not inclusive) to iterate to.

    No files are returned, as they are copied to some output directory.

    Arguments:
        in_list (list): Input name list
        out_list (list): Output name list
        max_val (int): Maximum value (not inclusive) to iterate to

    Returns:
        None
    '''

    for idx in range(0, max_val, 1):
        in_file = in_list[idx]
        out_file = out_list[idx]
        cp_img(file=in_file, dest_name=out_file)


def construct_args(list_opts, var_opts, dictionary=dict()):
    '''
    Constructs argument dictionary from an input options list and a corresponding variables list. Both lists must be
    of the same length. An input dictionary can be specified, in which case - a new dictionary is returned.

    Arguments:
        list_opts (list): Input options list
        var_opts (list): Corresponding variables list for each option
        dictionary (dict): Input dictionary that can be appended to

    Returns:
        new_dict (dict): Updated argument dictionary
    '''

    assert len(list_opts) == len(var_opts)

    new_dict = dictionary.copy()

    for idx in range(0, len(list_opts), 1):
        opt = list_opts[idx]
        var = var_opts[idx]
        if not var:
            var = ""
        tmp_dict = {opt: var}
        new_dict.update(tmp_dict)

    return new_dict

def which(program):
    '''
    Mimics UNIX 'which' command.

    Arguments:
        program (string): Input program/executable name as a string

    Returns:
        program, exe_file, or None: Searches system path and returns the requested executable/file. If it does not exist, None is returned.
    '''

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def make_ref_tar_img(nii_file, out_prefix, verbose=False):
    '''
    Wrapper function for fslmaths. Creates a mean image of an input 4D EPI timeseries to
    create a 'reference target image' - should one not be provided.

    Arugments:
        nii_file (file,string): Input 4D EPI timeseries.
        out_prefix (file,string): Output file prefix (no '.nii.gz')
        verbose (bool, optional): Enable verbose output [default: False]

    Returns:
        out_file (file,string): Mean image of input nifti file.
    '''

    # Create empty file and path to retrieve
    out_prefix = File.make_filename_path(out_prefix)

    # Initialize command
    mean_img = Command().init_cmd("fslmaths")

    # Add input/output arguments
    mean_img.append(nii_file)
    mean_img.append("-Tmean")
    mean_img.append(out_prefix)

    if verbose:
        print(' '.join(mean_img))

    # Perform/Execute command
    subprocess.call(mean_img)

    out_file = out_prefix + ".nii.gz"

    return out_file


def make_img_mask(nii_file, out_prefix, verbose=False):
    '''
    Wrapper function for fslmaths. Creates a binary mask for a nifti image file.

    Arugments:
        nii_file (file,string): Input 4D EPI timeseries.
        out_prefix (file,string): Output file prefix (no '.nii.gz')
        verbose (bool, optional): Enable verbose output [default: False]

    Returns:
        out_file (file,string): Binary mask image of input nifti file.
    '''

    # Get image dimensions
    [n_dim, x_dim, y_dim, z_dim, n_vol, x_vox, y_vox, z_vox, tr] = get_img_dims(nii_file=nii_file)

    # Create empty file and path to retrieve
    out_prefix = File.make_filename_path(out_prefix)

    # Initialize command
    mask_img = Command().init_cmd("fslmaths")

    # Add input/output arguments
    mask_img.append(nii_file)
    mask_img.append("-bin")
    mask_img.append(out_prefix)

    if verbose:
        print(' '.join(mask_img))

    # Perform/Execute command
    subprocess.call(mask_img)

    out_file = out_prefix + ".nii.gz"

    return out_file


if __name__ == '__main__':

    # Parse arguments
    args = docopt(__doc__, help=True, version='reg4D v0.0.1', options_first=False)
    # print(args)

    # Check for required arguments
    if not args['--in'] and not args['--out'] and not args['--ref']:
        print("")
        print("Usage:   reg4D.py --in IMAGE --out PREFIX --ref IMAGE   |   -h,-help,--help")
        print("")
        print("Please see help menu for details.")
        print("")

    # Check if FSL is installed
    fsl_install = which('applywarp')

    if not fsl_install:
        print("")
        print("FSL's applywarp is either not installed (properly) or added to the system path. Exiting.")
        sys.exit(1)

    # Check if MIRTK is installed (if image resampling options are selected)
    if args['--resamp-vox'] or args['--resamp-dim']:
        mirtk_install = which('mirtk')
        if not mirtk_install:
            print("")
            print("MIRTK is either not installed (properly) or added to the system path.")
            print("")
            print("Disabling any image resampling options")
            args['--resamp-vox'] = False
            args['--resamp-dim'] = False

    def apply_xfm_4D(nii_file, out_prefix, ref_xfm, ref_target="", num_jobs="infer", dim="-tr",
                     tr="infer", warp="", warp_app="relative", premat="", postmat="", resamp_vox=True,
                     resamp_dim=False, mask="", interp="", padding_size="", use_qform=False, data_type="",
                     super_sampling=False, super_level=2, verbose=False):
        '''
        Applies pre-computed (affine) linear and non-linear transforms to 4D EP (bold) image data in parallel.
        Additional options include the ability to resample the image back to its native voxel size and/or native
        space dimensions (in the case a native space image is provided).

        Arguments:
            nii_file (nifti file): Input 4D nifti file
            out_prefix (filename): Ouput prefix filename (and path, no '.nii.gz')
            ref_xfm (nifti file): Reference volume that corresponds to the transform(s)
            ref_target (nifti file, optional): Target nifti file to be resampled to (usually native space)
            num_jobs (int, optional): Number of jobs to run simulataneously [default: Number of availble cores]
            dim (string, optional): Dimension to separate and merge along [default: '-tr']
            tr (float, optional): Repetition Time (TR, in sec.) [default: Inferred from header or 2.000]
            warp (nifti file, optional): Warp file to reference
            warp_app (nifti file, optional): Warp field treatment. Valid options are: "relative" (default), and "absolute"
            premat (FSL mat file, optional): Pre-transform linear transformation matrix
            postmat (FSL mat file, optional): Post-transform linear transformation matrix
            resamp_vox (bool, optional): Resample voxel size to reference image [default: True, if ref_target is provided]
            resamp_dim (bool, optional): Resample image dimension to reference image [default: False]
            mask (nifti file, optional): Mask, in reference space to use for applying transform(s)
            interp (string, optional): Interpolation method, options include: "nn","trilinear","sinc","spline"
            padding_size (int, optional): Extrapolates outside original volume by n voxels
            use_qform (bool, optional): Use s/qforms of ref_vol and nii_file images - NOTE: no other transorms can be applied with this option
            data_type (string,optional): Force output data type ["char" "short" "int" "float" "double"]
            super_sampling (bool, optional): Intermediary supersampling of output [default: False]
            super_level (int, optional): Level of intermediary supersampling, a for 'automatic' or integer level [default: 2]
            verbose (bool, optional): Enable verbose output [default: False]
        Returns:
            out_file (nifti file): Transformed 4D nifti file
        '''

        # max number for random number generator
        n = 10000  # maximum N

        # Make temporay directories
        out_tmp = os.path.join(os.path.dirname(out_prefix), 'tmp_dir' + str(random.randint(0, n)))

        out_tmp_1 = os.path.join(os.path.abspath(out_tmp), 'tmp_dir' + str(random.randint(0, n)))
        out_tmp_2 = os.path.join(os.path.abspath(out_tmp), 'tmp_dir' + str(random.randint(0, n)))

        if not os.path.exists(out_tmp_1):
            os.makedirs(out_tmp_1)

        if not os.path.exists(out_tmp_2):
            os.makedirs(out_tmp_2)

        # Get absolute paths
        out_tmp = os.path.abspath(out_tmp)
        out_tmp_1 = os.path.abspath(out_tmp_1)
        out_tmp_2 = os.path.abspath(out_tmp_2)

        out_name = os.path.join(out_tmp, 'func_img')
        out_name_1 = os.path.join(out_tmp_1, 'func_split')
        out_name_2 = os.path.join(out_tmp_2, 'func_xfm')

        # Get image info
        [n_dim, x_dim, y_dim, z_dim, n_vol, x_vox, y_vox, z_vox, tr_img] = get_img_dims(nii_file=nii_file)

        # Check split dimension variable
        if dim == "-tr":
            split_dim = "-t"
        else:
            split_dim = dim

        # Get number of jobs for parallelization
        if num_jobs == 'infer':
            num_jobs = multiprocessing.cpu_count()
            if verbose:
                print("")
                print(f"Using {num_jobs} cores")

        # Get TR
        if tr == 'infer':
            tr = tr_img
        if verbose:
            print("")
            print(f"TR is: {tr} sec.")

        # Create initial images
        if not ref_target:
            ref_target = make_ref_tar_img(nii_file=nii_file, out_prefix=out_name + "_mean", verbose=verbose)

        ref_mask = make_img_mask(nii_file=ref_xfm, out_prefix=out_name + "_ref_mask", verbose=verbose)

        if resamp_vox or resamp_dim:
            ref_tar_rsm = File.make_filename_path(out_prefix=out_name + "_ref_tar_rsm")
            ref_tar_rsm_mask = File.make_filename_path(out_prefix=out_name + "_ref_tar_rsm_mask")
        else:
            ref_tar_rsm = ref_xfm
            ref_tar_rsm_mask = ref_mask

        # Make keyword argument dictionaries
        # Command 1 dictionary
        # Options list and variables
        list_opts_1 = ["warp", "warp_app", "premat", "postmat", "mask", "interp", "padding_size", "use_qform",
                       "data_type", "super_sampling", "super_level", "verbose"]
        var_opts_1 = [warp, warp_app, premat, postmat, ref_mask, interp, padding_size, use_qform, data_type,
                      super_sampling,
                      super_level, verbose]

        cmd_dict_1 = {"ref_vol": ref_xfm}  # Required argument(s) not in command list
        cmd_dict_1 = construct_args(list_opts=list_opts_1, var_opts=var_opts_1, dictionary=cmd_dict_1)

        # Command 2 dictionary
        # Options list and variables
        list_opts_2 = ["resamp_vox", "resamp_dim", "verbose"]
        var_opts_2 = [resamp_vox, resamp_dim, verbose]

        if resamp_vox or resamp_dim:
            cmd_dict_2 = {"ref": ref_target}  # Argument(s) not in command list
            cmd_dict_2 = construct_args(list_opts=list_opts_2, var_opts=var_opts_2, dictionary=cmd_dict_2)
        else:
            cmd_dict_2 = dict()

        # Command 3 dictionary
        # Options list and variables
        list_opts_3 = ["warp", "warp_app", "premat", "postmat", "mask", "interp", "padding_size", "use_qform",
                       "data_type", "super_sampling", "super_level", "verbose"]
        var_opts_3 = [warp, warp_app, premat, postmat, ref_tar_rsm_mask, interp, padding_size, use_qform, data_type,
                      super_sampling,
                      super_level, verbose]

        cmd_dict_3 = {"ref_vol": ref_tar_rsm}  # Required argument(s) not in command list
        cmd_dict_3 = construct_args(list_opts=list_opts_3, var_opts=var_opts_3, dictionary=cmd_dict_3)

        # Perform initial transform
        img_xfm = apply_xfm_3D(nii_file=ref_target, out_prefix=out_name + "_xfm", **cmd_dict_1)

        # Perform image resampling
        if resamp_vox or resamp_dim:
            img_rsm = img_res(nii_file=img_xfm, out_prefix=out_name + "_xfm_rsm", **cmd_dict_2)
        else:
            img_rsm = img_xfm

        if num_jobs != 1:
            # Split 4D EP image timeseries
            if verbose:
                print("")
                print("Splitting 4D timeseries")
            pre_xfm_list = split_4D(nii_4d=nii_file, out_prefix=out_name_1, dim=split_dim, verbose=verbose)

            # Parallelize process(es)
            out_tmp_names_list_2 = make_out_list(out_name=out_name_2, max_val=n_vol)
            cmd_tup_2 = multiproc_cmd_list_tuple(in_list=pre_xfm_list, out_list=out_tmp_names_list_2, **cmd_dict_3)

            # Apply transfrom(s) to 4D EP image (in parallel)
            if verbose:
                print("")
                print("Applying transform(s) to data in parallel")
            with multiprocessing.Pool(processes=num_jobs) as pool:
                proc = pool.starmap(apply_xfm_3D, cmd_tup_2)

            # Create list from output xfm-ed files, and sort, then merge
            if verbose:
                print("")
                print("Merging transformed data into 4D file")
            xfm_list = sorted(glob.glob(out_name_2 + "*.nii*"))
            out_file = merge_4D(nii_list=xfm_list, out_prefix=out_prefix, dim=dim, tr=tr, verbose=verbose)
        else:
            if verbose:
                print("")
                print("Applyting transform(s) to data - not in parallel - this may take some time...")
            apply_xfm_3D(nii_file=nii_file, out_prefix=out_prefix, **cmd_dict_3)

        # remove temporary directory and leftover files
        if verbose:
            print("")
            print("Cleaning up")
        shutil.rmtree(out_tmp)

        if verbose:
            print("")
            print("Done")

        return out_file

    # Apply transforms to 4D EP image
    if args['--in'] and args['--out'] and args['--ref']:
        # convert strings to int and/or float
        try:
            args['--TR'] = float(args['--TR'])
        except ValueError:
            pass
        try:
            args['--num-jobs'] = int(args['--num-jobs'])
        except ValueError:
            pass
        try:
            args['--super-level'] = int(args['--super-level'])
        except ValueError:
            pass
        if args['--no-parallel']:
            args['--num-jobs'] = 1
        # Run apply xfm 4D
        apply_xfm_4D(nii_file=args['--in'], out_prefix=args['--out'], ref_xfm=args['--ref'], ref_target=args['--ref-tar'], num_jobs=args['--num-jobs'], dim=args['--dim'],
                     tr=args['--TR'], warp=args['--warp'], warp_app=args['--warp-app'], premat=args['--premat'], postmat=args['--postmat'], resamp_vox=args['--resamp-vox'],
                     resamp_dim=args['--resamp-dim'], mask=args['--mask'], interp=args['--interp'], padding_size=args['--padding-size'], use_qform=args['--use-qform'], data_type=args['--data-type'],
                     super_sampling=args['--super-sampling'], super_level=args['--super-level'], verbose=args['--verbose'])
