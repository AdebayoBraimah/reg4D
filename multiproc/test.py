#!/usr/bin/env python3

import nibabel as nib
import os
import sys
import subprocess
import glob
import random
import shutil
import multiprocessing

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
            filename = os.path.abspath(file_name)
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
    Creates a command and an empty command list for UNIX command line programs/applications

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

    # Get image dimensions
    [n_dim, x_dim, y_dim, z_dim, n_vol, x_vox, y_vox, z_vox, tr] = get_img_dims(nii_file=ref)

    if n_dim > 3:
        print(f"Image {nii_file} is greater than 3 dimensions. Exiting")
        sys.exit(1)

    # Add required input/output arguments
    if not nii_file in '.nii.gz':
        # print(True)
        nii_file = nii_file + '.nii.gz'
    rsm.append(nii_file)
    # out_prefix = out_prefix + ".nii.gz"
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

    # out_file = out_prefix + ".nii.gz"
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

    if not file in '.nii.gz':
        file = file + '.nii.gz'

    if not dest_name in '.nii.gz':
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

# Test functions here

if __name__ == "__main__":

    def apply_xfm_4D(nii_file, out_prefix, ref_xfm, ref_target, num_jobs="infer", dim="-tr",
                     tr="infer", warp="", warp_app="relative", premat="", postmat="", resamp_vox=True,
                     resamp_dim=False, mask="", interp="", padding_size="", use_qform=False, data_type="",
                     super_sampling=False, super_level=2, verbose=False):
        '''working doc-string'''

        # random number for random number generator
        n = 10000  # maximum N

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

        # Get TR
        if tr == 'infer':
            tr = tr_img

        # Make keyword argument dictionaries
        # Command 1 dictionary
        cmd_dict_1 = {"ref_vol": ref_xfm}  # Required argument(s) not in command list

        if warp:
            tmp_dict = {"warp": warp}
            cmd_dict_1.update(tmp_dict)
        else:
            tmp_dict = {"warp": ""}
            cmd_dict_1.update(tmp_dict)

        if warp_app:
            tmp_dict = {"warp_app": warp_app}
            cmd_dict_1.update(tmp_dict)

        if premat:
            tmp_dict = {"premat": premat}
            cmd_dict_1.update(tmp_dict)
        else:
            tmp_dict = {"premat": ""}
            cmd_dict_1.update(tmp_dict)

        if postmat:
            tmp_dict = {"postmat": postmat}
            cmd_dict_1.update(tmp_dict)
        else:
            tmp_dict = {"postmat": ""}
            cmd_dict_1.update(tmp_dict)

        if mask:
            tmp_dict = {"mask": mask}
            cmd_dict_1.update(tmp_dict)
        else:
            tmp_dict = {"mask": ""}
            cmd_dict_1.update(tmp_dict)

        if interp:
            tmp_dict = {"interp": interp}
            cmd_dict_1.update(tmp_dict)
        else:
            tmp_dict = {"interp": ""}
            cmd_dict_1.update(tmp_dict)

        if padding_size:
            tmp_dict = {"padding_size": padding_size}
            cmd_dict_1.update(tmp_dict)
        else:
            tmp_dict = {"padding_size": ""}
            cmd_dict_1.update(tmp_dict)

        if use_qform:
            tmp_dict = {"use_qform": True}
            cmd_dict_1.update(tmp_dict)
        else:
            tmp_dict = {"use_qform": ""}
            cmd_dict_1.update(tmp_dict)

        if data_type:
            tmp_dict = {"data_type": data_type}
            cmd_dict_1.update(tmp_dict)
        else:
            tmp_dict = {"data_type": ""}
            cmd_dict_1.update(tmp_dict)

        if super_sampling:
            tmp_dict = {"super_sampling": True}
            cmd_dict_1.update(tmp_dict)
        else:
            tmp_dict = {"super_sampling": ""}
            cmd_dict_1.update(tmp_dict)

        if super_level:
            tmp_dict = {"super_level": super_level}
            cmd_dict_1.update(tmp_dict)
        else:
            tmp_dict = {"super_level": 2}
            cmd_dict_1.update(tmp_dict)

        if verbose:
            tmp_dict = {"verbose": True}
            cmd_dict_1.update(tmp_dict)
        else:
            tmp_dict = {"verbose": ""}
            cmd_dict_1.update(tmp_dict)

        # Command 2 dictionary
        cmd_dict_2 = {"ref": ref_target}  # Required argument(s) not in command list

        if not resamp_vox and not resamp_dim:
            cp_dir = True
        else:
            cp_dir = False
            if resamp_vox:
                tmp_dict = {"resamp_vox": True}
                cmd_dict_2.update(tmp_dict)
            else:
                tmp_dict = {"resamp_vox": ""}
                cmd_dict_2.update(tmp_dict)

            if resamp_dim:
                tmp_dict = {"resamp_dim": True}
                cmd_dict_2.update(tmp_dict)
            else:
                tmp_dict = {"resamp_dim": ""}
                cmd_dict_2.update(tmp_dict)

            if verbose:
                tmp_dict = {"verbose": True}
                cmd_dict_2.update(tmp_dict)
            else:
                tmp_dict = {"verbose": ""}
                cmd_dict_2.update(tmp_dict)

        # Make temporay directories
        out_tmp = os.path.join(os.path.dirname(out_prefix), 'tmp_dir' + str(random.randint(0, n)))
        # out_tmp = os.path.join(os.path.dirname(nii_file), 'tmp_dir' + str(random.randint(0, n)))

        out_tmp_1 = os.path.join(os.path.abspath(out_tmp), 'tmp_dir' + str(random.randint(0, n)))
        out_tmp_2 = os.path.join(os.path.abspath(out_tmp), 'tmp_dir' + str(random.randint(0, n)))
        out_tmp_3 = os.path.join(os.path.abspath(out_tmp), 'tmp_dir' + str(random.randint(0, n)))

        if not os.path.exists(out_tmp_1):
            os.makedirs(out_tmp_1)

        if not os.path.exists(out_tmp_2):
            os.makedirs(out_tmp_2)

        if not os.path.exists(out_tmp_3):
            os.makedirs(out_tmp_3)

        out_tmp_1 = os.path.abspath(out_tmp_1)
        out_tmp_2 = os.path.abspath(out_tmp_2)
        out_tmp_3 = os.path.abspath(out_tmp_3)

        out_name_1 = os.path.join(out_tmp_1, 'func_split')
        out_name_2 = os.path.join(out_tmp_2, 'func_xfm')
        out_name_3 = os.path.join(out_tmp_3, 'func_rsm')

        if verbose:
            print("")
            print("Splitting 4D timeseries")

        pre_xfm_list = split_4D(nii_4d=nii_file, out_prefix=out_name_1, dim=split_dim, verbose=verbose)

        ###### parallel starts here ######

        out_tmp_names_list_2 = make_out_list(out_name=out_name_2, max_val=n_vol)
        out_tmp_names_list_3 = make_out_list(out_name=out_name_3, max_val=n_vol)

        cmd_tup_2 = multiproc_cmd_list_tuple(in_list=pre_xfm_list, out_list=out_tmp_names_list_2, **cmd_dict_1)
        cmd_tup_3 = multiproc_cmd_list_tuple(in_list=out_tmp_names_list_2, out_list=out_tmp_names_list_3, **cmd_dict_2)

        if verbose:
            print("")
            print("Applying transform(s) to data")

        with multiprocessing.Pool(processes=num_jobs) as pool:
            proc_1 = pool.starmap(apply_xfm_3D, cmd_tup_2)

        if not cp_dir:
            with multiprocessing.Pool(processes=num_jobs) as pool:
                proc_2 = pool.starmap(img_res, cmd_tup_3)
        else:
            cp_dir_imgs(in_list=out_tmp_names_list_2,out_list=out_tmp_names_list_3,max_val=n_vol)

        # print(cmd_tup_2[0])
        # print(cmd_dict_2)

        # print(out_tmp)
        # print(out_tmp_names_list_2)
        # print(out_tmp_names_list_3)

        # Do stuff here

        ###### parallel ends here ######

        if verbose:
            print("")
            print("Merging transformed data into 4D file")

        xfm_list = sorted(glob.glob(out_name_3 + "*.nii*"))

        out_file = merge_4D(nii_list=xfm_list, out_prefix=out_prefix, dim=dim, tr=tr, verbose=False)

        # remove temporary directory and leftover files
        shutil.rmtree(out_tmp)

        return out_file

    # Test variables
    nii_file="/Users/brac4g/Desktop/reg4D/test_data/filtered_func_data.nii.gz"
    out="/Users/brac4g/Desktop/reg4D/test_data/img_test"
    ref_xfm="/Users/brac4g/Desktop/reg4D/test_data/highres.nii.gz"
    ref_target="/Users/brac4g/Desktop/reg4D/test_data/example_func.nii.gz"
    premat="/Users/brac4g/Desktop/reg4D/test_data/example_func2highres.mat"
    num_jobs=20
    resamp_vox = True
    resamp_dim = False
    apply_xfm_4D(nii_file=nii_file, out_prefix=out, ref_xfm=ref_xfm, ref_target=ref_target, num_jobs=num_jobs, dim="-tr",
                 tr="infer", warp="", warp_app="relative", premat=premat, postmat="", resamp_vox=resamp_vox,
                 resamp_dim=resamp_dim, mask="", interp="", padding_size="", use_qform=False, data_type="",
                 super_sampling=False, super_level=2, verbose=True)
