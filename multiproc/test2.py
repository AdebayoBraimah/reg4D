#!/usr/bin/env python3

# """Naval Fate.
#
# Usage:
#   naval_fate.py ship new <name>...
#   naval_fate.py ship <name> move <x> <y> [--speed=<kn>]
#   naval_fate.py ship shoot <x> <y>
#   naval_fate.py mine (set|remove) <x> <y> [--moored | --drifting]
#   naval_fate.py (-h | --help)
#   naval_fate.py --version
#
# Options:
#   -h --help     Show this screen.
#   --version     Show version.
#   --speed=<kn>  Speed in knots [default: 10].
#   --moored      Moored (anchored) mine.
#   --drifting    Drifting mine.
#
# """

# """test2.
#
# Usage:
#   test2.py [required options] [options]
#
# Required arguments
#     -i,--in INPUT      Input
#     -o,--out OUTPUT    Output
#
# Options:
#     -h --help       Show this screen.
#     --version       Show version.
#
# """

# """Usage: my_program.py [-hso FILE] [--quiet | --verbose] [INPUT ...]
#
# -h --help    show this
# -s --sorted  sorted output
# -o FILE      specify output file [default: ./test.txt]
# --quiet      print less text
# --verbose    print more text
#
# """

"""

Applies pre-computed (affine) linear and non-linear transforms to 4D EP (bold) image data in parallel.
Additional options include the ability to resample the image back to its native voxel size and/or native
space dimensions (in the case a native space image is provided).

This script requires that both FSL and MIRTK be installed and added to system path.

Usage:
  reg4D [options] [required arguments]

Required arguments
    -i,--in INPUT       Input 4D file
    -o,--out PREFIX     Output prefix for transformed 4D file
    -r,--ref IMAGE      Reference image that corresponds to the transform(s)

Options:
    -rt,--ref-tar IMAGE Reference (target) image to resample to post-transform
                        (usually the native space image) - behaves best with
                        linear transfroms.
    -n,--num-jobs INT   Number of jobs to run in parallel [default: Max number
                        of cores available]
    -d,--dim CMD        (Command) Dimension to split and concatenate along the 4D image
                        (valid options: "-x", "-y", "-z", "t", "-tr") [default: "-tr"]
    -w,--warp IMAGE     FSL-style non-linear warp field file to reference image
    --warp-app CMD      (Command) Warp field treatment and application
                        (valid options: "relative", "absolute") [default: "relative"]
    --premat MAT        FSL-style pre-transform linear transformation matrix
    --postmat MAT       FSL-style post-transform linear transformation matrix
    --resamp-vox        Resample voxel-size to reference target image (can only be
                        enabled when the '--ref-tar' option is used)
    --resamp-dim        Resample image dimension to reference target image (can only be
                        enabled when the '--ref-tar' option is used)
    -m,--mask IMAGE     Binary mask image file in reference space to use for applying transform(s)
    --interp CMD        Interpolation method, options include: "nn","trilinear","sinc","spline" [default: None]
    --padding-size INT  Extrapolates outside original volume by n voxels
    --use-qform         Use s/qforms of ref_vol and nii_file images - NOTE: no other transorms can be applied with
                        this option
    --data-type CMD     Force output data type (valid options: "char" "short" "int" "float" "double")
    --super-sampling    Intermediary supersampling of output
    --super-level CMD   Level of intermediary supersampling, a for 'automatic' or integer level.
                        Only used when '--super-sampling' option is enabled. [default: 2]
    -v,--verbose        Enable verbose output [default: False]

    -h,--help           Prints help message, then exits.
    --version           Prints version, then exits.
"""


from docopt import docopt


if __name__ == '__main__':
    args = docopt(__doc__, help=True, version='test2 v0.0.1', options_first=False)
    print(args)

    if not args['--in'] or not args['--out'] or not args['--ref']:
        print("")
        print("Usage:   test2.py --in IMAGE --out PREFIX --ref IMAGE   |   -h,-help,--help")
        print("")


        # def apply_xfm_4D(nii_file, out_prefix, ref_xfm, ref_target, num_jobs="infer", dim="-tr",
        #                  tr="infer", warp="", warp_app="relative", premat="", postmat="", resamp_vox=True,
        #                  resamp_dim=False, mask="", interp="", padding_size="", use_qform=False, data_type="",
        #                  super_sampling=False, super_level=2, verbose=False):
        #     '''
        #     Applies pre-computed (affine) linear and non-linear transforms to 4D EP (bold) image data in parallel.
        #     Additional options include the ability to resample the image back to its native voxel size and/or native
        #     space dimensions (in the case a native space image is provided).
        #
        #     Arguments:
        #         nii_file (nifti file): Input 4D nifti file
        #         out_prefix (filename): Ouput prefix filename (and path, no '.nii.gz')
        #         ref_xfm (nifti file): Reference volume that correspond to the transform(s)
        #         ref_target (nifti file, optional): Target nifti file to be resampled to (usually native space)
        #         num_jobs (int, optional): Number of jobs to run simulataneously [default: Number of availble cores]
        #         dim (string, optional): Dimension to separate and merge along [default: '-tr']
        #         tr (float, optional): Repetition Time (TR, in sec.) [default: Inferred from header or 2.000]
        #         warp (nifti file, optional): Warp file to reference
        #         warp_app (nifti file, optional): Warp field treatment. Valid options are: "relative" (default), and "absolute"
        #         premat (FSL mat file, optional): Pre-transform linear transformation matrix
        #         postmat (FSL mat file, optional): Post-transform linear transformation matrix
        #         resamp_vox (bool, optional): Resample voxel size to reference image [default: True, if ref_target is provided]
        #         resamp_dim (bool, optional): Resample image dimension to reference image [default: False]
        #         mask (nifti file, optional): Mask, in reference space to use for applying transform(s)
        #         interp (string, optional): Interpolation method, options include: "nn","trilinear","sinc","spline"
        #         padding_size (int, optional): Extrapolates outside original volume by n voxels
        #         se_qform (bool, optional): Use s/qforms of ref_vol and nii_file images - NOTE: no other transorms can be applied with this option
        #         data_type (string,optional): Force output data type ["char" "short" "int" "float" "double"]
        #         super_sampling (bool, optional): Intermediary supersampling of output [default: False]
        #         super_level (int, optional): Level of intermediary supersampling, a for 'automatic' or integer level [default: 2]
        #         verbose (bool, optional): Enable verbose output [default: False]
        #     Returns:
        #         out_file (nifti file): Transformed 4D nifti file
        #     '''
