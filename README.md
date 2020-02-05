# reg4D
-----

Applies pre-computed (affine) linear and non-linear transforms to 4D EP (bold) image data in parallel.
Additional options include the ability to resample the image back to its native voxel size and/or native
space dimensions (in the case a native space image is provided).

This script at minimum requires that `FSL` be installed and added to system path.
Image resampling options required that `MIRTK` be installed and added to the system path.

```
Usage:
  reg4D [options | expert options] [required arguments]

Required arguments
    -i,--in INPUT       Input 4D file
    -o,--out PREFIX     Output prefix for transformed 4D file
    -r,--ref IMAGE      Reference image that corresponds to the transform(s)

Options:
    --ref-tar IMAGE Reference (target) image to resample to post-transform
                        (usually the native space image) - behaves best with
                        linear transfroms.
    -TR,--TR FLOAT      Repetition time (TR) of the input 4D EPI. If this value is not
                        specified, it will then be read from the nifti header [default: infer]
    -n,--num-jobs INT   Number of jobs to run in parallel. Default behavior is to use the
                        max number of cores available [default: infer]

Expert Options:
    -d,--dim CMD        (Command) Dimension to split and concatenate along the 4D image
                        (valid options: -x, -y, -z, t, -tr) [default: -tr]
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
    --interp CMD        Interpolation method, options include: "nn","trilinear","sinc","spline"
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

NOTE: '--resamp-dim' and '--resamp-vox' options require MIRTK to be installed. If these options are specified
```
