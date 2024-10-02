# mrsupres function

The `mrsupres` function is a super-resolution algorithm designed to enhance the resolution of diffusion-weighted images (DWIs) while preserving the fidelity of the original measurements. It is compatible with the [MRtrix3](https://www.mrtrix.org/) toolbox.

<img width="518" alt="Screenshot 2024-10-02 at 1 09 19 PM" src="https://github.com/user-attachments/assets/e27438f9-3dde-4d1b-bf08-7f5cb3a1e8d0">

---

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Algorithm Details](#algorithm-details)
- [Testing Dataset](#testing-dataset)
- [References](#references)

---

## Installation

### Step 1: Add `mrsupres.cpp` to MRtrix3

Copy the `mrsupres.cpp` file into the `cmd` directory of your MRtrix3 installation:

```bash
cp mrsupres.cpp /path/to/mrtrix3/cmd/
```
### Step 2: Build MRtrix3
Navigate to your MRtrix3 directory and rebuild the toolbox to include the new command:
```bash
cd /path/to/mrtrix3
./build
```
---

## Usage
Run the mrsupres command with the following syntax:
```bash
/path/to/mrtrix3/bin/mrsupres input_dwi.nii output_dwi.nii -mask mask.nii -hires t1.nii -kernel 3 -iter_params 1,2,4,8,16 -force -nthread 8
```
- `input_dwi.nii`: Low-resolution diffusion-weighted image.
- `output_dwi.nii`: Super-resolved output image.
- `-mask mask.nii`: Mask image corresponding to the input DWI.
- `-hires t1.nii`: High-resolution T1-weighted image registered to the upsampled DWI.
- `-kernel 3`: Kernel size for weight computation (default is 3).
- `-iter_params 1,2,4,8,16`: Hyper-parameters for each iteration.
- `-force`: Overwrite output files if they exist.
- `-nthread 64`: Number of threads for parallel processing.

---

## Algorithm Details

The super-resolution algorithm enhances the resolution of DWIs by leveraging self-similarity and high-resolution anatomical information from a T1-weighted image.

**Preprocessing**

Upsampling: The low-resolution 2-mm DWIs are linearly interpolated to 1-mm resolution using the `mrgrid` function with nearest-neighbor interpolation, providing the initial input for the `mrsupres` function.

Registration: The high-resolution T1-weighted (T1w) image is registered to the upsampled DWIs using `bbregister` to ensure accurate spatial alignment. The aligned high-resolution T1w image is used as a reference.

**Iterative Refinement**

For each voxel in the DWI, weights are computed within a kernel of size 5×5×5. The algorithm iteratively refines the high-resolution estimate using the following steps:

1. Weight Computation: At each iteration t, weights are calculated based on the similarity between 3×3×3 patches of the DWI and the corresponding T1w image voxels. For a voxel at position (i,j,k), with (x,y,z) being one of the neighboring voxels in the 5×5×5 kernel, the specific weight is calculated using the following function:

$$
W_{[x, y, z]} = \frac{\sum_{patch} \left(DWI_{[\hat{x}, \hat{y}, \hat{z}]} - DWI_{[\hat{i}, \hat{j}, \hat{k}]} \right)^2}{3^3 \times h^2} \times \frac{\left( T1W_{[x, y, z]} - T1W_{[i, j, k]} \right)^2}{h^2}
$$

2. Intensity Update: The high-resolution DWI at the next iteration is updated using the computed weights:

$$
DWI_{[x, y, z]}^{t+1} = \frac{\sum_{kernel} W_{[i, j, k]}^{t} \times DWI_{[i, j, k]}^{t}}{\sum_{kernel} W_{[i, j, k]}^{t}}
$$

3. Data Consistency: Data consistency is maintained by enforcing a constraint that the down-sampled high-resolution estimate must be consistent with the original low-resolution data. This is achieved by:

i) Down-sampling the high-resolution estimate.
ii) Computing the residual difference with the original low-resolution image.
iii) Up-sampling the residual using nearest-neighbor interpolation.
iv) Updating the high-resolution estimate by adding the up-sampled residual.

By iteratively refining the high-resolution estimate and ensuring consistency with the original data, the algorithm effectively enhances the resolution of DWIs while preserving the fidelity of the original measurements.

**Parallel Computing**

The `mrsupres` function is designed to leverage parallel computing in `mrtrix3` framework, enhancing performance significantly. Key computational steps—including weight calculation, intensity update, and data consistency adjustment—are parallelized using multi-threading. This approach distributes the computational load across multiple CPU cores, allowing for efficient processing of large-scale datasets.

---

## Testing Dataset
A testing dataset is available for download:
[Testing Dataset on Google Drive](https://drive.google.com/drive/folders/14QLuT5WieGiJ1WWljUbR75PhSO_BPcdk)

---

## References
Tournier, J.-D., et al. (2019). MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. NeuroImage, 202, 116137.

Greve, D. N., & Fischl, B. (2009). Accurate and robust brain image alignment using boundary-based registration. NeuroImage, 48(1), 63–72.

Coupé, P., Manjón, J. V., Chamberland, M., Descoteaux, M., & Collins, D. L. (2013). Collaborative patch-based super-resolution for diffusion-weighted images. NeuroImage, 83, 245–261.

Manjón, J. V., Coupé, P., Buades, A., Collins, D. L., & Robles, M. (2010). MRI superresolution using self-similarity and image priors. International Journal of Biomedical Imaging, 2010, Article ID 425891.

Lee, H. H., Lin, Y. C., Lemberskiy, G., Ades-Aron, B., Baete, S., Boada, F. E., Fieremans, E., & Novikov, D. S. (n.d.). SUper-REsolution TRACTography (SURE-TRACT) pipeline using self-similarity between diffusional and anatomical images. www.cai2r.net









