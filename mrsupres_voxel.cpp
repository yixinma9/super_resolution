/* Copyright (c) 2008-2022 the MRtrix3 contributors.
 *
 *
 * For more details, see http://www.mrtrix.org/.
 */

#include <iostream>

#include "axes.h"
#include "command.h"
#include "image.h"
#include "progressbar.h"
#include "algo/threaded_loop.h"
#include <stdio.h>
#include <math.h>
#include <numeric>
#include <cmath>
#include <vector>
#include <limits.h>
#include <float.h>

#include <stdint.h>  // uint64_t

//#include "filter/resize.h"
//#include "filter/reslice.h"
//#include "interp/nearest.h"
//#include "interp/linear.h"
//#include "interp/cubic.h"
//#include "interp/sinc.h"
#include "progressbar.h"
#include "algo/copy.h"
#include "adapter/regrid.h"


#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


using namespace MR;
using namespace App;

void usage ()
{
  AUTHOR = "Yixin Ma (yma12@mgh.harvard.edu) & Hong-Hsi Lee (HLEE84@mgh.harvard.edu)";

  SYNOPSIS = "Self-similarity super resolution DWI";

  DESCRIPTION
    + "This application upsample the low-resolution DWI volume by a scale of 2 "
    "with the self-similarity info from additional high-resolution anatomical volume. (see reference below for details).";

  ARGUMENTS
    + Argument ("lowres", "the input low-res image.").type_image_in ()
    + Argument ("out", "the output image.").type_image_out ();
    //+ Argument ("out_DS", "the temporary output image.").type_image_out ();

  OPTIONS
  + Option ("mask", "Output ouside the mask will be set to 0.")
  +   Argument ("image").type_image_in ()
  + Option ("out_DS", "Temporary file to save output after downsampled to the original dimension as lowres.")
  +   Argument ("image").type_image_in ()
  + Option ("hires", "Reference hi-res either T1 or T2 weighted images.")
  +   Argument ("image").type_image_in ()
  + Option ("kernel", "Set the patch size of the super-resolution kernel filter. By default, kernal = 3, which means that we are using a patch size of (3*2+1)^3 = 7x7x7.")
  +   Argument ("kernel").type_integer (0, 128)
  + Option ("iter_params", "Determine how many iteration we want to use and the parameters at each iteration, e.g. : 2,4,8,16,32.")
  +   Argument ("values").type_sequence_float ()

  + DataType::options();

  REFERENCES
    + "Manjón JV, Coupé P, Buades A, Collins DL, Robles M. MRI superresolution using self-similarity and image priors. Int J Biomed Imaging."
    "2010;2010:425891. doi: 10.1155/2010/425891. Epub 2010 Dec 8. PMID: 21197094; PMCID: PMC3004412.";

}

class SuperResolution
{MEMALIGN (SuperResolution)  // assign memory

  public:
    SuperResolution (Image<float> lowres, Image<float> out, Image<float> out_temp, Image<float> hires, Image<bool>&mask, const int length, const int kernel, const int iter, const float param, const float mean_hires, const float mean_lowres)
    :mask (mask), out(out), out_temp(out_temp), lowres (lowres), hires (hires), kernel (kernel), length(length), iter(iter), param(param),  mean_lowres(mean_lowres), mean_hires(mean_hires)
    { }
    //template <typename ImageType>
    
    void operator () (Image<float> hires, Image<float> out)
    {
        assign_pos_of(hires).to(mask);
        if(mask.value() == 0){
            out.value() = 0;
        }else{
        load_image (lowres, out, hires, mask, kernel);
        float w12_sum = 0.0;
            float weighting = param*param;
        for (int vox = 0; vox < length; vox++){
            //if (X[vox] != -100){ // add if loop to tackle the mask beter, but ended up using jsut zero-filling
            w12[vox] = exp(-w1[vox]*weighting/mean_hires-w2[vox]/27/mean_lowres*weighting);
            //print("For voxel #" + str(vox) + ", weighting 1 from hi-res image is " + str(w1[vox]) + ", weighting 2 from low-res image is " + str(w2[vox]/27.0) + "\n");
            w12_sum += w12[vox];
            //}
        }
        //
        float vox_value = 0.0;
        for (int vox = 0; vox < length; vox++){
            w12[vox] = w12[vox]/w12_sum;
            vox_value = vox_value + w12[vox]*X[vox];
        }
        out.value() = vox_value;
        }
     //print("finish generating the output for voxel nums:" + str(out.index(0)) + "," + str(out.index(1)) + "," +  str(out.index(2)) + "\n");
    }
    
  private:
    Image<bool> mask;
    Image<float> out;
    Image<float> out_temp;
    Image<float> lowres;
    Image<float> hires;
    const int kernel;
    const int length;
    const int iter;
    const float param;
    float mean_lowres;
    float mean_hires;
    std::array<ssize_t, 3> pos;
    std::array<ssize_t, 3> pos1;
    std::array<ssize_t, 3> pos_inner;
    std::array<ssize_t, 3> pos_inner_center;
    //Following definition of dynamically sized array will cause some problem when parallel computing is turned on
    //float* w1 = new float[length];
    //float* w2 = new float[length];
    //float* X = new float[length];
    float w1[729];
    float w2[729];
    float w12[729];
    float X[729];
    float lowres_center_inner[27];
    float hires_center;
    //float lowres_center;
    
    //template <typename ImageType>
    void load_image (Image<float> lowres, Image<float> out, Image<float> hires, Image<bool> mask, int kernel) {
        //print("voxel nums are:" + str(out.index(0)) + "," + str(out.index(1)) + "," +  str(out.index(2)) + "\n");
        assign_pos_of(hires).to(out); // assign the index from out to hires
        hires_center = hires.value(); //get the high resolution center
        
        pos[0] = out.index(0); pos[1] = out.index(1); pos[2] = out.index(2); // pos: the position of the voxel for hires image at each Thread_Loop, ranging from [0,0,0] to [127,127,127];
        
        int q1 = 0;
        for (int zz = -1; zz <= 1; zz++) {
            pos_inner_center[2] = pos[2] + zz;
            for (int yy = -1; yy <= 1; yy++) {
                pos_inner_center[1] = pos[1] + yy;
                for (int xx = -1; xx <= 1; xx++) {
                    pos_inner_center[0] = pos[0] + xx;
                    mask.index(0) = pos_inner_center[0];
                    mask.index(1) = pos_inner_center[1];
                    mask.index(2) = pos_inner_center[2];
                    // pos_inner_center: the position of the voxel when moving around the 3x3x3 kernel, centering at pos, ranging from [-1,-1,-1] to [1,1,1];
                    //if (mask.value() == 0 or pos_inner_center[2] >= out.size(2) or pos_inner_center[2] < 0 or pos_inner_center[1] >= out.size(1) or pos_inner_center[1] < 0 or pos_inner_center[0] >= out.size(0) or pos_inner_center[0] < 0){
                    if (pos_inner_center[2] >= out.size(2) or pos_inner_center[2] < 0 or pos_inner_center[1] >= out.size(1) or pos_inner_center[1] < 0 or pos_inner_center[0] >= out.size(0) or pos_inner_center[0] < 0){
                        lowres_center_inner[q1] = 0;
                    }else{
                        if (iter == 0){
                            lowres.index(0) = floor(pos_inner_center[0]/2);
                            lowres.index(1) = floor(pos_inner_center[1]/2);
                            lowres.index(2) = floor(pos_inner_center[2]/2);
                            lowres_center_inner[q1] = lowres.value();
                        }else{
                            out_temp.index(0) = pos_inner_center[0];
                            out_temp.index(1) = pos_inner_center[1];
                            out_temp.index(2) = pos_inner_center[2];
                            lowres_center_inner[q1] = out_temp.value();
                        }
                    }
                    q1 += 1;
                }
            }
        }
        
        
        int k = 0;
        for (int z = -kernel; z <= kernel; z++) {
            pos1[2] = pos[2] + z;
            for (int y = -kernel; y <= kernel; y++) {
                pos1[1] = pos[1] + y;
                for (int x = -kernel; x <= kernel; x++) {
                    pos1[0] = pos[0] + x;
                    // pos1: the position of the voxel when moving around the 7x7x7 kernel, ranging from [-3,-3,-3] to [3,3,3];
                    mask.index(0) = pos1[0];
                    mask.index(1) = pos1[1];
                    mask.index(2) = pos1[2];
                    X[k] = -100;
                    // now the out.index, hires.index, and mask.index, they all carry the search grid in the 7x7x7 kernel
                    //if( pos1[2] >= out.size(2) or pos1[2] < 0 or pos1[1] >= out.size(1) or pos1[1] < 0 or pos1[0] >= out.size(0) or pos1[0] < 0 ){
                    if( pos1[2] >= out.size(2) or pos1[2] < 0 or pos1[1] >= out.size(1) or pos1[1] < 0 or pos1[0] >= out.size(0) or pos1[0] < 0 or mask.value() == 0){
                        w1[k] = DBL_MAX;
                        w2[k] = DBL_MAX;
                        X[k] = -100;
                    }else{
                    // start the inner loop 3x3x3 small kernel within the 7x7x7 kernel
                    w2[k] = 0;
                    int  q2 = 0;
                    for (int zz = -1; zz <= 1; zz++) {
                        pos_inner[2] = pos1[2] + zz;
                        for (int yy = -1; yy <= 1; yy++) {
                            pos_inner[1] = pos1[1] + yy;
                            for (int xx = -1; xx <= 1; xx++) {
                                pos_inner[0] = pos1[0] + xx;
                                // pos_inner: the position of the voxel when moving around the 3x3x3 kernel, centering at pos1, ranging from [-1,-1,-1] to [1,1,1];
                                mask.index(0) = pos_inner[0];
                                mask.index(1) = pos_inner[1];
                                mask.index(2) = pos_inner[2];
                                if ( pos_inner[2] >= out.size(2) or pos_inner[2] < 0 or pos_inner[1] >= out.size(1) or pos_inner[1] < 0 or pos_inner[0] >= out.size(0) or pos_inner[0] < 0)   {
                                //if (mask.value() == 0 or lowres_center_inner[q2] ==  0 or pos_inner[2] >= out.size(2) or pos_inner[2] < 0 or pos_inner[1] >= out.size(1) or pos_inner[1] < 0 or pos_inner[0] >= out.size(0) or pos_inner[0] < 0)   {
                                    w2[k] += lowres_center_inner[q2]*lowres_center_inner[q2];
                                } else {
                                    if (iter == 0){
                                        lowres.index(0) = floor(pos_inner[0]/2);
                                        lowres.index(1) = floor(pos_inner[1]/2);
                                        lowres.index(2) = floor(pos_inner[2]/2);
                                        if (q2 == 13){
                                            X[k] = lowres.value();
                                        }
                                        w2[k] += (lowres.value() - lowres_center_inner[q2])*(lowres.value() - lowres_center_inner[q2]);
                                    }else{
                                        out_temp.index(0) = pos_inner[0];
                                        out_temp.index(1) = pos_inner[1];
                                        out_temp.index(2) = pos_inner[2];
                                        if (q2 == 13){
                                            X[k] = out_temp.value();
                                        }
                                        w2[k] += (out_temp.value() - lowres_center_inner[q2])*(out_temp.value() - lowres_center_inner[q2]);
                                    }
                                    
                                }
                                //print(str(w2[k]) + "\n");
                                q2 += 1;
                            }
                        }
                    }
                    
                   // build the w1
                   hires.index(0) = pos1[0];
                   hires.index(1) = pos1[1];
                   hires.index(2) = pos1[2];
                   w1[k] = (hires.value() - hires_center)*(hires.value() - hires_center);
                   //print(str(hires_center) + "+" + str(lowres_center) + "\n");
                        }
                   k += 1;
                }
              }
            }
        }
};

void DownSample(Image<float>out, Image<float>lowres, Image<float>out_DS, Image<bool> mask){
for (auto l = Loop (lowres) (lowres); l; ++l) {
    float box_average = 0;
    float voxel_count = 0;
    for (int x = 0; x <=1; x++){
        for (int y = 0; y <=1; y++){
            for (int z = 0; z <=1; z++){
                out.index(0) = lowres.index(0)*2 + x;
                out.index(1) = lowres.index(1)*2 + y;
                out.index(2) = lowres.index(2)*2 + z;
                assign_pos_of(out).to(mask);
                //if (mask.value() != 0){
                box_average += out.value();
                voxel_count += 1;
                //}
            }
        }
    }
    assign_pos_of(lowres).to(out_DS);
    if (voxel_count == 0) {
        out_DS.value() = 0;
    }else{
        out_DS.value() = box_average/voxel_count - lowres.value();
    }
    }
}

void UpSample(Image<float>out, Image<float>lowres, Image<float>out_DS, Image<bool> mask, Image<float> out_temp){
for (auto l = Loop (out) (out); l; ++l) {
    assign_pos_of(out).to(mask);
    assign_pos_of(out).to(out_temp);
    //if (mask.value() != 0){
        out_DS.index(0) = floor(out.index(0)/2);
        out_DS.index(1) = floor(out.index(1)/2);
        out_DS.index(2) = floor(out.index(2)/2);
        out.value() = out.value() - out_DS.value();
        out_temp.value() = out.value();
    //    }
    }
}

 void CalMean(Image<float>hires, Image<float>lowres, Image<bool> mask, float& mean_hires, float& mean_lowres){
    float num_voxels = 0;
    float sum_lowres = 0;
    float sum_hires = 0;
    
    for (int z = 0; z < hires.size(2); z++) {
        for (int y = 0; y < hires.size(1); y++){
            for (int x = 0; x < hires.size(0); x++){
                hires.index(0) = x;
                hires.index(1) = y;
                hires.index(2) = z;
                
                mask.index(0) = x;
                mask.index(1) = y;
                mask.index(2) = z;
                
                lowres.index(0) = floor(x/2);
                lowres.index(1) = floor(y/2);
                lowres.index(2) = floor(z/2);
                
                if (mask.value() != 0){
                    num_voxels += 1;
                    sum_lowres += lowres.value()*lowres.value();
                    sum_hires += hires.value()*hires.value();
                }
            }
        }
    }
    mean_lowres = sum_lowres/num_voxels;
    mean_hires = sum_hires/num_voxels;
}

void run ()
{
    auto lowres = Image<float>::open (argument[0]);
    auto data2 = Header::open (argument[0]);
    Header header2 (data2);
    header2.datatype() = DataType::from_command_line (DataType::Float64);
    auto out_DS = Image<float>::create ("out_DS.nii", header2);
    
    
    Image<float> hires;
    auto opt = get_options ("hires");
    hires = Image<float>::open (opt[0][0]);
    auto data = Header::open (opt[0][0]);
    Header header (data);
    // The output will be same size as hires-T1 image
    header.datatype() = DataType::from_command_line (DataType::Float64);
    auto out = Image<float>::create (argument[1], header);
    auto out_temp = Image<float>::create ("out_temp.nii", header);
    
    Image<bool> mask;
    opt = get_options ("mask");
    mask = Image<bool>::open (opt[0][0]);
    
    opt = get_options("iter_params");
    auto iter_params = parse_floats(opt[0][0]);
    print("Number of iterations is: " + str(iter_params.size()) + "\n");
    print("weighting for each iteration is: " + str(iter_params) + "\n");

    
    const int kernel = App::get_option_value ("kernel", 3);
    const int length = std::pow(kernel*2+1, 3);
    
    print("selected patch size: " + str(kernel*2+1) + " x " + str(kernel*2+1) + " x " + str(kernel*2+1) + "\n");
    float mean_lowres = 0.0;
    float mean_hires = 0.0;
    CalMean(hires, lowres, mask, mean_hires, mean_lowres);
    print(" The mean value of the lowres image is: " + str(mean_lowres) + "\n");
    print(" The mean value of the hires image is: " + str(mean_hires) + "\n");
    for (int iter = 0; iter < iter_params.size(); iter++) {
        print("running " + str(iter) + " iteration, the parameter being used is " + str(iter_params[iter]) + "\n");
        ThreadedLoop ("running hi-res", out).run (SuperResolution (lowres, out, out_temp, hires, mask, length, kernel, iter, iter_params[iter], mean_hires, mean_lowres), hires, out);
        print("One iteration of SuperResolution finished \n");
	DownSample(out, lowres, out_DS, mask);
        print("Downsample finished \n");
        UpSample(out, lowres, out_DS, mask, out_temp);
	print("Upsample finished \n");
    }
}

