#!/bin/bash

set -euo pipefail
set -x

BEASTLIBRARY_DIR="${QUARANTINE_PATH}/resources/BEaST_libraries/combined"
RESAMPLEMODEL="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c.mnc"
RESAMPLEMASK="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc"

# REGISTRATIONMODEL="${QUARANTINE_PATH}/resources/adni_model_3d_v2/model_t1w.mnc"
# REGISTRATIONMODELMASK="${QUARANTINE_PATH}/resources/adni_model_3d_v2/model_t1w_mask.mnc"
# WMPRIOR="${QUARANTINE_PATH}/resources/adni_model_3d_v2/wm.mnc"
# GMPRIOR="${QUARANTINE_PATH}/resources/adni_model_3d_v2/gm.mnc"
# CSFPRIOR="${QUARANTINE_PATH}/resources/adni_model_3d_v2/csf.mnc"

REGISTRATIONMODEL="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c.mnc"
REGISTRATIONMODELMASK="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc"
WMPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_wm_tal_nlin_sym_09c.mnc"
GMPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_gm_tal_nlin_sym_09c.mnc"
CSFPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_csf_tal_nlin_sym_09c.mnc"
DEEPGMPRIOR=""

# 4 classes
#REGISTRATIONMODEL="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c.mnc"
#REGISTRATIONMODELMASK="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc"
#WMPRIOR="/home/gdevenyi/scratch/MNI-reclassify/recombine/final_prior3.mnc"
#GMPRIOR="/home/gdevenyi/scratch/MNI-reclassify/recombine/final_prior2.mnc"
#CSFPRIOR="/home/gdevenyi/scratch/MNI-reclassify/recombine/final_prior1.mnc"
#DEEPGMPRIOR="/home/gdevenyi/scratch/MNI-reclassify/recombine/final_prior4.mnc"

# 4 classes
# REGISTRATIONMODEL="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c.mnc"
# REGISTRATIONMODELMASK="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc"
# WMPRIOR="/gpfs/fs0/scratch/m/mchakrav/gdevenyi/MNI-reclassify/recombine/final_prior3.mnc"
# GMPRIOR="/gpfs/fs0/scratch/m/mchakrav/gdevenyi/MNI-reclassify/recombine/final_prior2.mnc"
# CSFPRIOR="/gpfs/fs0/scratch/m/mchakrav/gdevenyi/MNI-reclassify/recombine/final_prior1.mnc"
# DEEPGMPRIOR="/gpfs/fs0/scratch/m/mchakrav/gdevenyi/MNI-reclassify/recombine/final_prior4.mnc"

#Calculator for maths
calc () { awk "BEGIN{ print $* }" ;}

function make_qc() {
    #Generate a standardized view of the final correct brain in MNI space, with classification overlayed
    #Create animated version if img2webp is available
    mkdir -p ${tmpdir}/qc

    #Resample into MNI space for all the inputs
    antsApplyTransforms -d 3 ${MNI_XFM:+-t ${MNI_XFM}} -t ${tmpdir}/${n}/mni0_GenericAffine.xfm \
        -i ${tmpdir}/${n}/classify2.mnc -o ${tmpdir}/qc/classify.mnc -r ${RESAMPLEMODEL} -n GenericLabel
    antsApplyTransforms -d 3 ${MNI_XFM:+-t ${MNI_XFM}} -t ${tmpdir}/${n}/mni0_GenericAffine.xfm \
        -i $(dirname ${output})/$(basename ${output} .mnc).rescale.mnc -o ${tmpdir}/qc/corrected.mnc -r ${RESAMPLEMODEL} -n BSpline[5]
    antsApplyTransforms -d 3 ${MNI_XFM:+-t ${MNI_XFM}} -t ${tmpdir}/${n}/mni0_GenericAffine.xfm \
        -i ${tmpdir}/origqcref.mnc -o ${tmpdir}/qc/orig.mnc -r ${RESAMPLEMODEL} -n BSpline[5]
    mincmath -clobber -quiet ${N4_VERBOSE:+-verbose} -clamp -const2 0 65535 ${tmpdir}/qc/corrected.mnc ${tmpdir}/qc/corrected.clamp.mnc
    mv -f ${tmpdir}/qc/corrected.clamp.mnc ${tmpdir}/qc/corrected.mnc
    mincmath -clobber -quiet ${N4_VERBOSE:+-verbose} -clamp -const2 0 65535 ${tmpdir}/qc/orig.mnc ${tmpdir}/qc/orig.clamp.mnc
    mv -f ${tmpdir}/qc/orig.clamp.mnc ${tmpdir}/qc/orig.mnc

    #Create the bounding box for create_verify_image
    mincresample -clobber -quiet ${N4_VERBOSE:+-verbose} $(mincbbox -mincresample ${tmpdir}/qc/classify.mnc) ${tmpdir}/qc/classify.mnc ${tmpdir}/qc/label-crop.mnc
    minccalc -quiet ${N4_VERBOSE:+-verbose} -unsigned -byte -expression '1' ${tmpdir}/qc/label-crop.mnc ${tmpdir}/qc/bounding.mnc

    #Trasverse
    create_verify_image -range_floor 0 ${tmpdir}/qc/trans_classify.rgb \
        -width 1920 -autocols 10 -autocol_planes t \
        -bounding_volume ${tmpdir}/qc/bounding.mnc \
        -row ${tmpdir}/qc/corrected.mnc color:gray:0:65535 \
        volume_overlay:${tmpdir}/qc/classify.mnc:0.4

    create_verify_image -range_floor 0 ${tmpdir}/qc/trans_corrected.rgb \
        -width 1920 -autocols 10 -autocol_planes t \
        -bounding_volume ${tmpdir}/qc/bounding.mnc \
        -row ${tmpdir}/qc/corrected.mnc color:spect:0:65535

    create_verify_image -range_floor 0 ${tmpdir}/qc/trans_corrected_gray.rgb \
        -width 1920 -autocols 10 -autocol_planes t \
        -bounding_volume ${tmpdir}/qc/bounding.mnc \
        -row ${tmpdir}/qc/corrected.mnc color:gray:0:65535

    create_verify_image -range_floor 0 ${tmpdir}/qc/trans_orig.rgb \
        -width 1920 -autocols 10 -autocol_planes t \
        -bounding_volume ${tmpdir}/qc/bounding.mnc \
        -row ${tmpdir}/qc/orig.mnc color:spect:0:65535

    #Sagital
    create_verify_image -range_floor 0 ${tmpdir}/qc/sag_classify.rgb \
        -width 1920 -autocols 10 -autocol_planes s \
        -bounding_volume ${tmpdir}/qc/bounding.mnc \
        -row ${tmpdir}/qc/corrected.mnc color:gray:0:65535 \
        volume_overlay:${tmpdir}/qc/classify.mnc:0.4

    create_verify_image -range_floor 0 ${tmpdir}/qc/sag_corrected.rgb \
        -width 1920 -autocols 10 -autocol_planes s \
        -bounding_volume ${tmpdir}/qc/bounding.mnc \
        -row ${tmpdir}/qc/corrected.mnc color:spect:0:65535

    create_verify_image -range_floor 0 ${tmpdir}/qc/sag_corrected_gray.rgb \
        -width 1920 -autocols 10 -autocol_planes s \
        -bounding_volume ${tmpdir}/qc/bounding.mnc \
        -row ${tmpdir}/qc/corrected.mnc color:gray:0:65535

    create_verify_image -range_floor 0 ${tmpdir}/qc/sag_orig.rgb \
        -width 1920 -autocols 10 -autocol_planes s \
        -bounding_volume ${tmpdir}/qc/bounding.mnc \
        -row ${tmpdir}/qc/orig.mnc color:spect:0:65535

    #Coronal
    create_verify_image -range_floor 0 ${tmpdir}/qc/cor_classify.rgb \
        -width 1920 -autocols 10 -autocol_planes c \
        -bounding_volume ${tmpdir}/qc/bounding.mnc \
        -row ${tmpdir}/qc/corrected.mnc color:gray:0:65535 \
        volume_overlay:${tmpdir}/qc/classify.mnc:0.4

    create_verify_image -range_floor 0 ${tmpdir}/qc/cor_corrected.rgb \
        -width 1920 -autocols 10 -autocol_planes c \
        -bounding_volume ${tmpdir}/qc/bounding.mnc \
        -row ${tmpdir}/qc/corrected.mnc color:spect:0:65535

    create_verify_image -range_floor 0 ${tmpdir}/qc/cor_corrected_gray.rgb \
        -width 1920 -autocols 10 -autocol_planes c \
        -bounding_volume ${tmpdir}/qc/bounding.mnc \
        -row ${tmpdir}/qc/corrected.mnc color:gray:0:65535

    create_verify_image -range_floor 0 ${tmpdir}/qc/cor_orig.rgb \
        -width 1920 -autocols 10 -autocol_planes c \
        -bounding_volume ${tmpdir}/qc/bounding.mnc \
        -row ${tmpdir}/qc/orig.mnc color:spect:0:65535

    convert -background black -strip -append \
        ${tmpdir}/qc/cor_corrected.rgb \
        ${tmpdir}/qc/cor_classify.rgb \
        ${tmpdir}/qc/sag_corrected.rgb \
        ${tmpdir}/qc/sag_classify.rgb \
        ${tmpdir}/qc/trans_corrected.rgb \
        ${tmpdir}/qc/trans_classify.rgb \
        ${tmpdir}/qc/corrected.mpc

    convert -background black -strip -append \
        ${tmpdir}/qc/cor_orig.rgb \
        ${tmpdir}/qc/cor_corrected_gray.rgb \
        ${tmpdir}/qc/sag_orig.rgb \
        ${tmpdir}/qc/sag_corrected_gray.rgb \
        ${tmpdir}/qc/trans_orig.rgb \
        ${tmpdir}/qc/trans_corrected_gray.rgb \
        ${tmpdir}/qc/orig.mpc

    #Save static QC jpg
    convert -background black -strip -interlace Plane -sampling-factor 4:2:0 -quality "85%" \
        ${tmpdir}/qc/corrected.mpc $(dirname ${output})/$(basename ${output} .mnc).jpg

    #If webp software is available animate a before/after image
    if command -v img2webp; then
        convert -background black ${tmpdir}/qc/corrected.mpc ${tmpdir}/qc/corrected.png
        convert -background black ${tmpdir}/qc/orig.mpc ${tmpdir}/qc/orig.png
        img2webp -d 750 -lossy -min_size ${tmpdir}/qc/orig.png ${tmpdir}/qc/corrected.png -o $(dirname ${output})/$(basename ${output} .mnc).webp || true
    fi
}

isotropize() {
  # Need smoothing for downsampling to avoid aliasing
  # Ideas stolen from https://discourse.itk.org/t/resampling-to-isotropic-signal-processing-theory/1403
  inputres=$(python -c "print('\n'.join([str(abs(x)) for x in [float(x) for x in \"$(PrintHeader ${1} 1)\".split(\"x\")]]))")
  blurs=""

  for dim in ${inputres}; do
      if [[ $(python -c "print(${dim}>(${isostep}-1e-6))") == True ]]; then
         #Special casing for zero/negative blurs
          blurs+=1e-12x
      else
          blurs+=$(python -c "import math; print(math.sqrt((${isostep}**2.0 - ${dim}**2.0)/(2.0*math.sqrt(2.0*math.log(2.0)))**2.0))")x
      fi
  done

  mincmath -mult ${tmpdir}/vessels.mnc ${1} ${tmpdir}/${n}/presmooth_novessels.mnc

  SmoothImage 3 ${tmpdir}/${n}/presmooth_novessels.mnc "${blurs%?}" ${tmpdir}/${n}/smoothed.mnc 1 0
  ResampleImage 3 ${tmpdir}/${n}/smoothed.mnc ${tmpdir}/${n}/isotropized.mnc ${isostep}x${isostep}x${isostep} 0 4

  mincmath -quiet -clamp -const2 0 65535 ${tmpdir}/${n}/isotropized.mnc ${tmpdir}/${n}/downsample.mnc

}

#Iterative multi-scale N3 implementation
do_N3() {
    # This code does cycles of N3 at a given scale, and then halves the scale and repeats until the limit is reached
    # Input images are downsampled before N3
    isotropize ${n3input}
    n3input=${tmpdir}/${n}/downsample.mnc
    minccalc -unsigned -byte -expression 'A[0]>1?1:0' ${tmpdir}/${n}/downsample.mnc ${tmpdir}/${n}/nonzero.mnc
    antsApplyTransforms -d 3 -i ${tmpdir}/${n}/weight.mnc -o ${tmpdir}/${n}/tmpweight.mnc -r ${n3input} -n GenericLabel --verbose
    antsApplyTransforms -d 3 -i ${tmpdir}/bgmask.mnc -o ${tmpdir}/${n}/tmpbg.mnc -r ${n3input} -n GenericLabel --verbose
    ImageMath 3 ${tmpdir}/${n}/tmpweight.mnc m ${tmpdir}/${n}/tmpweight.mnc ${tmpdir}/${n}/nonzero.mnc
    #ImageMath 3 ${tmpdir}/${n}/tmpweight.mnc GetLargestComponent ${tmpdir}/${n}/tmpweight.mnc
    distance=${origdistance}
    j=0
    while (( j < levels )); do
        i=0
        while (( i < cycles )); do
            nu_correct -clobber -normalize_field \
                -stop ${stop} -distance ${distance} -iterations ${iters} -fwhm ${fwhm} -shrink 1 -lambda ${lambda} \
                -mask ${tmpdir}/${n}/tmpweight.mnc ${n3input} ${tmpdir}/${n}/corrected_${distance}_${i}.mnc

            evaluate_field -unsigned -double -clobber -like ${n3input} ${tmpdir}/${n}/corrected_${distance}_${i}.imp ${tmpdir}/${n}/corrected_${distance}_${i}_field.mnc

            biasmean=$(mincstats -mean -quiet -mask ${tmpdir}/${n}/tmpweight.mnc -mask_binvalue 1 ${tmpdir}/${n}/corrected_${distance}_${i}_field.mnc)
            origmean=$(mincstats -mean -quiet -mask ${tmpdir}/${n}/tmpweight.mnc -mask_binvalue 1 ${n3input})

            minccalc -clobber -unsigned -short \
              -expression "clamp((A[0]/${origmean})/(A[1]/${biasmean})*32767,0,65535)" \
              ${n3input} ${tmpdir}/${n}/corrected_${distance}_${i}_field.mnc ${tmpdir}/${n}/corrected_${distance}_${i}.mnc -clobber

            n3input=${tmpdir}/${n}/corrected_${distance}_${i}.mnc

            ((++i))
        done
        distance=$(calc "${distance} / 2")
        ((++j))
    done

  for file in ${tmpdir}/${n}/*imp; do
     echo evaluate_field -unsigned -double -clobber -like ${tmpdir}/originput.mnc ${file} ${tmpdir}/${n}/$(basename $file .imp)_field.mnc
  done | parallel
}

do_N3_reordered() {
    # This code does repeated cycles of the full step-down of levels of N3
    isotropize ${n3input}
    n3input=${tmpdir}/${n}/downsample.mnc
    minccalc -unsigned -byte -expression 'A[0]>1?1:0' ${tmpdir}/${n}/downsample.mnc ${tmpdir}/${n}/nonzero.mnc
    antsApplyTransforms -d 3 -i ${tmpdir}/${n}/weight.mnc -o ${tmpdir}/${n}/tmpweight.mnc -r ${n3input} -n GenericLabel --verbose
    antsApplyTransforms -d 3 -i ${tmpdir}/bgmask.mnc -o ${tmpdir}/${n}/tmpbg.mnc -r ${n3input} -n GenericLabel --verbose
    ImageMath 3 ${tmpdir}/${n}/tmpweight.mnc m ${tmpdir}/${n}/tmpweight.mnc ${tmpdir}/${n}/nonzero.mnc
    ImageMath 3 ${tmpdir}/${n}/tmpweight.mnc GetLargestComponent ${tmpdir}/${n}/tmpweight.mnc

    i=0
    while (( i < cycles )); do
        distance=${origdistance}
        j=0
        while (( j < levels )); do
            nu_correct -clobber -normalize_field \
                -stop ${stop} -distance ${distance} -iterations ${iters} -fwhm ${fwhm} -shrink 1 -lambda ${lambda} \
                -mask ${tmpdir}/${n}/tmpweight.mnc ${n3input} ${tmpdir}/${n}/corrected_${distance}_${i}.mnc

            nu_evaluate -verbose -clobber -mapping ${tmpdir}/${n}/corrected_${distance}_${i}.imp \
              -mask ${tmpdir}/${n}/tmpbg.mnc \
              ${n3input} ${tmpdir}/${n}/corrected_${distance}_${i}.mnc

            mincmath -clobber -clamp -const2 0 65535 ${tmpdir}/${n}/corrected_${distance}_${i}.mnc ${tmpdir}/${n}/clamp.mnc
            mv -f ${tmpdir}/${n}/clamp.mnc ${tmpdir}/${n}/corrected_${distance}_${i}.mnc
            n3input=${tmpdir}/${n}/corrected_${distance}_${i}.mnc

            ((++j))
            distance=$(calc "${distance} / 2")
        done

        ((++i))
    done

  for file in ${tmpdir}/${n}/*imp; do
     echo evaluate_field -double -clobber -like ${tmpdir}/originput.mnc ${file} ${tmpdir}/${n}/$(basename $file .imp)_field.mnc
  done | parallel
}

do_N4() {

  i=0
  while (( i < cycles )); do
    N4BiasFieldCorrection -d 3 --verbose  -i ${n3input} \
    -b [ ${distance} ] -s $(calc "int(${shrink} / (4.0/${isostep}))") \
    --histogram-sharpening [ 0.1,0.01,200 ] \
    -c [ $(printf "${iters}")$(printf "x${iters}%.0s" $(eval echo {1..$(calc "${levels} - 1")})),1e-5 ] \
    -x ${tmpdir}/bgmask.mnc -w ${tmpdir}/${n}/weight.mnc \
    -o [ ${tmpdir}/${n}/corrected_${i}.mnc, ${tmpdir}/${n}/corrected_${i}_field.mnc ]
    n3input=${tmpdir}/${n}/corrected_${i}.mnc
    ((++i))
  done

  for file in ${tmpdir}/${n}/*field.mnc; do
    ImageMath 3 ${file} / ${file} $(mincstats -quiet -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 -mean ${file} )
  done

}

make_outlier_map() {
  outlier_input=$1
  mask_input=$2
  outlier_output=$3

  median=$(mincstats -mask ${mask_input} -mask_binvalue 1 -median -quiet ${outlier_input})

  minccalc -clobber \
    -expression "abs(A[0]-${median})" \
    ${outlier_input} ${tmpdir}/${n}/madmap.mnc

  mad=$(mincstats -mask ${mask_input} -mask_binvalue 1 -median -quiet ${tmpdir}/${n}/madmap.mnc)

  minccalc -expression "((0.6745*(A[0] - ${median}))/${mad})>3.5?0:1" \
    ${outlier_input} ${outlier_output}
}


tmpdir=$(mktemp -p . -d $(basename $1)XXXXX)


function finish {
  if [[ ! -s ${tmpdir}/keep ]]; then
    #rm -rf "${tmpdir}"
    echo done
  fi
}
trap finish EXIT

#Add handler for failure to show where things went wrong
failure() {
    local lineno=$1
    local msg=$2
    echo "Failed at $lineno: $msg"
}
trap 'failure ${LINENO} "$BASH_COMMAND"' ERR

input=$1
output=$2

#Defaults
origdistance=400
distance=${origdistance}
levels=5
cycles=3
iters=100
lambda=2e-6
shrink=4.0
fwhm=0.1
stop=1e-5
isostep=4.0

# Calulate a scaling factor for mm
dx=$(mincinfo -attvalue xspace:step ${input})
dy=$(mincinfo -attvalue yspace:step ${input})
dz=$(mincinfo -attvalue zspace:step ${input})
shrink=$(python -c "print(${shrink} / ( ( abs(${dx}) + abs(${dy}) + abs(${dz}) ) / 3.0 ))")

# Forceably convert to MINC2, and clamp range to avoid negative numbers, rescale to 0-65535
mincconvert -2 ${input} ${tmpdir}/originput.mnc
#Rescale initial data into entirely positive range (fix for completely negative data)
ImageMath 3 ${tmpdir}/originput.mnc RescaleImage ${tmpdir}/originput.mnc 0 65535
cp -f ${tmpdir}/originput.mnc ${tmpdir}/origqcref.mnc

mincnorm -short -clamp -out_floor 0 -out_ceil 65535 ${tmpdir}/originput.mnc ${tmpdir}/originput.clamp.mnc
mv -f ${tmpdir}/originput.clamp.mnc ${tmpdir}/originput.mnc

# Pad image for processing
ImageMath 3 ${tmpdir}/originput.mnc PadImage ${tmpdir}/originput.mnc 20

n=0
mkdir -p ${tmpdir}/${n}

minc_anlm --clobber --mt $(nproc) ${tmpdir}/originput.mnc ${tmpdir}/${n}/denoise.mnc

#Masking of blood vessels
itk_vesselness --clobber --scales 8 --rescale ${tmpdir}/${n}/denoise.mnc ${tmpdir}/vessels.mnc
ThresholdImage 3 ${tmpdir}/vessels.mnc ${tmpdir}/vessels.mnc 35 Inf 0 1

#Round 1, Otsu mask of foreground
ThresholdImage 3 ${tmpdir}/originput.mnc ${tmpdir}/${n}/bgmask.mnc 1 Inf 1 0
cp -f ${tmpdir}/${n}/bgmask.mnc ${tmpdir}/bgmask.mnc
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/bgmask.mnc ${tmpdir}/vessels.mnc
ThresholdImage 3 ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/weight.mnc Otsu 4 ${tmpdir}/${n}/weight.mnc
ThresholdImage 3 ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/weight.mnc 2 Inf 1 0

iMath 3 ${tmpdir}/${n}/fgmask.mnc MC ${tmpdir}/${n}/weight.mnc 10 1 ball 1
ImageMath 3 ${tmpdir}/${n}/fgmask.mnc FillHoles ${tmpdir}/${n}/fgmask.mnc 2
ImageMath 3 ${tmpdir}/${n}/fgmask.mnc GetLargestComponent ${tmpdir}/${n}/fgmask.mnc

ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/vessels.mnc
minccalc -unsigned -byte -expression 'A[0]>1?1:0' ${tmpdir}/originput.mnc ${tmpdir}/${n}/nonzero.mnc
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/nonzero.mnc
ImageMath 3 ${tmpdir}/${n}/weight.mnc GetLargestComponent ${tmpdir}/${n}/weight.mnc

make_outlier_map ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/fgmask.mnc ${tmpdir}/${n}/hotmask.mnc
#ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/hotmask.mnc

n3input=${tmpdir}/originput.mnc


#Initial correction not used for the final outputs but rather just to enable good foreground-background masking
N4BiasFieldCorrection -d 3 --verbose  -i ${n3input} \
  -b [ ${distance} ] -s $(awk -v flt=${shrink} 'BEGIN { printf("%.0f", flt); }') \
  --histogram-sharpening [ 0.1,0.01,200 ] \
  -c [ 300x300x300x300x300,1e-5 ] \
  -x ${tmpdir}/bgmask.mnc -w ${tmpdir}/${n}/weight.mnc \
  -o [ ${tmpdir}/${n}/correct.mnc,${tmpdir}/${n}/bias.mnc ]

biasmean=$(mincstats -mean -quiet -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/${n}/bias.mnc)
origmean=$(mincstats -mean -quiet -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/originput.mnc)

minccalc -unsigned -short -expression "clamp((A[0]/${origmean})/(A[1]/${biasmean})*32767,0,65535)" ${tmpdir}/originput.mnc ${tmpdir}/${n}/bias.mnc ${tmpdir}/${n}/correct.mnc -clobber


((++n))
mkdir -p ${tmpdir}/${n}

#Redo normalization
antsApplyTransforms -d 3 -i ${tmpdir}/$(( n - 1 ))/fgmask.mnc -r ${tmpdir}/origqcref.mnc -o ${tmpdir}/fgmask_orig.mnc -n GenericLabel
mincnorm -short -clamp -mask ${tmpdir}/fgmask_orig.mnc -out_floor 0 -out_ceil 65535 ${tmpdir}/origqcref.mnc ${tmpdir}/originput.clamp.mnc
mv -f ${tmpdir}/originput.clamp.mnc ${tmpdir}/originput.mnc
ImageMath 3 ${tmpdir}/originput.mnc PadImage ${tmpdir}/originput.mnc 20

minc_anlm --clobber --mt $(nproc) ${tmpdir}/$(( n - 1 ))/correct.mnc ${tmpdir}/${n}/denoise.mnc

#Masking of blood vessels
itk_vesselness --clobber --scales 8 --rescale ${tmpdir}/${n}/denoise.mnc ${tmpdir}/vessels.mnc
ThresholdImage 3 ${tmpdir}/vessels.mnc ${tmpdir}/vessels.mnc 35 Inf 0 1

#Redo the Otsu mask a second time using the precorrected image
ThresholdImage 3 ${tmpdir}/$(( n - 1 ))/correct.mnc ${tmpdir}/${n}/bgmask.mnc 1 Inf 1 0
cp -f ${tmpdir}/${n}/bgmask.mnc ${tmpdir}/bgmask.mnc
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/bgmask.mnc ${tmpdir}/vessels.mnc
ThresholdImage 3 ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/weight.mnc Otsu 4 ${tmpdir}/${n}/weight.mnc
ThresholdImage 3 ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/weight.mnc 2 Inf 1 0
ImageMath 3 ${tmpdir}/${n}/weight.mnc GetLargestComponent ${tmpdir}/${n}/weight.mnc
iMath 3 ${tmpdir}/${n}/weight.mnc MC ${tmpdir}/${n}/weight.mnc 10 1 ball 1
ImageMath 3 ${tmpdir}/${n}/weight.mnc FillHoles ${tmpdir}/${n}/weight.mnc 2
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/vessels.mnc
ThresholdImage 3 ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/weight.mnc Otsu 4 ${tmpdir}/${n}/weight.mnc

ThresholdImage 3 ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/fgmask.mnc 2 Inf 1 0
ImageMath 3 ${tmpdir}/${n}/fgmask.mnc GetLargestComponent ${tmpdir}/${n}/fgmask.mnc
iMath 3 ${tmpdir}/${n}/fgmask.mnc MC ${tmpdir}/${n}/fgmask.mnc 10 1 ball 1
ImageMath 3 ${tmpdir}/${n}/fgmask.mnc FillHoles ${tmpdir}/${n}/fgmask.mnc 2

make_outlier_map ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/fgmask.mnc ${tmpdir}/${n}/hotmask.mnc

ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/fgmask.mnc ${tmpdir}/${n}/hotmask.mnc
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/vessels.mnc

ThresholdImage 3 ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/weight.mnc Otsu 4 ${tmpdir}/${n}/weight.mnc

ThresholdImage 3 ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/weight.mnc 3 Inf 1 0
ImageMath 3 ${tmpdir}/${n}/weight.mnc GetLargestComponent ${tmpdir}/${n}/weight.mnc

n3input=${tmpdir}/originput.mnc
do_N3

mincmath -clobber -mult ${tmpdir}/${n}/*field.mnc ${tmpdir}/${n}/field_combined.mnc

correct_field ${tmpdir}/${n}/field_combined.mnc ${tmpdir}/${n}/fgmask.mnc ${tmpdir}/${n}/field_combined_correct.mnc
mincmath -clobber -clamp -const2 0.1 1.79769e+308 ${tmpdir}/${n}/field_combined_correct.mnc ${tmpdir}/${n}/field_combined_correct_clamp.mnc

origmean=$(mincstats -mean -quiet -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/originput.mnc)
biasmean=$(mincstats -mean -quiet -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/${n}/field_combined_correct_clamp.mnc)

minccalc -unsigned -short -expression "clamp((A[0]/${origmean})/(A[1]/${biasmean})*32767,0,65535)" ${tmpdir}/originput.mnc ${tmpdir}/${n}/field_combined_correct_clamp.mnc ${tmpdir}/${n}/correct.mnc -clobber

((++n))
mkdir -p ${tmpdir}/${n}

minc_anlm --clobber --mt $(nproc) ${tmpdir}/$(( n - 1 ))/correct.mnc ${tmpdir}/${n}/denoise.mnc

antsRegistration_affine_SyN.sh --verbose --float --convergence 1e-7 \
    --skip-nonlinear --fixed-mask ${REGISTRATIONMODELMASK} \
    ${tmpdir}/${n}/denoise.mnc ${REGISTRATIONMODEL} ${tmpdir}/${n}/mni

minccalc -unsigned -byte -expression '1' ${RESAMPLEMODEL} ${tmpdir}/model_fov.mnc
mincresample -nofill -clobber -labels -near -like ${tmpdir}/originput.mnc -transform ${tmpdir}/${n}/mni0_GenericAffine.xfm ${tmpdir}/model_fov.mnc ${tmpdir}/subject_fov.mnc
ImageMath 3 ${tmpdir}/fgmask_fov.mnc m ${tmpdir}/$(( n - 1 ))/fgmask.mnc ${tmpdir}/subject_fov.mnc
cp -f ${tmpdir}/$(( n - 1 ))/fgmask.mnc ${tmpdir}/fgmask.mnc

antsApplyTransforms -d 3 -i ${REGISTRATIONMODELMASK} \
    -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] \
    -n GenericLabel --verbose -r ${tmpdir}/originput.mnc \
    -o ${tmpdir}/${n}/mnimask.mnc

iMath 3 ${tmpdir}/${n}/mnimask.mnc MD ${tmpdir}/${n}/mnimask.mnc 1 1 ball 1
ImageMath 3 ${tmpdir}/${n}/mnimask.mnc m ${tmpdir}/${n}/mnimask.mnc ${tmpdir}/vessels.mnc
ThresholdImage 3 ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/weight.mnc Otsu 4 ${tmpdir}/${n}/mnimask.mnc
ThresholdImage 3 ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/weight.mnc 2 Inf 1 0
iMath 3 ${tmpdir}/${n}/weight.mnc ME ${tmpdir}/${n}/weight.mnc 1 1 ball 1
ImageMath 3 ${tmpdir}/${n}/weight.mnc GetLargestComponent ${tmpdir}/${n}/weight.mnc
iMath 3 ${tmpdir}/${n}/weight.mnc MD ${tmpdir}/${n}/weight.mnc 1 1 ball 1
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/vessels.mnc

iMath 3 ${tmpdir}/${n}/mnimask.mnc MC ${tmpdir}/${n}/weight.mnc 5 1 ball 1
ImageMath 3 ${tmpdir}/${n}/mnimask.mnc FillHoles ${tmpdir}/${n}/mnimask.mnc 2

#New test
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/mnimask.mnc ${tmpdir}/vessels.mnc
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/$(( n - 1 ))/hotmask.mnc
ThresholdImage 3 ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/weight.mnc Kmeans 2 ${tmpdir}/${n}/weight.mnc
ThresholdImage 3 ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/weight.mnc 2 Inf 1 0

make_outlier_map ${tmpdir}/${n}/denoise.mnc ${tmpdir}/fgmask.mnc ${tmpdir}/${n}/hotmask.mnc
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/hotmask.mnc
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/vessels.mnc

n3input=${tmpdir}/$(( n - 1 ))/correct.mnc
do_N3

iMath 3 ${tmpdir}/${n}/correct_mask.mnc MD ${tmpdir}/${n}/mnimask.mnc 2 1 ball 1
ImageMath 3 ${tmpdir}/${n}/correct_mask.mnc FillHoles ${tmpdir}/${n}/correct_mask.mnc

mincmath -clobber -unsigned -double -mult ${tmpdir}/${n}/*field.mnc ${tmpdir}/${n}/field_combined.mnc

correct_field ${tmpdir}/${n}/field_combined.mnc ${tmpdir}/${n}/correct_mask.mnc ${tmpdir}/${n}/field_combined_correct.mnc
mincmath -clobber -clamp -const2 0.1 1.79769e+308 ${tmpdir}/${n}/field_combined_correct.mnc ${tmpdir}/${n}/field_combined_correct_clamp.mnc

mincmath -clobber -unsigned -double -mult ${tmpdir}/$(( n - 1 ))/field_combined_correct_clamp.mnc ${tmpdir}/${n}/field_combined_correct_clamp.mnc ${tmpdir}/${n}/field_combined_correct_clamp2.mnc
mv -f ${tmpdir}/${n}/field_combined_correct_clamp2.mnc ${tmpdir}/${n}/field_combined_correct_clamp.mnc

origmean=$(mincstats -mean -quiet -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/originput.mnc)
biasmean=$(mincstats -mean -quiet -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/${n}/field_combined_correct_clamp.mnc)

minccalc -unsigned -short -expression "clamp((A[0]/${origmean})/(A[1]/${biasmean})*32767,0,65535)" ${tmpdir}/originput.mnc ${tmpdir}/${n}/field_combined_correct_clamp.mnc ${tmpdir}/${n}/correct.mnc -clobber

((++n))
mkdir -p ${tmpdir}/${n}

minc_anlm --clobber --mt $(nproc) ${tmpdir}/$(( n - 1 ))/correct.mnc ${tmpdir}/${n}/denoise.mnc

 antsRegistration_affine_SyN.sh --clobber --verbose --close --convergence 1e-7 \
         --initial-transform ${tmpdir}/$(( n - 1 ))/mni0_GenericAffine.xfm \
         --skip-nonlinear \
         --fixed-mask ${REGISTRATIONMODELMASK} \
         --moving-mask ${tmpdir}/$(( n - 1 ))/mnimask.mnc \
         ${tmpdir}/${n}/denoise.mnc ${REGISTRATIONMODEL} ${tmpdir}/${n}/mni

antsApplyTransforms -d 3 --verbose -i ${tmpdir}/${n}/denoise.mnc -r ${RESAMPLEMODEL} \
    -t ${tmpdir}/${n}/mni0_GenericAffine.xfm -n BSpline[5] -o ${tmpdir}/${n}/mni.mnc

antsApplyTransforms -d 3 --verbose -i ${tmpdir}/$(( n - 1 ))/mnimask.mnc -r ${RESAMPLEMODEL} \
    -t ${tmpdir}/${n}/mni0_GenericAffine.xfm -n GenericLabel -o ${tmpdir}/${n}/mnimask_in_mni.mnc

mincmath -clamp -const2 0 65535 ${tmpdir}/${n}/mni.mnc ${tmpdir}/${n}/mni.clamp.mnc
mv -f ${tmpdir}/${n}/mni.clamp.mnc ${tmpdir}/${n}/mni.mnc

volume_pol --order 1 --min 0 --max 100 --noclamp ${tmpdir}/${n}/mni.mnc ${RESAMPLEMODEL} \
  --source_mask ${tmpdir}/${n}/mnimask_in_mni.mnc --target_mask ${RESAMPLEMASK} --clobber ${tmpdir}/${n}/mni.norm.mnc
mincbeast -verbose -fill -median -same_res -flip -v2 -conf ${BEASTLIBRARY_DIR}/default.1mm.conf ${BEASTLIBRARY_DIR} ${tmpdir}/${n}/mni.norm.mnc ${tmpdir}/${n}/beastmask.mnc

antsApplyTransforms -d 3 -i ${tmpdir}/${n}/beastmask.mnc -t [ ${tmpdir}/${n}/mni0_GenericAffine.xfm,1 ] -n GenericLabel --verbose \
    -r ${tmpdir}/originput.mnc -o ${tmpdir}/${n}/bmask.mnc

ImageMath 3 ${tmpdir}/${n}/bmask_temp.mnc m ${tmpdir}/${n}/bmask.mnc ${tmpdir}/vessels.mnc
ThresholdImage 3 ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/bmask_fix.mnc Otsu 4 ${tmpdir}/${n}/bmask_temp.mnc
ThresholdImage 3 ${tmpdir}/${n}/bmask_fix.mnc ${tmpdir}/${n}/bmask_fix.mnc 2 Inf 1 0
iMath 3 ${tmpdir}/${n}/bmask_fix.mnc ME ${tmpdir}/${n}/bmask_fix.mnc 1 1 ball 1
ImageMath 3 ${tmpdir}/${n}/bmask_fix.mnc GetLargestComponent ${tmpdir}/${n}/bmask_fix.mnc
iMath 3 ${tmpdir}/${n}/bmask_fix.mnc MD ${tmpdir}/${n}/bmask_fix.mnc 1 1 ball 1
iMath 3 ${tmpdir}/${n}/bmask_fix.mnc MC ${tmpdir}/${n}/bmask_fix.mnc 5 1 ball 1
ImageMath 3 ${tmpdir}/${n}/bmask_fix.mnc FillHoles ${tmpdir}/${n}/bmask_fix.mnc 2
iMath 3 ${tmpdir}/${n}/bmask_fix.mnc MD ${tmpdir}/${n}/bmask_fix.mnc 1 1 ball 1

cp -f ${tmpdir}/${n}/bmask_fix.mnc ${tmpdir}/bmask_fix.mnc

antsRegistration_affine_SyN.sh --clobber --verbose \
    --fast --mask-extract \
    --initial-transform ${tmpdir}/${n}/mni0_GenericAffine.xfm \
    --skip-linear --fixed-mask ${REGISTRATIONMODELMASK} --moving-mask ${tmpdir}/${n}/bmask_fix.mnc \
    ${tmpdir}/${n}/denoise.mnc ${REGISTRATIONMODEL} ${tmpdir}/${n}/mni

antsApplyTransforms -d 3 -i ${WMPRIOR} -t [ ${tmpdir}/${n}/mni0_GenericAffine.xfm,1 ] -t ${tmpdir}/${n}/mni1_inverse_NL.xfm \
    -n Linear --verbose \
    -r ${tmpdir}/originput.mnc -o ${tmpdir}/${n}/prior3.mnc
antsApplyTransforms -d 3 -i ${GMPRIOR} -t [ ${tmpdir}/${n}/mni0_GenericAffine.xfm,1 ] -t ${tmpdir}/${n}/mni1_inverse_NL.xfm \
    -n Linear --verbose \
    -r ${tmpdir}/originput.mnc -o ${tmpdir}/${n}/prior2.mnc
antsApplyTransforms -d 3 -i ${CSFPRIOR} -t [ ${tmpdir}/${n}/mni0_GenericAffine.xfm,1 ] -t ${tmpdir}/${n}/mni1_inverse_NL.xfm \
    -n Linear --verbose \
    -r ${tmpdir}/originput.mnc -o ${tmpdir}/${n}/prior1.mnc

if [[ ! -z ${DEEPGMPRIOR} ]]; then
  antsApplyTransforms -d 3 -i ${DEEPGMPRIOR} -t [ ${tmpdir}/${n}/mni0_GenericAffine.xfm,1 ] -t ${tmpdir}/${n}/mni1_inverse_NL.xfm \
      -n Linear --verbose \
      -r ${tmpdir}/originput.mnc -o ${tmpdir}/${n}/prior4.mnc
fi

iMath 3 ${tmpdir}/${n}/bmask_D.mnc MD ${tmpdir}/${n}/bmask_fix.mnc 1 1 ball 1

ImageMath 3 ${tmpdir}/${n}/atropos_mask.mnc m ${tmpdir}/${n}/bmask_D.mnc ${tmpdir}/vessels.mnc
ImageMath 3 ${tmpdir}/${n}/atropos_mask.mnc m ${tmpdir}/${n}/atropos_mask.mnc ${tmpdir}/$(( n - 1 ))/hotmask.mnc

if [[ ! -z ${DEEPGMPRIOR} ]]; then
  Atropos --verbose -d 3 -a ${tmpdir}/${n}/denoise.mnc -x ${tmpdir}/${n}/atropos_mask.mnc  -c [ 25, 0.005 ] \
      -m [ 0.1,1x1x1 ] --posterior-formulation Aristotle[ 0 ]  -s 1x2 -s 2x3 -s 1x3 -s 1x4 -s 3x4 \
      -l [ 0.69314718055994530942,1 ] \
      -i PriorProbabilityImages[ 4,${tmpdir}/${n}/prior%d.mnc,0.25 ] -o [ ${tmpdir}/${n}/classify1.mnc,${tmpdir}/${n}/posterior%d.mnc ] \
      --winsorize-outliers BoxPlot
else
  Atropos --verbose -d 3 -a ${tmpdir}/${n}/denoise.mnc -x ${tmpdir}/${n}/atropos_mask.mnc  -c [ 25, 0.005 ] \
      -m [ 0.1,1x1x1 ] --posterior-formulation Aristotle[ 0 ]  -s 1x2 -s 2x3 -s 1x3 \
      -l [ 0.69314718055994530942,1 ] \
      -i PriorProbabilityImages[ 3,${tmpdir}/${n}/prior%d.mnc,0.25 ] -o [ ${tmpdir}/${n}/classify1.mnc,${tmpdir}/${n}/posterior%d.mnc ] \
      --winsorize-outliers BoxPlot
fi

ThresholdImage 3 ${tmpdir}/${n}/classify1.mnc ${tmpdir}/${n}/weight.mnc 2 Inf 1 0
iMath 3 ${tmpdir}/${n}/weight.mnc ME ${tmpdir}/${n}/weight.mnc 1 1 ball 1
ImageMath 3 ${tmpdir}/${n}/weight.mnc GetLargestComponent ${tmpdir}/${n}/weight.mnc
iMath 3 ${tmpdir}/${n}/weight.mnc MD ${tmpdir}/${n}/weight.mnc 1 1 ball 1

ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/vessels.mnc
make_outlier_map ${tmpdir}/${n}/denoise.mnc ${tmpdir}/fgmask.mnc ${tmpdir}/${n}/hotmask.mnc
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/hotmask.mnc

isostep=2.0
n3input=${tmpdir}/$(( n - 1 ))/correct.mnc
do_N3

iMath 3 ${tmpdir}/${n}/correct_mask.mnc MD ${tmpdir}/${n}/bmask_fix.mnc 2 1 ball 1
ImageMath 3 ${tmpdir}/${n}/correct_mask.mnc FillHoles ${tmpdir}/${n}/correct_mask.mnc

mincmath -clobber -unsigned -double -mult ${tmpdir}/${n}/*field.mnc ${tmpdir}/${n}/field_combined.mnc

correct_field ${tmpdir}/${n}/field_combined.mnc ${tmpdir}/${n}/correct_mask.mnc ${tmpdir}/${n}/field_combined_correct.mnc
mincmath -clobber -clamp -const2 0.1 1.79769e+308 ${tmpdir}/${n}/field_combined_correct.mnc ${tmpdir}/${n}/field_combined_correct_clamp.mnc

mincmath -clobber -unsigned -double -mult ${tmpdir}/$(( n - 1 ))/field_combined_correct_clamp.mnc ${tmpdir}/${n}/field_combined_correct_clamp.mnc ${tmpdir}/${n}/field_combined_correct_clamp2.mnc
mv -f ${tmpdir}/${n}/field_combined_correct_clamp2.mnc ${tmpdir}/${n}/field_combined_correct_clamp.mnc


origmean=$(mincstats -mean -quiet -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/originput.mnc)
biasmean=$(mincstats -mean -quiet -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/${n}/field_combined_correct_clamp.mnc)

minccalc -unsigned -short -expression "clamp((A[0]/${origmean})/(A[1]/${biasmean})*32767,0,65535)" ${tmpdir}/originput.mnc ${tmpdir}/${n}/field_combined_correct_clamp.mnc ${tmpdir}/${n}/correct.mnc -clobber
cp -f ${tmpdir}/${n}/correct.mnc ${tmpdir}/corrected.mnc

minc_anlm --clobber --mt $(nproc) ${tmpdir}/corrected.mnc ${tmpdir}/denoise_corrected.mnc

if [[ ! -z ${DEEPGMPRIOR} ]]; then
  Atropos --verbose -d 3 -a ${tmpdir}/denoise_corrected.mnc -x ${tmpdir}/${n}/atropos_mask.mnc -c [ 25, 0.005 ] \
      -m [ 0.1,1x1x1 ] --posterior-formulation Aristotle[ 1 ] -s 1x2 -s 2x3 -s 1x3 -s 1x4 -s 3x4 \
      -l [ 0.69314718055994530942,1 ] \
      -i PriorProbabilityImages[ 4,${tmpdir}/${n}/posterior%d.mnc,0.25 ] -o ${tmpdir}/${n}/classify2.mnc \
      --winsorize-outliers BoxPlot
else
  Atropos --verbose -d 3 -a ${tmpdir}/denoise_corrected.mnc -x ${tmpdir}/${n}/atropos_mask.mnc -c [ 25, 0.005 ] \
      -m [ 0.1,1x1x1 ] --posterior-formulation Aristotle[ 1 ] -s 1x2 -s 2x3 -s 1x3 \
      -l [ 0.69314718055994530942,1 ] \
      -i PriorProbabilityImages[ 3,${tmpdir}/${n}/posterior%d.mnc,0.25 ] -o ${tmpdir}/${n}/classify2.mnc \
      --winsorize-outliers BoxPlot
fi

valuelow=$(mincstats -quiet -floor 1 -pctT 0.1 ${tmpdir}/corrected.mnc)
valuewm=$(mincstats -quiet -median -mask ${tmpdir}/${n}/classify2.mnc -mask_binvalue 3 ${tmpdir}/corrected.mnc)
valuegm=$(mincstats -quiet -median -mask ${tmpdir}/${n}/classify2.mnc -mask_binvalue 2 ${tmpdir}/corrected.mnc)
valuehigh=$(mincstats -quiet -floor 1 -pctT 99.9 ${tmpdir}/corrected.mnc)

mapping=($(python -c "import numpy as np; print(np.array2string(np.linalg.solve(np.array([[1, ${valuelow}, ${valuelow}**2], [1, ((${valuewm}+${valuegm})/2.0), ((${valuewm}+${valuegm})/2.0)**2], [1, ${valuehigh}, ${valuehigh}**2]]),np.array([0,32767,65535])),separator= ' ')[1:-1])"))

#Re pad final image using model FOV mask
ImageMath 3 ${tmpdir}/corrected.mnc PadImage ${tmpdir}/corrected.mnc 50

antsApplyTransforms -d 3 --verbose -i ${tmpdir}/fgmask_fov.mnc -n GenericLabel \
    -o ${tmpdir}/fgmask_fov.mnc \
    -r ${tmpdir}/corrected.mnc

ImageMath 3 ${tmpdir}/corrected.mnc m ${tmpdir}/corrected.mnc ${tmpdir}/fgmask_fov.mnc
ExtractRegionFromImageByMask 3 ${tmpdir}/corrected.mnc ${tmpdir}/repad.mnc ${tmpdir}/fgmask_fov.mnc 1 $(calc "int(10.0*4.0/${shrink})")
cp -f ${tmpdir}/repad.mnc ${tmpdir}/corrected.mnc

mincresample -unsigned -short ${tmpdir}/corrected.mnc ${output}

minccalc -clobber -quiet ${N4_VERBOSE:+-verbose} -short -unsigned -expression "clamp(A[0]^2*${mapping[2]} + A[0]*${mapping[1]} + ${mapping[0]},0,65535)" \
    ${tmpdir}/corrected.mnc $(dirname ${output})/$(basename ${output} .mnc).rescale.mnc

mincresample -like ${output} -keep -near -unsigned -byte -labels ${tmpdir}/bmask_fix.mnc $(dirname ${output})/$(basename ${output} .mnc).mask.mnc
mincresample -like ${output} -keep -near -unsigned -byte -labels ${tmpdir}/${n}/classify2.mnc $(dirname ${output})/$(basename ${output} .mnc).classify.mnc

make_qc

#Create LSQ6 version of affine transform
xfminvert ${tmpdir}/${n}/mni0_GenericAffine.xfm ${tmpdir}/mni0_GenericAffine_invert.xfm
param2xfm $(xfm2param ${tmpdir}/mni0_GenericAffine_invert.xfm | grep -E 'scale|shear') ${tmpdir}/scaleshear.xfm
xfminvert ${tmpdir}/scaleshear.xfm ${tmpdir}/unscaleshear.xfm
xfmconcat ${tmpdir}/mni0_GenericAffine_invert.xfm ${tmpdir}/unscaleshear.xfm ${tmpdir}/lsq6.xfm

mincresample -unsigned -short -tfm_input_sampling -transform ${tmpdir}/lsq6.xfm ${tmpdir}/corrected.mnc $(dirname ${output})/$(basename ${output} .mnc).lsq6.mnc
mincresample -transform ${tmpdir}/lsq6.xfm -like $(dirname ${output})/$(basename ${output} .mnc).lsq6.mnc \
  -keep -near -unsigned -byte -labels ${tmpdir}/bmask_fix.mnc $(dirname ${output})/$(basename ${output} .mnc).lsq6.mask.mnc
mincresample -transform ${tmpdir}/lsq6.xfm -like $(dirname ${output})/$(basename ${output} .mnc).lsq6.mnc \
  -keep -near -unsigned -byte -labels ${tmpdir}/${n}/classify2.mnc $(dirname ${output})/$(basename ${output} .mnc).lsq6.classify.mnc
