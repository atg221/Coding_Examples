inputPath=${PWD}/
outputPath=${PWD}/TemplateMultivariateBSplineSyN/

${ANTSPATH}/antsMultivariateTemplateConstruction2.sh \
  -d 2 \
  -o ${outputPath}T_ \
  -g 0.2 \
  -c 0 \
  -k 2 \
  -s 3x2x1x0 \
  -q 100x70x0x0 \
  -n 1 \
  -r 1 \
  -l 1 \
  -m CC[2] \
  -t BSplineSyN[0.1,26,0] \
  ${inputPath}/*_main.nii.gz

