#!/bin/bash

for directory in job*;
do
  echo "DIRECTORY: ${directory}"
  (cd "${directory}" ; mv *.data_pipeline.tar.gz inference_inputs/ )
done

exit 0
