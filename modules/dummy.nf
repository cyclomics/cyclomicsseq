#!/usr/bin/env nextflow
/*
========================================================================================
    Cyclomics/CycloSeq
========================================================================================
    Github : https://github.com/cyclomics/cycloseq
    Website: https://cyclomics.com
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// This module will just create an output file (with _2 to the basename eq test.tmp -> test_2.tmp) from the input file by copying it.

process DummyProcess{
  input:
    path something
  output:
    path "file.txt"
    
  script:
  // println 'Warning! Dummy process' 
  """
    echo dummy123 > file.txt
  """
}
