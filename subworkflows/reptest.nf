#!/usr/bin/env nextflow
/*
========================================================================================
    Cyclomics/CycloSeq subworkflow
========================================================================================
    Github : https://github.com/cyclomics/cycloseq
    Website: https://cyclomics.com
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

include {
    foo;
    bar;
    bara;
    baz
} from "../modules/dummy.nf"

workflow Testwf {
   main:
     foo()
     bar(foo.out)
   emit:
     my_data = bar.out
}

workflow Testwf2 {
   take:
     data
   main:
     bara(data)
     baz(bara.out)
   emit:
     my_data = baz.out
}

workflow SimpleTidehunter {
  take:
  main:
  emit:
  
}