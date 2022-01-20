#!/usr/bin/env nextflow
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

process foo {
    output:
      path 'foo.txt'
    script:
      """
      echo 'Defined in nf' >> foo.txt
      """
}

 process bar {
    input:
      path x
    output:
      path 'bar.txt'
    script:
      """
      cat $x > bar.txt
      """
}
 process bara {
    input:
      path x
    output:
      path 'bara.txt'
    script:
      """
      cat $x > bara.txt
      """
}

 process baz {
    input:
      path x
    output:
      path 'baz.txt'
    script:
      """
      cat $x >> baz.txt
      echo bonus_baz >> baz.txt
      """
}