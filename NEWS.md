#### DoddLabMetID 0.1.0
* Added a `NEWS.md` file to track changes to the package.

#### DoddLabMetID 0.1.1
* Modify some functions to automatically modify the peak table
* Add parameter set for calling
* Modify the MS2 plot function and add more details in the plot
* Fix the MS2 scoring bug

#### DoddLabMetID 0.1.11
* Fix the errors in the merge_one_mode function
* Update README file

#### DoddLabMetID 0.1.12
* Modify read_ms2 function to extract ms2_purity for mzml/mzxml data
* add the ms2_purify output

#### DoddLabMetID 0.1.13
* add peptide database & modify related functions
* modify initialize_annotation_parameter_class function to assign ms2_file in the ms2 folder

#### DoddLabMetID 0.1.14
* add export of combine_ms1_ms2 function

#### DoddLabMetID 0.1.15
* add gnps_bile_acid, gnps_acyl_amides, gnps_acyl_esters, all_public options in database
* add parameter export
* add object slot in the AnnotationParameterClass

#### DoddLabMetID 0.1.16
* adjust the merge_one_modes function
* add merge and export ms2 spec for manual check

#### DoddLabMetID 0.1.17 (20240424)
* add msdial_lipid option, and class_adduct_list parameter

#### DoddLabMetID 0.1.18 (20240626)
* modify the order function (dodd lab id, RT) for level 1 annotation in the merge_one_mode function

#### DoddLabMetID 0.1.19 (20240706)
* add functions: correct_annot_table, compare_previous_annotation, and add_class_info

#### DoddLabMetID 0.1.20 (20241015)
* add functions: generate_input_from_msp
* add a ms2_file variable in the AnnotationParameterClass to support defined ms2_file
* fix the bug that msp file can't be correctly aligned with MS1 table
