# nf-core/rnaseq pipeline parameters

RNA sequencing analysis pipeline for gene/isoform quantification and extensive quality control.

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row. See [usage docs](https://nf-co.re/rnaseq/usage#samplesheet-input).</small></details>| `string` |  | True |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |

## Reference genome options

Reference genome related files and options required for the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `genome` | Name of iGenomes reference. <details><summary>Help</summary><small>If using a reference genome configured in the pipeline using iGenomes (not recommended), use this parameter to give the ID for the reference. This is then used to build the full paths for all required reference genome files e.g. `--genome GRCh38`. <br><br>See the [nf-core website docs](https://nf-co.re/usage/reference_genomes) for more details.</small></details>| `string` |  |  |  |
| `fasta` | Path to FASTA genome file. <details><summary>Help</summary><small>This parameter is *mandatory* if `--genome` is not specified. If you don't have the appropriate alignment index available this will be generated for you automatically. Combine with `--save_reference` to save alignment index for future runs.</small></details>| `string` |  |  |  |
| `gtf` | Path to GTF annotation file. <details><summary>Help</summary><small>This parameter is *mandatory* if `--genome` is not specified.</small></details>| `string` |  |  |  |
| `gff` | Path to GFF3 annotation file. <details><summary>Help</summary><small>This parameter must be specified if `--genome` or `--gtf` are not specified.</small></details>| `string` |  |  |  |
| `gene_bed` | Path to BED file containing gene intervals. This will be created from the GTF file if not specified. | `string` |  |  |  |
| `transcript_fasta` | Path to FASTA transcriptome file. | `string` |  |  |  |
| `additional_fasta` | FASTA file to concatenate to genome FASTA file e.g. containing spike-in sequences. <details><summary>Help</summary><small>If provided, the sequences in this file will get concatenated to the existing genome FASTA file, a GTF file will be automatically created using the entire sequence as the gene, transcript, and exon features, and any alignment index will get created from the combined FASTA and GTF. It is recommended to save the reference with `--save_reference` to re-use the index for future runs so you do not need to create it again.</small></details>| `string` |  |  |  |
| `splicesites` | Splice sites file required for HISAT2. | `string` |  |  |  |
| `star_index` | Path to directory or tar.gz archive for pre-built STAR index. | `string` |  |  |  |
| `hisat2_index` | Path to directory or tar.gz archive for pre-built HISAT2 index. | `string` |  |  |  |
| `rsem_index` | Path to directory or tar.gz archive for pre-built RSEM index. | `string` |  |  |  |
| `salmon_index` | Path to directory or tar.gz archive for pre-built Salmon index. | `string` |  |  |  |
| `kallisto_index` | Path to directory or tar.gz archive for pre-built Kallisto index. | `string` |  |  |  |
| `hisat2_build_memory` | Minimum memory required to use splice sites and exons in the HiSAT2 index build process. <details><summary>Help</summary><small>HiSAT2 requires a huge amount of RAM to build a genome index for larger genomes, if including splice sites and exons e.g. the human genome might typically require 200GB. If you specify less than this threshold for the `HISAT2_BUILD` process then the splice sites and exons will be ignored, meaning that the process will require a lot less memory. If you are working with a small genome, set this parameter to a lower value to reduce the threshold for skipping this check. If using a larger genome, consider supplying more memory to the `HISAT2_BUILD` process.</small></details>| `string` | 200.GB |  |  |
| `gencode` | Specify if your GTF annotation is in GENCODE format. <details><summary>Help</summary><small>If your GTF file is in GENCODE format and you would like to run Salmon i.e. `--pseudo_aligner salmon`, you will need to provide this parameter in order to build the Salmon index appropriately.</small></details>| `boolean` |  |  |  |
| `gtf_extra_attributes` | By default, the pipeline uses the `gene_name` field to obtain additional gene identifiers from the input GTF file when running Salmon. <details><summary>Help</summary><small>This behaviour can be modified by specifying `--gtf_extra_attributes` when running the pipeline. Note that you can also specify more than one desired value, separated by a comma e.g. `--gtf_extra_attributes gene_id,...`.<br></small></details>| `string` | gene_name |  |  |
| `gtf_group_features` | Define the attribute type used to group features in the GTF file when running Salmon. | `string` | gene_id |  |  |
| `featurecounts_group_type` | The attribute type used to group feature types in the GTF file when generating the biotype plot with featureCounts. | `string` | gene_biotype |  |  |
| `featurecounts_feature_type` | By default, the pipeline assigns reads based on the 'exon' attribute within the GTF file. <details><summary>Help</summary><small>The feature type used from the GTF file when generating the biotype plot with featureCounts.</small></details>| `string` | exon |  |  |
| `igenomes_base` | Directory / URL base for iGenomes references. | `string` | s3://ngi-igenomes/igenomes |  | True |
| `igenomes_ignore` | Do not load the iGenomes reference config. <details><summary>Help</summary><small>Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.</small></details>| `boolean` |  |  | True |

## Read trimming options

Options to adjust read trimming criteria.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `trimmer` | Specifies the trimming tool to use - available options are 'trimgalore' and 'fastp'. | `string` | trimgalore |  |  |
| `extra_trimgalore_args` | Extra arguments to pass to Trim Galore! command in addition to defaults defined by the pipeline. | `string` |  |  |  |
| `extra_fastp_args` | Extra arguments to pass to fastp command in addition to defaults defined by the pipeline. | `string` |  |  |  |
| `min_trimmed_reads` | Minimum number of trimmed reads below which samples are removed from further processing. Some downstream steps in the pipeline will fail if this threshold is too low. | `integer` | 10000 |  |  |

## Read filtering options

Options for filtering reads prior to alignment

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `bbsplit_fasta_list` | Path to comma-separated file containing a list of reference genomes to filter reads against with BBSplit. You have to also explicitly set `--skip_bbsplit false` if you want to use BBSplit. <details><summary>Help</summary><small>The file should contain 2 columns: short name and full path to reference genome(s) e.g. <br>```<br>mm10,/path/to/mm10.fa<br>ecoli,/path/to/ecoli.fa<br>```</small></details>| `string` |  |  |  |
| `bbsplit_index` | Path to directory or tar.gz archive for pre-built BBSplit index. <details><summary>Help</summary><small>The BBSplit index will have to be built at least once with this pipeline (see `--save_reference` to save index). It can then be provided via `--bbsplit_index` for future runs.</small></details>| `string` |  |  |  |
| `remove_ribo_rna` | Enable the removal of reads derived from ribosomal RNA using SortMeRNA. <details><summary>Help</summary><small>Any patterns found in the sequences defined by the '--ribo_database_manifest' parameter will be used.</small></details>| `boolean` |  |  |  |
| `ribo_database_manifest` | Text file containing paths to fasta files (one per line) that will be used to create the database for SortMeRNA. <details><summary>Help</summary><small>By default, [rRNA databases](https://github.com/biocore/sortmerna/tree/master/data/rRNA_databases) defined in the SortMeRNA GitHub repo are used. You can see an example in the pipeline Github repository in `assets/rrna-default-dbs.txt`.<br>Please note that commercial/non-academic entities require [`licensing for SILVA`](https://www.arb-silva.de/silva-license-information) for these default databases.</small></details>| `string` | /SGRNJ06/randd/USER/zhouyiqi/2024/repo/rnaseq/assets/rrna-db-defaults.txt |  |  |

## UMI options

Options for processing reads with unique molecular identifiers

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_umi` | Enable UMI-based read deduplication. | `boolean` |  |  |  |
| `umitools_extract_method` | UMI pattern to use. Can be either 'string' (default) or 'regex'. <details><summary>Help</summary><small>More details can be found in the [UMI-tools documentation](https://umi-tools.readthedocs.io/en/latest/reference/extract.html#extract-method).<br></small></details>| `string` | string |  |  |
| `umitools_bc_pattern` | The UMI barcode pattern to use e.g. 'NNNNNN' indicates that the first 6 nucleotides of the read are from the UMI. <details><summary>Help</summary><small>More details can be found in the [UMI-tools documentation](https://umi-tools.readthedocs.io/en/latest/reference/extract.html#extract-method).</small></details>| `string` |  |  |  |
| `umitools_bc_pattern2` | The UMI barcode pattern to use if the UMI is located in read 2. | `string` |  |  |  |
| `umi_discard_read` | After UMI barcode extraction discard either R1 or R2 by setting this parameter to 1 or 2, respectively. | `integer` |  |  |  |
| `umitools_umi_separator` | The character that separates the UMI in the read name. Most likely a colon if you skipped the extraction with UMI-tools and used other software. | `string` |  |  |  |
| `umitools_grouping_method` | Method to use to determine read groups by subsuming those with similar UMIs. All methods start by identifying the reads with the same mapping position, but treat similar yet nonidentical UMIs differently. | `string` | directional |  |  |
| `umitools_dedup_stats` | Generate output stats when running "umi_tools dedup". <details><summary>Help</summary><small>It can be quite time consuming generating these output stats - see [#827](https://github.com/nf-core/rnaseq/issues/827).</small></details>| `boolean` |  |  |  |

## Alignment options

Options to adjust parameters and filtering criteria for read alignments.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `aligner` | Specifies the alignment algorithm to use - available options are 'star_salmon', 'star_rsem' and 'hisat2'. | `string` | star_salmon |  |  |
| `pseudo_aligner` | Specifies the pseudo aligner to use - available options are 'salmon'. Runs in addition to '--aligner'. | `string` |  |  |  |
| `pseudo_aligner_kmer_size` | Kmer length passed to indexing step of pseudoaligners <details><summary>Help</summary><small>Failure to set a good kmer size could cause issues with quantification with Kallisto or Salmon. This is mostly an issue for short reads (<50bp), where the default kmer size of 31 is an problem.</small></details>| `integer` | 31 |  |  |
| `bam_csi_index` | Create a CSI index for BAM files instead of the traditional BAI index. This will be required for genomes with larger chromosome sizes. | `boolean` |  |  |  |
| `star_ignore_sjdbgtf` | When using pre-built STAR indices do not re-extract and use splice junctions from the GTF file. | `boolean` |  |  |  |
| `salmon_quant_libtype` |  Override Salmon library type inferred based on strandedness defined in meta object. <details><summary>Help</summary><small>See [Salmon docs](https://salmon.readthedocs.io/en/latest/library_type.html).</small></details>| `string` |  |  |  |
| `min_mapped_reads` | Minimum percentage of uniquely mapped reads below which samples are removed from further processing. <details><summary>Help</summary><small>Some downstream steps in the pipeline will fail if this threshold is too low.</small></details>| `number` | 5.0 |  |  |
| `seq_center` | Sequencing center information to be added to read group of BAM files. | `string` |  |  |  |
| `stringtie_ignore_gtf` | Perform reference-guided de novo assembly of transcripts using StringTie i.e. dont restrict to those in GTF file. | `boolean` |  |  |  |
| `extra_star_align_args` | Extra arguments to pass to STAR alignment command in addition to defaults defined by the pipeline. Only available for the STAR-Salmon route. | `string` |  |  |  |
| `extra_salmon_quant_args` | Extra arguments to pass to Salmon quant command in addition to defaults defined by the pipeline. | `string` |  |  |  |
| `extra_kallisto_quant_args` | Extra arguments to pass to Kallisto quant command in addition to defaults defined by the pipeline. | `string` |  |  |  |
| `kallisto_quant_fraglen` | In single-end mode Kallisto requires an estimated fragment length. Specify a default value for that here. TODO: use existing RSeQC results to do this dynamically. | `integer` | 200 |  |  |
| `kallisto_quant_fraglen_sd` | In single-end mode, Kallisto requires an estimated standard error for fragment length. Specify a default value for that here. TODO: use existing RSeQC results to do this dynamically. | `integer` | 200 |  |  |

## Optional outputs

Additional output files produces as intermediates that can be saved

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `save_merged_fastq` | Save FastQ files after merging re-sequenced libraries in the results directory. | `boolean` |  |  |  |
| `save_umi_intermeds` | If this option is specified, intermediate FastQ and BAM files produced by UMI-tools are also saved in the results directory. | `boolean` |  |  |  |
| `save_non_ribo_reads` | If this option is specified, intermediate FastQ files containing non-rRNA reads will be saved in the results directory. | `boolean` |  |  |  |
| `save_bbsplit_reads` | If this option is specified, FastQ files split by reference will be saved in the results directory. | `boolean` |  |  |  |
| `save_reference` | If generated by the pipeline save the STAR index in the results directory. <details><summary>Help</summary><small>If an alignment index is generated by the pipeline use this parameter to save it to your results folder. These can then be used for future pipeline runs, reducing processing times.</small></details>| `boolean` |  |  |  |
| `save_trimmed` | Save the trimmed FastQ files in the results directory. <details><summary>Help</summary><small>By default, trimmed FastQ files will not be saved to the results directory. Specify this flag (or set to true in your config file) to copy these files to the results directory when complete.</small></details>| `boolean` |  |  |  |
| `save_align_intermeds` | Save the intermediate BAM files from the alignment step. <details><summary>Help</summary><small>By default, intermediate BAM files will not be saved. The final BAM files created after the appropriate filtering step are always saved to limit storage usage. Set this parameter to also save other intermediate BAM files.</small></details>| `boolean` |  |  |  |
| `save_unaligned` | Where possible, save unaligned reads from either STAR, HISAT2 or Salmon to the results directory. <details><summary>Help</summary><small>This may either be in the form of FastQ or BAM files depending on the options available for that particular tool.</small></details>| `boolean` |  |  |  |

## Quality Control

Additional quality control options.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `deseq2_vst` | Use vst transformation instead of rlog with DESeq2. <details><summary>Help</summary><small>See [DESeq2 docs](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization).</small></details>| `boolean` | True |  |  |
| `rseqc_modules` | Specify the RSeQC modules to run. | `string` | bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication |  |  |

## Process skipping options

Options to skip various steps within the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `skip_gtf_filter` | Skip filtering of GTF for valid scaffolds and/ or transcript IDs. <details><summary>Help</summary><small>If you're confident on the validity of the GTF with respect to the genome fasta file, or wish to disregard failures thriggered by the filtering module, activate this option.</small></details>| `boolean` |  |  |  |
| `skip_gtf_transcript_filter` | Skip the 'transcript_id' checking component of the GTF filtering script used in the pipeline. | `boolean` |  |  |  |
| `skip_bbsplit` | Skip BBSplit for removal of non-reference genome reads. | `boolean` | True |  |  |
| `skip_umi_extract` | Skip the UMI extraction from the read in case the UMIs have been moved to the headers in advance of the pipeline run. | `boolean` |  |  |  |
| `skip_trimming` | Skip the adapter trimming step. <details><summary>Help</summary><small>Use this if your input FastQ files have already been trimmed outside of the workflow or if you're very confident that there is no adapter contamination in your data.</small></details>| `boolean` |  |  |  |
| `skip_alignment` | Skip all of the alignment-based processes within the pipeline. | `boolean` |  |  |  |
| `skip_pseudo_alignment` | Skip all of the pseudoalignment-based processes within the pipeline. | `boolean` |  |  |  |
| `skip_markduplicates` | Skip picard MarkDuplicates step. | `boolean` |  |  |  |
| `skip_bigwig` | Skip bigWig file creation. | `boolean` |  |  |  |
| `skip_stringtie` | Skip StringTie. | `boolean` |  |  |  |
| `skip_fastqc` | Skip FastQC. | `boolean` |  |  |  |
| `skip_preseq` | Skip Preseq. | `boolean` | True |  |  |
| `skip_dupradar` | Skip dupRadar. | `boolean` |  |  |  |
| `skip_qualimap` | Skip Qualimap. | `boolean` |  |  |  |
| `skip_rseqc` | Skip RSeQC. | `boolean` |  |  |  |
| `skip_biotype_qc` | Skip additional featureCounts process for biotype QC. | `boolean` |  |  |  |
| `skip_deseq2_qc` | Skip DESeq2 PCA and heatmap plotting. | `boolean` |  |  |  |
| `skip_multiqc` | Skip MultiQC. | `boolean` |  |  |  |
| `skip_qc` | Skip all QC steps except for MultiQC. | `boolean` |  |  |  |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` |  |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |
| `test_data_base` | Base path / URL for data used in the test profiles <details><summary>Help</summary><small>Warning: The `-profile test` samplesheet file itself contains remote paths. Setting this parameter does not alter the contents of that file.</small></details>| `string` |  |  | True |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>| `integer` | 16 |  | True |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details>| `string` | 128.GB |  | True |
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>| `string` | 240.h |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `help` | Display help text. | `boolean` |  |  | True |
| `version` | Display version and exit. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `hook_url` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.</small></details>| `string` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | `string` |  |  | True |
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. | `string` |  |  |  |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `validationShowHiddenParams` | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.</small></details>| `boolean` |  |  | True |
| `validationFailUnrecognisedParams` | Validation of parameters fails when an unrecognised parameter is found. <details><summary>Help</summary><small>By default, when an unrecognised parameter is found, it returns a warinig.</small></details>| `boolean` |  |  | True |
| `validationLenientMode` | Validation of parameters in lenient more. <details><summary>Help</summary><small>Allows string values that are parseable as numbers or booleans. For further information see [JSONSchema docs](https://github.com/everit-org/json-schema#lenient-mode).</small></details>| `boolean` |  |  | True |

## Other parameters

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `protocol` |  | `string` | AccuraSCOPE-V1 |  |  |
| `well_sample` |  | `string` |  |  |  |
