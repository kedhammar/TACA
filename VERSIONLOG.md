# TACA Version Log

## 20240131.1

This is a test.

## 20240123.1

Exclude pod5 dir and files from being copied to metadata dir.

## 20240122.1

Adapt ONT analysis to new ONT JSON format (also backwards compatible).

## 20231204.1

Update ONT instrument transfer script to ignore runs started in the 3rd PromethION column, which will be used by Clinical Genomics.

## 20231201.1

Run Anglerfish v0.6.0 with --lenient and --ont_barcodes options.

## 20231130.1

Version 1.0.0
(1) Support BCL Convert
(2) Remove redundant codes and obsoleted platforms
(3) Switch MiSeq to Illumina V3 samplesheet
(4) Other refactors to improve performance

## 20231103.1

Fix bug with rsync permission issue cont.

## 20231031.1

Improve run_folder transfer

## 20231026.1

Fix bug with rsync permission issue

## 20231024.1

Fix functionality issues for Anglerfish launch by running via "conda run", fix erroneous file name reference and improve static typing.

## 20231023.1

Remove redundant cleanup_processing function in cleanup

## 20231016.1

Support multiple project ID for run folder transfer

## 20230927.1

Fix bug that NovaSeqXPlus date format cause error in writing pdc_archived timestamp

## 20230921.1

Remove the temp change of creating links

## 20230920.1

Supplement last PR, primary purpose is to differentiate user runs from QC runs in the instrument transfer script rather than the installed TACA.

## 20230915.1

Major overhaul of Nanopore stuff. Use generalized local script instead of installed TACA for both instruments and harmonize the way TACA handles Nanopore data from preprocessing. Implement automated handling of MinION QC runs.

## 20230913.1

Fix bugs for encrypting and archiving runs on ngi-preproc cont.

## 20230905.1

Fix bugs for encrypting and archiving runs on ngi-preproc

## 20230903.1

Adapt MinKNOW .json trimming to new format from Dorado update.

## 20230823.1

Allow manual database update of finished ONT runs

## 20230822.1

Add pandas to requiresments to accomodate last PR

## 20230821.1

Copy ONT metadata to ngi-nas-ns.

## 20230814.1

Update path to store Anglerfish results.

## 20230810.1

Implement logging for PromethION script.

## 20230809.2

Update handling of MinION QC runs to not run nanoseq, only anglerfish.

## 20230809.1

Update PromethION script to extend the scope of the log file parsing.

## 20230724.1

Enable TACA to retrieve error and warnings in bcl2fastq logs

## 20230718.1

Update PromethION script to run rsync w. -u flag and clarify archiving code.

## 20230713.1

Let PromethION script search through device logs to dump flow cell pore count history into the run dir.

## 20230711.1

Rework how PromethION script detects runs to catch mis-named ones, too.

## 20230621.1

Add support for NovaSeqXPlus and improve readability

## 20230609.1

Add functionality to update DB of specified run dirs.

## 20230607.1

Trim out unused data acquisition outputs from ONT report .json files before sending them to CouchDB.

## 20230510.1

Add storage_systems to server_status command to allow disk space surveillance of mounted virtual NAS:es. Also added Dockerfile and devcontainer setup.

## 20230503.1

Change how MinKNOW reports are synced to GenStat server to increase traceability and enable transfer of the reports of finished runs.

## 20230502.1

Enforce MinKNOW reports retain MinKNOW run ID upon transfer to ngi-internal. Improve logging.

## 20230428.1

Change offload location in promethion_transfer.py

## 20230419.1

Use a hidden file to indicate when the final rsync of ONT data to ngi-nas is done

## 20230331.1

Move MinKNOW reports to ngi-internal instead of embedding in StatusDB doc

## 20230307.1

Handle demux case that ordered read length is different from seq setup

##20230213.1
Further updates to ONT runs TACA <-> CouchDB interplay after local troubleshooting. Improve code readability, logging and exception handling.

##20230207.1
Add functionality for monitoring PromethION status

##20230117.2
More robust handling of ONT transfers

##20230117.1
Integrate ONT data flow with CouchDB nanopore_runs

##20221102.1
Include MinION data in transfer to HPC cluster

##20221028.1
Cleaner check for ONT transfer destination

##20221011.1
Add versioning for PromethION offload script

##20220830.1
Change promethion directory levels

## 20220811.1

Set short single index to 8nt

## 20220614.1

Updates to promethion offload script

## 20220613.1

Include promethion offload script

##20220610.1
Convert statusdb urls to https

## 20220524.1

Handle special demux for NoIndex cases

## 20220520.1

Include index list from undet of complex lanes

## 20220427.1

Support addtional 10X index types

## 20220414.1

Allow 0 mismatch in demux for short single index for MiSeq

## 20220412.3

Bug fix and refactor handling for 10X indexes

## 20220412.2

Handle cases with different UMI lengths

## 20220412.1

Add promethion to server_status

## 20220409.1

Small refactor for NoIndex

## 20220404.1

Change cleanup irma to cleanup miarka

## 20220330.1

Remove "Irma" from email/log

## 20220328.1

Re-define signal for sequencing done

## 20220314.1

Refactor cleanup_nas

## 20220302.1

Fix bug that samplesheet of MiSeq is overwritten

## 20220220.1

Change samplesheet location for MiSeq

## 20211216.1

Add option for syncing MinION delivery data to storage

## 20211208.1

Updates for compability with nanoseq v2.0.1

## 20211201.1

Only process MinION QC runs

## 20211011.1

Refactor of MinION code

## 20211005.1

Specify CopyComplete.txt as indicator for cleaning NovaSeq runs

## 20210819.1

Allow 0 mismatch in demux for short single index

## 20210716.1

Handle special case with index 2 only for NovaSeq

## 20210617.1

Support addtional 10X index types

## 20210604.1

Handle special case with index 2 only

## 20210330.1

Small adjust of 10X ATACseq masks

## 20210323.1

Support different FC ID pattern for NextSeq2000 in bioinfo_tab.py

## 20210318.1

Support different FC ID pattern for NextSeq2000

## 20210313.1

Support addtional 10X index types

## 20210303.1

Fix bug that runtype was overwritten

## 20210302.2

Fix FC name pattern for NextSeq2000

## 20210302.1

Setup VERSIONLOG.md
