# Batch sort

This is spike sorting software that's part of the standard data processing
pipeline for the lab. There are more detailed instructions in `SAC_doc.pdf`, 
but the main usage is calling `sac_batch` with arguments directing it towards
the NEV files to sort. Note that it will overwrite these files with the sorted
versions, so you should only ever use it on a copy of your data.
