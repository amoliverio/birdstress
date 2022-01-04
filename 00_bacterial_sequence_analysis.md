## Demultiplex reads
```
python prep_fastq_for_uparse_paired.py -o 16S_demult -m 16S_map.txt -i Madden_16s_NoIndex_L001_R1_001.fastq -r Madden_16s_NoIndex_L001_R3_001.fastq -b Madden_16s_NoIndex_L001_R2_001.fastq -c

#copy reads to a new folder for analysis
cp 16S_demult/demultiplexed_seqs_* /home/pattyjk/Bird_Stress_Project/Reads/mergedfastq/
```
## Fix headers on reads
```
cd ..
cd mergedfastq

sed -i 's/barcodelabel/sample/g' demultiplexed_seqs_1.fq
sed -i 's/barcodelabel/sample/g' demultiplexed_seqs_2.fq
sed -i 's/ 1:N:0:/1:N:0/g' demultiplexed_seqs_1.fq
sed -i 's/ 3:N:0:/3:N:0/g' demultiplexed_seqs_2.fq

cd ..
```

## Join paired ends
```
#usearch v 10.0.24
./usearch64 -fastq_mergepairs mergedfastq/demultiplexed_seqs_1.fq -reverse mergedfastq/demultiplexed_seqs_2.fq -fastqout mergedfastq/16S_merged.fq -fastq_maxdiffs 10 -fastq_merge_maxee 1.0
##88% could be joined, normal amount
```

## Dereplicate sequences
```
./usearch64 -fastx_uniques mergedfastq/16S_merged.fq -fastqout mergedfastq/16S_uniques_combined_merged.fastq -sizeout
```

## Remove Singletons
```
./usearch64 -sortbysize mergedfastq/16S_uniques_combined_merged.fastq -fastqout mergedfastq/16S_nosigs_uniques_combined_merged.fastq -minsize 2
```

## Precluster sequences
```
./usearch64 -cluster_fast mergedfastq/16S_nosigs_uniques_combined_merged.fastq -centroids_fastq mergedfastq/16S_denoised_nosigs_uniques_combined_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size
```

## Closed reference pick against Silva v. 1.23
```
./usearch64 -usearch_global mergedfastq/16S_denoised_nosigs_uniques_combined_merged.fastq -id 0.97 -db silva_132_97_16S.fna  -strand plus -uc mergedfastq/16S_ref_seqs.uc -dbmatched mergedfastq/16S_closed_reference.fasta -notmatchedfq mergedfastq/16S_failed_closed.fq
```

## Sort by size and de novo pick
```
./usearch64 -sortbysize mergedfastq/16S_failed_closed.fq -fastaout mergedfastq/16S_sorted_failed_closed.fq

./usearch64 -cluster_otus mergedfastq/16S_sorted_failed_closed.fq -minsize 2 -otus mergedfastq/16S_denovo_otus.fasta -relabel OTU_dn_ -uparseout mergedfastq/16S_denovo_out.up
```

## Combine the rep sets between de novo and reference-based OTU picking rep sets

``` 
cat mergedfastq/16S_closed_reference.fasta mergedfastq/16S_denovo_otus.fasta > mergedfastq/16S_full_rep_set.fna
```

## Map rep. set back to pre-dereplicated sequences and make OTU tables
```
./usearch64 -usearch_global mergedfastq/16S_merged.fq -db mergedfastq/16S_full_rep_set.fna  -strand plus -id 0.97 -uc mergedfastq/16S_OTU_map.uc -otutabout mergedfastq/16S_OTU_table.txt -biomout mergedfastq/16S_OTU_jsn.biom
```

## Assign taxonomy to the OTUs in QIIME (v. 1.9.1) with the RDP Classifier
```
source activate qiime1
assign_taxonomy.py -i mergedfastq/16S_full_rep_set.fna -o mergedfastq/16S_taxonomy -t '/home/pattyjk/SILVA_132_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_7_levels.txt' -r '/home/pattyjk/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna' -m rdp
```

## Add taxonomy to OTU table
```
biom convert -i mergedfastq/16S_OTU_table.txt -o mergedfastq/16S_OTU_table.biom --to-hdf5 --table-type='OTU table'

biom add-metadata -i mergedfastq/16S_OTU_table.biom -o mergedfastq/16S_table_tax.biom --observation-metadata-fp=mergedfastq/16S_taxonomy/16S_full_rep_set_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy
```

## Filter out unwanted taxa: Chloroplasts, mitochondria, and Archaea
```
#filter unwanted bacterial taxa
filter_taxa_from_otu_table.py -i  mergedfastq/16S_table_tax.biom -o mergedfastq/16S_table_tax_filt.biom -n D_4__Mitochondria,D_3__Chloroplast,D_0__Archaea

#summarize taxonomy
summarize_taxa.py -i mergedfastq/16S_table_tax_filt.biom -o mergedfastq/16S_taxa_sum
```

## Convert filtered OTU tables to text files for arrr
```
biom convert -i mergedfastq/16S_table_tax_filt.biom -o mergedfastq/16S_table_tax_filt.txt --table-type='OTU table' --header-key=taxonomy --to-tsv
```

## Split raw fastq reads by sampe for submission to NCBI
```
#QIIME 1.9.1

#demultiplex each read individually 
split_libraries_fastq.py -i Madden_16s_NoIndex_L001_R1_001.fastq -b Madden_16s_NoIndex_L001_R2_001.fastq -m 16S_map.txt --store_demultiplexed_fastq -o demult_r1 --barcode_type 12 --rev_comp_mapping_barcodes

split_libraries_fastq.py -i Madden_16s_NoIndex_L001_R3_001.fastq -b Madden_16s_NoIndex_L001_R2_001.fastq -m 16S_map.txt --store_demultiplexed_fastq -o demult_r2 --barcode_type 12 --rev_comp_mapping_barcodes

#split demultiplexed fastq files by sample ID
split_sequence_file_on_sample_ids.py -i demult_r1/seqs.fastq -o R1_by_sample --file_type fastq
split_sequence_file_on_sample_ids.py -i demult_r2/seqs.fastq -o R2_by_sample --file_type fastq
```

