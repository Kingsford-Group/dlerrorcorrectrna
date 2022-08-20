#!/bin/bash
#
# Usage: run.minimap2_alignment.sh <longreads_data_filename> <Organism>
#
# <Organism>: human or mouse

longreads_data=$1
organism=$2

curr_dir="$PWD"

bin_dir=/home/ltung/miniconda3/bin

if [ $organism == 'mouse' ]
then
    echo "mouse"
    ref_genome=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/GRCm38/GRCm38.fa
    ref_annotation=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/gtf/Mus_musculus.GRCm38.92.gtf
    ref_transcriptome=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/cDNA/Mus_musculus.GRCm38.cdna.all.fa
else
    echo "human"
    ref_genome=/mnt/disk40/user/ltung/transcriptomics/ensembl/human/GRCh38/GRCh38.fa
    ref_annotation=/home/mingfus/data/transcriptomics/ensembl/human/gtf/Homo_sapiens.GRCh38.90.gtf
    ref_transcriptome=/mnt/disk27/user/ltung/longreadscallop/data/longreads/mashmap/Homo_sapiens.GRCh38.cdna.all.fa
fi

minimap2_dir=minimap2_align_output
if [ ! -d $minimap2_dir ]
then
    mkdir $minimap2_dir
fi

# Run minimap2
echo "Running minimap2..."
{ time $bin_dir/minimap2 --cs -ax splice -t 85 $ref_genome $longreads_data -o $minimap2_dir/longreads.sam &> $minimap2_dir/run.minimap2.log; } 2> $minimap2_dir/run.minimap2.time.log
echo "Done minimap2"

