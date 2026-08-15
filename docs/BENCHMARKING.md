# Performance evaluation

- Use `-K 100000000` across all comparisons for reproducibility
- Use cram where possible
- Test aligner threads (12,16,24,32,48) and sort threads (2,4,8,12)

## Performance using recommended profile

- Uses [recommended](https://bwa-mem3.readthedocs.io/en/latest/best-practices/settings-profiles.html) settings `-y 0 --bam=0 --min-ext-len 30 --skip-contained-ext`

```bash
BWAMEM3INDEX=/mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/BWAmem3Index/Homo_sapiens_assembly38.fasta

FASTA=/mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta

TEST_R1=/home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_1.fastq.gz
TEST_R2=/home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_2.fastq.gz

OUTFILE=SRR7890943_WGS_cross-site_study.cram

# Align and sort
bwa-mem3 mem -t 30 -K 100000000 -m 10 -y 0 --min-ext-len 30 --bam=0 --skip-contained-ext $BWAMEM3INDEX $TEST_R1 $TEST_R2  | samtools sort -@ 8 --reference  -o $OUTFILE -

# Mark duplicates


# Calculate metrics
```


Time taken for main_mem function: 4534.06 sec


IO times (sec) :
Reading IO time (reads) avg: 885.09, (885.09, 885.09)
Writing IO time (SAM) avg: 2090.31, (2090.31, 2090.31)
Reading IO time (Reference Genome) avg: 0.00, (0.00, 0.00)
Index read time avg: 12.29, (12.29, 12.29)

Overall time (sec) (Excluding Index reading time):
PROCESS() (Total compute time + (read + SAM) IO time) : 4521.73
MEM_PROCESS_SEQ() (Total compute time (Kernel + SAM)), avg: 3444.41, (3444.41, 3444.41)

SAM Processing time (sec): --WORKER_SAM avg: 804.79, (804.79, 804.79)

Kernels' compute time (sec):
Total kernel (smem+sal+bsw) time avg: 2565.64, (2565.64, 2565.64)
SMEM compute avg: 719.75, (735.36, 714.29)
MEM_CHAIN avg: 547.07, (551.55, 537.14)
SAL compute avg: 545.05, (549.37, 535.05)
MEM_SA avg: 311.40, (313.22, 306.07)

BSW time, avg: 964.81, (965.68, 961.84)



## Performance using `--fast`

- `--fast` applies `-m 10 -y 0 --min-ext-len 30 --smem-dedup --skip-contained-ext --max-extend-chains 20 --adaptive-band --extend-mate-concordant`
- bwa-mem3’s best speed/accuracy trade-off

```bash
BWAMEM3INDEX=/mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/BWAmem3Index/Homo_sapiens_assembly38.fasta

FASTA=/mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta

TEST_R1=/home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_1.fastq.gz
TEST_R2=/home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_2.fastq.gz

OUTFILE=SRR7890943_WGS_cross-site_study.cram

# Align and sort
bwa-mem3 mem -t 16 -K 100000000 --fast $BWAMEM3INDEX $TEST_R1 $TEST_R2 | samtools sort -@ 16 --reference $FASTA -o $OUTFILE -

# Mark duplicates
samtools collate -O -u --threads 6 SRR7890943_WGS_cross-site_study.cram | \
samtools fixmate -m -u --threads 6 - - | \
samtools sort -u --threads 6 - | \
samtools markdup -f SRR7890943_WGS.md.metrics --threads 6 --reference /mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta --output-fmt cram - SRR7890943_WGS_cross-site_study.sorted.md.cram

# Calculate metrics
riker multi --tools basic isize alignment wgs
```

Time taken for main_mem function: 5081.57 sec

IO times (sec) :
Reading IO time (reads) avg: 846.12, (846.12, 846.12)
Writing IO time (SAM) avg: 1751.95, (1751.95, 1751.95)
Reading IO time (Reference Genome) avg: 0.00, (0.00, 0.00)
Index read time avg: 13.59, (13.59, 13.59)

Overall time (sec) (Excluding Index reading time):
PROCESS() (Total compute time + (read + SAM) IO time) : 5067.94
MEM_PROCESS_SEQ() (Total compute time (Kernel + SAM)), avg: 4088.36, (4088.36, 4088.36)

SAM Processing time (sec): --WORKER_SAM avg: 819.64, (819.64, 819.64)

Kernels' compute time (sec):
Total kernel (smem+sal+bsw) time avg: 3199.79, (3199.79, 3199.79)
SMEM compute avg: 1111.20, (1116.45, 1105.48)
MEM_CHAIN avg: 808.75, (813.07, 804.11)
SAL compute avg: 806.20, (810.50, 801.64)
MEM_SA avg: 453.07, (456.27, 449.36)
BSW time, avg: 577.40, (577.73, 576.24)


## Minibwa

```bash
MINIBWAINDEX=/mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/MinibwaIndex/Homo_sapiens_assembly38.fasta

FASTA=/mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta

TEST_R1=/home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_1.fastq.gz
TEST_R2=/home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_2.fastq.gz

OUTFILE=SRR7890943_WGS_cross-site_study.cram
OUTFILE2=SRR7890943_WGS_cross-site_study.sorted.md.cram

# Align and sort
minibwa map -t 16 -K 100m $MINIBWAINDEX $TEST_R1 $TEST_R2 | samtools sort -@ 8 --reference $FASTA -o $OUTFILE -

# Mark duplicates
samtools collate -O -u --threads 6 SRR7890943_WGS_cross-site_study.cram | \
samtools fixmate -m -u --threads 6 - - | \
samtools sort -u --threads 6 - | \
samtools markdup -f SRR7890943_WGS.md.metrics --threads 6 --reference $FASTA --output-fmt cram - $OUTFILE2

# Calculate metrics
riker multi --tools basic isize alignment wgs
```

[M::main] Version: 0.4-r400
[M::main] CMD: minibwa map -t 16 -K 100m /mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/MinibwaIndex/Homo_sapiens_assembly38.fasta /home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_1.fastq.gz /home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_2.fastq.gz
[M::main] Real time: 6297.376 sec; CPU: 88607.878 sec; Peak RSS: 10.228 GB
[bam_sort_core] merging from 17 files and 8 in-memory blocks...
(base) ubuntu@sarek-v5:~/test_data$ 


[M::main] Version: 0.4-r400
[M::main] CMD: minibwa map -t 30 -K 100m /mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/MinibwaIndex/Homo_sapiens_assembly38.fasta /home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_1.fastq.gz /home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_2.fastq.gz
[M::main] Real time: 4342.631 sec; CPU: 95693.064 sec; Peak RSS: 10.254 GB
[bam_sort_core] merging from 17 files and 8 in-memory blocks...


## Classic bwa mem


[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 30 -K 100000000 /mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta /home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_1.fastq.gz /home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_2.fastq.gz
[main] Real time: 21180.933 sec; CPU: 606038.558 sec
[bam_sort_core] merging from 21 files and 8 in-memory blocks...
(base) ubuntu@sarek-v5:~/test_data$ 



## GPU-accelerated bwa

```bash
REFERENCE_FILE=/home/debian/parabricks_bundle/Homo_sapiens_assembly38.fasta
KNOWN_SITES_FILE=/home/debian/parabricks_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

INPUT_FASTQ_1=/home/debian/test_data/SRR7890943_WGS_cross-site_study_1.fastq.gz
INPUT_FASTQ_2=/home/debian/test_data/SRQ7890943_WGS_cross-site_study_2.fastq.gz

OUTPUT_BAM=/home/debian/work/SRR7890943_WGS_cross-site_study.cram
OUTPUT_RECAL_FILE=/home/debian/work/SRR7890943_WGS_cross-site_study.recal.table

# Create the work directory if it doesn't exist yet
mkdir -p /home/debian/work

# Mount the home directory so all absolute paths map 1:1 inside the container
docker run --rm --gpus all \
    --volume /home/debian:/home/debian \
    --workdir /home/debian/work \
    nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1 \
    pbrun fq2bam \
        --ref ${REFERENCE_FILE} \
        --in-fq ${INPUT_FASTQ_1} ${INPUT_FASTQ_2} \
        --knownSites ${KNOWN_SITES_FILE} \
        --out-bam ${OUTPUT_BAM} \
        --out-recal-file ${OUTPUT_RECAL_FILE} \
        --align-only \
        --bwa-options "-K 10000000" \
        --bwa-cpu-thread-pool 32 \
        --verbose --monitor-usage
```


> Using 2x L4 GPUs

[PB Info 2026-Jul-24 13:07:10] Time spent monitoring (multiple of 10): 3630.426
[PB Info 2026-Jul-24 13:07:10] bwalib run finished in 3622.398 seconds
[PB Info 2026-Jul-24 13:07:10] ------------------------------------------------------------------------------
[PB Info 2026-Jul-24 13:07:10] ||        Program:                    GPU-PBBWA mem, Sorting Phase-I        ||
[PB Info 2026-Jul-24 13:07:10] ||        Version:                                           4.7.1-1        ||
[PB Info 2026-Jul-24 13:07:10] ||        Start Time:                       Fri Jul 24 12:06:40 2026        ||
[PB Info 2026-Jul-24 13:07:10] ||        End Time:                         Fri Jul 24 13:07:10 2026        ||
[PB Info 2026-Jul-24 13:07:10] ||        Total Time:                          60 minutes 30 seconds        ||
[PB Info 2026-Jul-24 13:07:10] ------------------------------------------------------------------------------




```bash
BWAMEM3INDEX=/home/stefan/references/Homo_sapiens/GATK/GRCh38/Sequence/BWAmem3Index/Homo_sapiens_assembly38.fasta

FASTA=/home/stefan/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta

TEST_R1=/home/stefan/references/test_data/SRR7890943_WGS_cross-site_study_1.fastq.gz
TEST_R2=/home/stefan/references/test_data/SRR7890943_WGS_cross-site_study_2.fastq.gz

OUTFILE=/home/stefan/references/test_data/SRR7890943_WGS_cross-site_study.cram
MARKDUPMETRICS=/home/stefan/references/test_data/SRR7890943_WGS.md.metrics
MARKDUPOUT=/home/stefan/references/test_data/SRR7890943_WGS_cross-site_study.sorted.md.cram


# Align and sort
bwa-mem3 mem -t 12 -K 100000000 -Y -m 10 -y 0 --min-ext-len 30 --bam=0 --skip-contained-ext $BWAMEM3INDEX $TEST_R1 $TEST_R2  | samtools sort -@ 3 --output-fmt cram --reference $FASTA -o $OUTFILE -

# Mark duplicates
samtools collate -O -u --threads 12 $OUTFILE | \
samtools fixmate -m -u --threads 12 - - | \
samtools sort -u --threads 12 - | \
samtools markdup -f $MARKDUPMETRICS --threads 12 --reference $FASTA --output-fmt cram - $MARKDUPOUT

# Calculate metrics
riker multi --threads 8 --tools basic isize alignment wgs --reference $FASTA --input $MARKDUPOUT --output SRR7890943_WGS_cross-site_study.sorted.md.cram_riker
```