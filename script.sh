INPUT_READ_FILE="./sample_files/test.ccs.bam"
INPUT_SUBREAD_FILE="./sample_files/test.subreads_to_ccs.bam"
OUTPUT_FILE="./sample_files/test.ccs.modified.bam"
NUMBER_THREADS=20
samtools view ${INPUT_SUBREAD_FILE} | cargo run --release -- ${NUMBER_THREADS} ${INPUT_READ_FILE}
python src/evaluator.py -i ${INPUT_READ_FILE} -o ${OUTPUT_FILE}