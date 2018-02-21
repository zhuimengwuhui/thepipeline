'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from stages import Stages


def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='dog_pipeline')
    # Get a list of paths to all the FASTQ files
    fastq_files = state.config.get_option('fastqs')
    # Stages are dependent on the state
    stages = Stages(state)

    # The original FASTQ files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.original_fastqs,
        name='original_fastqs',
        output=fastq_files)

    # Align paired end reads in FASTQ to the reference producing a BAM file
    pipeline.transform(
        task_func=stages.align_bwa,
        name='align_bwa',
        input=output_from('original_fastqs'),
        # Match the R1 (read 1) FASTQ file and grab the path and sample name.
        # This will be the first input to the stage.
        # We assume the sample name may consist of only alphanumeric
        # characters.
        # filter=formatter('(?P<path>.+)/(?P<readid>[a-zA-Z0-9-\.]+)_(?P<lib>[a-zA-Z0-9-]+)_(?P<lane>[a-zA-Z0-9]+)_(?P<sample>[a-zA-Z0-9]+)_1.fastq.gz'),
        # Ot9399_index28_CAAAAG_L001_R1.fastq.gz
        filter=formatter(
            '.+/(?P<sample>[a-zA-Z0-9-]+)_(?P<readid>[a-zA-Z0-9-:]+)_(?P<lib>[a-zA-Z0-9]+)_(?P<lane>[a-zA-Z0-9-]+)_R1.fastq.gz'),
        # Add one more inputs to the stage:
        #    1. The corresponding R2 FASTQ file
        # e.g. C2WPF.5_Solexa-201237_5_X4311_1.fastq.gz
        add_inputs=add_inputs(
            '{path[0]}/{sample[0]}_{readid[0]}_{lib[0]}_{lane[0]}_R2.fastq.gz'),
        # Add an "extra" argument to the state (beyond the inputs and outputs)
        # which is the sample name. This is needed within the stage for finding out
        # sample specific configuration options
        extras=['{readid[0]}', '{lib[0]}', '{lane[0]}', '{sample[0]}'],
        # extras=['{sample[0]}'],
        # The output file name is the sample name with a .bam extension.
        output='alignments/{sample[0]}/{readid[0]}_{lib[0]}_{lane[0]}_{sample[0]}.bam')

    # Sort the BAM file using Picard
    pipeline.transform(
        task_func=stages.sort_bam_picard,
        name='sort_bam_picard',
        input=output_from('align_bwa'),
        filter=suffix('.bam'),
        output='.sort.bam')

    return pipeline
