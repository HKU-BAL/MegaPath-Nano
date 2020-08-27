### Options
```
General options:
    --query
        Query file (fastq or fasta)
    --human
        Human genome set in config folder, default human.genome_set.
    --decoy
        Decoy genome set in config folder, default plasmid.genome_set.
    --species
        Genome set for species identification in config folder, default species_id.genome_set.
    --assembly
        Genome set for assembly identification in config folder, default assembly_id.genome_set.
   
Environment config:
    --python 
        Path to python, default python3.
    --temp_folder
        Temprary folder, default ''.
    --RAM_folder
        Temporary folder in RAM, default '/run/shm'.
    --taxonomy_db
        Taxonomy database, default 'db/ncbi_taxonomy.db'.
    --tool_folder
        Tool folder, default 'tools/'.
    --config_folder
        Config folder, default 'config/'.
    --assembly_folder
        Assembly folder, default 'genomes/'.
    --aligner
        Path to aligner, default 'minima2'.
    --read_simulator
        Path to read simulation program, default 'tools/nanosim/simulator.py'.
    --read_simulation_profiles
        Path to read simulation profiles, default 'tools/nanosim/nanosim_profiles'.
    --human_similar_filter_assembly_id
        Assembly ID for human similar region filter, default 'GCF_000001405.37'.
    --max_aligner_thread INT
        Maximum number of threads used by aligner, default 64.
    --max_qcat_thread INT
        Maximum number of threads used by qcat, default 64.
    --genus_height INT
        Height in taxonomy to be considered as genus, default 7.

Fastq filter options:
    --head_crop INT
        Crop INT bp from read head, default 0.
    --tail_crop INT
        Crop INT bp from read tail, default 0.
    --min_read_length INT
        Reads with length smaller than INT read length are filtered, default 500.
    --min_read_quality FLOAT
        Reads with average base quality smaller than FLOAT are filtered, default 7.0.

Alignment filter options:
    --adaptor_trimming, --no-adaptor_trimming
        Adapter trimming is on by default.
    --read_trimming, --no-read_trimming
        Read trimming is on by default.
    --read_filter, --no-read_filter
        Read quality filter is on by default.
    --human_filter, --no-human_filter
        Human filter is on by default.
    --decoy_filter, --no-decoy_filter
        Decoy filter is on by default.
    --variable_region_adjustment, --no-variable_region_adjustment
		Default false.
    --spike_filter, --no-spike_filter
        Spike filter is on by default.
    --closing_spike_filter, --no-closing_spike_filter
        Closing spike filter is on by default.
    --human_similar_filter, --no-human_similar_filter
        Human similar region filter is on by default.
    --microbe_similar_filter, --no-microbe_similar_filter
        Micorbe similar region filter is on by default.
    --short_alignment_filter, --no-short_alignment_filter
        Short alignment filter is on by default.
    --unique_alignment, --no-unique_alignment
        Unique alignment filter is on by default
    --noise_projection, --no-noise_projection
        off
    --similar_species_marker, --no-similar_species_marker
        on
    --mapping_only, --no-mapping_only
        off
    --reassign_read_id, --no-reassign_read_id
        off
    --all_steps, --filter_fq_only

    --min_alignment_score INT
        Minimal alignment score, default 0.
    --human_filter_alignment_score_threshold INT
        Alignment score threshold for flagging a read as a human read, default 1000.
    --human_filter_alignment_score_percent_threshold INT
        Alignment score (normalized by read length) threshold (in percent) for flagging a read as a human read, default 100.
    --decoy_filter_alignment_score_threshold INT
        Alignment score threshold for flagging a read as a decoy read, default 1000.
    --decoy_filter_alignment_score_percent_threshold INT
        Alignment score (normalized by read length) threshold (in percent) for flagging a read as a decoy read, default 100.
    --species_id_min_aligned_bp INT
        Minimal covered BP to include a species for analysis, default 0.
    --good_alignment_threshold INT
        Alignment score threshold in percentage of best alignment score, default 80.
    --assembly_id_min_average_depth FLOAT
        Minimal average depth to perform assembly selection, default 0.5.
    --variable_region_percent INT
        Maximum percentage of strands aligned for a region to be labeled as variable, default 50.
    --expected_max_depth_stdev INT
        Number of standard deviations for calculating expected max depth, default 6.
    --microbe_similar_filter_abundance_threshold_80 FLOAT
        Difference (no. of times) in apparent abundance to trigger similar region filter with 80% similarity, default 160.
    --microbe_similar_filter_abundance_threshold_90 FLOAT
        Difference (no. of times) in apparent abundance to trigger similar region filter with 90% similarity, default 80.
    --microbe_similar_filter_abundance_threshold_95 FLOAT
        Difference (no. of times) in apparent abundance to trigger similar region filter with 95% similarity, default 40.
    --microbe_similar_filter_abundance_threshold_98 FLOAT
        Difference (no. of times) in apparent abundance to trigger similar region filter with 98% similarity, default 16.
    --microbe_similar_filter_abundance_threshold_99 FLOAT
        Difference (no. of times) in apparent abundance to trigger similar region filter with 99% similarity, default 8.
    --microbe_similar_filter_abundance_threshold_99_2 FLOAT
        Difference (no. of times) in apparent abundance to trigger similar region filter with 99.2% similarity, default 6.4.
    --microbe_similar_filter_targeted_max_span_percent INT
        Maximum percent of regions (targeted) to be marked as similar region, default 90.
    --microbe_similar_filter_allowed_max_span_percent INT
        Maximum percent of regions (allowed) to be marked as similar region, default 97.
    --microbe_similar_filter_min_average_depth FLOAT
        Minimum average depth to be considered as source of noise, default 0.2.
    --microbe_similar_filter_max_span_percent_overall INT
        Maximum percent of regions to be marked as similar region (overall), default 97.
    --max_alignment_noise_overlap INT
        The maximum percent for an alignment to overlap with noise regions without being removed, default 50.
    --min_alignment_length INT
        Minimum alignment length to be considered as evidence, default 250.
    --closing_expected_max_depth_stdev INT
        Number of standard deviations for calculating expected max depth for closing spike filter, default 9.
    --unique_alignment_threshold INT
        Unique alignments shall have no alignments with alignment score within this percent, default 80.
    --number_of_genus_to_perform_noise_projection INT
        Number of genus to perform noise projection, default 3.
    --min_percent_abundance_to_perform_noise_projection INT
        Minimum percent of abundance relative to the most abundant species in a genus to perform noise projection, default 25
    --noise_projection_simulated_read_length_bin_size INT
        Read length bin size for generating simulated reads, default 1000
    --noise_projection_simulated_read_length_multiplier FLOAT
        Multiplier over average read length to obtain maximum read length, default 0.5.
    --noise_projection_simulated_read_error_profile 
        Error profile for generating simulated reads, default 'ecoli_R91D'.
    --noise_projection_num_read_to_simulate INT
        Number of simulated reads to generate, default 10000.
    --similar_species_marker_num_genus INT
        Number of top most abundant species (1 per genus) to be considered as possible source of noise, default 3.
    --similar_species_marker_alignment_similarity_1 INT
        Similarity cutoff (1) used for alignment, available choices are 99, 98, 95, 90, 80. default 98.
    --similar_species_marker_aligned_region_threshold_1 INT
        Percentage of aligned region (1) to be considered as highly similar, default 50.
    --similar_species_marker_alignment_similarity_2 INT
        Similarity cutoff (2) used for alignment, available choices are 99, 98, 95, 90, 80. default 95.
    --similar_species_marker_aligned_region_threshold_2 INT
        Percentage of aligned region (2) to be considered as highly similar, default 75.
    --similar_species_marker_similarity_combine_logic
        Logic for combining criteria 1 and 2, choices are 'and', 'or', default 'or'.

Output options:
    --output_prefix
        Output Prefix, query file name will be used for output prefix by default.
    --output_folder
        Output folder, default ./.
    --archive_format
        Format used for output archive file, choices are 'zip', 'tar', 'gztar', 'bztar', default 'gztar'.
    --quality_score_bin_size FLOAT
        Bin size for quality score histogram, default 0.2.
    --read_length_bin_size INT
        Bin size for read length histogram, default 100.
    --output_adaptor_trimmed_query, --no-output_adaptor_trimmed_query
        Fastq after adaptor trimming is not outputted by default.
    --output_trimmed_and_filtered_query, --no-output_trimmed_and_filtered_query
        Fastq after trimming and quality filtering is outputted by default.
    --output_human_and_decoy_filtered_query, --no-output_human_and_decoy_filtered_query
        Fastq passed human and decoy filter is outputted by default.
    --output_PAF, --no-output_PAF
        Alignment in PAF format is outputted by default.
    --output_noise_stat, --no-output_noise_stat
        Noise regions are outputted by default.
    --output_separate_noise_bed, --no-output_separate_noise_bed
        default yes
    --output_raw_signal, --no-output_raw_signal
        default yes
    --output_id_signal, --no-output_id_signal
        default yes
    --output_per_read_data, --no-output_per_read_data
        default yes
    --output_quality_score_histogram, --no-output_quality_score_histogram
        default yes
    --output_read_length_histogram, --no-output_read_length_histogram
        default yes
    --output_human_stat, --no-output_human_stat
        yes
    --output_decoy_stat, --no-output_decoy_stat
        yes
    --output_genome_set, --no-output_genome_set
        default yes
    --aligner_log
        Log for stderr output from aligner program, default 'minimap2.log'.
    --read_sim_log
        Log for stderr output from read simulator, default 'read_sim.log'.
    --adaptor_trimming_log
        Log for stdout output from adaptor trimming program, default 'adaptor_trimming.log'.
```
