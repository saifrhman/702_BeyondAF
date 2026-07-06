# Central path config for COMP702 / Beyond AlphaFold reproduction

export COMP702_ROOT="${COMP702_ROOT:-$HOME/COMP702_BeyondAF}"
export COMP390_ROOT="${COMP390_ROOT:-$COMP702_ROOT/code/COMP390_code}"
export COMP702_RESULT_ROOT="${COMP702_RESULT_ROOT:-$COMP390_ROOT/Result}"

export COMP702_LOG_DIR="${COMP702_LOG_DIR:-$COMP702_ROOT/logs/slurm}"
export COMP702_SCRATCH_ROOT="${COMP702_SCRATCH_ROOT:-$HOME/scratch}"
export COMP702_FASTSCRATCH_ROOT="${COMP702_FASTSCRATCH_ROOT:-$HOME/fastscratch}"

# OpenFold code/environment
export OPENFOLD_CODE_DIR="${OPENFOLD_CODE_DIR:-$COMP390_ROOT/openfold}"
export OPENFOLD_ENV="${OPENFOLD_ENV:-$COMP702_ROOT/envs/openfold_env}"
export COMP702_DATA_ENV="${COMP702_DATA_ENV:-$COMP702_ROOT/envs/comp702-data}"

# OpenFold large data roots
export OPENFOLD_DATA_ROOT="${OPENFOLD_DATA_ROOT:-$COMP702_SCRATCH_ROOT/openfold_data}"
export OPENFOLD_RUNS_ROOT="${OPENFOLD_RUNS_ROOT:-$COMP702_SCRATCH_ROOT/openfold_runs}"
export TRITON_CACHE_DIR="${TRITON_CACHE_DIR:-$COMP702_FASTSCRATCH_ROOT/triton_cache}"

# RODA/OpenFold alignment download
export RODA_S3_SRC="${RODA_S3_SRC:-s3://openfold/pdb}"
export RODA_ALIGNMENT_DST="${RODA_ALIGNMENT_DST:-$COMP702_SCRATCH_ROOT/alignment_data/alignment_dir_roda}"

# Alignment extraction
export ALIGNMENT_ARCHIVE="${ALIGNMENT_ARCHIVE:-$COMP702_SCRATCH_ROOT/alignment_data_node037_20260210_083640.tar.zst}"
export ALIGNMENT_EXTRACT_DIR="${ALIGNMENT_EXTRACT_DIR:-$OPENFOLD_DATA_ROOT/alignment_openfold_full_20260210}"

# OpenFold dataset/cache paths
export PDB_DATA_ROOT="${PDB_DATA_ROOT:-$OPENFOLD_DATA_ROOT/pdb_data}"
export PDB_DATA_FULL_ROOT="${PDB_DATA_FULL_ROOT:-$OPENFOLD_DATA_ROOT/pdb_data_full}"
export ALIGNMENT_SUBSET_ROOT="${ALIGNMENT_SUBSET_ROOT:-$OPENFOLD_DATA_ROOT/alignment_subset}"
export ALIGNMENT_TESTSET_ROOT="${ALIGNMENT_TESTSET_ROOT:-$OPENFOLD_DATA_ROOT/alignment_testset}"

# Repeated OpenFold inference outputs
export OPENFOLD_RUNS5_ROOT="${OPENFOLD_RUNS5_ROOT:-$COMP702_RESULT_ROOT/openfold_inference_runs5times}"

# pLDDT analysis outputs
export PLDDT_OUT_DIR="${PLDDT_OUT_DIR:-$COMP702_RESULT_ROOT/plddt_analysis_results}"
export ALL_PLDDT_CSV="${ALL_PLDDT_CSV:-$PLDDT_OUT_DIR/all_plddt_summary.csv}"
export UNRELAXED_PLDDT_CSV="${UNRELAXED_PLDDT_CSV:-$PLDDT_OUT_DIR/unrelaxed_plddt_summary.csv}"

# Relaxation analysis outputs
export RELAX_ANALYSIS_BY_FILES_DIR="${RELAX_ANALYSIS_BY_FILES_DIR:-$COMP702_RESULT_ROOT/relax_analysis_results_by_files}"
export RELAX_SAMPLE_CSV="${RELAX_SAMPLE_CSV:-$RELAX_ANALYSIS_BY_FILES_DIR/relax_sample_table.csv}"
export RELAX_CHAIN_CSV="${RELAX_CHAIN_CSV:-$RELAX_ANALYSIS_BY_FILES_DIR/summary_by_chain.csv}"
export RELAX_EPOCH_CSV="${RELAX_EPOCH_CSV:-$RELAX_ANALYSIS_BY_FILES_DIR/summary_by_epoch.csv}"
export RELAX_REPEAT_CSV="${RELAX_REPEAT_CSV:-$RELAX_ANALYSIS_BY_FILES_DIR/summary_by_repeat.csv}"

# Error-vs-pLDDT analysis outputs
export ERROR_PLDDT_OUT_DIR="${ERROR_PLDDT_OUT_DIR:-$COMP702_RESULT_ROOT/error_plddt_analysis_results}"

# Task 01: test-set cache generation
export TESTSET_MMCIF_DIR="${TESTSET_MMCIF_DIR:-$PDB_DATA_ROOT/mmcif/mmcif_files}"
export TESTSET_CACHE_DIR="${TESTSET_CACHE_DIR:-$PDB_DATA_ROOT/data_caches}"
export TESTSET_CLUSTER_FILE="${TESTSET_CLUSTER_FILE:-$ALIGNMENT_TESTSET_ROOT/all-seqs_clusters-40.txt}"
export TESTSET_MMCIF_CACHE="${TESTSET_MMCIF_CACHE:-$TESTSET_CACHE_DIR/mmcif_cache_testset.json}"
export TESTSET_CHAIN_CACHE="${TESTSET_CHAIN_CACHE:-$TESTSET_CACHE_DIR/chain_data_cache_testset.json}"

# Task 01: OpenFold subset training paths
export SUBSET_MMCIF_DIR="${SUBSET_MMCIF_DIR:-$PDB_DATA_ROOT/mmcif/mmcif_files}"
export SUBSET_ALIGN_DB_DIR="${SUBSET_ALIGN_DB_DIR:-$ALIGNMENT_SUBSET_ROOT/alignment_dbs}"
export SUBSET_ALIGN_INDEX="${SUBSET_ALIGN_INDEX:-$SUBSET_ALIGN_DB_DIR/alignment_db.index}"
export SUBSET_CHAIN_CACHE="${SUBSET_CHAIN_CACHE:-$PDB_DATA_ROOT/data_caches/chain_data_cache.json}"
export SUBSET_MMCIF_CACHE="${SUBSET_MMCIF_CACHE:-$PDB_DATA_ROOT/data_caches/mmcif_cache.json}"
export SUBSET_OBSOLETE_FILE="${SUBSET_OBSOLETE_FILE:-$PDB_DATA_ROOT/obsolete.dat}"

# Task 01: OpenFold resume-training checkpoints
export OPENFOLD_RESUME_CKPT_EP27="${OPENFOLD_RESUME_CKPT_EP27:-$OPENFOLD_RUNS_ROOT/run_subset_h100_ep27/csv_logs/version_0/checkpoints/26-270000.ckpt}"
export OPENFOLD_RESUME_CKPT_EP33="${OPENFOLD_RESUME_CKPT_EP33:-$OPENFOLD_RUNS_ROOT/run_subset_h100_ep33/csv_logs/version_0/checkpoints/32-330000.ckpt}"

# Task 02: optional GPU-local data staging
export OPENFOLD_GPU_LOCAL_BASE="${OPENFOLD_GPU_LOCAL_BASE:-/localscratch/$USER}"
export OPENFOLD_GPU_DATA_ROOT="${OPENFOLD_GPU_DATA_ROOT:-$OPENFOLD_GPU_LOCAL_BASE/openfold_data}"
export OPENFOLD_GPU_DATA_TAR="${OPENFOLD_GPU_DATA_TAR:-$OPENFOLD_GPU_LOCAL_BASE/openfold_data.tar.gz}"
export OPENFOLD_GPU_RUNS_ROOT="${OPENFOLD_GPU_RUNS_ROOT:-$OPENFOLD_GPU_LOCAL_BASE/openfold_runs}"

# Task 03: OpenFold inference inputs and outputs
export OPENFOLD_INFERENCE_FASTA_DIR="${OPENFOLD_INFERENCE_FASTA_DIR:-$COMP702_SCRATCH_ROOT/openfold_inference/fasta_dir}"
export OPENFOLD_INFERENCE_CIF_DIR="${OPENFOLD_INFERENCE_CIF_DIR:-$COMP702_SCRATCH_ROOT/taskcif}"
export OPENFOLD_INFERENCE_ALIGN_DIR="${OPENFOLD_INFERENCE_ALIGN_DIR:-$COMP702_SCRATCH_ROOT/openfold_inference/alignments_cf}"
export OPENFOLD_INFERENCE_OUT_ROOT="${OPENFOLD_INFERENCE_OUT_ROOT:-$COMP702_RESULT_ROOT/openfold_inference}"
export OPENFOLD_INFERENCE_RUNS_ROOT="${OPENFOLD_INFERENCE_RUNS_ROOT:-$COMP702_RESULT_ROOT/openfold_inference_runs}"

# OpenFold parameter files
export OPENFOLD_AF_MODEL1_PARAM="${OPENFOLD_AF_MODEL1_PARAM:-$OPENFOLD_CODE_DIR/openfold/resources/params/params_model_1.npz}"
export OPENFOLD_OFFICIAL_FINETUNING_CKPT="${OPENFOLD_OFFICIAL_FINETUNING_CKPT:-$OPENFOLD_CODE_DIR/openfold/resources/openfold_params/finetuning_2.pt}"

# Trained checkpoint root
export OPENFOLD_TRAINED_CKPT_ROOT="${OPENFOLD_TRAINED_CKPT_ROOT:-$OPENFOLD_RUNS_ROOT}"

# Task 03: copied/archived trained checkpoints available inside the COMP390 codebase
export OPENFOLD_ARCHIVE_CKPT_ROOT="${OPENFOLD_ARCHIVE_CKPT_ROOT:-$COMP390_ROOT/checkpoints}"
export OPENFOLD_EP32_CKPT="${OPENFOLD_EP32_CKPT:-$OPENFOLD_ARCHIVE_CKPT_ROOT/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/32-330000.ckpt}"
