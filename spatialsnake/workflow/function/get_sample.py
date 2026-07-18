import os
import sys
import re
import glob
from spatialsnake.workflow.function.stereoseq_spec import parse_stereoseq_input_spec
from spatialsnake.workflow.function.logging_utils import setup_logger

logger = setup_logger("sample_parser")

def build_compare_sample_key(sample_id, group_id):
  return f"{group_id}::{sample_id}"


def infer_subset_name(input_path):
    """Return a stable output name for a split SpatialData input.

    Splitting writes names such as ``celltype_selected_Tumor.zarr`` and
    ``clusters_selected_1.zarr``.  The technical prefix is removed so that
    downstream directories describe the selected population itself.
    """
    input_name = os.path.basename(os.path.normpath(str(input_path).strip()))
    if input_name.lower().endswith(".zarr"):
        input_name = input_name[:-5]
    input_name = re.sub(
        r"^(?:(?:celltype|clusters?|recluster)_selected_|cluster_)",
        "",
        input_name,
        flags=re.IGNORECASE,
    )
    subset_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", input_name).strip("_.")
    if not subset_name:
        raise ValueError(f"Cannot derive a cell-type subset name from input path: {input_path}")
    return subset_name


def build_subset_analysis_tasks(samples, input_paths):
    """Pair parent sample IDs with independently addressable subset inputs."""
    if len(samples) != len(input_paths):
        raise ValueError("sample IDs and subset input paths must have the same length")
    tasks = []
    used = set()
    for parent_sample, input_path in zip(samples, input_paths):
        parent_sample = str(parent_sample).strip()
        if not parent_sample:
            raise ValueError("parent sample_id cannot be empty")
        subset_name = infer_subset_name(input_path)
        task_key = (parent_sample, subset_name)
        if task_key in used:
            raise ValueError(
                "Duplicate subset output detected for parent sample "
                f"'{parent_sample}' and subset '{subset_name}'. Rename one input zarr."
            )
        used.add(task_key)
        tasks.append(
            {
                "parent_sample": parent_sample,
                "subset_name": subset_name,
                "input_path": str(input_path),
            }
        )
    return tasks


def read_compare_sample_table(sample_list_file):
    """Read the three-column compare_analysis sample sheet.

    Like the other Spatialsnake sample readers, columns may be separated by
    one or more whitespace characters (spaces or tabs).  Platform options and
    filtering thresholds belong in YAML and are intentionally rejected here.
    """
    if not os.path.isfile(sample_list_file):
        raise FileNotFoundError(f"compare_analysis sample table not found: {sample_list_file}")
    expected = ["sample_id", "input_path", "group"]
    rows = []
    header_seen = False
    with open(sample_list_file, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            fields = re.split(r"\s+", stripped)
            if not header_seen:
                header_seen = True
                if fields != expected:
                    raise RuntimeError(
                        "compare_analysis sample.txt must contain exactly three whitespace-separated "
                        "columns in this order: sample_id, input_path, group"
                    )
                continue
            if len(fields) != len(expected):
                raise RuntimeError(
                    f"sample.txt line {line_number} must contain exactly three whitespace-separated fields"
                )
            rows.append(dict(zip(expected, fields)))
    if not header_seen:
        raise RuntimeError("compare_analysis sample.txt is empty")
    if not rows:
        raise RuntimeError("compare_analysis sample.txt contains no sample rows")
    sample_ids = [row["sample_id"] for row in rows]
    seen = set()
    duplicated = set()
    for sample in sample_ids:
        if sample in seen:
            duplicated.add(sample)
        seen.add(sample)
    if duplicated:
        raise RuntimeError("Duplicate sample_id values in sample.txt: " + ", ".join(sorted(duplicated)))
    return rows

def check_file_exit(type,dir_path):
  invalid_samples=[]
  if type == "stereoseq":
    if os.path.isfile(dir_path):
      lower_name = dir_path.lower()
      if lower_name.endswith(".cellbin.gef") or lower_name.endswith(".gef") or ".gem" in lower_name:
        return "stereoseq"
      logger.error(f"Stereo-seq input file not supported: {dir_path}")
      return False
    if os.path.isdir(dir_path):
      feature_dir = os.path.join(dir_path, "feature_expression")
      if os.path.isdir(feature_dir):
        stereo_patterns = [
          "*.cellbin.gef",
          "*cellbin*.gef",
          "*.tissue.gef",
          "*.gef",
          "*.gem",
          "*.gem.gz",
        ]
        for pattern in stereo_patterns:
          if glob.glob(os.path.join(feature_dir, pattern)):
            return "stereoseq"
        logger.error(f"Stereo-seq feature_expression files not found under {feature_dir}")
        return False
      stereo_patterns = [
        "**/*.cellbin.gef",
        "**/*cellbin*.gef",
        "**/*.gef",
        "**/*.gem",
        "**/*.gem.gz",
        "**/*.tif",
        "**/*.tiff",
      ]
      for pattern in stereo_patterns:
        if glob.glob(os.path.join(dir_path, pattern), recursive=True):
          return "stereoseq"
      logger.error(f"Stereo-seq key files not found under {dir_path}. Expect *.cellbin.gef, *.gef, *.gem(.gz), or tif images")
      return False
  if os.path.isdir(dir_path):
    if type == "Merfish":
      has_cell_by_gene = len(glob.glob(os.path.join(dir_path, "**", "cell_by_gene.csv"), recursive=True)) > 0
      has_transcripts = len(glob.glob(os.path.join(dir_path, "**", "*transcripts*.csv*"), recursive=True)) > 0
      has_transcripts_parquet = len(glob.glob(os.path.join(dir_path, "**", "*transcripts*.parquet"), recursive=True)) > 0
      if not (has_cell_by_gene or has_transcripts or has_transcripts_parquet):
        logger.error(f"MERSCOPE/MERFISH key files not found under {dir_path}. Expect cell_by_gene.csv and/or detected_transcripts*.csv/parquet")
        return False
      return "merscope"
    required_files = {
        'visium': [
            'spatial/tissue_positions_list.csv',
            'spatial/scalefactors_json.json',
            'spatial/tissue_lowres_image.png',
            'spatial/tissue_hires_image.png'],
        'visium_HD': [
            'spatial/tissue_positions.parquet',
            'spatial/scalefactors_json.json',
            'spatial/tissue_lowres_image.png'],
        "visium_segment":[
            "spatial/tissue_hires_image.png",
            "spatial/scalefactors_json.json",
            "cell_segmentations.geojson"],
        'xenium': [
            'cells.parquet',
            'transcripts.parquet',
            'morphology.ome.tif',
            'experiment.xenium']}
    for file_pattern in required_files[type]:
      if os.path.isfile(os.path.join(dir_path,file_pattern)):
        continue
      else:
        invalid_samples.append(file_pattern)
        return False
    main_candidates = [
      'filtered_feature_bc_matrix.h5',
      'raw_feature_bc_matrix.h5',
      "cell_feature_matrix.h5",
      "filtered_feature_cell_matrix.h5",
      "raw_feature_cell_matrix.h5"
    ]
    for candidate in main_candidates:
      if os.path.isfile(os.path.join(dir_path, candidate)):
        return candidate
    dir_files = sorted(os.listdir(dir_path))
    for candidate in main_candidates:
      soft_matched = [
        file_name for file_name in dir_files
        if file_name.endswith(candidate) and os.path.isfile(os.path.join(dir_path, file_name))
      ]
      if soft_matched:
        return soft_matched[0]
  logger.error(f"Invalid data path detected: {dir_path}; invalid entries: {invalid_samples}")
  return False


def get_stereoseq_input_spec_map(
    sample_list_file,
    channel,
    sample_parameters=None,
    default_input_spec=None,
):
    input_spec_map = {}
    if not os.path.isfile(sample_list_file):
        return input_spec_map
    if channel == "compare_analysis":
        sample_parameters = sample_parameters or {}
        for row in read_compare_sample_table(sample_list_file):
            sample_id = row["sample_id"]
            parameters = sample_parameters.get(sample_id, {}) or {}
            input_spec = parameters.get("input_spec", default_input_spec)
            if input_spec in [None, "", False, "None", "False", "false", "NULL", "null"]:
                sys.exit(
                    f"\nStereo-seq sample '{sample_id}' requires input_spec in YAML "
                    "(global input_spec or sample_parameters override)"
                )
            try:
                parse_stereoseq_input_spec(input_spec)
            except ValueError as exc:
                sys.exit(f"\nInvalid Stereo-seq input_spec for sample '{sample_id}': {exc}")
            input_spec_map[build_compare_sample_key(sample_id, row["group"])] = input_spec
            input_spec_map[sample_id] = input_spec
        return input_spec_map

    with open(sample_list_file) as sample_list:
        next(sample_list, None)
        for line_number, line in enumerate(sample_list, start=2):
            line = line.strip()
            if not line:
                continue
            parts = re.split(r'\s+', line)
            min_columns = 4 if channel == "compare_analysis" else 3
            if len(parts) < min_columns:
                sys.exit("\nStereo-seq integrate requires sample.txt columns: sample_id input_dir input_spec [group]")
            sample_id = parts[0].strip()
            input_spec = parts[2].strip()
            group_id = parts[3].strip() if channel == "compare_analysis" and len(parts) > 3 else None
            try:
                parse_stereoseq_input_spec(input_spec)
            except ValueError as exc:
                sys.exit(f"\nInvalid Stereo-seq input_spec at line {line_number}: {exc}")
            if channel == "compare_analysis" and group_id is not None:
                input_spec_map[build_compare_sample_key(sample_id, group_id)] = input_spec
            input_spec_map[sample_id] = input_spec
    return input_spec_map
    


def get_sample_paths(
    sample_list_file,
    type,
    channel,
    option,
    require_non_empty=False,
    data_fold=None,
    sample_parameters=None,
    default_bin_size=None,
    default_input_spec=None,
    annotation_reference_required=True,
):
    valid_samples = []
    main_file=[]
    group=[]
    bin_size=[] 
    if not require_non_empty:
        return _downstream_analysis_samples(
            sample_list_file,
            channel,
            option,
            annotation_reference_required=annotation_reference_required,
        )
    if not os.path.isfile(sample_list_file):
        if require_non_empty:
            sys.exit(
                "\n未提供必要样本信息，请提供 sample.txt 或通过命令指定样本表路径"
                f"；当前样本表路径不可用: {sample_list_file}"
            )
    if channel == "compare_analysis":
        sample_parameters = sample_parameters or {}
        rows = read_compare_sample_table(sample_list_file)
        for row in rows:
            sample = row["sample_id"]
            dir_path = row["input_path"]
            parameters = sample_parameters.get(sample, {}) or {}
            valid_samples.append(sample)
            group.append(row["group"])
            if option == "compare_stage":
                main_file.append(dir_path)
                continue

            if type == "visium_HD":
                requested_bin = parameters.get("bin_size", default_bin_size)
                if requested_bin in [None, ""]:
                    sys.exit(f"\nVisium HD sample '{sample}' requires bin_size in YAML")
                bins = "{:03d}".format(int(requested_bin))
                bin_size.append(bins)
                dir_path = os.path.join(dir_path, "binned_outputs", f"square_{bins}um")
            elif type == "stereoseq":
                requested_spec = parameters.get("input_spec", default_input_spec)
                if requested_spec in [None, ""]:
                    sys.exit(f"\nStereo-seq sample '{sample}' requires input_spec in YAML")
                parse_stereoseq_input_spec(requested_spec)
                bin_size.append(str(requested_spec))
            elif type == "visium_segment":
                segment_root = data_fold if data_fold else dir_path
                dir_path = os.path.join(segment_root, sample, "segmented_outputs")

            detected_main = check_file_exit(type, dir_path)
            if not detected_main:
                sys.exit(f"\nExiting: dependent input for sample '{sample}' was not found under {dir_path}")
            main_file.append(detected_main)
        if type in ["visium_HD", "stereoseq"]:
            return valid_samples, main_file, bin_size, group
        return valid_samples, main_file, group

    with open(sample_list_file) as sample_list:
        next(sample_list, None)
        has_data_line = False
        for line_number, line in enumerate(sample_list, start=2):
            line = line.strip()
            if not line:
                continue
            has_data_line = True
            line = re.split(r'\s+', line)
            if channel=="single_analysis" and len(line)<2:
                sys.exit(f"\n please set up the correct parameter according to your type and analysis channel")
            elif channel=="compare_analysis" and len(line)<3:
                sys.exit(f"\n please set up the correct parameter according to your type and analysis channel")
            if type=="stereoseq":
                required_columns = 4 if channel=="compare_analysis" else 3
                if len(line) < required_columns:
                    sys.exit("\nStereo-seq requires sample.txt columns: sample_id input_dir input_spec [group]")
                try:
                    parse_stereoseq_input_spec(line[2].strip())
                except ValueError as exc:
                    sys.exit(f"\nInvalid Stereo-seq input_spec at line {line_number}: {exc}")
            sample = line[0].strip()
            dir_path = line[1].strip()
            if channel=="compare_analysis":
                if type == "stereoseq":
                    if len(line) < 4:
                        sys.exit("\nStereo-seq compare_analysis requires sample.txt columns: sample_id input_dir input_spec group")
                    bin_size.append(line[2].strip())
                    group.append(line[3].strip())
                elif type!="visium_HD":
                    group_index = 3 if type=="stereoseq" else 2
                    group.append(line[group_index].strip())
                elif len(line) > 3:
                    group.append(line[3].strip())
                else:
                    sys.exit(f"\n please set up the correct parameter according to your type and analysis channel")
            if type=="visium_HD":
                if len(line) < 3:
                    sys.exit(f"\n please set up the correct parameter according to your type and analysis channel")
                bins="{:03d}".format(int(line[2].strip()))
                bin_size.append(bins)
                dir_path=os.path.join(f"{line[1].strip()}","binned_outputs",f"square_{bins}um")
            if type=="visium_segment":
                # Visium Segment follows the same project-root convention as
                # the other raw-data branches.  ``data_fold`` is authoritative
                # and each sample is stored under its sample_id.
                segment_root = data_fold if data_fold else line[1].strip()
                dir_path=os.path.join(segment_root, sample, "segmented_outputs")
            detected_main = check_file_exit(type,dir_path)
            if detected_main:
                main_file.append(detected_main)
                if len(list(set(main_file)))==1:
                    valid_samples.append(sample)
                else:
                    logger.error(f"Line {line_number} in {sample} is missing a file or directory")
                    sys.exit(f"\nExiting:  the dir or file found.")
            else:
                sys.exit(f"\nExiting:  dependent file not found.")
        if require_non_empty and (not has_data_line or len(valid_samples) == 0):
            sys.exit(f"\nsample.txt 不能为空，请提供有效的样本参数")
    sample_list.close()
    if channel=="compare_analysis":
      if type in ["visium_HD", "stereoseq"]:
          return valid_samples,main_file,bin_size,group
      else:
          return valid_samples,main_file,group
    else:
      if type=="visium_HD":
          return valid_samples,main_file,bin_size
      else:
          return valid_samples,main_file

def _downstream_analysis_samples(
    sample_list_file,
    channel,
    option,
    annotation_reference_required=True,
):
    valid_samples = []
    main_file=[]
    reference=[]
    scale_factors=[]
    if not os.path.isfile(sample_list_file):
        sys.exit(
            "\n未提供必要样本信息，请提供 sample.txt 或通过命令指定样本表路径"
            f"；当前样本表路径不可用: {sample_list_file}"
        )
    with open(sample_list_file) as sample_list:
        next(sample_list, None)
        has_data_line = False
        for line_number, line in enumerate(sample_list, start=2):
            line = line.strip()
            if not line:
                continue
            has_data_line = True
            line = re.split(r'\s+', line)
            required_columns = 3 if option == "annotation" and annotation_reference_required else 2
            if len(line) < required_columns:
                expected = (
                    "sample_id input_path sc_reference"
                    if option == "annotation" and annotation_reference_required
                    else "sample_id input_path"
                )
                sys.exit(f"\n{sample_list_file} 第 {line_number} 行缺少必要列: {expected}")
            sample_id = line[0].strip()
            input_st = line[1].strip()
            if option=="annotation":
                input_sc = line[2].strip() if len(line) > 2 else ""
                reference.append(input_sc)
            else:
                scale_factors.append(line[2].strip() if len(line) > 2 else "")
            valid_samples.append(sample_id)
            main_file.append(input_st)
            # Reclustering fans out one independent Snakemake job per zarr,
            # including when the command is launched through compare_analysis.
            if channel=="compare_analysis" and option not in ["compare_stage", "advance_analysis", "reclustering"]:
                break
    if not has_data_line or len(valid_samples) == 0:
        sys.exit(f"\n{sample_list_file} 不能为空，请提供有效的样本参数")
    sample_list.close()
    if option == "advance_analysis" and channel == "compare_analysis" and len(main_file) > 0:
        unique_inputs = list(dict.fromkeys(main_file))
        if len(unique_inputs) > 1:
            sys.exit("\ncompare_analysis 下 advance_analysis 仅支持一个整合数据路径")
        main_file = [unique_inputs[0]]
    if option=="annotation":
        return valid_samples,main_file,reference
    else:
        return valid_samples,main_file,scale_factors

def seg_filter_sample(sample_list_file, required_samples=None):
    """Read per-sample filter params from sample.txt.

    Expected columns:
      - first column: sample_id
      - last three columns: min_cells, min_counts, mt_threshold

    For standard compare_analysis sample sheets the complete layout is:
      sample_id input_dir group min_cells min_counts mt_threshold
    Visium HD and Stereo-seq keep their bin/input_spec column before group.
    """
    sample_dict = {}
    try:
        with open(sample_list_file, "r") as sample_list:
            try:
                next(sample_list)
            except StopIteration:
                raise ValueError("sample.txt is empty")
            for line_number, line in enumerate(sample_list, start=2):
                line = line.strip()
                if not line:
                    continue
                line_parts = re.split(r"\s+", line)
                if len(line_parts) < 4:
                    raise RuntimeError(
                        f"{line_number}: sample.txt requires at least 4 columns: "
                        "sample_id ... min_cells min_counts mt_threshold"
                    )
                sample_id = line_parts[0]
                if sample_id in sample_dict:
                    raise RuntimeError(f"{line_number}: duplicate sample_id in sample.txt: '{sample_id}'")
                min_cells_raw = line_parts[-3]
                min_counts_raw = line_parts[-2]
                mt_threshold_raw = line_parts[-1]
                try:
                    min_cells = int(min_cells_raw)
                    min_counts = int(min_counts_raw)
                    mt_threshold = float(mt_threshold_raw)
                except ValueError as e:
                    raise RuntimeError(
                        f"{line_number}: invalid filter values for sample '{sample_id}': "
                        f"min_cells={min_cells_raw}, min_counts={min_counts_raw}, mt_threshold={mt_threshold_raw}"
                    ) from e
                if min_cells < 1 or min_counts < 1:
                    raise RuntimeError(
                        f"{line_number}: min_cells and min_counts must both be >= 1 for sample '{sample_id}'"
                    )
                if not 0 <= mt_threshold <= 100:
                    raise RuntimeError(
                        f"{line_number}: mt_threshold must be between 0 and 100 for sample '{sample_id}'"
                    )
                sample_dict[sample_id] = [min_cells, min_counts, mt_threshold]
    except FileNotFoundError as e:
        raise FileNotFoundError(f"'{sample_list_file}' not found") from e
    if len(sample_dict) == 0:
        raise RuntimeError("sample.txt has no valid sample lines for seg_filter_sample")
    if required_samples:
        missing = [str(sample) for sample in required_samples if str(sample) not in sample_dict]
        if missing:
            raise RuntimeError(
                "sample.txt is missing per-sample filter thresholds for: " + ", ".join(missing)
            )
    return sample_dict
    

def get_annotation(file_path, samples, channel, results_folder):
    """
    Reads annotation.txt, skipping the first line.
    The second line contains comma-separated cell type annotations.
    These annotations are mapped to cluster IDs (0, 1, 2...) based on their order.
    Returns a dictionary mapping cluster IDs (as strings) to cell type annotations.
    """
    if not os.path.isfile(file_path):
        sys.exit(f"the {file_path} file not found")
    
    sample_annotations = {}
    
    # Determine which samples to process
    target_samples = ["concatenated_sdata"] if channel == "compare_analysis" else samples

    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            # Skip the first line (header)
            next(file, None)
            
            # Read the second line
            line = next(file, None)
            
            if line:
                line = line.strip()
                if line:
                    # Split by comma
                    parts = line.split(',')
                    # Create dictionary: cluster ID (str) -> annotation
                    # Order matters: first item corresponds to cluster 0, second to cluster 1, etc.
                    anno_dict = {str(i): part.strip() for i, part in enumerate(parts)}
                    
                    # Assign this annotation dictionary to all target samples
                    # Assuming the annotation file provides a single set of annotations 
                    # that applies to the current analysis context
                    for sample_name in target_samples:
                        sample_annotations[sample_name] = anno_dict
                        logger.info(f"Loaded {len(anno_dict)} annotations for sample: {sample_name}")
                else:
                    logger.warning(f"Second line in {file_path} is empty")
            else:
                logger.warning(f"{file_path} has less than 2 lines")

    except Exception as e:
        sys.exit(f"Error reading annotation file {file_path}: {e}")

    return sample_annotations
