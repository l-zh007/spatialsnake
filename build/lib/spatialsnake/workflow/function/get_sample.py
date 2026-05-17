import os
import sys
import re
import glob
from spatialsnake.workflow.function.stereoseq_spec import parse_stereoseq_input_spec

def build_compare_sample_key(sample_id, group_id):
  return f"{group_id}::{sample_id}"

def check_file_exit(type,dir_path):
  invalid_samples=[]
  if type == "stereoseq":
    if os.path.isfile(dir_path):
      lower_name = dir_path.lower()
      if lower_name.endswith(".cellbin.gef") or lower_name.endswith(".gef") or ".gem" in lower_name:
        return "stereoseq"
      print(f"Stereo-seq input file not supported: {dir_path}")
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
        print(f"Stereo-seq feature_expression files not found under {feature_dir}")
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
      print(f"Stereo-seq key files not found under {dir_path}. Expect *.cellbin.gef, *.gef, *.gem(.gz), or tif images")
      return False
  if os.path.isdir(dir_path):
    if type == "Merfish":
      has_cell_by_gene = len(glob.glob(os.path.join(dir_path, "**", "cell_by_gene.csv"), recursive=True)) > 0
      has_transcripts = len(glob.glob(os.path.join(dir_path, "**", "*transcripts*.csv*"), recursive=True)) > 0
      has_transcripts_parquet = len(glob.glob(os.path.join(dir_path, "**", "*transcripts*.parquet"), recursive=True)) > 0
      if not (has_cell_by_gene or has_transcripts or has_transcripts_parquet):
        print(f"MERSCOPE/MERFISH key files not found under {dir_path}. Expect cell_by_gene.csv and/or detected_transcripts*.csv/parquet")
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
            'experiment.xenium'],
        'slide_seq': [
            'BeadLocationsForR.csv',
            'MappedDGEForR.csv']}
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
      "raw_feature_cell_matrix.h5",
      "MappedDGEForR.csv"
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
  print(f"no illigal {dir_path} {invalid_samples} file or dir")
  return False


def get_stereoseq_input_spec_map(sample_list_file, channel):
    input_spec_map = {}
    if not os.path.isfile(sample_list_file):
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
    


def get_sample_paths(sample_list_file, type,channel,option,require_non_empty=False):
    valid_samples = []
    main_file=[]
    group=[]
    bin_size=[] 
    if not require_non_empty:
        return _downstream_analysis_samples(sample_list_file, channel , option)
    if not os.path.isfile(sample_list_file):
        if require_non_empty:
            sys.exit(f"\nsample.txt is not a file")
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
                dir_path=os.path.join(f"{line[1].strip()}","segmented_outputs")
            detected_main = check_file_exit(type,dir_path)
            if detected_main:
                main_file.append(detected_main)
                if len(list(set(main_file)))==1:
                    valid_samples.append(sample)
                else:
                    print(f"!!!! line {line_number} in {sample} LOSS FILE OR DICTIONARY")
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

def _downstream_analysis_samples(sample_list_file, channel , option):
    valid_samples = []
    main_file=[]
    reference=[]
    scale_factors=[]
    if not os.path.isfile(sample_list_file):
        L.info(f"sample.txt 中没有样本id 随机设置样本输出ID custom_project")
        valid_samples.append("custom_project")
    with open(sample_list_file) as sample_list:
        next(sample_list, None)
        for line_number, line in enumerate(sample_list, start=2):
            line = line.strip()
            if not line:
                continue
            line = re.split(r'\s+', line)
            sample_id = line[0].strip()
            input_st = line[1].strip()
            if option=="annotation":
                input_sc = line[-1].strip()
                reference.append(input_sc)
            else:
                scale_factors.append(line[2].strip() if len(line) > 2 else "")
            valid_samples.append(sample_id)
            main_file.append(input_st)
            if channel=="compare_analysis" and option not in ["compare_stage", "advance_analysis"]:
                break
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

def seg_filter_sample(sample_list_file):
    """Read per-sample filter params from sample.txt.

    Expected columns:
      - first column: sample_id
      - last three columns: min_cells, min_genes, mt_threshold
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
                        "sample_id ... min_cells min_genes mt_threshold"
                    )
                sample_id = line_parts[0]
                min_cells_raw = line_parts[-3]
                min_genes_raw = line_parts[-2]
                mt_threshold_raw = line_parts[-1]
                try:
                    min_cells = float(min_cells_raw)
                    min_genes = float(min_genes_raw)
                    mt_threshold = float(mt_threshold_raw)
                except ValueError as e:
                    raise RuntimeError(
                        f"{line_number}: invalid filter values for sample '{sample_id}': "
                        f"min_cells={min_cells_raw}, min_genes={min_genes_raw}, mt_threshold={mt_threshold_raw}"
                    ) from e
                sample_dict[sample_id] = [min_cells, min_genes, mt_threshold]
    except FileNotFoundError as e:
        raise FileNotFoundError(f"'{sample_list_file}' not find") from e
    if len(sample_dict) == 0:
        raise RuntimeError("sample.txt has no valid sample lines for seg_filter_sample")
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
                        print(f"Loaded {len(anno_dict)} annotations for sample: {sample_name}")
                else:
                    print(f"Warning: Second line in {file_path} is empty.")
            else:
                print(f"Warning: {file_path} has less than 2 lines.")

    except Exception as e:
        sys.exit(f"Error reading annotation file {file_path}: {e}")

    return sample_annotations
