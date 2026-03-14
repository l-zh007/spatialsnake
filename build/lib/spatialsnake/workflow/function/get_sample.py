import os
import sys
import re

def check_file_exit(type,dir_path):
  invalid_samples=[]
  if os.path.isdir(dir_path):
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
    main_file=""
    if os.path.isfile(os.path.join(dir_path,'filtered_feature_bc_matrix.h5')):
      main_file='filtered_feature_bc_matrix.h5'
      return main_file
    elif os.path.isfile(os.path.join(dir_path,'raw_feature_bc_matrix.h5')):
      main_file='raw_feature_bc_matrix.h5'
      return main_file
    elif os.path.isfile(os.path.join(dir_path,"cell_feature_matrix.h5")):
      main_file="cell_feature_matrix.h5"
      return main_file
    elif os.path.isfile(os.path.join(dir_path,'filtered_feature_cell_matrix.h5')):
      main_file="filtered_feature_cell_matrix.h5"
      return main_file
    elif os.path.isfile(os.path.join(dir_path,'raw_feature_cell_matrix.h5')):
      main_file="raw_feature_cell_matrix.h5"
      return main_file
    elif os.path.isfile(os.path.join(dir_path,"MappedDGEForR.csv")):
      main_file="MappedDGEForR.csv"
      return main_file
  print(f"no illigal {dir_path} {invalid_samples} file or dir")
  return False
    


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
            sample = line[0].strip()
            dir_path = line[1].strip()
            if channel=="compare_analysis":
                if type!="visium_HD":
                    group.append(line[2].strip())
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
      if type=="visium_HD":
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
            if option=="annotion":
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
    if option=="annotion":
        return valid_samples,main_file,reference
    else:
        return valid_samples,main_file,scale_factors

def seg_filter_sample(filter_list):
    type = run_type
    sample_dict = {}
    try:
        with open(filter_list, 'r') as sample_list:
            try:
                next(sample_list)
            except StopIteration:
                raise ValueError(f"invalid file")
            for line_number, line in enumerate(sample_list, start=2):
                line = line.strip()
                if not line:
                    continue
                try:
                    line_parts = re.split(r'\s+', line)
                    if len(line_parts) < 3:
                        raise IndexError(f"without necessary params,just {len(line_parts)}")
                    sample_id = line_parts[0]
                    val1 = line_parts[1]
                    val2 = line_parts[2]
                    val3 = line_parts[3]
                    try:
                        val1_int = float(val1)
                        val2_int = float(val2)
                        val3_int = float(val3)
                    except ValueError as e:
                        raise ValueError(f"'{val1}','{val2}' not int") from e
                    sample_dict[sample_id] = [val1_int, val2_int,val3_int]
                except (IndexError, ValueError) as e:
                    raise RuntimeError(f"{line_number} : {str(e)}") from e
    except FileNotFoundError:
        raise FileNotFoundError(f"'{filter_list}' not find")
    sample_list.close()
    return sample_dict    
    

def get_annotion(file_path, samples, channel, results_folder):
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
