import os
import sys
import re
import pathlib
L = logging.getLogger("spatialsnak")
L.setLevel(logging.INFO)
L.propagate = False
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
if not L.handlers:
    L.addHandler(log_handler)

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
  L.info(f"no illigal {dir_path} {invalid_samples} file or dir")
  return False
    


def get_sample_paths(sample_list_file):
    valid_samples = []
    main_file=[]
    group=[]
    type = run_type
    with open(sample_list_file) as sample_list:
        next(sample_list)
        for line_number, line in enumerate(sample_list, start=2):
            line = re.split(r'\s+', line.strip())
            if channel=="single_analysis" and len(line)<2:
                L.info(f"\n please set up the correct parameter according to your type and analysis channel")
                sys.exit()
            elif channel=="compare_analysis" and len(line)<3:
                L.info(f"\n please set up the correct parameter according to your type and analysis channel")
                sys.exit()
            sample = line[0].strip()
            dir_path = line[1].strip()
            if channel=="compare_analysis":
                group.append(line[2].strip()) if type!="visium_HD" else group.append(line[3].strip())
            if type=="visium_HD":
              bins="{:03d}".format(int(line[2].strip()))
              bin_size.append(bins)
              dir_path=os.path.join(f"{line[1].strip()}","binned_outputs",f"square_{bins}um")
            if type=="visium_segment":
              dir_path=os.path.join(f"{line[1].strip()}","segmented_outputs")
            if check_file_exit(type,dir_path):
              main_file.append(check_file_exit(type,dir_path))
              if len(list(set(main_file)))==1:
                valid_samples.append(sample)
              else:
                L.info(f"!!!! line {line_number} in {sample} LOSS FILE OR DICTIONARY")
                sys.exit(f"\nExiting:  the dir or file found.")
            else:
              L.info(f"\ndependent file not found.")
              sys.exit()
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
    return sample_dict    
    

def get_annotion(file_path, samples):
    if not os.path.isfile(file_path):
      L.info(f"the {file_path} file not found")
      sys.exit()
    sample_annotations = {}
    annotation_counts = {}
    sample=["concatenated_sdata"] if channel=="compare_analysis" else samples
    for i in sample:
      target_path=os.path.join(results_folder, f"{i}", 'clustering') if channel=="single_analysis" else os.path.join(results_folder,"merge_data", 'clustering')
      numeric_folder_count = 0
      items = os.listdir(target_path)
      for item in items:
          item_full_path = os.path.join(target_path, item)
          if os.path.isdir(item_full_path) and item.isdigit():
              numeric_folder_count += 1
          else:
            continue
      annotation_counts[i]=int(numeric_folder_count)
    with open(file_path, 'r', encoding='utf-8') as file:
        next(file)
        for line_num, line in enumerate(file, start=2):
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            if len(parts) < 2:
                L.info(f"please check the line:{line_num} confirm it ligally")
                continue
            sample_name = parts[0].strip()
            annotations = [part.strip() for part in parts[1:]]
            anno = dict(zip([str(i) for i in range(len(annotations))],annotations))
            sample_annotations[sample_name] = anno
            annotation_counts[sample_name] = len(annotations)
            if len(annotations)!=int(annotation_counts[sample_name]) or sample_name not in sample:
              L.info(f"\nwrong!!!:some anno of clusters are mising or the sample_name wrong with the sample_list.txt")
              sys.exit()
    return sample_annotations

