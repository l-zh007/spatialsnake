import os
import sys
import re
import pathlib
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
            'morphology_focus.ome.tif',
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
      print(main_file)
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
      print(main_file)
      return main_file
  print(f"no illigal {dir_path} {invalid_samples} file or dir")
  return False
    


def get_sample_paths(sample_list_file):
    valid_samples = []
    main_file=[]
    group=[]
    type = run_type
    print(type)
    with open(sample_list_file) as sample_list:
        next(sample_list)
        for line_number, line in enumerate(sample_list, start=2):
            line = re.split(r'\s+', line.strip())
            print(len(line),channel)
            if channel=="single_analysis" and len(line)<2:
                sys.exit(f"\n please set up the correct parameter according to your type and analysis channel")
            elif channel=="compare_analysis" and len(line)<3:
                sys.exit(f"\n please set up the correct parameter according to your type and analysis channel")
            print('correct')
            sample = line[0].strip()
            dir_path = line[1].strip()
            if channel=="compare_analysis":
              group.append(line[2].strip())
            if type=="visium_HD":
              
              if channel=="compare_analysis":
                group.append(line[3].strip())
              bin_size.append(line[2].strip())
              dir_path=os.path.join(f"{line[1].strip()}","binned_outputs",f"square_{line[2].strip()}um")
            if type=="visium_segment":
              dir_path=os.path.join(f"{line[1].strip()}","segmented_outputs")
            if check_file_exit(type,dir_path):
              main_file.append(check_file_exit(type,dir_path))
              if len(list(set(main_file)))==1:
                valid_samples.append(sample)
              else:
                print(f"!!!! line {line_number} in {sample} LOSS FILE OR DICTIONARY")
                sys.exit(f"\nExiting:  the dir or file found.")
            else:
              sys.exit(f"\nExiting:  dependent file not found.")
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
                    try:
                        val1_int = int(val1)
                        val2_int = int(val2)
                    except ValueError as e:
                        raise ValueError(f"'{val1}','{val2}' not int") from e
                    sample_dict[sample_id] = [val1_int, val2_int]
                except (IndexError, ValueError) as e:
                    raise RuntimeError(f"{line_number} : {str(e)}") from e
    except FileNotFoundError:
        raise FileNotFoundError(f"'{filter_list}' not find")
    return sample_dict    
    

def get_annotion(file_path, samples):
    sample_annotations = {}
    annotation_counts = {}
    sample=["concatenated_sdata"] if channel=="compare_analysis" else samples
    print(sample)
    for i in sample:
      target_path=os.path.join(results_folder, f"{i}", 'clustering') if channel=="single_analysis" else os.path.join(results_folder,"merge_data", 'clustering')
      numeric_folder_count = 0
      items = os.listdir(target_path)
      for item in items:
          item_full_path = os.path.join(target_path, item)
          if os.path.isdir(item_full_path) and item.isdigit():
              print(item_full_path)
              numeric_folder_count += 1
          else:
            continue
      print(numeric_folder_count,i)
      annotation_counts[i]=int(numeric_folder_count)
    with open(file_path, 'r', encoding='utf-8') as file:
        next(file)
        for line_num, line in enumerate(file, start=2):
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            if len(parts) < 2:
                print(f"警告：第{line_num}行格式不正确，跳过该行")
                continue
            sample_name = parts[0].strip()
            annotations = [part.strip() for part in parts[1:]]
            anno = dict(zip([str(i) for i in range(len(annotations))],annotations))
            sample_annotations[sample_name] = anno
            annotation_counts[sample_name] = len(annotations)
            print(len(annotations),annotation_counts[sample_name])
            if len(annotations)!=int(annotation_counts[sample_name]) or sample_name not in sample:
              sys.exit(f"\nwrong!!!:some anno of clusters are mising or the sample_name wrong with the sample_list.txt")
    print(sample_annotations)
    return sample_annotations
# def get_annotion(anno_data_path):
#     anno_dict = {}
#     numeric=[]
#     flag=True
#     for i in samples:
#       target_path=os.path.join(results_folder, f"{i}", 'clustering') if channel=="single_analysis" else os.path.join(results_folder,"merge_data", 'clustering')
#       numeric_folder_count = 0
#       items = os.listdir(target_path)
#       for item in items:
#           item_full_path = os.path.join(target_path, item)
#           if os.path.isdir(item_full_path) and item.isdigit():
#               print(item_full_path)
#               numeric_folder_count += 1
#           else:
#             continue
#       print(numeric_folder_count,i)
#       numeric.append(numeric_folder_count)
#     if not numeric==list(set(numeric)):
#       flag=False
#     try:
#         with open(anno_data_path, 'r') as sample_list:
#             try:
#                 next(sample_list)
#             except StopIteration:
#                 raise ValueError(f"invalid file")
#             for line_number, line in enumerate(sample_list, start=2):
#                 line = line.strip()
#                 if not line:
#                     continue
#                 try:
#                     line_parts = re.split(",", line)
#                     cluster = line_parts[0]
#                     annotion = line_parts[1]
#                     anno_dict[cluster] = annotion
#                 except (IndexError, ValueError) as e:
#                     raise RuntimeError(f"{line_number} : {str(e)}") from e
#     except FileNotFoundError:
#         raise FileNotFoundError(f"'{filter_list}' not find")
#     if not numeric_folder_count==int(cluster):
#         print(numeric_folder_count,cluster)
#         sys.exit(f"\nwrong!!!:some anno of clusters are mising")
#     return anno_dict
