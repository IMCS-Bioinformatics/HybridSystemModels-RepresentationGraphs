from pathlib import Path
from decrypt_juris import read_juris_file

data = Path('Lambda_Core_blue/')
files = [x for x in data.iterdir() if '.txt' in str(x).lower()]

for file_name in files:
    file_contents = read_juris_file(file_name, False)
    group_sizes = []
    for group_id in file_contents.groups:
        group_sizes.append(len(file_contents.groups[group_id]['st_list']))

    print(file_name, "node_count:", len(file_contents.states), "group_sizes:", group_sizes)