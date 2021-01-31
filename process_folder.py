from pathlib import Path
from decrypt_juris import read_juris_file

# find states with no incoming links
def find_sources(file_contents):
    in_node_set = set()
    all_node_set = set()
    #print(file_contents.states)
    for state_id in file_contents.states:
        state = file_contents.states[state_id]
        group = state.group if state.group is not None else state_id
        all_node_set.add(group)
        for line_id in state.lines:
            line = state.lines[line_id]
            to_id = line['to_state']
            to_group = file_contents.states[to_id].group
            if to_group is None: to_group = to_id
            if to_group != group: in_node_set.add(to_group)

    without_incoming = all_node_set.difference(in_node_set)
    without_incoming = sorted(list(without_incoming))
    print(without_incoming)
    #print([file_contents.states[x].name for x in without_incoming])
    #print([x for x in without_incoming if 'G' in x]) #print only groups

def print_groups(file_contents):
    group_sizes = []
    for group_id in file_contents.groups:
        #print(file_contents.groups[group_id])
        group_type = int(file_contents.groups[group_id]['main'])
        if group_type>0: group_sizes.append(len(file_contents.groups[group_id]['st_list']))
    print(file_name, "node_count:", len(file_contents.states), "group_sizes:", group_sizes)

data = Path('Lambda_Core_blue/')
files = [x for x in data.iterdir() if '.txt' in str(x).lower()]

for file_name in files:
    file_contents = read_juris_file(file_name, False)
    #print_groups(file_contents)
    find_sources(file_contents)
