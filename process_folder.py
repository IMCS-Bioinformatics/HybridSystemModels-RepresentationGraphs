from pathlib import Path
from decrypt_juris import read_juris_file
import networkx as nx


def create_state_graph(file_contents):
    """ Creates the state graph from the Juris data holder
    """
    G = nx.Graph()
    for state_id in file_contents.states:
        state = file_contents.states[state_id]
        for line_id in state.lines:
            line = state.lines[line_id]
            to_id = line['to_state']
            G.add_edge(state_id, to_id)
    return G


def find_descendant_attractors_for_one(G, state_id):
    """ Creates the state graph from the Juris data holder
        :returns the list of descendant groups
    """
    # find descendant groups
    groups = [file_contents.states[x].group for x in nx.dfs_preorder_nodes(G, state_id) if file_contents.states[x].group is not None]
    # filter only main groups
    main_groups = [group for group in groups if int(file_contents.groups[group]["main"]) >= 1]
    main_groups = list(set(main_groups))
    return main_groups

def find_descendant_attractors(file_contents, state_id_list):
    """ Creates the map of descendant stated for each given state
    """
    G = create_state_graph(file_contents)
    return {state:find_descendant_attractors_for_one(G, state) for state in state_id_list}


def find_sources(file_contents):
    """ finds states with no incoming links
        prints descendant groups reachable from those states
    """
    in_node_set = set()
    all_node_set = set()
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
    #print(without_incoming)
    descendants = find_descendant_attractors(file_contents, without_incoming)
    print(descendants)
    #print([x for x in without_incoming if 'G' in x]) #print only groups
    return descendants

# print simple group statistics
def print_groups(file_contents):
    group_sizes = []
    for group_id in file_contents.groups:
        group_type = int(file_contents.groups[group_id]['main'])
        if group_type>0: group_sizes.append(len(file_contents.groups[group_id]['st_list']))
    print(file_name, "node_count:", len(file_contents.states), "group_sizes:", group_sizes)

########################## program start #############################

data = Path('Lambda_Core_blue/')
files = [x for x in data.iterdir() if '.txt' in str(x).lower()]

for file_name in files:
    file_contents = read_juris_file(file_name, False)
    #print_groups(file_contents)
    find_sources(file_contents)
