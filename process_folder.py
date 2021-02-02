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
    main_groups = tuple(sorted(set(main_groups)))
    return main_groups

def find_descendant_attractors(file_contents, state_id_list):
    """ Creates the map of descendant states for each given state
        group by attractor sets
    """
    G = create_state_graph(file_contents)
    state2groups = {state:find_descendant_attractors_for_one(G, state) for state in state_id_list} # this is the map from states to attractor set reachable from that state

    # group together states with the same attractors
    unique_attractors = set(state2groups.values())
    d = {}
    for n in unique_attractors:
        state_list = [k for k in state2groups.keys() if state2groups[k] == n]
        d[n] = len(state_list)# store only state count
    return d

# replace group ids to group names in the group_list
def replace_group_names(file_contents, group_list):
    value= tuple(file_contents.groups[x]['name'] for x in group_list)
    return tuple(sorted(value))

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
    # print([x for x in without_incoming if 'G' in x]) #print only groups
    descendants = find_descendant_attractors(file_contents, without_incoming)

    #replace group ids with names
    renamed = {replace_group_names(file_contents, x):descendants[x] for x in descendants}
    print(renamed)

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
#data = Path('Lambda_Complete/')
files = [x for x in data.iterdir() if '.txt' in str(x).lower()]

for file_name in files:
    file_contents = read_juris_file(file_name, False)
    #print_groups(file_contents)
    find_sources(file_contents)
