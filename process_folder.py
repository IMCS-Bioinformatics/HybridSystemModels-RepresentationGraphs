from pathlib import Path
from decrypt_juris import read_juris_file
import networkx as nx


def create_state_graph(file_contents, reverse=False):
    """ Creates the state graph from the Juris data holder
    """
    G = nx.DiGraph()
    for state_id in file_contents.states:
        state = file_contents.states[state_id]
        for line_id in state.lines:
            line = state.lines[line_id]
            to_id = line['to_state']
            if reverse:
                G.add_edge(to_id, state_id)
            else:
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

############### find decision nodes ##################
def analyze_decision_nodes(file_contents):
    G_forward = create_state_graph(file_contents)
    G = create_state_graph(file_contents, reverse=True)# find reverse graph
    pred_list = []

    for group_id in file_contents.groups:
        group_type = int(file_contents.groups[group_id]['main'])
        if group_type>0:
            state_list = file_contents.groups[group_id]['st_list']
            repr_state = state_list[0] # assume that group is strongly connected so we can take any its node for finding predecessors
            predecessor_nodes = nx.dfs_preorder_nodes(G, repr_state)
            pred_list.append(list(predecessor_nodes))

    if len(pred_list)!=2: raise Exception("there should be 2 attractors")
    nodes1 = set(pred_list[0])
    nodes2 = set(pred_list[1])
    common = nodes1.intersection(nodes2)
    #decision nodes are those from which both attractors are reachable but at least for one descendant only one attractor is reachable
    decision_nodes = []
    for node in common:
        for child in G_forward.neighbors(node):
            if child not in nodes1 or child not in nodes2:
                decision_nodes.append(node)
                break

    print(len(file_contents.states), len(decision_nodes))

# print simple group statistics
def print_groups(file_contents):
    group_sizes = []
    group_genes = []
    for group_id in file_contents.groups:
        group_type = int(file_contents.groups[group_id]['main'])
        if group_type>0:
            genes = set()
            for state in file_contents.groups[group_id]['st_list']:
                lines = file_contents.states[state].lines
                for key in lines:
                    gene = lines[key]['gene']
                    genes.add(gene)
            #print(genes)
            group_genes.append(sorted(list(genes)))
            group_sizes.append(len(file_contents.groups[group_id]['st_list']))
    #print(file_name, "node_count:", len(file_contents.states), "group_sizes:", group_sizes)
    #print(file_name, "node_count:", len(file_contents.states), "group_genes:", group_genes)
    print(group_genes)

########################## program start #############################

#data = Path('Lambda_Core_blue/')
#data = Path('Lambda_Complete/')
#data = Path('HK22_Complete/')
data = Path('Lambda_Oppenheim/')
files = [x for x in data.iterdir() if '.txt' in str(x).lower()]

for file_name in files:
    file_contents = read_juris_file(file_name, False)
    print_groups(file_contents)
    #find_sources(file_contents)
    #analyze_decision_nodes(file_contents)
    #break
