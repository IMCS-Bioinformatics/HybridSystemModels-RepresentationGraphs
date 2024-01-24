# This program reads a state space file
# version 0.9

class SpaceDataHolder:
    def __init__(self) -> None:
        self.genes = {}
        self.bs = {'bin': 0, 'tern': 0}
        self.orders = {}
        self.vn = {}
        self.groups = {}
        self.states = {}
        super().__init__()

class SpaceState:
    def __init__(self, st_id = None, genes = None, name=None, from0=0, group=None) -> None:
        self.id = st_id
        self.genes = genes
        self.name = name
        self.from0 = from0
        self.group = group
        self.lines = {}
        super().__init__()


def split_string_with_delimiter(line, delimiter):
    return line.split(delimiter)

def toNumbers(number_list):
    return [int(x) for x in number_list]

def remove_empty_lines(data):
    return [line.rstrip() for line in data if len(line.rstrip())>0]

def read_space_file(filename, only_0_reachable):
    paz = 1 if only_0_reachable else 0
    def get_ordering(info, data):
        cI_poz = []
        cro_poz = []
        if info.bs["tern"] == 6:
            cI_poz = [x-1 for x in [1, 3, 5, 7, 9, 11]]
            cro_poz = [x-1 for x in [25, 27, 29, 31, 33, 35]]

        if info.bs["tern"] == 3:
            cI_poz = [x-1 for x in [1, 3, 5]]
            cro_poz = [x-1 for x in [13, 15, 17]]

        orders = toNumbers(split_string_with_delimiter(data[1], " "))
        info.orders['all'] = orders
        cI_table = []
        for o in cI_poz:
            cI_table.append(orders[o])

        info.orders['cI_order'] = toNumbers(cI_table)
        cro_table = []
        for o in cro_poz:
            cro_table.append(orders[o])

        info.orders['cro_order'] = toNumbers(cro_table)

    def nosaukums(data, bs):
        sep = []
        n = 0
        if bs['bin'] == 1 and bs['tern'] == 6:
            sep = ["-", "", ".", "", ".", "", "-", "", ".", "", ".", "", ""]
            n = 13

        if bs['bin'] == 1 and bs['tern'] == 3:
            sep = ["-", "", ".", "", ".", "", ""]
            n = 7

        if bs['bin'] == 4 and bs['tern'] == 6:
            sep = ["", "", "", "-", "", ".", "", ".", "", "-", "", ".", "", ".", "", ""]
            n = 16

        if bs['bin'] == 5 and bs['tern'] == 6:
            sep = ["", "", "", "", "-", "", ".", "", ".", "", "-", "", ".", "", ".", "", ""]
            n = 17

        data_list = split_string_with_delimiter(data, " ")
        rez = []

        for i in range(n):
            rez.append(data_list[i])
            rez.append(sep[i])

        return "".join(rez)

    with open(filename,"r") as f:
        data = f.readlines()

    data = remove_empty_lines(data)

    info = SpaceDataHolder()
    iii = 2
    g_count = int(data[iii])
    plus = 2 * g_count + 1
    info.bs['bin'] = int(split_string_with_delimiter(data[iii + 2], " ")[0])
    info.bs['tern'] = int(split_string_with_delimiter(data[iii + 2], " ")[1])
    info.vn['from0'] = int(split_string_with_delimiter(data[iii + 5+plus], " ")[0])
    info.vn['count'] = int(split_string_with_delimiter(data[iii + 5+plus], " ")[1])
    info.vn['all'] = int(split_string_with_delimiter(data[iii + 5+plus], " ")[2])
    get_ordering(info, data)
    ii = 0
    for g in split_string_with_delimiter(data[iii + 1], " "):
        info.genes[ii] = g
        ii += 1

    iii = 3+plus
    gr_count = int(data[iii + 5])
    gr_info = split_string_with_delimiter(data[iii + 6], " ")
    ii = 0
    for i in range(gr_count):
        gid = "G" + str(gr_info[ii])
        name = gr_info[ii + 2]+ "("+gr_info[ii + 1]+")" + "_"+gid
        if gr_info[ii + 2] == "1":
            info.groups[gid] = {'id': gid, 'st_count': int(gr_info[ii + 1]), 'main': gr_info[ii + 2],
                                                   'st_list': split_string_with_delimiter(data[iii + 7 + i], " "), 'name':name}
        elif paz == 0:
            info.groups[gid] = {'id': gid, 'st_count': int(gr_info[ii + 1]), 'main': gr_info[ii + 2],
                                                   'st_list': split_string_with_delimiter(data[iii + 7 + i], " "), 'name':name}

        ii += 3

    ii = gr_count + 10+plus
    for i in range(info.vn['count']):
        st_id = data[ii]
        g_list = split_string_with_delimiter(data[ii + 2], " ")
        g_string = ""
        G = {}
        for j in range(g_count):
            t = info.genes[j]
            if g_list[j] == "1":
                G[t] = "+"
            else:
                G[t] = "-"

        from0 = 0
        if int(split_string_with_delimiter(data[ii + 1], " ")[1]) == 0:
            from0 = 1

        stG = "G" + (split_string_with_delimiter(data[ii + 1], " ")[0])
        if not stG in info.groups:
            stG = None
        elif info.groups[stG]['main'] == "0":
            stG = None

        if paz == 0 or (paz == 1 and from0 == 1):
            info.states[st_id] = SpaceState(st_id=st_id, genes=G, name = nosaukums(data[ii + 3], info.bs), from0=from0, group=stG)
            #info.states[st_id] = {'id': st_id, 'genes': G, 'name': nosaukums(data[ii + 3], info.bs), 'from0': from0, 'group': stG, 'lines': {}}
            l_sk = int(data[ii + 4])
            l_info = split_string_with_delimiter(data[ii + 5], " ")
            for j in range(l_sk):
                l_id = st_id + "." + l_info[j * 2]
                gene_id = l_info[j * 2+1]
                gene = info.genes[int(gene_id)]
                info.states[st_id].lines[l_id] = {'id': l_id, 'to_state': l_info[j * 2], 'gene': gene, 'from_state':st_id}
        ii += 6

    info.genes = split_string_with_delimiter(data[3], " ")
    return info
