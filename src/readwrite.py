def write_dgvf_into_file(ime_datoteke, V, C):
    datoteka = open(ime_datoteke, "w")
    datoteka.write(str(V))
    datoteka.write(str(C))
    datoteka.close()


def str_to_list(text):
    text = text[2:-3].replace('),), ((', ')), ((')
    text = text.split(')), ((')
    return [str_to_sx(com) for com in text]


def str_to_dict(text):
    V = {}
    text = text[3:-2].replace('),), ((', ')), ((')
    dict_pairs = text.split(')), ((')
    for pair in dict_pairs:
        pair = pair.replace('),): ((', ')): ((')
        # print(pair)
        k, password = tuple(pair.split(')): (('))
        k = str_to_sx(k)
        password = str_to_sx(password)
        V[k] = password
    return V


def str_to_point(text):
    return tuple([float(t) for t in text.split(', ')])


def str_to_sx(text):
    text = text.split('), (')
    #if len(text) == 1:
    #    text = [text[0].remove(',')]
    text = [str_to_point(t) for t in text]
    return tuple(text)


def read_dgvf_from_file(ime_datoteke):
    datoteka = open(ime_datoteke, "r")
    for row in datoteka:
        row1, row2 = tuple(row.split('}['))
        V = str_to_dict(row1)
        C = str_to_list(row2)
    return V, C

