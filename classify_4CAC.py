
import argparse
import numpy as np

def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'classify contigs into categories'
    )
    parser.add_argument('--assembler',
                        help='contigs to classify were assembled from SPAdes (metaSPAdes) or Flye (metaFlye)',
                        required=True, type=str
                        )
    parser.add_argument('--asmdir',
                        help='path for assembly graph file, path file: scaffolds.fasta (SPAdes) and assembly_info.txt (Flye)',
                        required=True, type=str
                        )
    parser.add_argument('--classdir',
                        help='path for the XBGoost 4-class classification result',
                        required=False, type=str
                        )
    parser.add_argument('--outdir',
                        help='output path for the 4CAC classification result',
                        required=False, type=str
                        )
    parser.add_argument('--phage',
                        help='phage score threshold',
                        required=False, type=float, default=0.95
                        )
    parser.add_argument('--plasmid',
                        help='plasmid score threshold',
                        required=False, type=float, default=0.95
                        )
    return parser.parse_args()


def getEdgeID(edge_label, flyeAssembler):
    if flyeAssembler:
        return int(edge_label[5:])
    else:
        return int(edge_label)

def readFlyePathFile(file_path, contigs, contigs_length, contigs_nodeList):
    headline = True
    contig_count = 0
    with open (file_path) as file:
        for row in file:
            if headline:
                headline = False
            else:
                row_array = row.split()
                path = row_array[7].split(',')
                node_list = []
                for node in path:
                    if node != '*' and node != '??':
                        node_ID = abs(int(node))
                        node_list.append(node_ID)

                contigs[row_array[0]] = contig_count
                contigs_length.append(int(row_array[1]))
                contigs_nodeList.append(node_list)
                contig_count += 1
    file.close()
    print(str(len(contigs))+" contigs contained in the path file: "+file_path)

def readSPAdesPathFile(file_path, contigs, contigs_length, contigs_nodeList):
    node_list = []
    readNode = False
    contig_count = 0

    with open(file_path) as file:
        for row in file:
            row = row[0:-1]
            if len(row) > 4 and row[0:4] == 'NODE':
                if row[-1] == '\'': # each node has both forward and reverse version, we read only one of them
                    readNode = False
                else:
                    readNode = True
                    if len(node_list) > 0:
                        cur_list = node_list[:]
                        contigs_nodeList.append(cur_list)
                        node_list.clear()

                    contigs[row] = contig_count
                    node_label = row.split('_')
                    contigs_length.append(int(node_label[3]))
                    contig_count += 1
            else:
                if readNode:
                    if row[-1] == ';':
                        row = row[0:-1]
                    path = row.split(',')
                    for node in path:
                        node_ID = int(node[0:-1])
                        node_list.append(node_ID)

    file.close()
    if len(node_list) > 0:
        cur_list = node_list[:]
        contigs_nodeList.append(cur_list)
        node_list.clear()

    print(str(len(contigs)) + " contigs contained in the path file: "+file_path)


def contigCategoryToNodes(file_xgb, phage_score, plas_score, contigs_map, contigs_length, contigs_nodeList, nodes_map):
    nodes_cate = [0] * len(nodes_map)
    initial_class = [0] * len(contigs_map)

    with open(file_xgb) as file:
        for row in file:
            row_array = row.split()
            contig_index = contigs_map[row_array[0]]
            contig_cate = 0
            initial_cate = 0
            if float(row_array[1]) >= phage_score:
                contig_cate = 1
                initial_cate = 1
            elif float(row_array[2]) >= plas_score:
                contig_cate = 2
                initial_cate = 2
            elif float(row_array[4]) >= float(row_array[3]) and float(row_array[4]) > float(row_array[1]) and float(row_array[4]) > float(row_array[2]):
                contig_cate = 4
                initial_cate = 4
            elif float(row_array[3]) >= float(row_array[4]) and float(row_array[3]) > float(row_array[1]) and float(row_array[3]) > float(row_array[2]):
                contig_cate = 3
                initial_cate = 3
            elif float(row_array[3]) >= 0.1:
                initial_cate = 3

            initial_class[contig_index] = initial_cate
            if contig_cate != 0 and contigs_length[contig_index] >= 1000: # if this is a long contig and it has been classified
                for node in contigs_nodeList[contig_index]: # path of nodes for this contig
                    node_ID = nodes_map[node]
                if nodes_cate[node_ID] == 0:      # if the edge isn't classified yet
                    nodes_cate[node_ID] = contig_cate
                elif nodes_cate[node_ID] != contig_cate:   # edge with conflict classification as it was contained in multiple contigs with different classes
                    nodes_cate[node_ID] = -1
    file.close()

    for index in range(len(nodes_cate)):
        if nodes_cate[index] == -1:  # set conflict edges to be unclassified
            nodes_cate[index] = 0

    print("edges contained in contig paths are assigned with contig classification! ")
    return nodes_cate, initial_class


def readAssemblyGraph(file_graph, flyeAssembler):
    nodes_map = {}
    nodes_length = []
    nodes_adj = []
    node_count = 0
    start_edge = False

    with open(file_graph) as file:
        for row in file:
            row_array = row.split()
            if row_array[0] == 'S':
                edge_ID = getEdgeID(row_array[1], flyeAssembler)
                nodes_map[edge_ID] = node_count
                nodes_length.append(len(row_array[2]))
                node_count += 1

            elif row_array[0] == 'L':
                if start_edge == False:
                    start_edge = True
                    nodes_adj = [[] for _ in range(node_count)] # initialize the adjacent list for nodes

                if row_array[1] != row_array[3]:
                    preEdge_ID = getEdgeID(row_array[1], flyeAssembler)
                    surEdge_ID = getEdgeID(row_array[3], flyeAssembler)

                    preEdge_ID = nodes_map[preEdge_ID]
                    surEdge_ID = nodes_map[surEdge_ID]
                    if surEdge_ID not in nodes_adj[preEdge_ID]:
                        nodes_adj[preEdge_ID].append(surEdge_ID)
                        nodes_adj[surEdge_ID].append(preEdge_ID)

            elif row_array[0] == 'P':
                break
    file.close()
    print("complete reading assembly graph! node size = %.3f" % len(nodes_map))
    return nodes_map, nodes_length, nodes_adj

def correctionMaxAdj(nodes_adj, nodes_cate):
    correct_list = np.array([[0, 0, 0]])
    correct_list = np.delete(correct_list, 0, axis=0)

    for index in range(len(nodes_adj)):
        if nodes_cate[index] != 0 and len(nodes_adj[index]) >= 2: # read a classified edge with at least two adjacent neighbors
            multi_adjClass = False
            support_num = 0
            correct_cate = 0

            for adjNode in nodes_adj[index]:
                adjCate = nodes_cate[adjNode]
                if adjCate != 0:   # find a classified neighbor
                    if correct_cate == 0:  # the first classified neighbor
                        correct_cate = adjCate
                        support_num = 1
                    elif correct_cate == adjCate: # classified neighbors with same class
                        support_num += 1
                    else:    # classified neighbors with different classes
                        multi_adjClass = True

                if multi_adjClass:
                    break

            if correct_cate != nodes_cate[index] and multi_adjClass == False and support_num >=2:
                correct_list = binaryInsertNP(correct_list, support_num, index, correct_cate)


    if len(correct_list) > 0:
        print('%.3f classified edges has classified neighbors with differenct classes, the maximum degree is: %.3f' % (len(correct_list), correct_list[-1][0]))
    else:
        print('0 edge need to be corrected! ')


    correctEdges = 0
    num_multiSup = 0
    while len(correct_list) > 0:
        curNode = correct_list[-1][1]
        nodes_cate[curNode] = correct_list[-1][2]  # correct one node
        correctEdges += 1
        if correct_list[-1][0] > 2:
            num_multiSup += 1
        correct_list = np.delete(correct_list, -1, axis=0)

        for adjNode in nodes_adj[curNode]:
            if nodes_cate[adjNode] != 0:
                pos = np.where(correct_list[:,1] == adjNode)
                if len(pos[0]) == 1:
                    index = pos[0][0]
                    correct_list = np.delete(correct_list, index, axis=0)

    print('Complete the correction step. %.3f classified edges were corrected, and %.3f with >2 supporting neighbors!' % (correctEdges, num_multiSup))


def propagationMaxAdj(nodes_adj, nodes_cate):
    propagate_list = np.array([[0, 0, 0]])
    propagate_list = np.delete(propagate_list, 0, axis=0)

    for index in range(len(nodes_adj)):
        if nodes_cate[index] == 0 and len(nodes_adj[index]) >= 1: # read an unclassified edge with at least one adjacent neighbors
            multi_adjClass = False
            support_num = 0
            propagate_cate = 0

            for adjNode in nodes_adj[index]:
                adjCate = nodes_cate[adjNode]
                if adjCate != 0:   # find a classified neighbor
                    if propagate_cate == 0:  # the first classified neighbor
                        propagate_cate = adjCate
                        support_num = 1
                    elif propagate_cate == adjCate: # classified neighbors with same class
                        support_num += 1
                    else:    # classified neighbors with different classes
                        multi_adjClass = True
                if multi_adjClass:
                    break

            if multi_adjClass == False and propagate_cate != 0:
                propagate_list = binaryInsertNP(propagate_list, support_num, index, propagate_cate)


    if len(propagate_list) > 0:
        print('%.3f unclassified edges has classified neighbors, the maximum degree is: %.3f' % (len(propagate_list), propagate_list[-1][0]))
    else:
        print('0 edge need to be propagated! ')


    propagateEdges = 0
    num_multiSup = 0
    while len(propagate_list) > 0:
        curNode = propagate_list[-1][1]
        nodes_cate[curNode] = propagate_list[-1][2]
        propagateEdges += 1
        if propagate_list[-1][0] > 1:
            num_multiSup += 1
        propagate_list = np.delete(propagate_list, -1, axis=0)

        for adjNode in nodes_adj[curNode]:
            if nodes_cate[adjNode] == 0: # update unclassified neighbors
                pos = np.where(propagate_list[:,1] == adjNode)
                if len(pos[0]) == 1:
                    index = pos[0][0]
                    if nodes_cate[curNode] == propagate_list[index][2]:
                        support_num = propagate_list[index][0]+1
                        propagate_list = np.delete(propagate_list, index, axis=0)
                        propagate_list = binaryInsertNP(propagate_list, support_num, adjNode, nodes_cate[curNode])
                    else:
                        propagate_list = np.delete(propagate_list, index, axis=0)
                else:
                    support_num = 0
                    for node in nodes_adj[adjNode]:
                        if nodes_cate[node] != 0:
                            support_num += 1
                        if support_num > 1:
                            break

                    if support_num == 1:
                        propagate_list = binaryInsertNP(propagate_list, support_num, adjNode, nodes_cate[curNode])

    print('Complete the propagation step. %.3f uncertain edges were classified, and %.3f with >1 supporting neighbors!' % (propagateEdges, num_multiSup))


def binaryInsertNP(list, sup, edge, cate):
    i=list[:,0].searchsorted(sup, side='left')
    return np.vstack((list[:i,], np.array([sup, edge, cate]), list[i:,]))


def finalContigClassifyFromEdges(file_out, contigs_map, contigs_nodeList, nodes_map, nodes_length, nodes_cate, initial_class):
    with open(file_out, 'w') as o:
        for contig in contigs_map:
            cate_length = [0] * 5
            maxLength = 0
            sumLength = 0
            finalCate = -1
            contig_ID = contigs_map[contig]

            for node in contigs_nodeList[contig_ID]:
                node_ID = nodes_map[node]
                nodeCate = nodes_cate[node_ID]
                cate_length[nodeCate] += nodes_length[node_ID]
                sumLength += nodes_length[node_ID]

                if maxLength < cate_length[nodeCate]:
                    maxLength = cate_length[nodeCate]
                    finalCate = nodeCate

            if maxLength/float(sumLength) < 0.8: # if the class with longest length doesn't account for 80% of the contig length, this contig will be classified as uncertain
                finalCate = 0

            if finalCate == 0:
                finalCate = initial_class[contig_ID]

            contigCate = 'uncertain'
            if finalCate == 1:
                contigCate = 'phage'
            elif finalCate == 2:
                contigCate = 'plasmid'
            elif finalCate == 3:
                contigCate = 'prokarya'
            elif finalCate == 4:
                contigCate = 'eukarya'
            o.write(contig+","+contigCate+"\n")
    o.close()
    print('The final classification of 4CAC could be found in: '+file_out)


def main(args):

    if args.classdir:
        classdir = args.classdir
    else:
        classdir = args.asmdir
    if args.outdir:
        outdir = args.outdir
    else:
        outdir = args.asmdir

    # input files of Flye assembly
    file_graph = args.asmdir+"assembly_graph.gfa"
    file_path = args.asmdir+"assembly_info.txt"
    file_xgb = classdir+"assembly.fasta.probs_xgb_4class.out"
    flyeAssembler = True

    # input files of SPAdes assembly
    if args.assembler.upper() == 'SPADES' or args.assembler.upper() == 'METASPADES':
        file_graph = args.asmdir+"assembly_graph_with_scaffolds.gfa"
        file_path = args.asmdir+"scaffolds.paths"
        file_xgb = classdir+"scaffolds.fasta.probs_xgb_4class.out"
        flyeAssembler = False


    # read nodes and edges in the assembly graph
    nodes_map, nodes_length, nodes_adj = readAssemblyGraph(file_graph, flyeAssembler)


    # each contig corresponds to a path consisting of one or more edges in the graph.
    contigs_map = {}
    contigs_length = []
    contigs_nodeList = []
    if flyeAssembler:
        readFlyePathFile(file_path, contigs_map, contigs_length, contigs_nodeList)
    else:
        readSPAdesPathFile(file_path, contigs_map, contigs_length, contigs_nodeList)


    # initial classification of contigs using XGBoost 4-class classifier
    # assign contigs' classification to its nodes in the assembly graph
    score_phage = args.phage
    score_plasmid = args.plasmid
    nodes_cate, initial_class = contigCategoryToNodes(file_xgb, score_phage, score_plasmid, contigs_map, contigs_length, contigs_nodeList, nodes_map)


    # correction and propagation step
    correctionMaxAdj(nodes_adj, nodes_cate)
    propagationMaxAdj(nodes_adj, nodes_cate)

    # generate final classification of contigs by its edges' classification in the assembly graph
    file_out = outdir+"4CAC_classification.fasta"
    finalContigClassifyFromEdges(file_out, contigs_map, contigs_nodeList, nodes_map, nodes_length, nodes_cate, initial_class)

    print("Thanks for using the 4CAC classification algorithm! ")


if __name__=='__main__':
    args = parse_user_input()
    main(args)
