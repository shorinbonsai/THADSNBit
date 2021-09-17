import os
import sys
import networkx as nx
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyvis.network import Network



pathy = "./edgegraphmut2.txt"

graphy = nx.MultiGraph()
# graphy = nx.read_adjlist(pathy)
# print(graphy)
# print("Test")

edges = []

with open(pathy) as f:
        lines = f.readlines()
        county = 0
        for line in lines:
            line = line.strip()
            spl = line.split(' ')
            listA = []
            listC = []
            for j in range(len(spl)):
                spl[j] = int(spl[j])
            for i in range(len(spl)-1):
                listA.append(county)
                tempdict = {"weight":0.5}
                listC.append(tempdict)
            
            listB = spl[1:]
            zipped = zip(listA,listB, listC)
            edges.extend(zipped)
            county += 1
            # raise ValueError("ValueError exception thrown")

graphy.add_edges_from(edges)
G = nx.Graph()
for u,v,data in graphy.edges(data=True):
    w = data["weight"] if "weight" in data else 1.0
    if G.has_edge(u,v):
        G[u][v]["weight"] += w
    else:
        G.add_edge(u,v,weight=w)
for u,v,data in G.edges(data=True):
    if G.has_edge(u,v) and G.has_edge(v,u):
        print("duplicate")
        

print(G.edges(data=True))
# raise ValueError("ValueError exception thrown")


nt = Network('1000px', '1000px')
nt.from_nx(graphy)
# print(nt.get_edges())
# nt.show_buttons()
nt.show('nx.html')

# pos = nx.spring_layout(graphy)
# nx.draw_networkx_nodes(graphy, pos, node_color = 'r', node_size = 100, alpha = 1)
# ax = plt.gca()
# for e in graphy.edges:
#     ax.annotate("",
#                 xy=pos[e[0]], xycoords='data',
#                 xytext=pos[e[1]], textcoords='data',
#                 arrowprops=dict(arrowstyle="->", color="0.5",
#                                 shrinkA=5, shrinkB=5,
#                                 patchA=None, patchB=None,
#                                 connectionstyle="arc3,rad=rrr".replace('rrr',str(0.3*e[2])
#                                 ),
#                                 ),
#                 )
# plt.axis('off')
# nx.draw_networkx_labels(graphy, pos, font_size=8)
# labels = nx.get_edge_attributes(graphy, 'weight')

# nx.draw_networkx_edge_labels(graphy,pos,edge_labels=labels, font_size=6)
# nx.draw_networkx_labels(graphy, pos, font_size=8)
# plt.savefig("testvis.png")
# plt.show()
plt.close()
