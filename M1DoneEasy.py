from networkx import *
from matplotlib.pyplot import *
options = {
    'node_color': 'black',
    'node_size': 100,
    'width': 3,
}
def printGraph(G):
    subplot(221)
    draw_random(G, **options)
    subplot(222)
    draw_circular(G, **options)
    subplot(223)
    draw_spectral(G, **options)
    subplot(224)
    draw_shell(G, nlist=[range(5,10), range(5)], **options)
    show()
