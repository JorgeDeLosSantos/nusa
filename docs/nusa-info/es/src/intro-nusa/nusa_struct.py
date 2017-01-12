from graphviz import Digraph

dot = Digraph(comment='The Round Table', format="png", filename="nusa_structure")
dot.edge("NuSA","core")
dot.edge("NuSA","element")
dot.edge("NuSA","model")
dot.edge("NuSA","lib")
dot.edge("NuSA","io")
dot.render()
