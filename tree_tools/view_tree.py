import dendropy
import sys
tree = dendropy.Tree.get_from_path(sys.argv[1], 'newick')
print(tree.as_ascii_plot())
