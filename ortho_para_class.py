""" Module information:
    dendropy - Provides tools for phylogenetic tree manipulation
    itertools - provides a way to iterate through the combination of two list"""

import itertools

class Event:
    """ Class representing an event (speciation or duplication) in a phylogentic tree."""
    def __init__(self):
        self.type = None
        # Sets of OTU that have been through speciation or not at a one event/node in the tree
        self.left_speciation = None
        self.left_no_speciation = None
        self.right_speciation = None
        self.right_no_speciation = None

class RelTree:
    """ Class representing the tree of species relationships """
    def __init__(self, nwk, sep, *, targets=None, id_first=False):
        self.tree = nwk
        self.separator = sep
        self.targets = targets
        self.id_first = id_first

    def display_tree(self):
        """ Prints a tree using ascii characters. """
        print(self.tree.as_ascii_plot())

    def get_species(self, label):
        """ Gets the species when interpreting a Newick Tree. """
        i = 0
        if self.id_first:
            i = 1

        if " " in label:
            return label.split(" ")[i]

        elif self.separator not in label:
            print("Separator does not exist in all OTU labels\n")
            exit()

        return label.split(self.separator)[i]

    def label_tree_events(self, tree):
        """ Labels pairs of taxons paralogs or orthologs \n
            Orthologs -> Result of a speciation event \n
            Paralogs -> Result of a duplication event"""
        relationship_dict = {}
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                node.species = [self.get_species(str(node.taxon))]
                node.taxon.sp_occurred = "NO"
                node.clade = [node.taxon]
                relationship_dict[node.taxon.label] = {}
            else:
                self.assign_relationships(node, relationship_dict)
        return relationship_dict

    def assign_orthologous(self, clade1, clade2, relationship_dict):
        """ Method used when speciation occurs. \n 
            Labels combination of species in each clade orthologous \n
            (Ex. Child 1 has species {man-a, man-b} and child 2 has species {mouse-a,mouse-b},
            so the combinations {man-a,mouse-a}, {man-a,mouse-b}, {man-b,mouse-a}, {man-b,mouse-b}
            would be orthologs)"""
        event = "orthologous"
        for tax_1, tax_2 in itertools.product(clade1, clade2):
            tax_1.sp_occurred = "YES"
            tax_2.sp_occurred = "YES"
            relationship_dict[tax_1.label][tax_2.label] = event
            relationship_dict[tax_2.label][tax_1.label] = event

    def assign_paralogous(self, clade1, clade2, relationship_dict):
        """ Method used when duplication occurs \n
            Labels combination of species in each clade paralogous"""
        for tax_1, tax_2 in itertools.product(clade1, clade2):
            if tax_1.sp_occurred == "YES" or tax_2.sp_occurred == "YES":
                event = "out-paralogous"
            elif (
                tax_1.sp_occurred == "NO" and tax_2.sp_occurred == "NO"
            ):
                event = "in-paralogous"
            else:
                if self.get_species(tax_1.label) == self.get_species(tax_2.label):
                    event = "paralogous"
                else:
                    event = "out-paralogous"

            relationship_dict[tax_1.label][tax_2.label] = event
            relationship_dict[tax_2.label][tax_1.label] = event

    def assign_ambigious(self, child1, child2, relationship_dict):
        """ Method used for polytomies when the event (duplication or speciation)
            cannot be determined \n
            Labels combination of species in each clade paralogous if there is a speciation event
            more recently in time between the OTU and the current event node. \n
            Otherwise, labels the relationship ambiguous."""
        if set(child1.species) & set(child2.species):
            special_event = "paralogous"
        else:
            special_event = "ambiguous"

        for tax_1, tax_2 in itertools.product(
            child1.clade, child2.clade
        ):
            if tax_1.sp_occurred == "NO":
                tax_1.sp_occurred = "UNKNOWN"
            if tax_2.sp_occurred == "NO":
                tax_2.sp_occurred = "UNKNOWN"

            if special_event == "paralogous" and (
                (
                    tax_1.sp_occurred == "YES"
                    or tax_2.sp_occurred == "YES"
                )
                or not (
                    self.get_species(tax_1.label)
                    == self.get_species(tax_2.label)
                )
            ):
                event = "out-paralogous"
            elif special_event == "paralogous":
                event = "paralogous"
            else:
                event = "ambigious"

            relationship_dict[tax_1.label][tax_2.label] = event
            relationship_dict[tax_2.label][tax_1.label] = event

    def get_poly_event(self,num_children, children):
        """ Determines if all relationships in a polytomy
            are orthologous, paralogous, or if it could be both (ambiguous)"""
        speciation = False
        duplication = False

        for i in range(0, num_children):
            for k in range(i, num_children):
                if i != k:
                    if speciation and duplication:
                        break
                    if not set(children[i].species) & set(children[k].species):
                        speciation = True
                    else:
                        duplication = True

        if speciation & duplication:
            return "ambiguous"
        elif speciation:
            return "orthologous"
        else:
            return "paralogous"

    def assign_relationships(self, node, relationship_dict):
        """ Assign paralogous or orthologous relationship(s)
            at a given node and store it in a dictionary \n
            Input:\n
            node -> The node to assign relationships at. 
            (Each children's clade will be assigned a relationship to the other clades).\n
            relationship_dict -> The dictionary to store the relationship information in.\n
            Relationships can be orthologous, in-paralogous, out-paralogous, and ambigious.
            Ambigious cases only happen after/during polytomies 
            when speciation events are unknown between two OTU."""
        children = node.child_nodes()
        num_children = node.num_child_nodes()

        if num_children == 2:
            node.species = children[0].species + children[1].species
            node.clade = children[0].clade + children[1].clade

            if not set(children[0].species) & set(children[1].species):
                self.assign_orthologous(children[0].clade,
                                        children[1].clade,
                                        relationship_dict)
            else:
                self.assign_paralogous(children[0].clade,
                                       children[1].clade,
                                       relationship_dict)
        else:
            # There are more than 2 children (polytomy)
            poly_event = self.get_poly_event(num_children, children)

            node.clade = []
            node.species = []

            for child in children:
                node.clade += child.clade
                node.species += child.species

            for i in range(0, num_children):
                for k in range(i, num_children):
                    if i != k:
                        if poly_event == "ambiguous":
                            self.assign_ambigious(children[i], children[k], relationship_dict)
                        elif poly_event == "orthologous":
                            self.assign_orthologous(children[i].clade,
                                                    children[k].clade,
                                                    relationship_dict)
                        else:
                            self.assign_paralogous(children[i].clade,
                                                   children[k].clade,
                                                   relationship_dict)

    #Compact version of labeling a tree
    def label_tree_events_compact(self):
        """ Labels by sets of clades at each tree event
            as opposed to a matrix of each individual relationship"""
        tree = self.tree
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                node.species = [self.get_species(str(node.taxon))]
                node.sp_occur = set()
                node.no_sp_occur = set([node.taxon.label])
            else:
                self.assign_relationships_compact(node)

    def assign_relationships_compact(self, node):
        """ Function to assign relationships to sets of OTUs"""
        children = node.child_nodes()
        num_children = node.num_child_nodes()

        if num_children == 2:
            node.species = children[0].species + children[1].species

            if not set(children[0].species) & set(children[1].species):
                node.event = "S" #-> Speciation
                node.sp_occur = children[0].sp_occur \
                              | children[1].sp_occur \
                              | children[0].no_sp_occur \
                              | children [1].no_sp_occur
                node.no_sp_occur = set()
            else:
                node.event = "D" #-> Duplication
                node.sp_occur = children[0].sp_occur | children[1].sp_occur
                node.no_sp_occur = children[0].no_sp_occur | children [1].no_sp_occur
        elif num_children == 1:
            print(num_children)

    def print_compact_relationship(self):
        """ Print a compact version of relationships in the following format:
            (ORTHOLOGOUS/PARALOGOUS: OUT1,OUT2 <=====> OTU3,OTU4)
            Ex: ORTHOLOGOUS mouse_a, mouse_b <===> mouse_c
            Indicates that that pairs (mouse_a, mouse_c) and (mouse_b,mouse_c) are orthologous"""
        self.label_tree_events_compact()

        tree = self.tree

        for node in tree.postorder_internal_node_iter():
            children = node.child_nodes()
            if node.event == "S":
                print('   ORTHOLOGY RELATIONSHIP:',
                       ','.join((children[0].sp_occur | children[0].no_sp_occur)),
                       "<====>", ','.join(children[1].sp_occur | children[1].no_sp_occur))
            else:
                if(children[0].no_sp_occur and children[1].no_sp_occur):
                    print('   IN-PARALOGOUS RELATIONSHIP:',
                           ','.join(children[0].no_sp_occur),
                           "<====>", ','.join(children[1].no_sp_occur))

                if(children[0].sp_occur and (children[0].no_sp_occur or children[1].no_sp_occur)):
                    print('   OUT-PARALOGOUS RELATIONSHIP:',
                           ','.join(children[0].sp_occur),
                           "<====>",','.join(children[1].sp_occur | children[1].no_sp_occur))
                if children[1].sp_occur:
                    print('   OUT-PARALOGOUS RELATIONSHIP:',
                           ','.join(children[0].sp_occur | children[0].no_sp_occur),
                           "<====>", ','.join(children[1].sp_occur))



    def print_relationships_of_child(self, target_child, relationship_dict):
        """Print relationships of a target OTU to every other OTU"""
        child_list = relationship_dict.keys()

        print(target_child + ":")

        for child in child_list:
            if target_child != child:
                print("\t" + child + ": " + relationship_dict[target_child][child])

    def print_all_relationships(self, child_list, relationship_dict):
        """ Print all relationships of each combination of OTUs in the given
            child list and all OTUs in the tree"""
        for child in child_list:
            self.print_relationships_of_child(child, relationship_dict)
            print("")

    def write_relationships_of_child(self, target_child, relationship_dict, directory):
        """ Given a 2d dictionary of relationships,
            output all relationship for the target taxon to a text file"""
        child_list = relationship_dict.keys()

        fname = directory + target_child.replace(self.separator, "") + ".csv"

        with open(fname,"w",encoding="utf-8") as f:
            f.write("taxon,relationship,extra_info\n")

            for child in child_list:
                if target_child != child:
                    relationship = relationship_dict[target_child][child]
                    if relationship in {"in-paralogous","out-paralogous"}:
                        f.write(child + "," + "paralogous" + "," + relationship)
                    else:
                        f.write(child + "," + relationship)
                    f.write("\n")

    def write_all_relationships(self, child_list, relationship_dict, directory):
        """ Write all relationships of each combination of OTUs in the given
            child list and all OTUs in the tree \n
            Writes file to given directory"""
        for child in child_list:
            self.write_relationships_of_child(child, relationship_dict, directory)

    def get_relationship_dict(self):
        """ Returns a dictionary of all relationships between every pair of OTUs"""
        return self.label_tree_events(self.tree)
