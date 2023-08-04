import numpy as np
from lxml import etree
from itertools import combinations

molecules = ["Ca", "CaM"]
kdiff = {
    "Ca": "100",
    "CaER": "10",
    "CaM":"4",
}
ca = "Ca"

kfl = 1e-03
krl = 1e-04

kfcam = 2.2222222222222222e-07 #   10.1161/circresaha.114.302857
krcam = 3.3333333333333333e-06

kfca = 12e-10#2.4e-05   #   kd for Ca 15uM (Shifferaw?)
krca = 0.0025
kfC1O = 2
krC1O = 1  #  G&S 2022
kfOC2 = 0.06
krOC2 = 0.6  #  G&S 2022

include_C2 = False
def add_reaction(root1, specie1, specie2, product1, kf1, kr1, rxn_idx,
                 power="1"):
    my_reac = etree.SubElement(root1, "Reaction",
                               name="RyRCaM%d" % rxn_idx,
                               id="RyRCaM%d" % rxn_idx)
    etree.SubElement(my_reac, "Reactant", specieID=specie1)
    etree.SubElement(my_reac, "Reactant", specieID=specie2, power=power)
    etree.SubElement(my_reac, "Product", specieID=product1)
    kf = etree.SubElement(my_reac, "forwardRate")
    kf.text = str(kf1)
    kr = etree.SubElement(my_reac, "reverseRate")
    kr.text = str(kr1)
    return rxn_idx+1


def add_transition(root1, specie1, product1, kf1, kr1, rxn_idx):
    my_reac = etree.SubElement(root1, "Reaction",
                               name="RyRCaM%d" % rxn_idx,
                               id="RyRCaM%d" % rxn_idx)
    etree.SubElement(my_reac, "Reactant", specieID=specie1)
    etree.SubElement(my_reac, "Product", specieID=product1)
    kf = etree.SubElement(my_reac, "forwardRate")
    kf.text = str(kf1)
    kr = etree.SubElement(my_reac, "reverseRate")
    kr.text = str(kr1)
    return rxn_idx+1


def all_biding_states(my_states, starting_state, ending_state):
    possible_opening = []
    for state in my_states:
        if starting_state not in state:
            continue
        n = int(state.split(starting_state)[0][-1])
        new_base = state.split(str(n)+starting_state)[0][:-1]
        for i in range(1, n+1):
            new_l = n - i
            if new_l:
                new_name = "%s_%d%s_%d%s" % (new_base, new_l,
                                             starting_state, i,
                                             ending_state)
            else:
                new_name = "%s_%d%s" % (new_base, i,
                                        ending_state)
            possible_opening.append(new_name)
    return possible_opening



if __name__ == "__main__":
    cam_interactions = False
    if cam_interactions:
        my_file = open("Rxn_4_states.xml", "w")
    else:
        my_file = open("Rxn_4_states_no_cam.xml", "w")
    
    #RyR binding calmodulin
    RyR_states = ["RyR"]
    root = etree.Element("ReactionScheme")

    for molecule in molecules:
        if molecule in kdiff:
            m_kdiff = kdiff[molecule]
        else:
            m_kdiff = "0"
        child = etree.SubElement(root, "Specie", name=molecule,
                                 id=molecule, kdiff=m_kdiff, kdiffunit="mu2/s")


    lcam_comp = []
    #RyR binding CaM states:
    for i in range(1, 5):
        new_name = "RyR_%dL" % i
        RyR_states.append(new_name)
        lcam_comp.append(new_name)
        if cam_interactions:
            new_name = "RyR_%dCaM" % i
            RyR_states.append(new_name)
            lcam_comp.append(new_name)
            for j in range(1, 5-i):
                new_name = "RyR_%d%s_%d%s" % (i, "CaM", j, "L")
                RyR_states.append(new_name)
                lcam_comp.append(new_name)

    #Ca binding:
    without_Ca = RyR_states[:]
    possible_opening = all_biding_states(without_Ca, "L", "C1")
    RyR_states += possible_opening

    open_states = all_biding_states(possible_opening, "C1", "O")
    RyR_states += open_states
    if include_C2:
        other_closed = all_biding_states(open_states, "O", "C2")
        RyR_states += other_closed
    for state in sorted(list(set(RyR_states))):
        child = etree.SubElement(root, "Specie", name=state,
                                 id=state, kdiff="0",
                                 kdiffunit="mu2/s")
        
    rxn_idx = 0

    #L transitions (actually binding, but whatever it's on the same protein)
    for specie in lcam_comp:
        if "L" not in specie:
            continue
        n_l = int(specie.split("L")[0][-1])
        new_n_l = n_l - 1
        if new_n_l == 0:
            rxn_idx = add_transition(root, specie[:-3], specie, kfl, krl,
                                     rxn_idx)
        else:
            rxn_idx = add_transition(root, specie[:-2]+str(new_n_l)+"L", specie,
                                     kfl, krl, rxn_idx)

    for specie in RyR_states:
        if "CaM" not in specie:
            continue
        #start with the product
        n_cam = int(specie.split("CaM")[0][-1])
        rest = specie.split("CaM")[-1]
        new_n_cam = n_cam - 1
        if new_n_cam:
            reactant = specie.split(str(n_cam)+"CaM")[0]+str(new_n_cam)+"CaM"+rest
        else:
            reactant = specie.split(str(n_cam)+"CaM")[0][:-1] + rest
        rxn_idx = add_reaction(root, reactant, "CaM", specie, kfcam, krcam, rxn_idx)
    
    #Ca binding to L
    for specie in RyR_states:
        if "C1" not in specie:
            continue
        n_C1 = int(specie.split("C1")[0][-1])
        n_L = 0
        if "L" in specie:
            n_L = int(specie.split("L")[0][-1])
            new_base = specie.split("L")[0][:-2]
        else:
            new_base = specie.split("C1")[0][:-2]
        new_L = n_L + 1
        new_C1 = n_C1 - 1
        new_reactant = "%s_%dL" % (new_base, new_L)
        if new_C1:
            new_reactant = "%s_%dC1" % (new_reactant, new_C1)
        new_reactant += specie.split("C1")[-1]

        rxn_idx = add_reaction(root, new_reactant, "Ca", specie, kfca, krca, rxn_idx, power="2")
    #C1 to O transition

    for specie in RyR_states:
        if "O" not in specie:
            continue
        n_O = int(specie.split("O")[0][-1])
        n_C1 = 0
        if "C1" in specie:
            n_C1 = int(specie.split("C1")[0][-1])
            new_base = specie.split("C1")[0][:-2]
        else:
            new_base = specie.split("O")[0][:-2]
        new_C1 = n_C1 + 1
        new_O = n_O - 1
        new_reactant = "%s_%dC1" % (new_base, new_C1)
        if new_O:
            new_reactant = "%s_%dO" % (new_reactant, new_O)
        new_reactant += specie.split("O")[-1]
        rxn_idx = add_transition(root, new_reactant, specie, kfC1O, krC1O,
                                 rxn_idx)


    for specie in RyR_states:
        if "C2" not in specie:
            continue
        n_C2 = int(specie.split("C2")[0][-1])
        n_O = 0
        if "O" in specie:
            n_O = int(specie.split("O")[0][-1])
            new_base = specie.split("O")[0][:-2]
        else:
            new_base = specie.split("C2")[0][:-2]
        new_O = n_O + 1
        new_C2 = n_C2 - 1
        new_reactant = "%s_%dO" % (new_base, new_O)
        if new_C2:
            new_reactant = "%s_%dC2" % (new_reactant, new_C2)
        new_reactant += specie.split("C2")[-1]

        rxn_idx = add_transition(root, new_reactant, specie, kfOC2, krOC2,
                                 rxn_idx)
    

        
    my_file.write('<?xml version="1.0"?>\n')
    my_file.write(etree.tostring(root, pretty_print=True).decode('utf-8'))
    my_file.close()
    if cam_interactions:
        output_name = "IO_RyR.xml"
    else:
        output_name = "IO_RyR_no_cam.xml"
    with open(output_name, "w") as new_f:

        root = etree.Element("OutputScheme")
        dataset = etree.SubElement(root, "OutputSet",
                                   filename="all", dt="0.01")
        for specie in RyR_states+molecules:
            etree.SubElement(dataset, "OutputSpecie", name=specie)
        new_f.write('<?xml version="1.0"?>\n')
        new_f.write(etree.tostring(root, pretty_print=True).decode('utf-8'))
