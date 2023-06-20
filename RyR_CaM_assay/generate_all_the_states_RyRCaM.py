import numpy as np
from lxml import etree
from itertools import combinations
my_file = open("Rxn_RyRCaM.xml", "w")
molecules = ["Ca"]
kdiff = {
    "Ca": "100",
    "CaER": "10",
}
ryr_name = "RyR_4CaM"

ca = "Ca"
#calculated from EdSchutter's models
kf = {"2C": 17e-10, "2N": 1.46e-9}
kr = {"2C": 35e-4, "2N": 60e-3}
ryr_cam_binding = ["CaM",  "CaMCa2C", "CaMCa2N", "CaMCa4"]
kf_rev = {"2N": 17e-10, "2C": 1.46e-9}
kr_rev = {"2N": 35e-4, "2C": 60e-3}

def make_product_name_2(no_specie1, specie1, no_specie2, specie2, RyR_states):
    if no_specie1 == 4 and no_specie2 == 0:
        return "RyR_4%s" % specie1
    new_name = "RyR_%d%s_%d%s" % (no_specie1, specie1, no_specie2, specie2)
    return check_specie_2(new_name, RyR_states)


def add_reaction(root1, specie1, product1, kf1, kr1, rxn_idx):
    my_reac = etree.SubElement(root1, "Reaction",
                               name="RyRCaM%d" % rxn_idx,
                               id="RyRCaM%d" % rxn_idx)
    etree.SubElement(my_reac, "Reactant", specieID=specie1)
    etree.SubElement(my_reac, "Reactant", specieID="Ca", power="2")
    etree.SubElement(my_reac, "Product", specieID=product1)
    kf = etree.SubElement(my_reac, "forwardRate")
    kf.text = str(kf1)
    kr = etree.SubElement(my_reac, "reverseRate")
    kr.text = str(kr1)
    print(etree.tostring(my_reac))
    return my_reac

def check_specie_2(specie, RyR_states):
    if specie in RyR_states:
        return specie
    splicik = specie.split("_")
    if len(splicik) == 3:
        return "RyR_%s_%s" % (splicik[2], splicik[1])
    return

def check_specie_3(specie, RyR_states):
    if specie in RyR_states:
        return specie
    splicik = specie.split("_")
    alt = "RyR_%s_%s_%s" %( splicik[2], splicik[0], splicik[1])
    if alt in RyR_states:
        return alt
    alt = "RyR_%s_%s_%s" %( splicik[2], splicik[1], splicik[0])
    if alt in RyR_states:
        return alt
    return

    
#RyR binding calmodulin
RyR_states = []
root = etree.Element("ReactionScheme")

for molecule in molecules:
    if molecule in kdiff:
        m_kdiff = kdiff[molecule]
    else:
        m_kdiff = "0"
    child = etree.SubElement(root, "Specie", name=molecule, id=molecule, kdiff=m_kdiff, kdiffunit="mu2/s")

#RyR binding CaM states:
for specie in ryr_cam_binding:
    new_name = "RyR_4%s" % specie
    child = etree.SubElement(root, "Specie", name=new_name,
                                     id=new_name, kdiff="0",
                                     kdiffunit="mu2/s")
    RyR_states.append(new_name)

for comb1 in combinations(range(4), 2):
    specie1 = ryr_cam_binding[comb1[0]]
    specie2 = ryr_cam_binding[comb1[-1]]
    for i in range(1, 4):
        j = 4 - i
        new_name = "RyR_%d%s_%d%s" % (i, specie1, j, specie2)
        child = etree.SubElement(root, "Specie", name=new_name,
                                 id=new_name, kdiff="0",
                                 kdiffunit="mu2/s")
        RyR_states.append(new_name)
        if i != j:
            new_name = "RyR_%d%s_%d%s" % (j, specie2, i, specie1)
            child = etree.SubElement(root, "Specie", name=new_name,
                                     id=new_name, kdiff="0",
                                     kdiffunit="mu2/s")
            RyR_states.append(new_name)
            

numbers = [[1, 1, 2], [1, 2, 1], [2, 1, 1]]    
for comb2 in combinations(range(4), 3):
    specie1 = ryr_cam_binding[comb2[0]]
    specie2 = ryr_cam_binding[comb2[1]]
    specie3 = ryr_cam_binding[comb2[2]]
    for number in numbers:
        new_name = "RyR_%d%s_%d%s_%d%s" % (number[0], specie1,
                                           number[1], specie2,
                                           number[2], specie3)
        child = etree.SubElement(root, "Specie", name=new_name,
                                 id=new_name, kdiff="0",
                                 kdiffunit="mu2/s")
        RyR_states.append(new_name)

rxn_idx = 0

#start with RyR_4CaM:

for ryr_specie in RyR_states:
    components = ryr_specie.split("_")

    if len(components) == 2:
        if ryr_specie == "RyR_4CaMCa4":
            continue
        elif ryr_specie == "RyR_4CaMCa2C":
            product = "RyR_3CaMCa2C_1CaMCa4"
            product = check_specie_2(product, RyR_states)
            add_reaction(root, ryr_specie, product, kf["2N"], kr["2N"], rxn_idx)
            rxn_idx += 1
        elif ryr_specie == "RyR_4CaMCa2N":
            product = "RyR_3CaMCa2N_1CaMCa4"
            product = check_specie_2(product, RyR_states)
            add_reaction(root, ryr_specie, product, kf["2C"], kr["2C"], rxn_idx)
            rxn_idx += 1
        else:
            product = "RyR_3CaM_1CaMCa2C"
            product = check_specie_2(product, RyR_states)
            add_reaction(root, ryr_specie, product, kf["2C"], kr["2C"], rxn_idx)
            rxn_idx += 1
            product = "RyR_3CaM_1CaMCa2N"
            product = check_specie_2(product, RyR_states)
            add_reaction(root, ryr_specie, product, kf["2N"], kr["2N"], rxn_idx)
            rxn_idx += 1

    elif len(components) == 3:

        no_s = [int(components[1][0]), int(components[2][0])]
        species = [components[1][1:], components[2][1:]]
        
        if "CaMCa4" in species:
            camca4_idx = species.index("CaMCa4")
        else:
            camca4_idx = None

        if camca4_idx is not None:
            print(components, no_s, species, no_s[camca4_idx], no_s[not camca4_idx])
            no_camca4 = no_s[camca4_idx]
            specie1 = species[not camca4_idx]
            no_specie = no_s[not camca4_idx]
            endswith = specie1[-2:]
            print(endswith)
            #specie binding reactions:
            print(no_camca4, no_specie, ryr_specie)
            if endswith in ["2C", "2N"]:
                product = make_product_name_2(no_camca4+1, "CaMCa4", no_specie-1, specie1, RyR_states)
                add_reaction(root, ryr_specie, product, kf_rev[endswith], kr_rev[endswith], rxn_idx)
                rxn_idx += 1
        
                
my_file.write('<?xml version="1.0"?>\n')
my_file.write(etree.tostring(root, pretty_print=True).decode('utf-8'))
