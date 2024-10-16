import numpy as np
from lxml import etree
from itertools import combinations
my_file = open("Rxn_RyRCaM.xml", "w")
molecules = ["Ca"]
kdiff = {
    "Ca": "100",
    "CaER": "10",
}

other_indices = {
    0: [1, 2],
    1: [0, 2],
    2: [0, 1]
}

other_indices_4 = {
    0: [1, 2, 3],
    1: [0, 2, 3],
    2: [0, 1, 3],
    3: [0, 1, 2],
}

ca = "Ca"
#calculated from EdSchutter's models
kf = {"2C": 17e-13, "2N": 1.46e-12}
kr = {"2C": 35e-4, "2N": 60e-3}
ryr_cam_binding = ["CaM",  "CaMCa2C", "CaMCa2N", "CaMCa4"]
kf_rev = {"2N": 17e-13, "2C": 1.46e-12}
kr_rev = {"2N": 35e-4, "2C": 60e-3}

def make_product_name_2(no_specie1, specie1, no_specie2, specie2, RyR_states):
    if no_specie1 == 4 and no_specie2 == 0:
        return "RyR_4%s" % specie1
    if no_specie1 == 0 and no_specie2 == 4:
        return "RyR_4%s" % specie2
    new_name = "RyR_%d%s_%d%s" % (no_specie1, specie1, no_specie2, specie2)
    return check_specie_2(new_name, RyR_states)


def make_product_name_3(no_specie1, specie1, no_specie2, specie2, no_specie3,
                        specie3, RyR_states):
    if no_specie1 == 0:
        new_name = "RyR_%d%s_%d%s" % (no_specie2, specie2, no_specie3, specie3)
        return check_specie_2(new_name, RyR_states)
    if no_specie2 == 0:
        new_name = "RyR_%d%s_%d%s" % (no_specie1, specie1, no_specie3, specie3)
        return check_specie_2(new_name, RyR_states)
    if no_specie3 == 0:
        new_name = "RyR_%d%s_%d%s" % (no_specie1, specie1, no_specie2, specie2)
        return check_specie_2(new_name, RyR_states)
    
    if specie1 == specie2:
        new_name = "RyR_%d%s_%d%s" % (no_specie1+no_specie2, specie1,
                                      no_specie3, specie3)
        return check_specie_2(new_name, RyR_states)

    if specie1 == specie3:
        new_name = "RyR_%d%s_%d%s" % (no_specie1+no_specie3, specie1,
                                      no_specie2, specie2)
        return check_specie_2(new_name, RyR_states)
    
    if specie3 == specie2:
        new_name = "RyR_%d%s_%d%s" % (no_specie1, specie1,
                                      no_specie3+no_specie2, specie3)
        return check_specie_2(new_name, RyR_states)

    
    new_name = "RyR_%d%s_%d%s_%d%s" % (no_specie1, specie1, no_specie2,
                                      specie2, no_specie3, specie3)

    return check_specie_3(new_name, RyR_states)


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
    return rxn_idx+1


def check_specie_2(specie, RyR_states):
    if specie in RyR_states:
        return specie
    splicik = specie.split("_")
    if len(splicik) == 3:
        return "RyR_%s_%s" % (splicik[2], splicik[1])


def check_specie_3(specie, RyR_states):
    if specie in RyR_states:
        return specie
    splicik = specie.split("_")
    alt = "RyR_%s_%s_%s" %( splicik[3], splicik[1], splicik[2])
    if alt in RyR_states:
        return alt
    alt = "RyR_%s_%s_%s" %( splicik[3], splicik[2], splicik[1])
    if alt in RyR_states:
        return alt
    alt = "RyR_%s_%s_%s" %( splicik[1], splicik[3], splicik[2])
    if alt in RyR_states:
        return alt
    alt = "RyR_%s_%s_%s" %( splicik[2], splicik[1], splicik[3])
    if alt in RyR_states:
        return alt
    alt = "RyR_%s_%s_%s" %( splicik[2], splicik[3], splicik[1])
    if alt in RyR_states:
        return alt

    
#RyR binding calmodulin
RyR_states = []
root = etree.Element("ReactionScheme")

for molecule in molecules:
    if molecule in kdiff:
        m_kdiff = kdiff[molecule]
    else:
        m_kdiff = "0"
    child = etree.SubElement(root, "Specie", name=molecule,
                             id=molecule, kdiff=m_kdiff, kdiffunit="mu2/s")

#RyR binding CaM states:
for specie in ryr_cam_binding:
    new_name = "RyR_4%s" % specie
    RyR_states.append(new_name)

for comb1 in combinations(range(4), 2):
    specie1 = ryr_cam_binding[comb1[0]]
    specie2 = ryr_cam_binding[comb1[-1]]
    for i in range(1, 4):
        j = 4 - i
        new_name = "RyR_%d%s_%d%s" % (i, specie1, j, specie2)
        RyR_states.append(new_name)
        new_name = "RyR_%d%s_%d%s" % (j, specie1, i, specie2)
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
        RyR_states.append(new_name)

RyR_states.append("RyR_1CaM_1CaMCa2C_1CaMCa2N_1CaMCa4")

for state in sorted(list(set(RyR_states))):
    child = etree.SubElement(root, "Specie", name=state,
                             id=state, kdiff="0",
                             kdiffunit="mu2/s")
        
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
            rxn_idx = add_reaction(root, ryr_specie, product, kf["2N"], kr["2N"],
                                   rxn_idx)
        elif ryr_specie == "RyR_4CaMCa2N":
            product = "RyR_3CaMCa2N_1CaMCa4"
            product = check_specie_2(product, RyR_states)
            rxn_idx = add_reaction(root, ryr_specie, product, kf["2C"], kr["2C"],
                                   rxn_idx)

        else:
            product = "RyR_3CaM_1CaMCa2C"
            product = check_specie_2(product, RyR_states)
            rxn_idx = add_reaction(root, ryr_specie, product, kf["2C"], kr["2C"], rxn_idx)
            
            product = "RyR_3CaM_1CaMCa2N"
            product = check_specie_2(product, RyR_states)
            rxn_idx = add_reaction(root, ryr_specie, product, kf["2N"], kr["2N"], rxn_idx)
            

    elif len(components) == 3:

        no_s = [int(components[1][0]), int(components[2][0])]
        species = [components[1][1:], components[2][1:]]
        
        if "CaMCa4" in species:
            camca4_idx = species.index("CaMCa4")
        else:
            camca4_idx = None

        if camca4_idx is not None:
            no_camca4 = no_s[camca4_idx]
            specie1 = species[not camca4_idx]
            no_specie = no_s[not camca4_idx]
            endswith = specie1[-2:]
            #specie binding reactions:
            if endswith in ["2C", "2N"]:
                product = make_product_name_2(no_camca4+1, "CaMCa4",
                                              no_specie-1, specie1, RyR_states)
                rxn_idx = add_reaction(root, ryr_specie, product,
                                       kf_rev[endswith], kr_rev[endswith], rxn_idx)
                
            else:
                product = make_product_name_3(no_camca4, "CaMCa4", no_specie-1,
                                              specie1, 1, "CaMCa2N", RyR_states)
                rxn_idx = add_reaction(root, ryr_specie, product, kf["2N"],
                                       kr["2N"], rxn_idx)
                
                product = make_product_name_3(no_camca4, "CaMCa4", no_specie-1,
                                              specie1, 1, "CaMCa2C", RyR_states)
                rxn_idx = add_reaction(root, ryr_specie, product, kf["2C"],
                                       kr["2C"], rxn_idx)
                
            continue
        if "CaM" in species:
            cam_idx = species.index("CaM")
        else:
            cam_idx = None

        if cam_idx is not None:
            no_cam = no_s[cam_idx]
            specie1 = species[not cam_idx]
            no_specie = no_s[not cam_idx]
            endswith = specie1[-2:]
            if endswith in ["2C", "2N"]:
                product = make_product_name_3(1, "CaMCa4", no_cam, "CaM",
                                              no_specie-1,
                                              specie1, RyR_states)
                
                rxn_idx = add_reaction(root, ryr_specie, product,
                                       kf_rev[endswith], kr_rev[endswith], rxn_idx)
                
                product = make_product_name_2(no_cam-1, 'CaM', no_specie+1,
                                              specie1, RyR_states)
                rxn_idx = add_reaction(root, ryr_specie, product,
                                       kf[endswith], kr[endswith], rxn_idx)
                
            continue
        #   RyR_n1CaMCa2C_n2CaMCa2N
        specie1 = species[0]
        no_1 = no_s[0]
        specie2 = species[1]
        no_2 = no_s[1]

        #  1: RyR_n1CaMCa2C_n2CaMCa2N + 2 Ca <-> RyR_(n1-1)CaMCa2C_n2CaMCa2N_1CaMCa4 kf["2N"], kr["2N"]
        product = make_product_name_3(1, "CaMCa4", no_1-1, specie1,
                                      no_2,
                                      specie2, RyR_states)
        rxn_idx = add_reaction(root, ryr_specie, product, kf_rev[specie1[-2:]],
                     kr_rev[specie1[-2:]], rxn_idx)
        
        #  2: RyR_n1CaMCa2C_n2CaMCa2N + 2 Ca <-> RyR_n1CaMCa2C_(n2-1)CaMCa2N_1CaMCa4 kf["2C"], kr["2C"]
        product = make_product_name_3(1, "CaMCa4", no_1, specie1,
                                      no_2-1,
                                      specie2, RyR_states)
        rxn_idx = add_reaction(root, ryr_specie, product,
                     kf_rev[specie2[-2:]], kr_rev[specie2[-2:]],
                     rxn_idx)
        
    elif len(components) == 4:
        no_s = [int(components[1][0]), int(components[2][0]),
                int(components[3][0])]
        species = [components[1][1:], components[2][1:],
                   components[3][1:]]
        for i, cur_specie in enumerate(species):
            cur_no = no_s[i]
            if cur_specie == "CaMCa4":
                continue
            elif cur_specie == "CaM":
                #1 2C
                if "CaMCa2C" in species:
                    idx = species.index("CaMCa2C")
                    no_camca2c = no_s[idx]
                    if idx in [0, 1] and i in [0,1]:
                        other = 2
                    elif idx in [0, 2] and i in [0,2]:
                        other = 1
                    else:
                        other = 0
                    other_no = no_s[other]
                    other_sp = species[other]
                    product = make_product_name_3(no_camca2c+1,
                                                  "CaMCa2C",
                                                  cur_no-1, "CaM",
                                                  other_no,
                                                  other_sp,
                                                  RyR_states)
                    rxn_idx = add_reaction(root, ryr_specie, product, kf["2C"],
                                 kr["2C"], rxn_idx)
                    
                else:
                    if cur_no == 2:
                        product = "RyR_1CaM_1CaMCa2C_1CaMCa2N_1CaMCa4"
                        rxn_idx = add_reaction(root, ryr_specie, product,
                                     kf["2C"],
                                     kr["2C"], rxn_idx)
                        
                    else:
                        inds = other_indices[i]
                        product = make_product_name_3(1,
                                                      "CaMCa2C",
                                                      no_s[inds[0]],
                                                      species[inds[0]],
                                                      no_s[inds[1]],
                                                      species[inds[1]],
                                                      RyR_states)
                        rxn_idx = add_reaction(root, ryr_specie, product, kf["2C"],
                                               kr["2C"], rxn_idx)
                         
                #1 2N
                if "CaMCa2N" in species:
                    idx = species.index("CaMCa2N")
                    no_camca2c = no_s[idx]
                    if idx in [0, 1] and i in [0,1]:
                        other = 2
                    elif idx in [0, 2] and i in [0,2]:
                        other = 1
                    else:
                        other = 0
                    other_no = no_s[other]
                    other_sp = species[other]
                    product = make_product_name_3(no_camca2c+1,
                                                  "CaMCa2N",
                                                  cur_no-1, "CaM",
                                                  other_no,
                                                  other_sp,
                                                  RyR_states)
                    rxn_idx = add_reaction(root, ryr_specie, product, kf["2N"],
                                 kr["2N"], rxn_idx)
                    
                else:
                    if cur_no == 2:
                        product = "RyR_1CaM_1CaMCa2C_1CaMCa2N_1CaMCa4"
                        rxn_idx = add_reaction(root, ryr_specie, product,
                                     kf["2N"],
                                     kr["2N"], rxn_idx)
                        
                    else:
                        inds = other_indices[i]
                      
                        product = make_product_name_3(1,
                                                      "CaMCa2N",
                                                      no_s[inds[0]],
                                                      species[inds[0]],
                                                      no_s[inds[1]],
                                                      species[inds[1]],
                                                      RyR_states)
                        rxn_idx = add_reaction(root, ryr_specie, product, kf["2N"],
                                               kr["2N"], rxn_idx)
            else:
                endswith = cur_specie[-2:]
                if "CaMCa4" in species:
                    idx = species.index("CaMCa4")
                    no_camca4 = no_s[idx]
                    if idx in [0, 1] and i in [0,1]:
                        other = 2
                    elif idx in [0, 2] and i in [0,2]:
                        other = 1
                    else:
                        other = 0
                    other_no = no_s[other]
                    other_sp = species[other]
                    product = make_product_name_3(cur_no-1, cur_specie,
                                                  no_camca4+1, "CaMCa4",
                                                  other_no,
                                                  other_sp,
                                                  RyR_states)
                    rxn_idx = add_reaction(root, ryr_specie, product,
                                           kf_rev[endswith], kr_rev[endswith],
                                           rxn_idx)
                else:
                    if cur_no == 2:
                        product = "RyR_1CaM_1CaMCa2C_1CaMCa2N_1CaMCa4"
                    else:

                        inds = other_indices[i]
                        product = make_product_name_3(1, "CaMCa4",
                                                      no_s[inds[0]],
                                                      species[inds[0]],
                                                      no_s[inds[1]],
                                                      species[inds[1]],
                                                      RyR_states)
                    rxn_idx = add_reaction(root, ryr_specie, product,
                                           kf_rev[endswith], kr_rev[endswith],
                                           rxn_idx)
                        
                    
    else:
        # RyR_1CaM_1CaMCa2C_1CaMCa2N_1CaMCa4

        #  1 RyR_2CaMCa2C_1CaMCa2N_1CaMCa4                        
        product = make_product_name_3(1, "CaMCa4",
                                      1, "CaMCa2N",
                                      2, "CaMCa2C",
                                      RyR_states)
        rxn_idx = add_reaction(root, ryr_specie, product,
                               kf["2C"], kr["2C"],
                               rxn_idx)

        #  2 RyR_1CaMCa2C_2CaMCa2N_1CaMCa4
        product = make_product_name_3(1, "CaMCa4",
                                      2, "CaMCa2N",
                                      1, "CaMCa2C",
                                      RyR_states)
        rxn_idx = add_reaction(root, ryr_specie, product,
                               kf["2N"], kr["2N"],
                               rxn_idx)

        #  3 RyR_1CaM_1CaMCa2N_2CaMCa4
        product = make_product_name_3(2, "CaMCa4",
                                      1, "CaMCa2N",
                                      1, "CaM",
                                      RyR_states)
        rxn_idx = add_reaction(root, ryr_specie, product,
                               kf["2N"], kr["2N"],
                               rxn_idx)
        #  4 RyR_1CaM_1CaMCa2C_2CaMCa4
        product = make_product_name_3(2, "CaMCa4",
                                      1, "CaMCa2C",
                                      1, "CaM",
                                      RyR_states)
        rxn_idx = add_reaction(root, ryr_specie, product,
                               kf["2C"], kr["2C"],
                               rxn_idx)
        
        
my_file.write('<?xml version="1.0"?>\n')
my_file.write(etree.tostring(root, pretty_print=True).decode('utf-8'))
