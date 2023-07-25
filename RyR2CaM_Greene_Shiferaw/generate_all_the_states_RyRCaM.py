import numpy as np
from lxml import etree
from itertools import combinations
my_file = open("Rxn_4_states.xml", "w")
molecules = ["Ca"]
kdiff = {
    "Ca": "100",
    "CaER": "10",
}
ca = "Ca"

kfl = 2.4e-06
krl = 7e-05

def add_Ca_reaction(root1, specie1, product1, kf1, kr1, rxn_idx):
    my_reac = etree.SubElement(root1, "Reaction",
                               name="4RyRCaM%d" % rxn_idx,
                               id="4RyRCaM%d" % rxn_idx)
    etree.SubElement(my_reac, "Reactant", specieID=specie1)
    etree.SubElement(my_reac, "Reactant", specieID="Ca")
    etree.SubElement(my_reac, "Product", specieID=product1)
    kf = etree.SubElement(my_reac, "forwardRate")
    kf.text = str(kf1)
    kr = etree.SubElement(my_reac, "reverseRate")
    kr.text = str(kr1)
    return rxn_idx+1


def add_transition(root1, specie1, product1, kf1, kr1, rxn_idx):
    my_reac = etree.SubElement(root1, "Reaction",
                               name="4RyRCaM%d" % rxn_idx,
                               id="4RyRCaM%d" % rxn_idx)
    etree.SubElement(my_reac, "Reactant", specieID=specie1)
    etree.SubElement(my_reac, "Product", specieID=product1)
    kf = etree.SubElement(my_reac, "forwardRate")
    kf.text = str(kf1)
    kr = etree.SubElement(my_reac, "reverseRate")
    kr.text = str(kr1)
    return rxn_idx+1



    
#RyR binding calmodulin
RyR_states = ["4RyR"]
root = etree.Element("ReactionScheme")

for molecule in molecules:
    if molecule in kdiff:
        m_kdiff = kdiff[molecule]
    else:
        m_kdiff = "0"
    child = etree.SubElement(root, "Specie", name=molecule,
                             id=molecule, kdiff=m_kdiff, kdiffunit="mu2/s")

#RyR binding CaM states:
L_dict = {}
lcam_comp = []
for i in range(1, 5):
    new_name = "4RyR_%dCaM" % i
    RyR_states.append(new_name)
    lcam_comp.append(new_name)
    new_name = "4RyR_%dL" % i
    RyR_states.append(new_name)
    lcam_comp.append(new_name)
    for j in range(1, 5-i):
        new_name = "4RyR_%d%s_%d%s" % (i, "CaM", j, "L")
        RyR_states.append(new_name)
        lcam_comp.append(new_name)

#Ca binding:
without_Ca = RyR_states[:]
possible_opening = []
for state in without_Ca:
    if "L" not in state:
        continue
    n = int(state.split("L")[0][-1])
    new_base = state.split(str(n)+"L")[0][:-1]
    for i in range(1, n+1):
        new_l = n - i
        if new_l:
            new_name = "%s_%dL_%dC1" %(new_base, new_l, i)
        else:
            new_name = "%s_%dC1" %(new_base, i)
        RyR_states.append(new_name)
        possible_opening.append(new_name)

open_states = []

for specie in possible_opening:
    n_c1 = int(specie.split("C1")[0][-1])
    for k in range(1, n_c1+1):
        base = specie.split("C1")[0][:-2]
        new_n_c1 = n_c1 - k
        new_specie = base
        if new_n_c1:
            new_specie = "%s_%dC1_%dO" % (new_specie, new_n_c1, k)
        else:
            new_specie = "%s_%dO" % (new_specie, k)

        open_states.append(new_specie)
        RyR_states.append(new_specie)

#opening -- LCa is closed (C1)
other_closed = []
for specie in open_states:
    n_o = int(specie[-2])
    base = specie[:-3]
    for k in range(1, n_o + 1):
        new_o = n_o - k
        if not new_o:
            new_specie = "%s_%dC2" %(base, k)
        else:
            new_specie = "%s_%dO_%dC2" %(base, new_o, k)
        other_closed.append(new_specie)
        RyR_states.append(new_specie)

print(len(RyR_states))

for state in sorted(list(set(RyR_states))):
    child = etree.SubElement(root, "Specie", name=state,
                             id=state, kdiff="0",
                             kdiffunit="mu2/s")
        
rxn_idx = 0

#L transitions
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

    

# #start with 4RyR:
# for ryr_specie in RyR_states:
#     components = ryr_specie.split("_")
#     if len(components) == 2:
#         if ryr_specie == "4RyR_4CaMCa4":
#             continue
#         elif ryr_specie == "4RyR_4CaMCa2C":
#             product = "4RyR_3CaMCa2C_1CaMCa4"
#             product = check_specie_2(product, RyR_states)
#             rxn_idx = add_reaction(root, ryr_specie, product, kf["2N"], kr["2N"],
#                                    rxn_idx)
#         elif ryr_specie == "4RyR_4CaMCa2N":
#             product = "4RyR_3CaMCa2N_1CaMCa4"
#             product = check_specie_2(product, RyR_states)
#             rxn_idx = add_reaction(root, ryr_specie, product, kf["2C"], kr["2C"],
#                                    rxn_idx)

#         else:
#             product = "4RyR_3CaM_1CaMCa2C"
#             product = check_specie_2(product, RyR_states)
#             rxn_idx = add_reaction(root, ryr_specie, product, kf["2C"], kr["2C"], rxn_idx)
            
#             product = "4RyR_3CaM_1CaMCa2N"
#             product = check_specie_2(product, RyR_states)
#             rxn_idx = add_reaction(root, ryr_specie, product, kf["2N"], kr["2N"], rxn_idx)
            

#     elif len(components) == 3:

#         no_s = [int(components[1][0]), int(components[2][0])]
#         species = [components[1][1:], components[2][1:]]
        
#         if "CaMCa4" in species:
#             camca4_idx = species.index("CaMCa4")
#         else:
#             camca4_idx = None

#         if camca4_idx is not None:
#             no_camca4 = no_s[camca4_idx]
#             specie1 = species[not camca4_idx]
#             no_specie = no_s[not camca4_idx]
#             endswith = specie1[-2:]
#             #specie binding reactions:
#             if endswith in ["2C", "2N"]:
#                 product = make_product_name_2(no_camca4+1, "CaMCa4",
#                                               no_specie-1, specie1, RyR_states)
#                 rxn_idx = add_reaction(root, ryr_specie, product,
#                                        kf_rev[endswith], kr_rev[endswith], rxn_idx)
                
#             else:
#                 product = make_product_name_3(no_camca4, "CaMCa4", no_specie-1,
#                                               specie1, 1, "CaMCa2N", RyR_states)
#                 rxn_idx = add_reaction(root, ryr_specie, product, kf["2N"],
#                                        kr["2N"], rxn_idx)
                
#                 product = make_product_name_3(no_camca4, "CaMCa4", no_specie-1,
#                                               specie1, 1, "CaMCa2C", RyR_states)
#                 rxn_idx = add_reaction(root, ryr_specie, product, kf["2C"],
#                                        kr["2C"], rxn_idx)
                
#             continue
#         if "CaM" in species:
#             cam_idx = species.index("CaM")
#         else:
#             cam_idx = None

#         if cam_idx is not None:
#             no_cam = no_s[cam_idx]
#             specie1 = species[not cam_idx]
#             no_specie = no_s[not cam_idx]
#             endswith = specie1[-2:]
#             if endswith in ["2C", "2N"]:
#                 product = make_product_name_3(1, "CaMCa4", no_cam, "CaM",
#                                               no_specie-1,
#                                               specie1, RyR_states)
                
#                 rxn_idx = add_reaction(root, ryr_specie, product,
#                                        kf_rev[endswith], kr_rev[endswith], rxn_idx)
                
#                 product = make_product_name_2(no_cam-1, 'CaM', no_specie+1,
#                                               specie1, RyR_states)
#                 rxn_idx = add_reaction(root, ryr_specie, product,
#                                        kf[endswith], kr[endswith], rxn_idx)
                
#             continue
#         #   RyR_n1CaMCa2C_n2CaMCa2N
#         specie1 = species[0]
#         no_1 = no_s[0]
#         specie2 = species[1]
#         no_2 = no_s[1]

#         #  1: RyR_n1CaMCa2C_n2CaMCa2N + 2 Ca <-> RyR_(n1-1)CaMCa2C_n2CaMCa2N_1CaMCa4 kf["2N"], kr["2N"]
#         product = make_product_name_3(1, "CaMCa4", no_1-1, specie1,
#                                       no_2,
#                                       specie2, RyR_states)
#         rxn_idx = add_reaction(root, ryr_specie, product, kf_rev[specie1[-2:]],
#                      kr_rev[specie1[-2:]], rxn_idx)
        
#         #  2: RyR_n1CaMCa2C_n2CaMCa2N + 2 Ca <-> RyR_n1CaMCa2C_(n2-1)CaMCa2N_1CaMCa4 kf["2C"], kr["2C"]
#         product = make_product_name_3(1, "CaMCa4", no_1, specie1,
#                                       no_2-1,
#                                       specie2, RyR_states)
#         rxn_idx = add_reaction(root, ryr_specie, product,
#                      kf_rev[specie2[-2:]], kr_rev[specie2[-2:]],
#                      rxn_idx)
        
#     elif len(components) == 4:
#         no_s = [int(components[1][0]), int(components[2][0]),
#                 int(components[3][0])]
#         species = [components[1][1:], components[2][1:],
#                    components[3][1:]]
#         for i, cur_specie in enumerate(species):
#             cur_no = no_s[i]
#             if cur_specie == "CaMCa4":
#                 continue
#             elif cur_specie == "CaM":
#                 #1 2C
#                 if "CaMCa2C" in species:
#                     idx = species.index("CaMCa2C")
#                     no_camca2c = no_s[idx]
#                     if idx in [0, 1] and i in [0,1]:
#                         other = 2
#                     elif idx in [0, 2] and i in [0,2]:
#                         other = 1
#                     else:
#                         other = 0
#                     other_no = no_s[other]
#                     other_sp = species[other]
#                     product = make_product_name_3(no_camca2c+1,
#                                                   "CaMCa2C",
#                                                   cur_no-1, "CaM",
#                                                   other_no,
#                                                   other_sp,
#                                                   RyR_states)
#                     rxn_idx = add_reaction(root, ryr_specie, product, kf["2C"],
#                                  kr["2C"], rxn_idx)
                    
#                 else:
#                     if cur_no == 2:
#                         product = "4RyR_1CaM_1CaMCa2C_1CaMCa2N_1CaMCa4"
#                         rxn_idx = add_reaction(root, ryr_specie, product,
#                                      kf["2C"],
#                                      kr["2C"], rxn_idx)
                        
#                     else:
#                         inds = other_indices[i]
#                         product = make_product_name_3(1,
#                                                       "CaMCa2C",
#                                                       no_s[inds[0]],
#                                                       species[inds[0]],
#                                                       no_s[inds[1]],
#                                                       species[inds[1]],
#                                                       RyR_states)
#                         rxn_idx = add_reaction(root, ryr_specie, product, kf["2C"],
#                                                kr["2C"], rxn_idx)
                         
#                 #1 2N
#                 if "CaMCa2N" in species:
#                     idx = species.index("CaMCa2N")
#                     no_camca2c = no_s[idx]
#                     if idx in [0, 1] and i in [0,1]:
#                         other = 2
#                     elif idx in [0, 2] and i in [0,2]:
#                         other = 1
#                     else:
#                         other = 0
#                     other_no = no_s[other]
#                     other_sp = species[other]
#                     product = make_product_name_3(no_camca2c+1,
#                                                   "CaMCa2N",
#                                                   cur_no-1, "CaM",
#                                                   other_no,
#                                                   other_sp,
#                                                   RyR_states)
#                     rxn_idx = add_reaction(root, ryr_specie, product, kf["2N"],
#                                  kr["2N"], rxn_idx)
                    
#                 else:
#                     if cur_no == 2:
#                         product = "4RyR_1CaM_1CaMCa2C_1CaMCa2N_1CaMCa4"
#                         rxn_idx = add_reaction(root, ryr_specie, product,
#                                      kf["2N"],
#                                      kr["2N"], rxn_idx)
                        
#                     else:
#                         inds = other_indices[i]
                      
#                         product = make_product_name_3(1,
#                                                       "CaMCa2N",
#                                                       no_s[inds[0]],
#                                                       species[inds[0]],
#                                                       no_s[inds[1]],
#                                                       species[inds[1]],
#                                                       RyR_states)
#                         rxn_idx = add_reaction(root, ryr_specie, product, kf["2N"],
#                                                kr["2N"], rxn_idx)
#             else:
#                 endswith = cur_specie[-2:]
#                 if "CaMCa4" in species:
#                     idx = species.index("CaMCa4")
#                     no_camca4 = no_s[idx]
#                     if idx in [0, 1] and i in [0,1]:
#                         other = 2
#                     elif idx in [0, 2] and i in [0,2]:
#                         other = 1
#                     else:
#                         other = 0
#                     other_no = no_s[other]
#                     other_sp = species[other]
#                     product = make_product_name_3(cur_no-1, cur_specie,
#                                                   no_camca4+1, "CaMCa4",
#                                                   other_no,
#                                                   other_sp,
#                                                   RyR_states)
#                     rxn_idx = add_reaction(root, ryr_specie, product,
#                                            kf_rev[endswith], kr_rev[endswith],
#                                            rxn_idx)
#                 else:
#                     if cur_no == 2:
#                         product = "4RyR_1CaM_1CaMCa2C_1CaMCa2N_1CaMCa4"
#                     else:

#                         inds = other_indices[i]
#                         product = make_product_name_3(1, "CaMCa4",
#                                                       no_s[inds[0]],
#                                                       species[inds[0]],
#                                                       no_s[inds[1]],
#                                                       species[inds[1]],
#                                                       RyR_states)
#                     rxn_idx = add_reaction(root, ryr_specie, product,
#                                            kf_rev[endswith], kr_rev[endswith],
#                                            rxn_idx)
                        
                    
#     else:
#         # RyR_1CaM_1CaMCa2C_1CaMCa2N_1CaMCa4

#         #  1 RyR_2CaMCa2C_1CaMCa2N_1CaMCa4                        
#         product = make_product_name_3(1, "CaMCa4",
#                                       1, "CaMCa2N",
#                                       2, "CaMCa2C",
#                                       RyR_states)
#         rxn_idx = add_reaction(root, ryr_specie, product,
#                                kf["2C"], kr["2C"],
#                                rxn_idx)

#         #  2 RyR_1CaMCa2C_2CaMCa2N_1CaMCa4
#         product = make_product_name_3(1, "CaMCa4",
#                                       2, "CaMCa2N",
#                                       1, "CaMCa2C",
#                                       RyR_states)
#         rxn_idx = add_reaction(root, ryr_specie, product,
#                                kf["2N"], kr["2N"],
#                                rxn_idx)

#         #  3 RyR_1CaM_1CaMCa2N_2CaMCa4
#         product = make_product_name_3(2, "CaMCa4",
#                                       1, "CaMCa2N",
#                                       1, "CaM",
#                                       RyR_states)
#         rxn_idx = add_reaction(root, ryr_specie, product,
#                                kf["2N"], kr["2N"],
#                                rxn_idx)
#         #  4 RyR_1CaM_1CaMCa2C_2CaMCa4
#         product = make_product_name_3(2, "CaMCa4",
#                                       1, "CaMCa2C",
#                                       1, "CaM",
#                                       RyR_states)
#         rxn_idx = add_reaction(root, ryr_specie, product,
#                                kf["2C"], kr["2C"],
#                                rxn_idx)
        
        
my_file.write('<?xml version="1.0"?>\n')
my_file.write(etree.tostring(root, pretty_print=True).decode('utf-8'))
