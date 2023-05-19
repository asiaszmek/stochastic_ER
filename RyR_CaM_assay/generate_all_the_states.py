from lxml import etree
my_file = open("Rxn.xml", "w")
molecules = ["Ca", "RyR", "CaM", "CaMCa2C", "CaMCa2C1N",
             "CaMCa1N", "CaMCa2N", "CaMCa4", "CaER"]
kdiff = {"Ca": "100",
         "CaER": "10",
         "CaM": "4",
         "CaMCa2C": "4",
         "CaMCa1N": "4",
         "CaMCa2N": "4",
         "CaMCa4": "4"}
ryr_name = "RyR"
ryr_cam_binding = ["CaM",  "CaMCa2C", "CaMCa2N", "CaMCa4"]
ca = "Ca"
kf_ryr_cam = 6.66e-5
kr_ryr_cam = 3.6e-03
kf_ryr_camca = 6.66e-4
kr_ryr_camca = 3.6e-03
kf_ryr_2c = 17e-8
kr_ryr_2c = 35e-4
kf_ryr_2n = 6e-9
kr_ryr_2n = 6e-3

def check_specie(specie, RyR_states):
    if specie in RyR_states:
        return specie
    splicik = specie.split("_")
    if len(splicik) > 2:
        return "RyR_%s_%s" % (splicik[2], splicik[1])
    return specie

def fix_cam(specie_list):
    if specie_list[0].endswith("CaM"):
        return specie_list
    if specie_list[1].endswith("CaM"):
        return [specie_list[1], specie_list[0]]

def fix_camca4(specie_list):
    if specie_list[0].endswith("CaMCa4"):
        return specie_list
    if specie_list[1].endswith("CaMCa4"):
        return [specie_list[1], specie_list[0]]



def get_components_ca_reac(specie):
    #  e.g.   "RyR_1CaMCa2N_1CaMCa4"
    components = [specie.split("_")[1], specie.split("_")[2]]
    no_1, no_2 = int(components[0][0]), int(components[1][0])
    if not (no_2 == 1 or no_1 == 1):
        return 
    ifcam = fix_cam(components)
    if ifcam is not None:

        if ifcam[1][0] != "1":
            return
        if ifcam[-1].endswith("2C"):
            return ["RyR_%dCaM"%(int(ifcam[0][0])+1)], [kf_ryr_2c], [kr_ryr_2c]
        elif ifcam[-1].endswith("2N"):
            return ["RyR_%dCaM"%(int(ifcam[0][0])+1)], [kf_ryr_2n], [kr_ryr_2n]
        else:
            return ["RyR_%sCaM_1CaMCa2C"%(ifcam[0][0]),
                    "RyR_%sCaM_1CaMCa2N"%(ifcam[0][0])],\
                    [kf_ryr_2n, kf_ryr_2c], [kr_ryr_2n, kr_ryr_2c]
    
    ifcamca4 = fix_camca4(components)
    if ifcamca4 is None:
        if no_1 == 1 and components[0].endswith("2C"):
            outspecie = ["RyR_1CaM_%s" % components[1]]
            kf_rate = [kf_ryr_2c]
            kr_rate = [kr_ryr_2c]            
            if no_2 == 1:
                outspecie.append("RyR_1CaM_1CaMCa2C")
                kf_rate.append(kf_ryr_2n)
                kr_rate.append(kr_ryr_2n)
            return outspecie, kf_rate, kr_rate
        elif no_1 == 1 and components[0].endswith("2N"):
            outspecie = ["RyR_1CaM_%s" % components[1]]
            kf_rate = [kf_ryr_2n]
            kr_rate = [kr_ryr_2n]            
            if no_2 == 1:
                outspecie.append("RyR_1CaM_1CaMCa2N")
                kf_rate.append(kf_ryr_2c)
                kr_rate.append(kr_ryr_2c)
            return outspecie, kf_rate, kr_rate
        elif no_2 == 1 and components[1].endswith("2C"):
            outspecie = ["RyR_1CaM_%s" % components[0]]
            kf_rate = [kf_ryr_2c]
            kr_rate = [kr_ryr_2c]            
            if no_1 == 1:
                outspecie.append("RyR_1CaM_1CaMCa2C")
                kf_rate.append(kf_ryr_2n)
                kr_rate.append(kr_ryr_2n)
            return outspecie, kf_rate, kr_rate
        elif no_2 == 1 and components[1].endswith("2N"):
            outspecie = ["RyR_1CaM_%s" % components[0]]
            kf_rate = [kf_ryr_2n]
            kr_rate = [kr_ryr_2n]            
            if no_1 == 1:
                outspecie.append("RyR_1CaM_1CaMCa2N")
                kf_rate.append(kf_ryr_2c)
                kr_rate.append(kr_ryr_2c)
            return outspecie, kf_rate, kr_rate
        return
    if ifcamca4[0][0] == "1":
        non_camca4 = int(ifcamca4[1][0])
        non_camca4_sp = ifcamca4[1][1:]
        if non_camca4_sp.endswith("2C"):
            outspecie = ["RyR_%d%s_1CaMCa2N" % (non_camca4,
                                                    non_camca4_sp),
                             "RyR_%d%s" % (non_camca4+1,
                                           non_camca4_sp)]
            kf_rate = [kf_ryr_2c, kf_ryr_2n]
            kr_rate = [kr_ryr_2c, kr_ryr_2n]
            return outspecie, kf_rate, kr_rate
        elif non_camca4_sp.endswith("2N"):
            outspecie = ["RyR_%d%s_1CaMCa2C" % (non_camca4,
                                                non_camca4_sp),
                             "RyR_%d%s" % (non_camca4+1,
                                           non_camca4_sp)]
            kf_rate = [kf_ryr_2n, kf_ryr_2c]
            kr_rate = [kr_ryr_2n, kr_ryr_2c]
            return outspecie, kf_rate, kr_rate

    if ifcamca4[1][0] == "1":
        camca4_sp = ifcamca4[0]
        if ifcamca4[1].endswith("2C"):
            return ["RyR_1CaM_%s"%camca4_sp], [kf_ryr_2c], [kr_ryr_2c]
        elif ifcamca4[1].endswith("2N"):
            return ["RyR_1CaM_%s"%camca4_sp], [kf_ryr_2n], [kr_ryr_2n]
    
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
for i in [1, 2, 3, 4]: #how many CaM molecules can RyR bind
    for j in range(1, i+1):
        k = i - j
        for l, cam_state in enumerate(ryr_cam_binding):
            if k == 0:
                new_name = "RyR_%d%s" % (i, cam_state)
                RyR_states.append(new_name)
                child = etree.SubElement(root, "Specie", name=new_name,
                                         id=new_name, kdiff="0",
                                         kdiffunit="mu2/s")
                    
            else:
                for cam_state2 in ryr_cam_binding[l+1:]:
                    new_name = "RyR_%d%s_%d%s" % (j, cam_state, k, cam_state2)
                    
                    child = etree.SubElement(root, "Specie", name=new_name,
                                             id=new_name, kdiff="0",
                                             kdiffunit="mu2/s")
                    RyR_states.append(new_name)                    
    

rxn_idx = 0

for ryr_state in RyR_states:
    if len(ryr_state.split("_")) == 2:
        my_reac = etree.SubElement(root, "Reaction",
                                       name="RyRCaM%d" % rxn_idx,
                                       id="RyRCaM%d" % rxn_idx)
        rxn_idx += 1

        i = int(ryr_state.split("RyR_")[-1][0])
        specie2 = ryr_state.split("RyR_")[-1][1:]
        if i == 1:
                specie1 = "RyR"
        else:
                specie1 = "RyR_%d%s" % (i-1, specie2)
        etree.SubElement(my_reac, "Reactant", specieID=specie1)
        etree.SubElement(my_reac, "Reactant", specieID=specie2)
        etree.SubElement(my_reac, "Product", specieID=ryr_state)

        if specie2 == "CaM":
            kf = etree.SubElement(my_reac, "forwardRate")
            kf.text = str(kf_ryr_cam)
            kr = etree.SubElement(my_reac, "reverseRate")
            kr.text = str(kr_ryr_cam)
        else:
            kf = etree.SubElement(my_reac, "forwardRate")
            kf.text = str(kf_ryr_camca)
            kr = etree.SubElement(my_reac, "reverseRate")
            kr.text = str(kr_ryr_camca)
    else:
        n1 = int(ryr_state.split("RyR_")[-1][0])
        n2 = int(ryr_state.split("_")[-1][0])
        reac = [ryr_state.split("_")[1][1:] , ryr_state.split("_")[-1][1:]]
        #rxn 1
        for k, pair in enumerate([(n1-1, n2), (n1, n2-1)]):
            my_reac = etree.SubElement(root, "Reaction",
                                       name="RyRCaM%d" % rxn_idx,
                                       id="RyRCaM%d" % rxn_idx)
            rxn_idx += 1
            specie2 = reac[k]
            if pair[0] and pair[1]:
                specie1 = "RyR_%d%s_%d%s" % (pair[0], reac[0], pair[1], reac[1])
                if specie1 not in RyR_states:
                    specie1 = "RyR_%d%s_%d%s" % (pair[1], reac[1],
                                                 pair[0], reac[0])
            elif pair[1]:
                specie1 = "RyR_%d%s" % (pair[1], reac[1])
            else:
                specie1 = "RyR_%d%s" % (pair[0], reac[0])

            etree.SubElement(my_reac, "Reactant", specieID=specie1)
            etree.SubElement(my_reac, "Reactant", specieID=specie2)
            etree.SubElement(my_reac, "Product", specieID=ryr_state)
            if specie2 == "CaM":
                kf = etree.SubElement(my_reac, "forwardRate")
                kf.text = str(kf_ryr_cam)
                kr = etree.SubElement(my_reac, "reverseRate")
                kr.text = str(kr_ryr_cam)
            else:
                kf = etree.SubElement(my_reac, "forwardRate")
                kf.text = str(kf_ryr_camca)
                kr = etree.SubElement(my_reac, "reverseRate")
                kr.text = str(kr_ryr_camca)
 


for ryr_state in RyR_states:
    if len(ryr_state.split("_")) == 2:
      
        end_ca = ryr_state.split("RyR_")[-1][1:]
        if end_ca == "CaM":
            continue #  we take the  endstate and try to reproduce reactions
        n = int(ryr_state.split("RyR_")[-1][0]) - 1
        ending = ryr_state.split("_")[1][1:]
        if end_ca.endswith("2C"):  
            if n:
                species1 = ["RyR_1CaM_%d%s" % (n, ending)]
            else:
                species1 = ["RyR_1CaM"]   
            kf_rate = [kf_ryr_2c]
            kr_rate = [kr_ryr_2c]
        elif end_ca.endswith("2N"):
            if n:
                species1 = ["RyR_1CaM_%d%s" % (n, ending)]
            else:
                species1 = ["RyR_1CaM"]   
            
            kf_rate = [kf_ryr_2n]
            kr_rate = [kr_ryr_2n]
        else:
            if n:
                species1 = ["RyR_1CaMCa2C_%d%s" %(n, ending),
                            "RyR_1CaMCa2N_%d%s" %(n, ending)]
            else:
                species1 = [ryr_state[:-1]+"2C", ryr_state[:-1]+"2N"]
            kf_rate = [kf_ryr_2n, kf_ryr_2c]
            kr_rate = [kr_ryr_2n, kr_ryr_2c]
            
        for k, specie in enumerate(species1):
            my_reac = etree.SubElement(root, "Reaction",
                                       name="RyRCaM%d" % rxn_idx,
                                       id="RyRCaM%d" % rxn_idx)
            rxn_idx += 1
            new_specie = check_specie(specie, RyR_states)
            etree.SubElement(my_reac, "Reactant", specieID=new_specie)
            etree.SubElement(my_reac, "Reactant", specieID="Ca", power="2")
            etree.SubElement(my_reac, "Product", specieID=ryr_state)
            kf = etree.SubElement(my_reac, "forwardRate")
            kf.text = str(kf_rate[k])
            kr = etree.SubElement(my_reac, "reverseRate")
            kr.text = str(kr_rate[k])
    else:
        #  e.g.   "RyR_1CaMCa2N_1CaMCa4"
        out = get_components_ca_reac(ryr_state)
        if out is not None:
            out_species = out[0]
            kf_rate = out[1]
            kr_rate = out[2]
            for k, specie in enumerate(out_species):
                my_reac = etree.SubElement(root, "Reaction",
                                           name="RyRCaM%d" % rxn_idx,
                                           id="RyRCaM%d" % rxn_idx)
                rxn_idx += 1
                
                splicik = specie.split("_")
                if len(splicik) > 2:
                    if splicik[1][1:] == splicik[2][1:]:
                        n = int(splicik[1][0])+int(splicik[2][0])
                        specie = "RyR_%d%s" % (n, splicik[1][1:])
                new_specie = check_specie(specie, RyR_states)
                etree.SubElement(my_reac, "Reactant", specieID=new_specie)
                    
                etree.SubElement(my_reac, "Reactant", specieID="Ca", power="2")
                etree.SubElement(my_reac, "Product", specieID=ryr_state)
                kf = etree.SubElement(my_reac, "forwardRate")
                kf.text = str(kf_rate[k])
                kr = etree.SubElement(my_reac, "reverseRate")
                kr.text = str(kr_rate[k])
     
#last reactions:
cam = etree.SubElement(root,"Reaction", name="CaMC_bind", id="CaMC_bind")
etree.SubElement(cam, "Reactant", specieID="CaM")
etree.SubElement(cam, "Reactant", specieID="Ca", power="2")
etree.SubElement(cam, "Product", specieID="CaMCa2C",)
kf = etree.SubElement(cam, "forwardRate")
kf.text = str(17e-10)
kr = etree.SubElement(cam, "reverseRate")
kr.text = str(35e-4)

cam = etree.SubElement(root,"Reaction", name="CaM2C_bind", id="CaM2C_bind")
etree.SubElement(cam, "Reactant", specieID="CaMCa2C")
etree.SubElement(cam, "Reactant", specieID="Ca")
etree.SubElement(cam, "Product", specieID="CaMCa2C1N",)
kf = etree.SubElement(cam, "forwardRate")
kf.text = str(13.5e-6)
kr = etree.SubElement(cam, "reverseRate")
kr.text = str(228.3e-3)

cam = etree.SubElement(root,"Reaction", name="CaM2C_bind", id="CaM2C_bind2")
etree.SubElement(cam, "Reactant", specieID="CaMCa2C1N")
etree.SubElement(cam, "Reactant", specieID="Ca")
etree.SubElement(cam, "Product", specieID="CaMCa4",)
kf = etree.SubElement(cam, "forwardRate")
kf.text = str(26.3e-6)
kr = etree.SubElement(cam, "reverseRate")
kr.text = str(64e-3)

cam = etree.SubElement(root,"Reaction", name="CaM_bind", id="CaM_bind")
etree.SubElement(cam, "Reactant", specieID="CaM")
etree.SubElement(cam, "Reactant", specieID="Ca")
etree.SubElement(cam, "Product", specieID="CaMCa1N",)
kf = etree.SubElement(cam, "forwardRate")
kf.text = str(13.5e-6)
kr = etree.SubElement(cam, "reverseRate")
kr.text = str(228.3e-3)

cam = etree.SubElement(root,"Reaction", name="CaM_bind2", id="CaM_bind2")
etree.SubElement(cam, "Reactant", specieID="CaMCa1N")
etree.SubElement(cam, "Reactant", specieID="Ca")
etree.SubElement(cam, "Product", specieID="CaMCa2N",)
kf = etree.SubElement(cam, "forwardRate")
kf.text = str(26.3e-6)
kr = etree.SubElement(cam, "reverseRate")
kr.text = str(64e-3)

cam = etree.SubElement(root,"Reaction", name="CaM2N_bind", id="CaM2N_bind")
etree.SubElement(cam, "Reactant", specieID="CaMCa2N")
etree.SubElement(cam, "Reactant", specieID="Ca", power="2")
etree.SubElement(cam, "Product", specieID="CaMCa4",)
kf = etree.SubElement(cam, "forwardRate")
kf.text = str(17e-10)
kr = etree.SubElement(cam, "reverseRate")
kr.text = str(35e-4)
my_file.write('<?xml version="1.0"?>\n')
my_file.write(etree.tostring(root, pretty_print=True).decode('utf-8'))
