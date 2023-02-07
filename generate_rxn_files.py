from lxml import etree
#import argparse
import sys

file_list_full_ER = [                 
    "Rxn_module_Ca.xml",
    "Rxn_module_IP3R.xml",
    "Rxn_module_mGLuR.xml",
    "Rxn_module_RyR.xml",
    "Rxn_module_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",
    "Rxn_module_TG.xml"
]
file_tg =  "Rxn_module_TG.xml"
file_cabuf = "Rxn_module_CaOutBuf.xml"
file_fluo = "Rxn_module_Fluo4.xml"
file_fura ="Rxn_module_Fura2.xml"

file_list_ER_no_mGLur = [
    "Rxn_module_Ca.xml",
    "Rxn_module_IP3R.xml",
    "Rxn_module_IP3.xml",
    "Rxn_module_RyR.xml",
    "Rxn_module_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",

]
file_list_ER_no_mGLur_no_RyR = [
    "Rxn_module_Ca.xml",
    "Rxn_module_IP3R.xml",
    "Rxn_module_IP3.xml",
    "Rxn_module_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",

]
#"For checking if adding SOCE help Ca wave propagation
file_list_ER_no_IP3R = [
    "Rxn_module_Ca.xml",
    "Rxn_module_RyR.xml",
    "Rxn_module_SERCA_ER.xml",
    "Rxn_module_SOCE.xml",
    ]

#for reproducing Ca wave paper
file_list_ER_no_IP3R_no_SOCE = [
    "Rxn_module_Ca.xml",
    "Rxn_module_RyR.xml",
    "Rxn_module_SERCA_ER.xml",
    ]

def read_in_files(flist):
    roots = []
    species = set()
    specie_diff = {}
    rxn_ids = set()
    my_rxn_file = etree.Element("ReactionScheme")

    for fname_xml in flist:
        print(fname_xml)
        tree = etree.parse(fname_xml)
        roots.append(tree.getroot())
        for son in roots[-1]:
            if son.tag == "Specie":
                specie_id = son.attrib["id"]
                specie_kdiff = son.attrib["kdiff"]
                species.add(specie_id)
                if specie_id in specie_diff:
                    if specie_diff[specie_id]!= specie_kdiff:
                        print(specie_id, specie_diff[specie_id], specie_kdiff)
                else:
                    specie_diff[specie_id] = specie_kdiff
                    my_rxn_file.append(son)
            elif son.tag == "Reaction":
                if son.attrib["id"] in rxn_ids:
                    print(son.attrib["id"], " already exists", fname_xml)
                else:
                    rxn_ids.add(son.attrib["id"])
                    my_rxn_file.append(son)
    return my_rxn_file

if __name__ == "__main__":
    
    # 1 no mGluR no RyR
    my_rxn_f = read_in_files(file_list_ER_no_mGLur_no_RyR)
    f = open("Rxn_no_mGLuR_no_RyR.xml", "w")
    f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(file_list_ER_no_mGLur)
    f = open("Rxn_no_mGLuR.xml", "w")
    f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(file_list_ER_no_IP3R)
    f = open("Rxn_no_IP3R.xml", "w")
    f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(file_list_ER_no_IP3R_no_SOCE)
    f = open("Rxn_no_IP3R_no_SOCE.xml", "w")
    f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    
                    
